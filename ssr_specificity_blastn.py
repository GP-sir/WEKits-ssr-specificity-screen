import re
import csv
import argparse
from collections import defaultdict

try:
    import openpyxl
    HAS_XLSX = True
except ImportError:
    HAS_XLSX = False

TIER_RANK = {'Top_Tier': 0, 'High_Tier': 1, 'excluded': 2}

CHR_RE = re.compile(r'Chr(\d{1,2}[ABD])', re.I)
DIRECTION_RE = re.compile(r'[-_]([FR])(\d*)$', re.I)
LOCUS_TAIL_RE = re.compile(r'(?:[-_][pPfrFR]\d*)+$')


def parse_line(line):
    tok = line.split()
    if len(tok) < 11:
        return None
    qseqid, sseqid = tok[0], tok[1]
    try:
        sstart, send = int(tok[9]), int(tok[10])
    except ValueError:
        big = sorted((int(t) for t in tok[2:]
                      if t.lstrip('-').isdigit() and abs(int(t)) > 1000),
                     reverse=True)[:2]
        if len(big) < 2:
            return None
        sstart, send = big[0], big[1]
    m = CHR_RE.search(sseqid)
    chrom = m.group(1).upper() if m else sseqid
    return dict(qseqid=qseqid, chrom=chrom,
                sstart=min(sstart, send), send=max(sstart, send))


def classify_primer(n_loci, unique_target):
    if n_loci == 1 and unique_target:
        return 'Gold'
    if n_loci == 2:
        return 'Silver'
    if n_loci == 1:
        return 'unique_off_target'
    return 'multi'


def split_primer(q):
    m = DIRECTION_RE.search(q)
    if m:
        return q[:m.start()], m.group(1).upper()
    return q, '?'


def locus_of(base):
    return LOCUS_TAIL_RE.sub('', base)


def fmt(st):
    return '; '.join(f'{c}:{p}' for c, p in st['loci'])


def main():
    ap = argparse.ArgumentParser(
        description='Assign primer-pair tiers per SSR locus from BLAST output.')
    ap.add_argument('inputs', nargs='+', help='BLAST tabular output (one or more files)')
    ap.add_argument('-o', '--output', default='final_result.tsv',
                    help='output file for the 8-column summary (.tsv or .xlsx)')
    ap.add_argument('--chr', default='4A', help='target chromosome (default 4A)')
    ap.add_argument('--interval', nargs=2, type=int,
                    default=[520000000, 530000000], metavar=('START', 'END'),
                    help='target interval in bp (default 520000000 530000000)')
    ap.add_argument('--scope', choices=['interval', 'chromosome'],
                    default='interval',
                    help='interval: Gold requires a single hit inside the target '
                         'interval; chromosome: inside the target chromosome')
    args = ap.parse_args()

    tgt_chr = args.chr.upper()
    iv_start, iv_end = args.interval

    hits = defaultdict(set)
    n_records = 0
    for path in args.inputs:
        with open(path, encoding='utf-8') as fh:
            for line in fh:
                if not line.strip():
                    continue
                rec = parse_line(line)
                if rec is None:
                    continue
                n_records += 1
                hits[rec['qseqid']].add((rec['chrom'], rec['sstart']))

    primer_stat = {}
    for q, loci in hits.items():
        n = len(loci)
        if args.scope == 'chromosome':
            tgt = sorted(p for c, p in loci if c == tgt_chr)
        else:
            tgt = sorted(p for c, p in loci if c == tgt_chr and iv_start <= p <= iv_end)
        ut = n == 1 and bool(tgt)
        primer_stat[q] = dict(n_loci=n, cls=classify_primer(n, ut),
                              loci=sorted(loci), tgt=tgt)

    loci_F = defaultdict(list)
    loci_R = defaultdict(list)
    n_q = 0
    for q, st in primer_stat.items():
        base, side = split_primer(q)
        loc = locus_of(base)
        n_q += 1
        if side == 'F':
            loci_F[loc].append((q, st))
        elif side == 'R':
            loci_R[loc].append((q, st))

    rows = []
    n_top = n_high = 0
    for loc in sorted(set(loci_F) | set(loci_R)):
        Fs, Rs = loci_F.get(loc, []), loci_R.get(loc, [])
        best = None
        for fq, fst in Fs:
            for rq, rst in Rs:
                fc, rc = fst['cls'], rst['cls']
                if fc == 'Gold' and rc == 'Gold':
                    tier = 'Top_Tier'
                elif {fc, rc} == {'Gold', 'Silver'}:
                    tier = 'High_Tier'
                else:
                    tier = 'excluded'
                size = None
                if fst['tgt'] and rst['tgt']:
                    size = min(abs(r - f) for f in fst['tgt'] for r in rst['tgt'])
                size_pen = 0 if (size is not None and 100 <= size <= 500) else 1
                key = (TIER_RANK[tier],
                       fst['n_loci'] + rst['n_loci'],
                       size_pen,
                       abs(size - 300) if size is not None else 10 ** 9,
                       fq, rq)
                if best is None or key < best[0]:
                    best = (key, dict(
                        Locus=loc, F_primer=fq, R_primer=rq,
                        F_class=fc, R_class=rc,
                        F_loci=fmt(fst), R_loci=fmt(rst),
                        Tier=tier,
                        Amplicon_bp=size if size is not None else '',
                        n_F=len(Fs), n_R=len(Rs)))
        if best:
            rows.append(best[1])
            if best[1]['Tier'] == 'Top_Tier':
                n_top += 1
            elif best[1]['Tier'] == 'High_Tier':
                n_high += 1
    rows.sort(key=lambda r: (TIER_RANK[r['Tier']], r['Locus']))

    n_F_total = sum(len(v) for v in loci_F.values())
    n_R_total = sum(len(v) for v in loci_R.values())
    n_design = len(set(split_primer(q)[0] for q in primer_stat))
    summary = [
        ('Target', f'{tgt_chr}:{iv_start}-{iv_end}'),
        ('Scope', args.scope),
        ('BLAST records parsed', n_records),
        ('Primers (F+R)', n_q),
        ('  - F primers', n_F_total),
        ('  - R primers', n_R_total),
        ('Primer pairs as designed', n_design),
        ('SSR loci', len(rows)),
        ('Top_Tier loci', n_top),
        ('High_Tier loci', n_high),
        ('Specific loci (Top+High)', n_top + n_high),
        ('Loci without qualifying combo', len(rows) - n_top - n_high),
    ]
    print('\n===== SUMMARY =====')
    for k, v in summary:
        print(f'{k:34s} {v}')

    result_header = ['Locus', 'Primer_F', 'Primer_R', 'F_class', 'R_class',
                     'F_loci', 'R_loci', 'Tier']
    result_rows = [[r['Locus'], r['F_primer'], r['R_primer'],
                    r['F_class'], r['R_class'], r['F_loci'],
                    r['R_loci'], r['Tier']]
                   for r in rows if r['Tier'] in ('Top_Tier', 'High_Tier')]

    out_path = args.output
    if out_path.lower().endswith('.xlsx') and HAS_XLSX:
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = 'final_result'
        ws.append(result_header)
        for r in result_rows:
            ws.append(r)
        wb.save(out_path)
    else:
        if out_path.lower().endswith('.xlsx'):
            out_path = out_path[:-5] + '.tsv'
        with open(out_path, 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(result_header)
            w.writerows(result_rows)
    print(f'\nSaved: {out_path}')


if __name__ == '__main__':
    main()
