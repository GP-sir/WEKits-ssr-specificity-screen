import re
import argparse
from collections import defaultdict

try:
    import openpyxl
    HAS_XLSX = True
except ImportError:
    HAS_XLSX = False

TIER_RANK = {'Top_Tier': 0, 'High_Tier': 1, 'excluded': 2}

NAME_RE = re.compile(r'^(?P<locus>.+?)_p\d+(?:_P\d+)?$')
CHR_RE = re.compile(r'Chr(\d{1,2}[ABD])', re.I)


def parse_line(line):
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 10:
        fields = line.split()
    if len(fields) < 10:
        return None
    qseqid = fields[0].strip()
    sseqid = fields[1].strip()
    nums = []
    for f in fields[2:]:
        f = f.strip()
        try:
            nums.append(float(f) if '.' in f else int(f))
        except ValueError:
            nums.append(None)
    big = [x for x in nums if isinstance(x, int) and x > 1000]
    if len(big) < 2:
        return None
    sstart, send = min(big[0], big[1]), max(big[0], big[1])
    m = CHR_RE.search(sseqid)
    chrom = m.group(1).upper() if m else sseqid
    return dict(qseqid=qseqid, chrom=chrom, sstart=sstart, send=send)


def classify_primer(n_loci, unique_target):
    if n_loci == 1 and unique_target:
        return 'Gold'
    if n_loci == 2:
        return 'Silver'
    if n_loci == 1:
        return 'unique_off_target'
    return 'multi'


def split_primer(q):
    if q.endswith('_F'):
        return q[:-2], 'F'
    if q.endswith('_R'):
        return q[:-2], 'R'
    return q, '?'


def locus_of(base):
    m = NAME_RE.match(base)
    return m.group('locus') if m else base


def fmt(st):
    return '; '.join(f'{c}:{p}' for c, p in st['loci'])


def main():
    ap = argparse.ArgumentParser(
        description='Assign primer-pair tiers per SSR locus from BLAST output.')
    ap.add_argument('inputs', nargs='+', help='BLAST tabular output (one or more files)')
    ap.add_argument('-o', '--output', default='locus_screened.xlsx',
                    help='output workbook (.xlsx) or csv prefix')
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

    rows, top_combos = [], []
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
                if tier in ('Top_Tier', 'High_Tier'):
                    top_combos.append(dict(
                        Locus=loc, F_primer=fq, R_primer=rq,
                        F_class=fc, R_class=rc,
                        F_loci=fmt(fst), R_loci=fmt(rst), Tier=tier,
                        Amplicon_bp=size if size is not None else ''))
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

    cols = ['Locus', 'F_primer', 'R_primer', 'F_class', 'R_class',
            'F_loci', 'R_loci', 'Tier', 'Amplicon_bp', 'n_F', 'n_R']
    tc_cols = ['Locus', 'F_primer', 'R_primer', 'F_class', 'R_class',
               'F_loci', 'R_loci', 'Tier', 'Amplicon_bp']
    if HAS_XLSX and args.output.endswith('.xlsx'):
        wb = openpyxl.Workbook()
        ws1 = wb.active
        ws1.title = 'locus_result'
        ws1.append(cols)
        for r in rows:
            ws1.append([r.get(c, '') for c in cols])
        ws2 = wb.create_sheet('top_combos')
        ws2.append(tc_cols)
        for r in sorted(top_combos, key=lambda r: (TIER_RANK[r['Tier']], r['Locus'])):
            ws2.append([r.get(c, '') for c in tc_cols])
        ws3 = wb.create_sheet('primer_detail')
        ws3.append(['Primer_ID', 'n_loci', 'class', 'loci'])
        for q, st in sorted(primer_stat.items()):
            ws3.append([q, st['n_loci'], st['cls'], fmt(st)])
        ws4 = wb.create_sheet('summary')
        for k, v in summary:
            ws4.append([k, v])
        wb.save(args.output)
        print(f'\nSaved: {args.output}')
    else:
        base = args.output.rsplit('.', 1)[0]
        import csv
        with open(base + '_locus_result.csv', 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for r in rows:
                w.writerow([r.get(c, '') for c in cols])
        print(f'\nSaved: {base}_locus_result.csv (output written as csv)')


if __name__ == '__main__':
    main()
