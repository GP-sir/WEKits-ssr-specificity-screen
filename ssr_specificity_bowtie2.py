import re
import csv
import sys
import argparse
import subprocess
from collections import defaultdict
from pathlib import Path

try:
    import openpyxl
    HAS_XLSX = True
except ImportError:
    HAS_XLSX = False

DEFAULT_INDEX = '/opt/ref/wheat/cs1.1/database/wheat_bt2'
DEFAULT_MAX_PRODUCT_SIZE = 300
CHR_HEAD_RE = re.compile(r'^(?:[Cc]hr)?(?P<chr>[1-7][ABD])_')
DIRECTION_RE = re.compile(r'[-_](?P<dir>[FR])(?P<num>\d*)$', re.I)
DESIGN_TAIL_RE = re.compile(r'[-_][pP]\d+$')
RANGE_RE = re.compile(r'_(?P<start>\d+)-(?P<end>\d+)_(?P<motif>[^_]+)$')
POS_RE = re.compile(r'_(?P<motif>[^_]+)_(?P<pos>\d+)$')
QID_RE = re.compile(r'^Q([0-9]+)$')


def parse_header(header):
    dm = DIRECTION_RE.search(header)
    if not dm:
        raise SystemExit('Unrecognized FASTA header (no F/R direction): ' + header)
    direction = dm.group('dir').upper()
    design = dm.group('num')
    core = header[:dm.start()]

    pm = DESIGN_TAIL_RE.search(core)
    if pm:
        if not design:
            design = pm.group()[2:]
        while True:
            m = DESIGN_TAIL_RE.search(core)
            if not m:
                break
            core = core[:m.start()]
    if not design:
        design = '0'

    cm = CHR_HEAD_RE.match(core)
    if not cm:
        raise SystemExit('Unrecognized FASTA header (no chromosome): ' + header)
    chrom = cm.group('chr').upper()

    motif, pos = 'NA', None
    rm = RANGE_RE.search(core)
    if rm:
        motif = rm.group('motif')
        pos = int(rm.group('start'))
    else:
        pm2 = POS_RE.search(core)
        if pm2:
            motif = pm2.group('motif')
            pos = int(pm2.group('pos'))
    return dict(chrom=chrom, site=core, motif=motif, pos=pos,
                direction=direction, design=int(design))


def parse_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                if header is None:
                    raise SystemExit('Sequence before FASTA header in ' + str(path))
                chunks.append(line)
    if header is not None:
        yield header, ''.join(chunks)


def collect_inputs(items):
    files = []
    for item in items:
        p = Path(item)
        if p.is_dir():
            files.extend(sorted(p.glob('*.fasta')))
        elif p.is_file():
            files.append(p)
        else:
            raise SystemExit('Input not found: ' + item)
    if not files:
        raise SystemExit('No FASTA files found.')
    return files


CHR_TOKEN_RE = re.compile(r'(\d[ABD])', re.I)


def norm_chr(value):
    m = CHR_TOKEN_RE.search(str(value or ''))
    return m.group(1).upper() if m else str(value or '').strip()


def grade(count):
    return 'Unmapped' if count == 0 else ('Gold' if count == 1 else 'Silver')


def main():
    ap = argparse.ArgumentParser(
        description='Whole-genome SSR primer specificity screen (perfect Bowtie2 matches).')
    ap.add_argument('fasta', nargs='+', help='FASTA files or directories containing *.fasta')
    ap.add_argument('-x', '--index', default=DEFAULT_INDEX,
                    help='bowtie2 index prefix (default: %(default)s)')
    ap.add_argument('--sam', help='use an existing SAM file instead of running bowtie2')
    ap.add_argument('-o', '--output', default='final_result.tsv',
                    help='output file for the 8-column summary (.tsv or .xlsx)')
    ap.add_argument('-p', '--threads', type=int, default=8,
                    help='bowtie2 threads (default 8)')
    ap.add_argument('-k', type=int, default=2,
                    help='max alignments reported per read (default 2)')
    ap.add_argument('--interval', nargs=3, metavar=('CHR', 'START', 'END'),
                    help='count only loci whose motif position is in CHR:START-END')
    ap.add_argument('--max-product-size', type=int, default=DEFAULT_MAX_PRODUCT_SIZE,
                    help='max amplicon size for selecting one best pair per locus '
                         '(default 300)')
    ap.add_argument('--keep-sam', action='store_true',
                    help='save intermediate SAM when bowtie2 is run')
    args = ap.parse_args()

    primers = []
    seq_to_id = {}
    id_to_seq = {}
    min_len = None
    for path in collect_inputs(args.fasta):
        source = Path(path).name
        for header, seq in parse_fasta(path):
            info = parse_header(header)
            seq = seq.upper()
            if not seq or re.fullmatch(r'[ACGT]+', seq) is None:
                raise SystemExit('Invalid primer sequence: ' + header)
            if seq not in seq_to_id:
                seq_to_id[seq] = len(seq_to_id) + 1
            sid = seq_to_id[seq]
            id_to_seq[sid] = seq
            min_len = len(seq) if min_len is None else min(min_len, len(seq))
            primers.append(dict(
                id=header, site=info['site'], chrom=info['chrom'],
                motif=info['motif'], pos=info['pos'],
                direction=info['direction'], design=info['design'],
                seq=seq, sid=sid, source=source))
    if not primers:
        raise SystemExit('No primer records found.')

    site_source = {}
    for p in primers:
        site_source.setdefault(p['site'], p['source'])

    seed_len = min(15, min_len)
    if args.sam:
        with open(args.sam, encoding='utf-8') as fh:
            sam_text = fh.read()
    else:
        import tempfile
        with tempfile.NamedTemporaryFile('w', suffix='.fasta', delete=False) as tf:
            for sid in sorted(id_to_seq):
                tf.write(f'>Q{sid}\n{id_to_seq[sid]}\n')
            uniq_fasta = tf.name
        cmd = ['bowtie2', '-f', '--end-to-end', '-N', '0', '-L', str(seed_len),
               '-D', '25', '-R', '5', '-i', 'S,1,0.50', '--score-min', 'L,0,0',
               '--mp', '6,6', '--np', '1', '--rdg', '5,3', '--rfg', '5,3',
               '-k', str(args.k), '--no-unal', '--threads', str(args.threads),
               '-x', args.index, '-U', uniq_fasta]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
        sys.stderr.write(proc.stderr)
        if proc.returncode != 0:
            raise SystemExit('bowtie2 failed with exit code ' + str(proc.returncode))
        sam_text = proc.stdout
        if args.keep_sam:
            with open('primers.sam', 'w', encoding='utf-8') as fh:
                fh.write(sam_text)
        import os
        os.remove(uniq_fasta)

    alignments = defaultdict(set)
    for line in sam_text.splitlines():
        if not line or line.startswith('@'):
            continue
        f = line.split('\t')
        if len(f) < 11:
            continue
        qm = QID_RE.fullmatch(f[0])
        if not qm:
            raise SystemExit('Unexpected QNAME in SAM: ' + f[0])
        sid = int(qm.group(1))
        flag = int(f[1])
        if flag & 4:
            continue
        seq = id_to_seq.get(sid)
        if seq is None:
            continue
        chrom = f[2]
        start = int(f[3])
        strand = '-' if flag & 16 else '+'
        cigar = f[5]
        nm = as_ = None
        for tag in f[11:]:
            if tag.startswith('NM:i:'):
                nm = int(tag[5:])
            elif tag.startswith('AS:i:'):
                as_ = int(tag[5:])
        if nm != 0 or as_ != 0 or start < 1 or cigar != f'{len(seq)}M':
            continue
        end = start + len(seq) - 1
        alignments[sid].add((chrom, start, end, strand, cigar))

    hit_count = {sid: len(v) for sid, v in alignments.items()}
    tier_of = {sid: grade(hit_count.get(sid, 0)) for sid in id_to_seq}

    by_site = defaultdict(lambda: {'F': [], 'R': []})
    for p in primers:
        sid = p['sid']
        by_site[p['site']][p['direction']].append(dict(
            **p, count=hit_count.get(sid, 0), tier=tier_of[sid],
            hits=sorted(alignments.get(sid, []))))

    primer_class_rows = []
    alignment_rows = []
    for p in sorted(primers, key=lambda x: (x['chrom'], x['site'], x['design'], x['direction'])):
        sid = p['sid']
        cnt = hit_count.get(sid, 0)
        primer_class_rows.append([p['id'], p['site'], p['chrom'], p['design'],
                                  p['direction'], p['seq'], len(p['seq']), cnt,
                                  tier_of[sid]])
        hits = sorted(alignments.get(sid, []))
        if not hits:
            alignment_rows.append([f'Q{sid}', p['id'], p['site'], p['chrom'],
                                   p['design'], p['direction'], p['seq'], len(p['seq']),
                                   cnt, tier_of[sid], 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'])
        else:
            for h in hits:
                alignment_rows.append([f'Q{sid}', p['id'], p['site'], p['chrom'],
                                       p['design'], p['direction'], p['seq'], len(p['seq']),
                                       cnt, tier_of[sid], h[0], h[1], h[2], h[3],
                                       h[4], 0, 0])

    pair_rows = []
    warning_rows = []
    final_pairs = set()
    pair_tier_counts = defaultdict(int)
    combination_count = 0
    warning_count = 0
    site_stats = []
    sid_by_id = {p['id']: p['sid'] for p in primers}

    for site in sorted(by_site):
        group = by_site[site]
        fl, rl = group['F'], group['R']
        site_counts = defaultdict(int)
        if not fl or not rl:
            warning_count += 1
            warning_rows.append([site, 'NA', 'NA', 'NA',
                                 'SSR site has no F candidates or no R candidates'])
            site_stats.append([site, 0, 0, 0, 0, 0, 'Unmapped_or_incomplete', 0])
            continue
        motif_pos = next((p['pos'] for p in fl + rl if p['pos'] is not None), None)
        design_chr = fl[0]['chrom'] if fl else rl[0]['chrom']
        for f in fl:
            for r in rl:
                combination_count += 1
                if f['tier'] == 'Gold' and r['tier'] == 'Gold':
                    pt = 'Top_Tier'
                elif (f['tier'], r['tier']) in {('Gold', 'Silver'), ('Silver', 'Gold')}:
                    pt = 'High_Tier'
                elif f['tier'] == 'Silver' and r['tier'] == 'Silver':
                    pt = 'Low_Tier'
                else:
                    pt = 'Unmapped_or_incomplete'
                pair_tier_counts[pt] += 1
                site_counts[pt] += 1
                pair_rows.append([site, f"F{f['design']}+R{r['design']}", f['id'], r['id'],
                                  f['tier'], r['tier'], f['count'], r['count'], pt])
                if pt not in {'Top_Tier', 'High_Tier'}:
                    continue
                f_hits = [h for h in f['hits'] if norm_chr(h[0]) == design_chr and h[3] == '+']
                r_hits = [h for h in r['hits'] if norm_chr(h[0]) == design_chr and h[3] == '-']
                amps = set()
                for fh in f_hits:
                    for rh in r_hits:
                        if fh[2] >= rh[1]:
                            continue
                        if motif_pos is not None and not (fh[1] <= motif_pos <= rh[2]):
                            continue
                        product = rh[2] - fh[1] + 1
                        if product <= 0:
                            continue
                        amps.add((fh[1], rh[1], product))
                if not amps:
                    warning_count += 1
                    warning_rows.append([site, f['id'], r['id'], pt,
                                         'No exact inward-facing F(+)/R(-) amplicon on the '
                                         'target chromosome spanning the SSR coordinate'])
                    continue
                for fp, rp, size in sorted(amps):
                    final_pairs.add((site, motif_pos, design_chr, pt,
                                     f['id'], f['seq'], r['id'], r['seq'],
                                     size, fp, rp, f['design'], r['design']))
        total = sum(site_counts.values())
        best = next((x for x in ('Top_Tier', 'High_Tier', 'Low_Tier',
                                 'Unmapped_or_incomplete') if site_counts[x]),
                    'Unmapped_or_incomplete')
        retained = sum(1 for row in final_pairs if row[0] == site)
        site_stats.append([site, total, site_counts['Top_Tier'], site_counts['High_Tier'],
                           site_counts['Low_Tier'], site_counts['Unmapped_or_incomplete'],
                           best, retained])

    def tm(seq):
        return f'{2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C")):.2f}'

    def gc(seq):
        return f'{100.0 * (seq.count("G") + seq.count("C")) / len(seq):.2f}'

    final_sorted = sorted(final_pairs, key=lambda x: (
        x[2], x[1] if x[1] is not None else 0,
        0 if x[3] == 'Top_Tier' else 1, x[11], x[12], x[9], x[10]))
    final_rows = []
    for row in final_sorted:
        (site, pos, chrom, pt, fid, fseq, rid, rseq,
         size, fp, rp, fd, rd) = row
        final_rows.append([site, pos, chrom, pt, fid, fseq, tm(fseq), gc(fseq),
                           tier_of[sid_by_id[fid]], rid, rseq, tm(rseq), gc(rseq),
                           tier_of[sid_by_id[rid]], size, fp, rp])

    amp_by_site = defaultdict(lambda: {'Top_Tier': 0, 'High_Tier': 0, 'best': None})
    for row in final_sorted:
        s, pt = row[0], row[3]
        amp_by_site[s][pt] += 1
        if amp_by_site[s]['best'] is None:
            amp_by_site[s]['best'] = pt

    tier_rank = {'Top_Tier': 0, 'High_Tier': 1, 'Low_Tier': 2,
                 'Unmapped_or_incomplete': 3}
    best_amp_first = {}
    for row in final_sorted:
        if row[8] > args.max_product_size:
            continue
        best_amp_first.setdefault(row[0], row)
    over_size_skipped = sum(1 for row in final_sorted if row[8] > args.max_product_size)

    best_combo_by_site = {}
    for row in pair_rows:
        site = row[0]
        rank = tier_rank.get(row[8], 9)
        cur = best_combo_by_site.get(site)
        if cur is None or rank < cur[0]:
            best_combo_by_site[site] = (rank, row)

    site_meta = {}
    for p in primers:
        site_meta.setdefault(p['site'], (p['chrom'], p['pos']))
    site_stats_by_site = {r[0]: r for r in site_stats}

    def in_interval(site):
        if not args.interval:
            return True
        ichr, istart, iend = args.interval[0], int(args.interval[1]), int(args.interval[2])
        chrom, pos = site_meta.get(site, ('NA', None))
        if norm_chr(chrom) != norm_chr(ichr):
            return False
        if pos is None:
            return True
        return istart <= pos <= iend

    locus_best_rows = []
    for site in sorted(by_site):
        if not in_interval(site):
            continue
        st = site_stats_by_site.get(site)
        amp = amp_by_site.get(site)
        chrom, pos = site_meta.get(site, ('NA', None))
        total = st[1] if st else 0
        comb_best = st[6] if st else 'Unmapped_or_incomplete'
        amp_top = amp['Top_Tier'] if amp else 0
        amp_high = amp['High_Tier'] if amp else 0
        sel = best_amp_first.get(site)
        if sel is not None:
            amp_best = sel[3]
            fid, rid, size = sel[4], sel[6], sel[8]
        else:
            amp_best = 'No_amplicon'
            fid = rid = size = 'NA'
        locus_best_rows.append([
            site_source.get(site, 'NA'), site, chrom,
            pos if pos is not None else 'NA', total, comb_best,
            amp_top, amp_high, amp_best, fid, rid, size])

    amp_rows_by_fasta = defaultdict(lambda: {'Top_Tier': 0, 'High_Tier': 0})
    for r in final_rows:
        if not in_interval(r[0]):
            continue
        src = site_source.get(r[0], 'NA')
        if r[3] in amp_rows_by_fasta[src]:
            amp_rows_by_fasta[src][r[3]] += 1

    per_fasta_rows = []
    for src in sorted({p['source'] for p in primers}):
        rs = [r for r in locus_best_rows if r[0] == src]
        n = len(rs)
        n_top = sum(1 for r in rs if r[8] == 'Top_Tier')
        n_high = sum(1 for r in rs if r[8] == 'High_Tier')
        n_none = n - n_top - n_high
        c_top = sum(1 for r in rs if r[5] == 'Top_Tier')
        c_high = sum(1 for r in rs if r[5] == 'High_Tier')
        per_fasta_rows.append([
            src, n, n_top, n_high, n_top + n_high, n_none,
            f'{100.0 * n_top / n:.2f}' if n else '0.00',
            f'{100.0 * n_high / n:.2f}' if n else '0.00',
            f'{100.0 * (n_top + n_high) / n:.2f}' if n else '0.00',
            c_top, c_high,
            amp_rows_by_fasta[src]['Top_Tier'],
            amp_rows_by_fasta[src]['High_Tier']])

    best_pair_header = [
        'FASTA', 'SSR_ID', 'Chr', 'Motif_Pos', 'Status', 'Best_Combination_Tier',
        'Selected_F_ID', 'F_Sequence', 'F_Length', 'F_Tm', 'F_GC', 'F_Grade', 'F_Hit_Count',
        'Selected_R_ID', 'R_Sequence', 'R_Length', 'R_Tm', 'R_GC', 'R_Grade', 'R_Hit_Count',
        'Product_Size', 'F_Pos', 'R_Pos', 'Design',
        'Total_F_R_Combinations', 'Total_Valid_Amplicons']
    best_pair_rows = []
    for site in sorted(by_site):
        if not in_interval(site):
            continue
        st = site_stats_by_site.get(site)
        amp = amp_by_site.get(site)
        total_combo = st[1] if st else 0
        n_amp = (amp['Top_Tier'] + amp['High_Tier']) if amp else 0
        chrom, pos = site_meta.get(site, ('NA', None))
        pos_out = pos if pos is not None else 'NA'
        brow = best_amp_first.get(site)
        if brow is None:
            continue
        (s, p, c, pt, fid, fseq, rid, rseq,
         size, fp, rp, fd, rd) = brow
        status = pt
        comb_tier = pt
        design = f'F{fd}+R{rd}'
        fg = tier_of[sid_by_id[fid]]
        rg = tier_of[sid_by_id[rid]]
        fcnt = hit_count.get(sid_by_id[fid], 0)
        rcnt = hit_count.get(sid_by_id[rid], 0)
        best_pair_rows.append([
            site_source.get(site, 'NA'), s, chrom, pos_out, status, comb_tier,
            fid, fseq, len(fseq), tm(fseq), gc(fseq), fg, fcnt,
            rid, rseq, len(rseq), tm(rseq), gc(rseq), rg, rcnt,
            size, fp, rp, design, total_combo, n_amp])

    final_top = sum(1 for r in final_rows if r[3] == 'Top_Tier')
    final_high = sum(1 for r in final_rows if r[3] == 'High_Tier')
    gold_n = sum(1 for p in primers if tier_of[p['sid']] == 'Gold')
    silver_n = sum(1 for p in primers if tier_of[p['sid']] == 'Silver')
    unmapped_n = sum(1 for p in primers if tier_of[p['sid']] == 'Unmapped')
    n_top_loci = sum(1 for r in best_pair_rows if r[4] == 'Top_Tier')
    n_high_loci = sum(1 for r in best_pair_rows if r[4] == 'High_Tier')

    summary = [
        ('Target', 'ALL' if not args.interval else
         f'{args.interval[0]}:{args.interval[1]}-{args.interval[2]}'),
        ('Max product size', args.max_product_size),
        ('Original primer records', len(primers)),
        ('Unique primer sequences', len(id_to_seq)),
        ('Gold primers', gold_n),
        ('Silver primers', silver_n),
        ('Unmapped primers', unmapped_n),
        ('SSR loci with a best pair', len(best_pair_rows)),
        ('Top_Tier loci', n_top_loci),
        ('High_Tier loci', n_high_loci),
        ('Specific loci (Top+High)', n_top_loci + n_high_loci),
    ]
    print('\n===== SUMMARY =====')
    for k, v in summary:
        print(f'{k:34s} {v}')

    rank_out = {'Top_Tier': 0, 'High_Tier': 1}

    def hits_str(pid):
        sid = sid_by_id[pid]
        return '; '.join(f'{h[0]}:{h[1]}' for h in sorted(alignments.get(sid, [])))

    result_header = ['Locus', 'Primer_F', 'Primer_R', 'F_class', 'R_class',
                     'F_loci', 'R_loci', 'Tier']
    result_rows = [[b[1], b[6], b[13], b[11], b[18],
                    hits_str(b[6]), hits_str(b[13]), b[4]]
                   for b in best_pair_rows]
    result_rows.sort(key=lambda r: (rank_out.get(r[7], 9), r[0]))

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
