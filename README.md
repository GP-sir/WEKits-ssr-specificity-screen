# SSR primer specificity screening

Assigns Top_Tier / High_Tier classes to SSR primer pairs from WheatOmics BLAST
tabular output (outfmt 6).

## Usage

```
python screen_ssr_locus.py blast_hits.tsv -o screened.xlsx
```

Options:

```
--chr 4A                          target chromosome
--interval 520000000 530000000    target interval in bp
--scope interval|chromosome       interval (default): a Gold primer must have its
                                  single hit inside the target interval;
                                  chromosome: inside the target chromosome
-o, --output                      .xlsx workbook, or a prefix for csv output
```

Multiple input files are accepted; hits are merged before classification.

## Input

BLAST tabular output without header, one hit per line:

```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

Primer names are expected as `<locus>_p<N>_P<M>_F` / `_R` (TriticeaeSSRdb) or
`<locus>_p<N>_F` / `_R` (MISA + Primer3). All forward and reverse primers of the
same locus are combined pairwise, so alternative primer pairs within a locus are
evaluated together.

## Classification

- Gold: the primer has exactly one reported alignment, inside the target region.
- Silver: the primer has exactly two reported alignments.
- multi: more than two alignments; off_target: a single alignment elsewhere.
- Pair: Gold + Gold = Top_Tier; Gold + Silver = High_Tier; otherwise excluded.
- One combination per locus is kept, in the order Top_Tier > High_Tier > excluded,
  then fewer total alignments, then amplicon size closest to 300 bp.

## Output

`locus_result` (best combination per locus), `top_combos` (all qualifying
combinations), `primer_detail` (per primer), `summary` (counts).

Requires Python 3.8+ and openpyxl (`pip install openpyxl`); csv output is written
if openpyxl is unavailable.
