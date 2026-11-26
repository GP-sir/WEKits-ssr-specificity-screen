# SSR primer specificity screening

Two scripts for classifying SSR primer pairs by genome-wide specificity. Both use the same rules — a primer is Gold with one perfect hit, Silver with two; a pair is Top_Tier (Gold × Gold) or High_Tier (Gold × Silver) — and both write a single summary table.

| Script | Input | Alignment |
|---|---|---|
| `ssr_specificity_blastn.py` | BLAST tabular output, outfmt 6 | BLAST hits from WheatOmics |
| `ssr_specificity_bowtie2.py` | Primer FASTA file(s) | Bowtie2 against the wheat reference |

## BLAST

```bash
python ssr_specificity_blastn.py blast_hits.tsv -o final_result.tsv
```

Several hit files may be given and are merged before classification. Primer names end in `_F` / `_R`, with an optional pair number. By default the target region is 4A:520000000–530000000; change it with `--chr` and `--interval CHR START END`, or use `--scope chromosome` to accept Gold hits anywhere on the target chromosome.

The input is standard outfmt 6 (`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`), without a header line.

## Bowtie2

```bash
python ssr_specificity_bowtie2.py primers.fasta -o final_result.tsv
python ssr_specificity_bowtie2.py /path/to/fastas/ -o final_result.tsv --interval 4A 520000000 530000000 -p 16
python ssr_specificity_bowtie2.py primers.fasta --sam alignments.sam -o final_result.tsv
```

The first argument is a FASTA file or a directory of `*.fasta`. Primer names must end in `-F` / `-R` or `_F` / `_R`, optionally followed by a number; the leading token gives the chromosome. `--sam` reuses an existing alignment file instead of running Bowtie2.

Other options: `-x` sets the Bowtie2 index (default `/opt/ref/wheat/cs1.1/database/wheat_bt2`), `-p` threads (8), `-k` alignments reported per read (2), `--max-product-size` amplicon cutoff (300 bp), and `--keep-sam` keeps the intermediate SAM. Edit `DEFAULT_INDEX` at the top of the script to use a different reference.

## Output

`-o` is the output file, `.tsv` by default; `.xlsx` writes one `final_result` worksheet. Both scripts produce the same eight columns:

```
Locus  Primer_F  Primer_R  F_class  R_class  F_loci  R_loci  Tier
```

Only Top_Tier and High_Tier loci are written, one best pair per locus. `F_loci` and `R_loci` contain all perfect hits of each primer as `chr:position`, separated by `; `.

Python 3 is required. `openpyxl` is needed only for `.xlsx` output; Bowtie2 on `PATH` is needed only for the Bowtie2 script.
