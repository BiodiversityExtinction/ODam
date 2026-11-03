# ODam

`ODam` profiles selected transversion mismatch patterns directly from an aligned BAM and the corresponding reference FASTA. It is intended to screen for G>T and C>A excesses that may be consistent with oxidative damage, while retaining G>C and C>G as comparison transversions.

It computes:
- Whole-read enrichment ratios (`G>T/G>C`, `C>A/C>G`, combined)
- Whole-read mismatch proportions (`G>T/G`, `G>C/G`, `C>A/C`, `C>G/C`)
- A signed combined excess mismatch rate
- End-specific summaries (`3p`, `5p`)
- Terminal-vs-interior fold change
- Per-position mismatch proportions along reads
- Optional PDF plots and TSV outputs

## Requirements

- Python 3.9+
- `pysam`
- A readable reference FASTA and its `.fai` index
- R (`Rscript`) only when PDF output is requested

Install `pysam` if needed:

```bash
pip install pysam
```

## Files

- Main script: `ODam.py`

## Quick Start

Single sample:

```bash
python3 ODam.py \
  --bam /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --sample-label sample1 \
  --output-dir odam_out \
  --plot-log-y \
```

Batch mode:

```bash
python3 ODam.py \
  --sample-list sample_bams.tsv \
  --output-dir odam_batch_out \
  --plot-log-y \
```

## Sample List Format

`--sample-list` expects a tab-separated text file with either two or three columns.

Recommended format, with one reference FASTA per sample:

```tsv
sample	bam	reference
sample1	/path/to/sample1.bam	/path/to/reference1.fa
sample2	/path/to/sample2.bam	/path/to/reference2.fa
```

Alternative format, when all BAMs use the same reference:

```tsv
sample	bam
sample1	/path/to/sample1.bam
sample2	/path/to/sample2.bam
```

For the two-column format, provide the shared reference on the command line:

```bash
python3 ODam.py \
  --sample-list sample_bams.tsv \
  --reference /path/to/reference.fa \
  --output-dir odam_batch_out
```

The header row is optional. If included, it must be the first physical line in the file. Comment lines starting with `#` and blank lines are ignored after header detection.

Important details:

- Columns must be separated by tabs, not spaces.
- Column 1 is the sample label used in output filenames.
- Column 2 is the BAM path.
- Column 3 is optional and is the reference FASTA path.
- If column 3 is omitted, `--reference` must be supplied.
- If column 3 is present, different samples can use different references.

## Parameters

Run `python3 ODam.py --help` for full CLI help.

Core inputs:

- `--bam`: input BAM for single-sample mode
- `--sample-list`: 2- or 3-column TSV for batch mode
- `--reference`:
  - required in single-sample mode
  - optional in batch mode if sample list provides per-row reference (3rd column)
  - can be used as a global fallback for 2-column batch lists
- `--output-dir`:
  - optional base output directory
  - in single mode, fills in any unspecified summary, position, and PDF paths with default files in that directory
  - in batch mode, fills in any unspecified outputs with `ODam_batch_summary.tsv`, `plots/`, and `positions/` under that directory
  - because a default PDF path is created, using `--output-dir` requires `Rscript`; to run without R, omit `--output-dir` and request only TSV outputs explicitly

Output options:

- `--pos-tsv-out`: per-position TSV in single-sample mode
- `--plot-pdf-out`: PDF plot in single-sample mode
- `--summary-tsv-out`: one-row summary TSV in single-sample mode
- `--summary-tsv-append`: append row to `--summary-tsv-out` (header written if file is new)
- `--batch-summary-out`: summary TSV for batch mode
- `--batch-plot-dir`: write one PDF per sample in batch mode
- `--batch-pos-dir`: write one per-position TSV per sample in batch mode
- `--reuse-existing-pos`: reuse the expected per-position TSV for end-specific summaries and plotting; the BAM is still rescanned with the current filters and sampling settings to obtain valid independent whole-read statistics and corrected terminal/interior folds

Read selection / filtering:

- `--max-reads` (default `1000000`):
  - randomly samples eligible reads without replacement, with every eligible read having an equal probability of selection
  - processes up to the specified number of reads per sample
  - sampling uses two BAM passes: one to count eligible reads and one to process the selected reads
  - `<=0` means use all reads
- `--random-seed` (default `1`):
  - RNG seed for reproducible random sampling
- `--min-mapq` (default `30`)
- `--min-baseq` (default `30`)
- `--region`:
  - optional region string (`chr:start-end`) for testing/debug
  - requires an indexed BAM

Eligible alignments must be mapped, primary, non-supplementary, not marked duplicate, not marked QC-failed, have `MAPQ >= --min-mapq`, and contain sequence and base qualities. Duplicate removal is based on the BAM duplicate flag; ODam does not identify duplicates itself. At the base level, ODam counts only aligned query/reference pairs with `baseQ >= --min-baseq`, a reference base of G or C, and an observed read base of A, C, G, or T. Soft clips, insertions, deletions, ambiguous bases, and unaligned bases do not contribute to counts.

Profiling settings:

- `--window` (default `10`):
  - terminal window size used for terminal vs interior fold
  - interior bases must be more than `window` bp from both read ends
- `--max-pos` (default `70`):
  - maximum distance from each read end retained in the 5′ and 3′ positional profiles
  - does not limit or change whole-read counts and overall statistics
- `--plot-max-pos` (default `30`):
  - maximum distance shown in plot and used for mean/median/max mismatch summaries
- `--pseudocount` (default `0.5`):
  - stabilizes ratios when denominator counts are low/zero

Orientation/end handling:

- `--normalize-ends` (default):
  - recommended for robust 3p/5p interpretation
- `--no-normalize-ends`:
  - raw query-orientation end assignment

Distances from read ends use BAM query-sequence coordinates. Soft-clipped bases are not counted as observations, but because they remain part of the query sequence they can increase the recorded distance of an aligned base from a read end.

Performance:

- `--threads` (default `1`):
  - HTSlib BAM decompression threads

Plot axis options:

- `--plot-log-y`: log10 y-axis for mismatch proportion plots
- `--plot-y-max`: fixed y-axis maximum

## Output Meaning

### Console summary

For each sample, `ODam` prints:

- `Reads processed`: number of reads used after filter checks and read cap (`--max-reads`)
- `Read sampling mode`:
  - random eligible-read sampling when `--max-reads` is positive
  - all eligible reads when `--max-reads <= 0`
- `BAM decompression threads`: value of `--threads`
- `End binning mode`:
  - `strand-normalized` (default): reverse-mapped reads are flipped for 5p/3p assignment
  - `raw query-orientation`: end assignment follows query orientation directly
- `Overall G>T enrichment ratio [G>T / G>C]`:
  - relative abundance of G>T compared with G>C
  - `>1` means G>T occurs more often than G>C among the counted observations
- `Overall C>A enrichment ratio [C>A / C>G]`:
  - relative abundance of C>A compared with C>G
- `Overall combined enrichment [(G>T + C>A) / (G>C + C>G)]`:
  - pooled G>T and C>A mismatches relative to pooled G>C and C>G mismatches
- `Overall G>T/G`, `G>C/G`, `C>A/C`, and `C>G/C mismatch proportions`:
  - whole-read mismatch counts divided by the corresponding whole-read reference G or C observations
- `Overall combined excess mismatch rate`:
  - signed value calculated as `0.5 * [((G>T - G>C) / G) + ((C>A - C>G) / C)]`
  - calculated independently for each sample, without modern-sample correction or a pseudocount
  - negative values are retained
- `3p/5p end ... enrichment ratio`:
  - same ratios computed per end
- `3p/5p end terminal-vs-interior fold (combined ratio)`:
  - terminal is positions `1..window` from the focal end
  - interior is positions `window+1..max_pos` from the focal end, excluding any base within `window` bp of the opposite end
  - with defaults, this compares positions 1–10 against positions 11–70 that are also more than 10 bp from the opposite end
  - reads of length `<= 2 * window` contribute no true-interior bases
  - `>1` means terminal enrichment
- `G>T/G mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `G>T / G` for plotted range (`--plot-max-pos`)
- `G>C/G mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `G>C / G` for plotted range
- `C>A/C mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `C>A / C` for plotted range
- `C>G/C mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `C>G / C` for plotted range

### Per-position TSV (`--pos-tsv-out`, `--batch-pos-dir`)

Each row aggregates observations by distance from a read end across the selected alignments. It is not a genomic site-by-site table and does not retain chromosome or genomic coordinate.

Columns:

- `End`: `3p` or `5p`
- `Pos`: distance from end (1-based)
- `G`, `C`: reference-base observations at that position
- `G>T`, `G>C`, `C>A`, `C>G`: mismatch counts
- `gt_over_gc`: `(G>T + pseudocount) / (G>C + pseudocount)`
- `ca_over_cg`: `(C>A + pseudocount) / (C>G + pseudocount)`
- `combined_over_other`: `(G>T + C>A + pseudocount) / (G>C + C>G + pseudocount)`
- `p_gt_over_g`: `G>T / G`
- `p_gc_over_g`: `G>C / G`
- `p_ca_over_c`: `C>A / C`
- `p_cg_over_c`: `C>G / C`
- `excess_gt_rate`: `max(0, (G>T - G>C) / G)`; control-adjusted per-position rate for adding `G>T` damage
- `excess_ca_rate`: `max(0, (C>A - C>G) / C)`; control-adjusted per-position rate for adding `C>A` damage

The positional excess-rate columns do not use pseudocounts. A negative signal-minus-control estimate is written as zero because these columns are intended as non-negative input rates for simulation. If no relevant G or C observations are available at a position, the corresponding value is `NA`. This clipping applies only to the positional `excess_gt_rate` and `excess_ca_rate`; the whole-read `overall_combined_excess_rate` retains negative values.

### Summary TSV (`--summary-tsv-out`, `--batch-summary-out`)

One row per sample with:

- sample/path/reference
- `status`: `ok`, `ok_reused_pos`, `plot_error`, or `error` (`plot_error` is used in batch mode, where a failed plot does not stop the batch)
- `error`: error message for failed samples or failed PDF generation
- `reads_processed`
- `missing_reference_contigs`: number of BAM contigs absent from the reference FASTA
- `reads_skipped_missing_reference`: reads skipped because their contig was absent from the reference FASTA
- `g_bases_observed`, `c_bases_observed`, `gc_bases_observed`: reference G/C observations across entire aligned reads, with every eligible aligned base counted once
- `gt_count`, `gc_count`, `ca_count`, `cg_count`: whole-read transversion counts, with every eligible mismatch counted once
- `oxidative_transversion_count`: `G>T + C>A`
- `control_transversion_count`: `G>C + C>G`
- `total_tracked_transversion_count`: `G>T + C>A + G>C + C>G`
- `overall_gt_over_gc`: `(G>T + 0.5) / (G>C + 0.5)` by default
- `overall_ca_over_cg`: `(C>A + 0.5) / (C>G + 0.5)` by default
- `overall_combined`: `(G>T + C>A + 0.5) / (G>C + C>G + 0.5)` by default
- `overall_p_gt_over_g`, `overall_p_gc_over_g`, `overall_p_ca_over_c`, `overall_p_cg_over_c`: whole-read mismatch proportions
- `overall_combined_excess_rate`: signed whole-read excess mismatch rate, calculated independently within each sample as `0.5 * [((G>T - G>C) / G) + ((C>A - C>G) / C)]`. It uses no pseudocount, is not corrected using modern samples, and retains negative values. Negative values are retained because this summary field is intended for plotting and statistical comparisons; unlike the positional SCAR rates, it is not a mutation probability.
- mean/median/max per-position mismatch proportions over `--plot-max-pos`, including all four proportions: `G>T/G`, `G>C/G`, `C>A/C`, and `C>G/C`
- 3p and 5p enrichment ratios
- 3p and 5p terminal/interior fold values
- `interior_definition`: records that interior bases are more than `window` bp from both read ends

In batch mode, the summary header is created at startup and each sample row is appended and flushed as soon as that sample finishes. This allows the table to be inspected while ODam is still running and preserves completed rows if the job is interrupted. One failed sample does not stop the whole run. Failed samples are retained in the summary table with `status=error`. If only PDF plotting fails in batch mode, the computed metrics are still written with `status=plot_error`. In single-sample mode, a plotting failure is reported as an error and stops the command.

Whole-read, end-profile, and terminal/interior accounting are deliberately independent. Overall statistics scan every eligible aligned base across the entire selected read exactly once. The 5′ and 3′ profiles are calculated separately within `--max-pos`; a base can appear in both profiles when the two windows overlap on a short read. This overlap is correct for end profiling and is never used to derive overall counts. Terminal/interior folds use dedicated scan-time buckets: an interior base must be more than `--window` bp from both read ends, so a terminal signal at the opposite end cannot contaminate the focal-end interior estimate.

Suggested manuscript wording: “For each read end, terminal enrichment was calculated by comparing positions 1–10 with positions 11–70, excluding from the latter any bases situated within 10 bp of the opposite read end.”

If a previous batch run created per-position TSVs but failed before writing the batch summary, rerun with `--reuse-existing-pos` and the same `--batch-pos-dir`. ODam reuses existing position files for end-specific metrics and plots, but always rescans each BAM to calculate independent whole-read counts and corrected terminal/interior folds. Existing positional TSVs cannot reconstruct either result because they do not retain distance from the opposite read end. Reused rows are marked with `status=ok_reused_pos`; if the BAM or reference is unavailable, the sample fails rather than silently calculating the old fold from positional buckets. Use positional files generated with the same `--max-pos`, `--normalize-ends`, and end-profile settings. The BAM rescan uses the current `--max-reads`, `--random-seed`, `--min-mapq`, `--min-baseq`, `--window`, and region settings.

### PDF plot (`--plot-pdf-out`, `--batch-plot-dir`)

Two panels (`3p`, `5p`) showing per-position mismatch proportions:

- `G>T / G` (orange)
- `G>C / G` (blue)
- `C>A / C` (green)
- `C>G / C` (magenta)

Defaults to positions `1..30` (`--plot-max-pos 30`).

## Tests

Run the synthetic BAM regression suite with a Python environment containing
`pysam`:

```bash
python3 -m unittest discover -s tests -v
```

The tests cover whole-read counting, 50/100/150/200 bp reads, reverse-mapped
reads, CIGAR operations, filtering, `--max-pos` and `--plot-max-pos`
independence, signed combined excess rates, true-interior accounting,
reproducible random sampling, single/batch consistency, and safe reuse of
positional TSV files.

## Recommended defaults for screening

For broad sample screening:

- Keep defaults for `--max-reads` and `--plot-max-pos`
- Use `--plot-log-y`
- Use `--normalize-ends` (default)
- Keep filters consistent across all samples (`--min-mapq`, `--min-baseq`)

Example:

```bash
python3 ODam.py \
  --sample-list sample_bams.tsv \
  --reference /path/to/reference.fa \
  --plot-log-y \
  --batch-summary-out batch_summary.tsv \
  --batch-plot-dir batch_plots
```

## Notes / Caveats

- Positive `--max-reads` values randomly sample eligible reads without replacement. Every eligible read has an equal probability of selection. The default `--random-seed 1` makes repeated runs with the same input and settings reproducible.
- Random sampling requires two passes through each BAM and may therefore take longer than processing the first reads encountered.
- ODam samples and reports alignments, not read pairs or sequencing templates. Paired-end mates are separate eligible alignments.
- Very low-count positions can produce unstable peaks; interpret per-position maxima with denominator counts in mind (`G`, `C` columns in TSV).
- Whole-read summaries generally use more observations than individual positions, but enrichment ratios can still be affected by low control counts.
- Overall statistics count each eligible aligned base once across the full read and are independent of `--max-pos` and `--plot-max-pos`; end-specific profiles intentionally may count a short-read base in both end windows. Terminal/interior folds exclude opposite-terminal bases from the interior bucket.
- Reads aligned to BAM contigs absent from the reference FASTA are skipped and counted in the summary diagnostics.
- ODam measures reference/read mismatches; it does not by itself distinguish molecular damage from genuine biological differences, reference divergence, mapping error, sequencing error, or library-specific artifacts.
- G>C and C>G are comparison transversions, not guaranteed estimates of the background expected specifically for G>T and C>A. Enrichment and excess metrics should therefore be interpreted alongside per-position profiles, denominator counts, consistent processing, and appropriate controls.
- ODam does not call genotypes, mask known variants, model sequence context, or prove that an observed excess was caused by oxidation.
- `overall_combined_excess_rate` is a descriptive contrast that gives equal weight to the G-based and C-based channels. Inspect G>T and C>A separately when the biological question does not assume that both occur together.
