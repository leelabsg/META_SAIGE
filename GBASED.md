# Gbased SAIGE-META

A single-cohort variant of the MetaSAIGE rare-variant gene-based test that
uses **only** GWAS summary statistics and the SAIGE step-3 LD matrix file
(`(row, col, value)` triplet). No access to individual-level genotypes is
required at meta-analysis time.

## What changed vs the upstream `Run_MetaSAIGE`

| Area | Upstream `Run_MetaSAIGE` | This branch `Run_GbasedSAIGE` |
|---|---|---|
| Entry point | `Run_MetaSAIGE(...)` | `Run_GbasedSAIGE(...)` (asserts `n.cohorts == 1`) |
| Number of cohorts | multi-cohort meta-analysis | single cohort |
| Ancestry stratification | `--ancestry` flag + ancestry-specific aggregation | removed entirely |
| URV (ultra-rare variant) LD treatment | URV pseudo-variant variance was approximated as `sum(MAC_URV)`, URV–RV cross-products treated as 0 | URV pseudo-variant variance and URV–RV cross-products are computed exactly by streaming the SAIGE step-3 LD file once and routing every triplet into the collapsed `k × k` matrix via a `gid` lookup |
| Annotation × MAF grid | one-at-a-time invocation | `expand.grid(annotation, mafcutoff)` in a single invocation |
| Binary trait `load_cohort` | broken (referenced symbols that only existed in the continuous branch) | parity with continuous, including the streaming LD aggregation |

The streaming aggregation is the load-bearing fix: each `(row, col, value)`
entry in the LD file contributes to **exactly one cell** of the collapsed
`M = C'ΣC` matrix, so Σ is never materialized and the URV pseudo-variant
gets a correct variance (and a correct URV–RV cross-row).

## Inputs and outputs

Inputs (all white-space-delimited text from SAIGE outputs):

- GWAS summary statistics from SAIGE step 2 (`--gwas_path`)
- LD-matrix marker information from SAIGE step 3 (`--info_file_path`)
- LD-matrix triplet file prefix from SAIGE step 3 (`--gene_file_prefix`),
  one file per gene named `<prefix><gene>.txt`
- Group file with per-variant functional annotations (`--groupfile`),
  same format as SAIGE-GENE+

Output: a tab-separated table with one row per `(gene, annotation_mafcutoff)`
mask plus one Cauchy-combined row per gene:

```
CHR  GENE  Group                          Pval     MAC    #Rare Variants  #Ultra Rare Variants  P-value of Collapsed Ultra Rare
2    TTN   missense_lof_synonymous_0.01   0.34731  88004  1024            0                     NA
2    TTN   Cauchy                         0.34731  NA     NA              NA                    NA
```

## Running it

Single mask, continuous trait:

```bash
Rscript SAIGE_META/inst/scripts/RV_meta_GC.R \
    --num_cohorts 1 \
    --trait_type continuous \
    --chr 2 \
    --col_co 10 \
    --info_file_path WB_case1_group1_step3_TTN.marker_info.txt \
    --gene_file_prefix WB_case1_group1_step3_TTN_ \
    --gwas_path WB_case1_group1_step2_LDL_chr2.txt \
    --output_prefix WBonly_LDL_case1_chr2_TTN_meta_res.txt \
    --groupfile UKBexomeOQFE_chr2.gene.anno.hg38_PlinkMatch_v2_TTN.txt \
    --annotation missense_lof_synonymous \
    --mafcutoff 0.01
```

Multi-mask (3 annotations × 3 MAF cutoffs in one call):

```bash
... \
    --annotation lof missense_lof missense_lof_synonymous \
    --mafcutoff 0.01 0.001 0.0001
```

Passing `--num_cohorts 2` is rejected with a pointer to the upstream
multi-cohort `Run_MetaSAIGE`.

## Verification

Three end-to-end tests are included under `test/`:

| Test | Asserts |
|---|---|
| `test_script.sh` (project root) | Single-mask continuous run reproduces `p = 0.347307034421392` for TTN/LDL (matches `test/baseline_continuous_single.txt`) |
| `test/test_multi_mask.sh` | Multi-mask run produces 10 rows (`3 × 3 + Cauchy`); the shared `(missense_lof_synonymous, 0.01)` cell matches the baseline byte-for-byte |
| `test/test_binary.sh` | Binary `load_cohort` runs end-to-end on a synthesized fixture; produces 10 rows with all p-values in `(0, 1]` |
