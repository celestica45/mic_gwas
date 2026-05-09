# MIC GWAS pyseer pipeline

We are building this step by step.

Current metadata cleaning:

- remove `phenotype-` from column names
- filter to the configured antibiotic and organism
- drop samples without parseable `gen_measurement`
- add numeric `mic_value` and `mic_log2` columns
- write MIC log2 statistics, counts, and a labeled pie chart
- add `antibiotic_name` as the first column
- write a cleaning summary
- write a separate metadata CSV with only samples that have `assembly_ID`
- fetch NCBI Assembly metadata from `assembly_ID`

Current assembly download stage:

- write unique `assembly_ID` values to an accession list
- download genome FASTAs using NCBI Datasets
- unpack each FASTA to `results/oxacillin/assemblies/`
- name each FASTA with the matching `assembly_ID`, for example `GCA_003073635.1.fa`

Current assembly QC stage:

- run QUAST on downloaded assembly FASTAs
- optionally use a user-provided reference genome and matching annotation
- produce QUAST assembly reports only
- do not filter assemblies yet

Sequencing metadata/read-pair checking is available but off by default:

```yaml
skip_sequence_metadata: true
```

Pipeline stages can be skipped in `config/config.yaml`:

```yaml
skip_metadata: false
skip_assembly_metadata: true
skip_sequence_metadata: true
skip_assembly_download: false
skip_quast: false
skip_panaroo_qc: false
skip_phenotype_data: true
skip_gwas_runs: true
```

`skip_phenotype_data` controls GitHub-style phenotype/covariate files, and
`skip_gwas_runs` is reserved for later pyseer association tests.

QUAST settings:

```yaml
quast:
  threads: 4
  min_contig_length: 500
  extension: "fa"
  reference: null
  annotation: null
```

When `quast.reference` is `null`, QUAST runs reference-free. If you provide a
reference FASTA path, QUAST runs with `-r`. If you also provide a matching
annotation path, QUAST runs with `-g`. The workflow does not download the
reference genome or annotation.

Run:

```bash
snakemake --cores 4 --use-conda
```

Run only the QUAST stage:

```bash
snakemake --cores 4 --use-conda quast_stage
```

Run only the Prokka annotation stage:

```bash
snakemake --cores 4 --use-conda prokka_stage
```

Output:

```text
results/oxacillin/metadata/oxacillin_clean_metadata.csv
results/oxacillin/metadata/oxacillin_metadata_summary.yaml
results/oxacillin/metadata/oxacillin_clean_metadata_with_assembly.csv
results/oxacillin/metadata/oxacillin_mic_log2_stats.csv
results/oxacillin/metadata/oxacillin_mic_log2_counts.csv
results/oxacillin/metadata/oxacillin_mic_log2_pie_chart.png
results/oxacillin/metadata/oxacillin_assembly_metadata.csv
results/oxacillin/downloads/oxacillin_assembly_accessions.txt
results/oxacillin/downloads/oxacillin_assembly_manifest.csv
results/oxacillin/assemblies/
results/oxacillin/qc/quast/report.tsv
results/oxacillin/qc/quast/report.txt
results/oxacillin/qc/quast/transposed_report.tsv
```

The NCBI Datasets zip is a Snakemake `temp()` intermediate, so the final
download/collection marker is `oxacillin_assembly_manifest.csv`.

Current GWAS input stage:

- write MIC log2 phenotype and covariate metadata
- run Prokka on downloaded assemblies
- run Panaroo QC on Prokka GFFs
- run Panaroo on Prokka GFFs
- use Panaroo's native gene presence/absence `.Rtab`
- require Panaroo's core genome alignment
- build Mash MDS population-structure files and scree plot
- run fixed-effects gene GWAS with pyseer when requested
- print the active GWAS input configuration in the terminal
- keep pyseer association tests out of default runs while `skip_gwas_runs: true`

GWAS settings:

```yaml
gwas:
  active_inputs:
    - genes

  population_structure:
    mash_mds:
      sketch_size: 10000
      max_scree_dimensions: 30

  phenotype:
    sample_column: "assembly_ID"
    log2_column: "mic_log2"

  covariates:
    columns:
      - country
      - collection_year
      - platform

  inputs:
    genes:
      population_structure:
        methods:
          - none
          - mash_mds
        mds_dimensions: 8
        covariates: []

    snps:
      population_structure:
        methods:
          - mash_mds
        mds_dimensions: 8
        covariates: []

    kmers:
      population_structure:
        methods:
          - phylogeny_lmm
        covariates: []
```

Run only the GWAS input stage:

```bash
snakemake --cores 4 --use-conda gwas_inputs_stage
```

Run only the GWAS population-structure stage:

```bash
snakemake --cores 4 --use-conda gwas_population_structure_stage
```

Run only fixed-effects gene GWAS:

```bash
snakemake --cores 4 --use-conda gwas_gene_fixed_effects_stage
```

The fixed-effects gene GWAS output folder changes automatically from the
per-input covariate config. With `covariates: []`, outputs go to:

```text
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/
```

With `covariates: [country]`, outputs go to:

```text
results/oxacillin/gwas/fixed_effects/gene/mash_mds/covariates_country/
```

Run only phenotype/covariate data:

```bash
snakemake --cores 4 --use-conda phenotype_data_stage
```

Run only Panaroo QC:

```bash
snakemake --cores 4 --use-conda panaroo_qc_stage
```

`skip_panaroo_qc: false` requests Panaroo QC even when GWAS inputs and pangenome
construction are skipped. If `skip_prokka: true`, Panaroo QC assumes existing
Prokka GFFs are already present in `results/oxacillin/annotation/prokka/`.

GWAS input outputs:

```text
results/oxacillin/gwas/phenotypes/oxacillin_phenotype.tsv
results/oxacillin/gwas/phenotypes/oxacillin_metadata.tsv
results/oxacillin/gwas/population_structure/mash_mds/mash_sketch.msh
results/oxacillin/gwas/population_structure/mash_mds/mash_raw.tsv
results/oxacillin/gwas/population_structure/mash_mds/mash.tsv
results/oxacillin/gwas/population_structure/mash_mds/scree_plot.png
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/oxacillin_gene_pyseer.txt
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/oxacillin_gene_patterns.txt
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/oxacillin_gene_top_hits.tsv
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/oxacillin_gene_qq_plot.png
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/oxacillin_gene_significance_thresholds.txt
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/oxacillin_gene_notes_summary.tsv
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/oxacillin_gene_manhattan_plot.png
results/oxacillin/gwas/fixed_effects/gene/mash_mds/no_covariates/mash_mds.pkl
results/oxacillin/qc/panaroo_qc/.qc_done
results/oxacillin/pangenome/panaroo/gene_presence_absence.csv
results/oxacillin/pangenome/panaroo/gene_presence_absence.Rtab
results/oxacillin/pangenome/panaroo/core_gene_alignment.aln
```

`mash_raw.tsv` contains labels as emitted by Mash. `mash.tsv` is normalized
back to the versioned assembly IDs from the assembly manifest so pyseer sample
matching is consistent with phenotype and Panaroo `.Rtab` inputs.
