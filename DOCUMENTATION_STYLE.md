# Documentation Writing Style Guide

Use this guide to revise documentation for software, bioinformatics pipelines, command-line tools, R/Python packages, workflow repositories, and analysis code. It is informed by the documentation style of PLINK, bcftools, GENESIS/Bioconductor, Diataxis, Google developer documentation, MDN, The Turing Way, and peer-reviewed guidance on scientific software.

## 1. Core Standard

Good documentation lets a competent user do the correct thing without guessing.

Write documentation that is:

- Clear: each sentence has one purpose.
- Concise: remove filler, repetition, and vague promises.
- Reproducible: commands, inputs, outputs, versions, and assumptions are explicit.
- Task-oriented: users can complete real work from the page.
- Reference-ready: options, parameters, files, and outputs are easy to scan.
- Honest: document limitations, failure modes, and expected edge cases.

Prefer precise, direct language over promotional or conversational language.

## 2. Audience First

Before writing, define the primary reader.

Examples:

- New user: needs installation, minimal example, expected output.
- Returning user: needs exact command syntax and option meanings.
- Pipeline maintainer: needs architecture, dependencies, configuration, and testing.
- Reviewer or collaborator: needs reproducibility details, versions, and rationale.

Write for the least informed competent user, not for yourself.

## 3. Recommended Documentation Structure

Use a stable structure across projects.

```text
README.md
docs/
  installation.md
  quickstart.md
  usage.md
  configuration.md
  inputs.md
  outputs.md
  reference.md
  troubleshooting.md
  changelog.md
  contributing.md
```

For small projects, keep everything in `README.md`. For larger projects, keep the README short and move detail into `docs/`.

## 4. README Template

Use this order.

```markdown
# Project Name

One-sentence description of what the tool does and for whom.

## Purpose

Briefly state the problem solved, expected use case, and scope.

## Requirements

- Tool/package versions
- Operating system assumptions
- Required input files
- Required reference files

## Installation

Commands needed to install or prepare the software.

## Quick Start

Minimal working example using test or example data.

## Inputs

Required and optional inputs, formats, naming conventions, and constraints.

## Usage

Common commands or workflow steps.

## Outputs

Generated files, formats, interpretation, and where they are written.

## Troubleshooting

Common errors, likely causes, and fixes.

## Citation

How to cite the software, dependencies, and reference datasets.

## License

License name and link.
```

## 5. Use the Diataxis Model

Separate documentation by purpose.

| Type | Purpose | Example |
|---|---|---|
| Tutorial | Teach by doing | Run the pipeline on example data |
| How-to guide | Solve a specific task | Impute chromosome 22 only |
| Reference | Provide exact facts | Command-line options |
| Explanation | Clarify concepts | Why relatedness correction is needed |

Do not mix all four modes in one long page unless the project is small.

## 6. Style Principles

Use short sentences.

Prefer:

> Run sample quality control before imputation.

Avoid:

> It is generally recommended that users should consider running sample quality control before proceeding with imputation.

Use active voice when possible.

Prefer:

> The pipeline writes filtered variants to `results/qc/`.

Avoid:

> Filtered variants are written to `results/qc/` by the pipeline.

Use imperative verbs for instructions.

Prefer:

> Create the output directory.

Avoid:

> The output directory should be created.

Use present tense for tool behavior.

Prefer:

> `bcftools view` filters records by region.

Avoid:

> `bcftools view` will filter records by region.

## 7. Concision Rules

Delete words that do not change meaning.

Avoid:

- basically
- simply
- obviously
- very
- quite
- robust
- seamless
- easy
- user-friendly
- powerful
- state-of-the-art, unless justified

Replace vague claims with concrete information.

Prefer:

> Processes one chromosome per job and writes one compressed VCF per chromosome.

Avoid:

> Provides a powerful and scalable workflow for genomic data.

## 8. Command-Line Documentation

Follow the style used by tools such as PLINK and bcftools: syntax first, then options, then examples.

Recommended structure:

````markdown
## Command

```bash
tool subcommand --input INPUT --output OUTPUT [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `--input` | Input VCF, BCF, BED, or table file. |
| `--output` | Output prefix or directory. |

### Optional arguments

| Argument | Default | Description |
|---|---:|---|
| `--threads` | `1` | Number of CPU threads. |
| `--maf` | none | Minimum minor allele frequency. |

### Example

```bash
tool subcommand \
  --input data/chr22.vcf.gz \
  --output results/chr22 \
  --threads 8
```

### Output

| File | Description |
|---|---|
| `results/chr22.filtered.vcf.gz` | Filtered variant file. |
| `results/chr22.log` | Runtime log and filtering summary. |
````

Use placeholders in uppercase: `INPUT`, `OUTPUT`, `PREFIX`, `THREADS`.

Use literal file names in lowercase: `chr22.vcf.gz`, `samples.tsv`.

## 9. Bioinformatics Pipeline Documentation

For pipelines, document each stage in a fixed pattern.

````markdown
## Step 1: Variant Quality Control

### Purpose

Remove variants that fail missingness, allele frequency, or Hardy-Weinberg filters.

### Inputs

| File | Description |
|---|---|
| `data/raw.bed` | PLINK binary genotype file. |
| `data/raw.bim` | Variant metadata. |
| `data/raw.fam` | Sample metadata. |

### Command

```bash
plink2 \
  --bfile data/raw \
  --geno 0.02 \
  --maf 0.01 \
  --hwe 1e-6 \
  --make-bed \
  --out results/qc/variants_qc
```

### Outputs

| File | Description |
|---|---|
| `variants_qc.bed/.bim/.fam` | Filtered genotype dataset. |
| `variants_qc.log` | PLINK log file. |

### Notes

- Adjust thresholds before final analysis if the study design requires it.
- Record all threshold changes in the methods section or changelog.
````

Each pipeline step should state purpose, inputs, command, outputs, and notes.

## 10. Inputs and Outputs

Document file formats explicitly.

For each input, include:

- File name or pattern.
- Required or optional status.
- Format.
- Required columns or fields.
- Genome build or reference version, when relevant.
- Compression and indexing requirements.

Example:

| Input | Required | Format | Notes |
|---|---|---|---|
| `genotypes.vcf.gz` | Yes | bgzipped VCF | Must be indexed with `tabix`. |
| `samples.tsv` | Yes | TSV | Requires `sample_id`, `sex`, and `batch`. |
| `regions.bed` | No | BED | Coordinates must match the VCF genome build. |

For outputs, include:

- File path.
- Format.
- Meaning.
- Whether the file is final, intermediate, or diagnostic.
- Whether it can be regenerated.

## 11. Parameter Documentation

Every parameter should answer five questions.

| Question | Example |
|---|---|
| What does it do? | Filters variants by missing genotype rate. |
| What type is it? | Float between `0` and `1`. |
| What is the default? | `0.02`. |
| When should it change? | Use stricter thresholds for final GWAS. |
| What happens if omitted? | No missingness filter is applied. |

Avoid documenting only the name of the parameter.

Poor:

> `--geno`: genotype missingness.

Better:

> `--geno FLOAT`: Exclude variants with missing genotype rate greater than `FLOAT`.

## 12. Examples

Examples should be complete enough to run.

Good examples include:

- Realistic file names.
- Required indexes or companion files.
- Output location.
- Expected result.
- Small test dataset when possible.

Avoid examples that depend on hidden local paths.

Poor:

```bash
run_pipeline.sh input output
```

Better:

```bash
bash run_pipeline.sh \
  --vcf examples/chr22.vcf.gz \
  --samples examples/samples.tsv \
  --out results/example_chr22
```

Then state expected outputs.

## 13. Reproducibility Requirements

For scientific software, always document:

- Software name and version.
- Reference genome build.
- Reference panel version.
- Database release.
- Random seed, if stochastic methods are used.
- Container image or environment file.
- Hardware assumptions for large jobs.
- Date-sensitive external resources.

Example:

```markdown
This workflow was tested with PLINK 2.0 alpha, bcftools 1.20, R 4.4, and GENESIS 2.34. Genomic coordinates use GRCh38. TOPMed imputation was run with the r3 reference panel.
```

## 14. Error and Troubleshooting Style

Troubleshooting entries should be diagnostic, not generic.

Use this pattern:

````markdown
### Error: `Failed to open index file`

**Likely cause:** The VCF is compressed but not indexed.

**Fix:**

```bash
tabix -p vcf genotypes.vcf.gz
```

**Check:**

```bash
ls genotypes.vcf.gz.tbi
```
````

Avoid:

> Make sure your files are correct.

## 15. Scientific Tone

Use precise, neutral scientific language.

Prefer:

> The workflow estimates kinship coefficients before association testing to reduce inflation due to relatedness.

Avoid:

> The workflow smartly handles family structure.

Do not overclaim.

Prefer:

> This step reduces technical artifacts caused by low-quality genotypes.

Avoid:

> This step guarantees high-quality results.

## 16. Code Comments

Comments should explain why, not restate what.

Poor:

```python
# Set x to 5
x = 5
```

Better:

```python
# TOPMed rejects chromosome labels with the `chr` prefix.
vcf = strip_chr_prefix(vcf)
```

Use comments sparingly. If the code needs many comments to be understandable, improve names or structure first.

## 17. Markdown Style

Use Markdown consistently.

- Use one `#` title per file.
- Use sentence-case headings.
- Use fenced code blocks with language labels.
- Use tables for options, files, and parameters.
- Use bullets for short lists.
- Use numbered lists only for ordered procedures.
- Put commands in `bash` fences.
- Put paths, flags, variables, and file names in backticks.
- Keep lines reasonably short, but prioritize readability.

Example:

````markdown
## Filter variants

Run variant-level quality control with PLINK.

```bash
plink2 \
  --bfile data/raw \
  --geno 0.02 \
  --make-bed \
  --out results/variants_qc
```
````

## 18. Naming Conventions

Use predictable names.

| Item | Convention | Example |
|---|---|---|
| Directories | lowercase, descriptive | `results/qc/` |
| Files | lowercase with underscores | `sample_metadata.tsv` |
| Shell variables | uppercase | `VCF`, `OUTDIR` |
| Command options | tool-defined | `--maf`, `--geno` |
| Workflow steps | numbered when ordered | `01_variant_qc` |

Avoid spaces in file and directory names.

## 19. Versioning and Changelog

Maintain `CHANGELOG.md` for user-visible changes.

Use categories:

- Added
- Changed
- Fixed
- Deprecated
- Removed

Document changes that affect commands, outputs, parameters, dependencies, or interpretation.

Example:

```markdown
## 1.2.0

### Changed

- Updated default genotype missingness threshold from `0.05` to `0.02`.
- Renamed `results/final.vcf.gz` to `results/imputed_filtered.vcf.gz`.
```

## 20. Review Checklist

Before publishing documentation, check:

- The first paragraph states what the project does.
- Installation commands are current.
- The quick start runs from a clean checkout.
- Every input file is documented.
- Every output file is documented.
- All required parameters have defaults or clear instructions.
- Commands use realistic paths.
- Software versions are stated.
- Genome build and reference datasets are stated.
- Troubleshooting covers common failures.
- No sentence makes unsupported claims.
- No section duplicates another section.
- Links work.

## 21. Before and After Examples

### Example 1

Before:

> This script performs QC on genotype data and outputs cleaned files.

After:

> `run_variant_qc.sh` filters PLINK genotype files by variant missingness, minor allele frequency, and Hardy-Weinberg equilibrium. It writes GWAS-ready binary PLINK files to `results/qc/`.

### Example 2

Before:

> Use this option if needed.

After:

> Use `--keep RELATED_SAMPLES` to retain related individuals during kinship estimation.

### Example 3

Before:

> The pipeline is easy to use and fast.

After:

> On the test dataset of 1,000 samples and 50,000 variants, the pipeline completes variant QC in approximately 3 minutes using 8 threads.

## 22. Project Update Workflow

When revising documentation in any project:

1. Read the README and identify the primary user.
2. List all commands, inputs, outputs, and parameters.
3. Check whether the quick start runs.
4. Rewrite the README around purpose, installation, quick start, usage, outputs, and troubleshooting.
5. Move long reference material into separate files.
6. Standardize command examples.
7. Add missing reproducibility details.
8. Remove filler and unsupported claims.
9. Validate all commands and links.
10. Record user-visible changes in `CHANGELOG.md`.

## 23. Source Guidance

This guide follows principles from:

- PLINK 2.0 documentation: compact command syntax, explicit flags, file-format reference, and log-oriented behavior.
- bcftools manual: command-first reference documentation, concise option tables, and complete examples.
- GENESIS and Bioconductor documentation: vignette-based workflows, reproducible R examples, package metadata, and citation expectations.
- Diataxis: separation of tutorials, how-to guides, reference, and explanation.
- Google developer documentation style guide: active voice, present tense, direct instructions, and precise code formatting.
- MDN writing guidelines: clear page purpose, scannable structure, and practical examples.
- The Turing Way: reproducibility, transparent methods, computational environment capture, and research software openness.
- PLOS Computational Biology guidance on scientific software: usability, documentation, testing, versioning, and reproducibility.

## 24. Useful References

- PLINK 2.0 documentation: https://www.cog-genomics.org/plink/2.0/
- bcftools manual: https://samtools.github.io/bcftools/bcftools.html
- GENESIS Bioconductor package: https://bioconductor.org/packages/release/bioc/html/GENESIS.html
- Bioconductor package documentation guidance: https://contributions.bioconductor.org/docs.html
- Diataxis documentation framework: https://diataxis.fr/
- Google developer documentation style guide: https://developers.google.com/style
- MDN writing guidelines: https://developer.mozilla.org/en-US/docs/MDN/Writing_guidelines
- The Turing Way: https://the-turing-way.netlify.app/
- Lee BD. Ten simple rules for documenting scientific software. PLOS Computational Biology. 2018;14(12):e1006561. https://doi.org/10.1371/journal.pcbi.1006561
- Wilson G, Bryan J, Cranston K, et al. Good enough practices in scientific computing. PLOS Computational Biology. 2017;13(6):e1005510. https://doi.org/10.1371/journal.pcbi.1005510
