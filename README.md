# somatiq_re

Somatic repeat-expansion index from short-read Illumina WGS/WES CRAM.

**Approach (P0, prancSTR-first):** at a catalog of known loci, run **HipSTR** (germline
genotype + per-locus stutter model + per-read `MALLREADS`) then **prancSTR** (single-sample
mosaic fraction `f`, mosaic allele, p-value) to produce a per-sample × locus **somatic
score**, with level-1 QC and optional annotation from pre-computed ExpansionHunter results.

See `DESIGN.md` for the full design and tool comparison, and `docs/P0_SPEC.md` for the P0
spec. Later phases add the cohort tool + level-2 QC, a native pysam sizing engine (full
expansion-index ladder, large alleles), a PCR-exome mode, and a light GeneMapper/`.fsa` CE
index.

## Status
P0 in progress. The **core logic** (catalog, HipSTR/prancSTR/EH parsing, MALLREADS,
interruption detection, level-1 QC, score-row assembly) is implemented and unit-tested
(stdlib-only). HipSTR/prancSTR execution and CRAM read access are thin wrappers run inside
the container where the data lives.

## v1 loci
TCF4 (CTG18.1), AR (SBMA, X-linked), DMPK (DM1), ATXN7 (SCA7) — GRCh38, `catalogs/somatiq_grch38.json`.

## Install
```bash
pip install -e .            # core (stdlib only)
pip install -e '.[cram]'    # + pysam (read-length detection)
pip install -e '.[parquet]' # + pandas/pyarrow (Parquet output)
pip install -e '.[dev]'     # + pytest
```

## Usage
```bash
# validate the catalog / preview HipSTR regions
somatiq catalog-check --catalog catalogs/somatiq_grch38.json
somatiq regions --catalog catalogs/somatiq_grch38.json --out regions.bed

# score one sample (runs HipSTR -> prancSTR; needs those tools on PATH -> use the container)
somatiq call \
  --cram SAMPLE.cram --ref GRCh38.fa \
  --catalog catalogs/somatiq_grch38.json \
  --sample SAMPLE --sex M \
  --out SAMPLE.somatiq.tsv

# optionally annotate with a pre-computed ExpansionHunter VCF (default: off)
somatiq call ... --eh-annotate --eh-vcf SAMPLE.eh.vcf
```
Access is **targeted**: only the CRAM header, `.crai` index, and the locus slices are read —
never the full CRAM. Ensure the `.crai` is co-located and the reference is reachable.

## Container
```bash
docker build -t somatiq_re .
docker run --rm -v /data:/data somatiq_re call --cram /data/S.cram --ref /data/GRCh38.fa \
  --catalog catalogs/somatiq_grch38.json --sample S --sex M --out /data/S.tsv
```

## Test
```bash
pytest -q
```

## Score row (key columns)
`prancstr_f` (the P0 score), `prancstr_mosaic_allele`, `prancstr_pval`, germline genotype,
`informative_reads`, `detectable_f_floor`, `read_length`, `allele_exceeds_readlen`,
`interruption_flag`, `stutter_model_source`, `qc_level1` / `qc_reasons`.
