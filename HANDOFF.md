# somatiq_re — project handoff & notes

Context for resuming development in a fresh Claude Code CLI (another machine). Read this
first, then `DESIGN.md` (full design) and `docs/P0_SPEC.md` (P0 detail).

## What this project is
A lightweight Python tool + runner to quantify a **somatic repeat-expansion index** at
known tandem-repeat loci from standard Illumina **WGS/WES CRAM**, at cohort scale
(UK Biobank / All of Us / internal cohort). Primary target: reproduce the All of Us
**TCF4** blood-mosaicism signal, then apply to an internal cohort.

## Status @ v0.1.0
- **Design complete** (`DESIGN.md`, `docs/P0_SPEC.md`).
- **P0 implemented + unit-tested** (37 tests, stdlib-only core). CLI: `catalog-check`,
  `regions`, `call`.
- **Not yet run end-to-end** — `somatiq call` needs HipSTR + prancSTR + a real CRAM,
  which run in the owner's cloud (data is not accessible from the dev sandbox).
- **Reference build: GRCh38 only.** hg19/GRCh37 not supported or tested.

## Architecture (one line)
Pluggable front-ends → common per-read table → shared metric back-end. v1 engine =
**HipSTR-batch (germline + stutter) → prancSTR (mosaic fraction `f`)**. Native pysam
sizing engine + GeneMapper CE index come later. See `DESIGN.md` §6, §10.

## Repo map
- `src/somatiq/` — package. Modules: `catalog`, `mallreads`, `hipstr`, `prancstr`,
  `interruptions`, `eh` (optional annotation, default off), `readlen`, `qc`, `schema`,
  `call`, `cli`. Dependency flow: `cli→call→{hipstr,prancstr,eh,readlen}→schema→{catalog,
  mallreads,interruptions,qc}`. Core (up to `schema`) is stdlib-only.
- `catalogs/somatiq_grch38.json` — TCF4/AR/DMPK/ATXN7 (GRCh38, EH coords; ref- vs
  disease-strand motifs; interruption model; derived HipSTR regions to validate).
- `tests/` — 37 tests over the pure logic (no external tools needed).
- `Dockerfile`, `environment.yml`, `DEPLOY.md` — deployment (incl. air-gapped).
- `DESIGN.md`, `docs/P0_SPEC.md` — design + P0 spec.

## Decisions log (why things are the way they are)
1. **prancSTR-first**: fastest path to a peer-reviewed somatic score; it is **HipSTR-only**
   (needs `MALLREADS` + `INFRAME_*` stutter), so the pipeline runs HipSTR first.
2. **Native pysam engine deferred to P4**: recovers the full expansion-index *ladder*
   that prancSTR's single `f` collapses, and handles alleles beyond the read length.
3. **Read-length ceiling** (~40-44 units @150bp) is a fundamental short-read limit;
   flagged, never silently truncated (`allele_exceeds_readlen`).
4. **EH is not run by us** — it is the owner's upstream sample selector; we only optionally
   *annotate* from a pre-computed EH VCF (`--eh-annotate`, default off).
5. **Interruptions** (HTT/ATXN1/FMR1; DMPK variant repeats): P0 annotates from HipSTR ALT;
   pure-tract sizing + anchor-extension are native-engine (P4); correlation with `f` is a
   cohort (P1) job.
6. **Targeted access only** — header + `.crai` + locus slices; never a full-CRAM scan.
   `read_length` derived from the locus reads.
7. **Single-sample runs** allowed with HipSTR `--def-stutter-model`; `f` marked provisional
   (WARN) until a batch stutter model exists.

## Roadmap (next work, in order)
- **Finish P0 in the owner's cloud**: build the container, run `somatiq call` on a real
  CRAM, confirm exact HipSTR/prancSTR flag names + the prancSTR `.tab` columns match the
  installed versions (record in `docs/tool_flags.md`), validate the HipSTR region defs
  (ATXN7 compound `(GCA)*(GCC)+`, AR strand) against known genotypes.
- **simTR validation gate**: implement `somatiq validate-sim` (planned CLI) — plant known
  `f`, recover within tolerance; build the sensitivity-vs-depth curve.
- **QC hardening (level-1)**: add `dumpSTR`-style HipSTR gates (`DSTUTTER`/`DFLANKINDEL`
  fractions, `FILTER`) + `prancstr_pval`/`mosaic_support` thresholds. (Currently level-1
  uses depth, germline `Q`, stutter-model source, ceiling, interruption.)
- **P1 — `somatiq cohort` + level-2 QC**: aggregate → cohort FDR on `pval`, control-sample
  drift (NA06075), batch/stutter consistency, replicate concordance, `index ~ length+age`
  regression → reproduce AoU-TCF4.
- **P2** exome mode + batch runner; **P3** GeneMapper/`.fsa` CE index; **P4** native engine.

## How to run / verify
```bash
pip install -e '.[dev]'          # or pip install -e .  (core only)
pytest -q                        # 37 tests
somatiq catalog-check --catalog catalogs/somatiq_grch38.json
```
Full pipeline (HipSTR/prancSTR present): see `README.md` / `DEPLOY.md`.

## Gotchas
- prancSTR/HipSTR CLI flags and the `.tab`/VCF fields can drift between versions — verify
  once against the installed build (the wrappers isolate command construction + parsing).
- prancSTR false-positive rate rises under PCR/exome and with a mismatched stutter model —
  the exome arm needs a HipSTR-learned/validated stutter model.
- Never commit PHI, sample IDs, cohort data, real CRAM paths, or credentials.
