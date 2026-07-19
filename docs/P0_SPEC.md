# P0 Implementation Spec — prancSTR‑first somatic pipeline (WGS)

**Status:** Spec for sign‑off (no code yet). Scope = roadmap **P0** in `DESIGN.md` §9.
**Goal:** produce a per‑sample × locus **somatic score** (prancSTR mosaic fraction `f`)
for the four v1 loci from PCR‑free WGS CRAM, via **HipSTR‑batch → prancSTR**, with an
optional **ExpansionHunter triage**, validated on **simTR** truth and the **NA06075**
control.

## 0. Done‑criteria (what "P0 complete" means)
1. `somatiq call` runs end‑to‑end on a CRAM/BAM and emits the score row (§5) for all 4 loci.
2. Runs the HipSTR‑batch → prancSTR path with the shipped catalog (§3).
3. Recovers a known mosaic fraction on **simTR**‑simulated data within tolerance (§6.1).
4. On **NA06075** (real control), calls the DMPK short germline allele and **flags the
   read‑length ceiling** for the 56/70 CTG alleles rather than mis‑calling them (§6.2).
5. Minimal **level‑1 QC** flags attached to every row (§4).

**Explicitly NOT in P0** (later phases): exome mode (P2), `somatiq cohort` + level‑2 QC
(P1), native pysam engine (P4), `somatiq genemapper` (P3). P0 is the thin, correct
vertical slice on WGS.

## 1. Environment & dependencies
External tools (pinned versions recorded in `env.lock` at build time):
- **HipSTR** (≥ v0.7) — germline STR genotyping + per‑locus stutter EM + `MALLREADS`.
- **prancSTR** via **TRTools** (`pip install trtools`, ≥ v6) — the somatic caller; also
  ships **simTR** (our simulator).
- **ExpansionHunter** (≥ v5) — optional Stage‑0 triage.
- **samtools/bcftools**, **pysam** (QC read stats), Python ≥ 3.10, `numpy`/`pandas`/`pyarrow`, `click`.

**Build step 1 (before any wrapper code):** run `HipSTR --help`, `prancSTR --help`,
`simTR --help` against the pinned versions and record exact flag names in
`docs/tool_flags.md` — the commands in §3–§4 are skeletons to confirm, not guesses.

Reference build: **GRCh38** (matches the EH catalog coords below). CRAM decode requires
the exact reference FASTA used for alignment.

## 2. Repo layout (created when P0 is approved)
```
somatiq_re/
├── DESIGN.md
├── docs/P0_SPEC.md            # this file
├── pyproject.toml
├── catalogs/somatiq_grch38.json
├── src/somatiq/
│   ├── catalog.py             # load/validate catalog
│   ├── triage_eh.py           # Stage 0: run/parse ExpansionHunter
│   ├── hipstr.py              # Stage 1a: HipSTR batch wrapper + VCF field extract
│   ├── prancstr.py            # Stage 1b: prancSTR wrapper + output parse
│   ├── call.py                # orchestration → score row
│   ├── qc.py                  # level-1 QC
│   ├── schema.py              # output row dataclass + Parquet/TSV writer
│   └── cli.py                 # `somatiq` entry points
└── tests/
    ├── test_catalog.py
    ├── sim/                   # simTR configs + recovery test
    └── validation/na06075/    # real-control harness
```

## 3. Locus catalog (`catalogs/somatiq_grch38.json`)
Coordinates/motifs from the Illumina ExpansionHunter GRCh38 catalog (authoritative).
**Strand note:** EH stores the *reference‑strand* motif, which differs from the
conventional *disease* motif — we record both. `ref_region` is 1‑based inclusive (EH
convention); the HipSTR region file is derived (0‑based, §3.1).

| locus_id | disease | chrom | ref_region (GRCh38, 1‑based) | ref_motif | disease_motif | ploidy | period |
|---|---|---|---|---|---|---|---|
| TCF4 | Fuchs dystrophy (CTG18.1) | chr18 | 55586155‑55586227 | CAG | CTG | diploid | 3 |
| AR | SBMA (Kennedy) | chrX | 67545316‑67545385 | GCA | CAG | **sex‑dependent** (hemizygous male) | 3 |
| DMPK | Myotonic dystrophy 1 | chr19 | 45770204‑45770264 | CAG | CTG | diploid | 3 |
| ATXN7 | SCA7 | chr3 | 63912684‑63912714 | GCA | CAG | diploid | 3 |

Per‑locus JSON fields:
```json
{
  "locus_id": "DMPK",
  "disease": "DM1",
  "chrom": "chr19",
  "ref_region": {"start": 45770204, "end": 45770264, "base": 1},
  "ref_motif": "CAG",
  "disease_motif": "CTG",
  "period": 3,
  "ploidy": "diploid",
  "eh_variant_id": "DMPK",
  "hipstr_region": {"start": 45770203, "end": 45770264, "period": 3, "ref_copies": 20, "name": "DMPK"},
  "notes": "NA06075 control = 12/56/70 CTG; 56/70 exceed a 150bp read"
}
```
- **AR ploidy** resolved per sample from provided sex (male → hemizygous, single germline
  mode; female → diploid). Sex source = sample manifest column (P0) — no inference.
- **ATXN7** has a compound `(GCA)*(GCC)+` structure in EH; for P0 HipSTR we treat the
  pure `(GCA)*` tract (validate modal call — see §3.1 caveat).

### 3.1 HipSTR region file (derived, must be validated)
HipSTR consumes a BED‑like region file: `chrom  start  end  period  ref_copies  name`
(0‑based start). We generate it from the catalog, then **validate** that HipSTR's modal
genotype on a panel of normal samples matches expected reference copy numbers before
trusting any locus. Compound/interrupted loci (ATXN7) and strand are the known risks —
each locus is signed off individually against known genotypes.

## 4. Pipeline (`somatiq call`)
```
Stage 0 (optional, --triage):  ExpansionHunter  → per-locus carrier flag + germline size
Stage 1a:  HipSTR (batch)      → per-sample VCF: GT, Q, DP, MALLREADS, INFRAME_UP/DOWN/PGEOM
Stage 1b:  prancSTR (per sample) → f, mosaic allele, p-value
Merge + level-1 QC             → score row (§5)
```

### Stage 0 — ExpansionHunter triage (optional)
- Cmd skeleton: `ExpansionHunter --reads CRAM --reference REF --variant-catalog cat.json --output-prefix PFX [--analysis-mode streaming]`
- Parse the EH VCF/JSON for each locus → `eh_gt_short`, `eh_gt_long`, `carrier` (long
  allele ≥ locus threshold) and, if requested, retain the **realigned BAM** for the
  future FE‑2 cross‑check (not used by prancSTR).
- Purpose in P0: prioritize which sample×locus to quantify and detect large alleles that
  the spanning‑only Stage 1 cannot size. Triage never gates Stage 1 unless `--triage-gate`.

### Stage 1a — HipSTR batch (learns stutter)
- **Batch = many samples pooled per locus so HipSTR's EM learns the PCR‑stutter model.**
  Too few samples ⇒ `--def-stutter-model` (fixed) with a WARN flag.
- Cmd skeleton: `HipSTR --bams B1,B2,... --fasta REF --regions somatiq.hipstr.bed --str-vcf out.vcf.gz [--def-stutter-model]`
- Extract per sample×locus: `GT`, `Q`, `DP`, `DSTUTTER`, `DFLANKINDEL`, **`MALLREADS`**
  (per‑read allele counts), and INFO **`INFRAME_UP` / `INFRAME_DOWN` / `INFRAME_PGEOM`**
  (the stutter model prancSTR needs).
- P0 default: single‑sample `somatiq call` runs HipSTR on that one CRAM with
  `--def-stutter-model`; the **batch/joint** mode (recommended) is invoked by
  `somatiq call --hipstr-batch manifest.tsv` and is how cohorts should run.

### Stage 1b — prancSTR (the somatic caller)
- prancSTR is **HipSTR‑only** and reads the fields above from the HipSTR VCF.
- Cmd skeleton: `prancSTR --vcf out.vcf.gz --vcftype hipstr --out PFX --samples SAMPLE [--region-file loci.bed]`
  (confirm exact flags incl. the read field, e.g. `--readfield MALLREADS`, via `--help`).
- Output per sample×locus: mosaic fraction **`f`**, **mosaic allele** (size/Δ), **p‑value**
  for `H0: f=0`, and any posterior/confidence field.

## 5. Output score row (Parquet + TSV), one per sample × locus
| column | source | note |
|---|---|---|
| sample_id, locus_id, disease | manifest/catalog | |
| sex, ploidy_used | manifest/catalog | AR resolved here |
| data_mode | fixed `wgs` (P0) | |
| germline_gt_copies, germline_gt_bp, germline_q | HipSTR GT/Q | inherited allele(s) |
| dp, dstutter, dflankindel | HipSTR | read support |
| stutter_up, stutter_down, stutter_pgeom | HipSTR INFO | learned model |
| stutter_model_source | `learned` \| `default` | `default`→WARN |
| **prancstr_f** | prancSTR | **the P0 score** |
| prancstr_mosaic_allele, prancstr_pval | prancSTR | |
| informative_reads, detectable_f_floor | MALLREADS / `1/n` | power |
| allele_exceeds_readlen | derived | ceiling flag |
| qc_level1, qc_reasons | §4 qc.py | PASS/WARN/FAIL + reasons |
| eh_carrier, eh_gt_long | Stage 0 (if run) | optional |

## 4bis. Level‑1 QC computed in P0 (`qc.py`)
- **min spanning/informative depth gate** (default ≥ 20 informative reads) → below = FAIL
  (low power) with `detectable_f_floor ≈ 1/informative_reads` reported.
- **stutter_model_source == default** → WARN (single‑sample run; batch recommended).
- **germline Q** below threshold → WARN/FAIL.
- **allele_exceeds_readlen**: germline or mosaic allele bp > detected read length → flag
  (prancSTR `f` unreliable for that allele; the honest NA06075 outcome).
- **base‑quality binning detection** (read the CRAM quality histogram): binned → note that
  the future McCarroll artifact filter is weakened (informational in P0).

## 6. Validation harness
### 6.1 simTR synthetic truth (primary, always available)
- Use `simTR` to simulate reads at each locus with a planted mosaic allele at known
  `f ∈ {0, 0.05, 0.10, 0.20, 0.40}` and depth `∈ {30×, 50×}`, PCR‑free error model.
- Run the full pipeline; **acceptance:** recovered `f` within ±0.05 (or ±25% relative) for
  `f ≥ 0.10` at 30–50×; `f=0` yields non‑significant p‑values at the target FPR. Produces
  the per‑locus sensitivity‑vs‑depth curve.

### 6.2 NA06075 real control (DMPK 12/56/70 CTG)
- **Prerequisite:** obtain a PCR‑free WGS CRAM of NA06075 (confirm availability; if none,
  6.2 is deferred and 6.1 stands as the P0 gate).
- 12 CTG = 36 bp (spannable); **56/70 CTG = 168/210 bp exceed a 150 bp read**.
- **Acceptance:** pipeline calls the ~12 germline allele from spanning reads and raises
  `allele_exceeds_readlen` for the expanded alleles instead of a false clean call —
  i.e. correct, honest ceiling behavior. If 250 bp‑read data exists, expect partial
  recovery of the larger alleles.

## 7. CLI surface (P0)
- `somatiq call --cram X.cram --ref ref.fa --catalog cat.json --sex M --out X.parquet [--triage] [--hipstr-batch manifest.tsv]`
- `somatiq validate-sim --catalog cat.json --out simreport/` (runs 6.1)
- `somatiq catalog-check --catalog cat.json --ref ref.fa` (validates HipSTR regions vs known genotypes)

## 8. Open items / assumptions to confirm at sign‑off
1. **HipSTR region defs** for the 4 loci must be generated and validated (§3.1) — ATXN7
   compound structure and AR strand are the flagged risks.
2. **NA06075 WGS CRAM** availability (else 6.2 deferred; 6.1 is the gate).
3. **Batch policy:** single‑sample `--def-stutter-model` for the first smoke test vs
   requiring a ≥N‑sample batch for any trusted `f`. Proposed default: allow single‑sample
   runs but **WARN** and treat `f` as provisional until a batch stutter model exists.
4. **Exact prancSTR/HipSTR flag names** — locked in build step 1 (§1).
5. Read length is auto‑detected per CRAM for the ceiling flag.
