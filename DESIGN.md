# somatiq_re — Somatic Repeat‑Expansion Index from Short‑Read Illumina WGS/WES

**Status:** Research & design (no implementation yet)
**Goal:** A lightweight, scalable Python tool + runner that quantifies a *somatic
repeat‑expansion index* at known tandem‑repeat loci directly from standard
Illumina WGS/WES CRAM, suitable for large genomic datasets (UK Biobank /
All of Us / Regeneron scale).

---

## 1. Problem statement & scope

Repeat‑expansion disorders (HD, the SCAs, DM1, C9orf72, FXS, …) are driven by
unstable tandem repeats that not only are inherited at a given length but keep
mutating *somatically* throughout life. In HD, somatic CAG expansion in blood —
after correcting for inherited CAG length and age — correlates with age‑at‑onset
and, in recent work, with neurodegeneration biomarkers *decades before* motor
onset. "Somatic instability" is therefore both a biomarker and a therapeutic
readout.

The classical way to measure it is **PCR + fragment analysis (GeneMapper) or
MiSeq amplicon**, from which an **"expansion index" / "instability index"** is
computed. Those methods are PCR‑based and locus‑by‑locus. This project asks a
different question:

> Can we recover a comparable *somatic expansion index* from the reads that are
> **already present** in standard PCR‑free WGS (and, more cautiously, PCR‑positive
> WES) CRAMs, at population scale, without any new wet‑lab assay?

Recent large studies (Mukamel & McCarroll, *Nature* 2025, 900k UKB + All of Us;
the McCarroll single‑cell HTT work, *Cell*/*Nat Genet* 2025) show the answer is
"yes, for small‑to‑moderate somatic increments" — this project builds a focused,
reusable, auditable implementation of that idea.

**In scope**
- Targeted analysis at a **user‑supplied catalog of known loci** (motif + coordinates).
- Per‑read repeat‑length estimation from CRAM/BAM.
- A **somatic expansion index** (fragment‑analysis analog) plus a **mosaic‑fraction**
  estimate, with rigorous QC/artifact filtering.
- A **batch runner** that scales to tens/hundreds of thousands of CRAMs.
- PCR‑free WGS as the primary target; a PCR/stutter‑aware mode for WES/amplicon.

**Out of scope (at least v1)**
- De novo discovery of unknown expanded loci (that is EHdn/STRling territory).
- Sizing expansions much larger than the read length (short reads fundamentally
  cannot count individual units beyond a spanning read; only coarse IRR‑based
  bounds are possible — and somatic *index* work targets the sub‑read‑length regime).
- Long‑read (that is a separate, more direct measurement).

---

## 2. What a "somatic repeat‑expansion index" actually is

### 2.1 The classical PCR/fragment‑analysis expansion index (the thing to emulate)
From HD literature (Lee et al. 2011/2017; Ciosi et al. 2019, *eLife* 64674):

1. Run capillary fragment analysis / MiSeq → a distribution of peaks (one per
   repeat length) with heights ∝ molecule abundance.
2. Define the **main/modal peak** = the inherited (germline) allele (the tallest
   peak in the expanded "tail" allele).
3. Apply a **relative peak‑height threshold** (commonly **10%**, sometimes 5% or a
   conservative 20%) of the modal peak; peaks below it are discarded as noise.
4. **Normalize** each surviving peak height to the sum of surviving heights.
5. **Weight** each normalized height by its **CAG offset** from the modal peak
   (Δrepeats = length − modal).
6. **Sum** → the **expansion index**. `0` = no expansion beyond the inherited
   allele; `>0` = somatic mosaicism, larger = more expansion.

$$ \text{ExpIndex} = \sum_{i:\,h_i \ge \tau\,h_{\max}} \frac{h_i}{\sum_j h_j}\;(\ell_i - \ell_{\text{mode}}) $$

The whole method is deliberately simple and robust — it just needs a clean
length‑vs‑abundance distribution. **Our core insight: a per‑read repeat‑length
histogram from WGS is exactly such a distribution**, with read counts replacing
peak heights. Everything else (thresholding, normalization, offset weighting)
transfers directly.

Two related metrics come from the same CE peak table and both transfer to the
read histogram — worth emitting both:
- **Expansion index** — sum of normalized heights × *positive* offsets (expansion tail only).
- **Instability index** — mean length change vs the inherited allele (signed;
  decomposable into expansion vs contraction components).

**Existing CE/GeneMapper peak‑based software (the reference implementations to
mirror):** yes, this is already tooled on the PCR side. Most directly **TRACE**
(*Tandem Repeats Analysis by Capillary Electrophoresis*; R package + `traceshiny`
web app, MGH, bioRxiv 2026) takes raw fragment‑analysis files → expansion/
instability indices end‑to‑end. Older/related: **TraceTrack** (batch CE
processing) and the classic Monckton‑lab peak‑table "instability index" spreadsheet
(Lee et al. 2010). **We treat TRACE's index definitions as the spec** for
somatiq_re's histogram metrics, so a WGS index and a CE index are computed the same
way and are directly comparable on matched samples.

### 2.2 The mosaic‑fraction / mixture view (prancSTR analog)
Instead of (or in addition to) a summary index, model the reads as a **mixture**:
a germline component (one or two alleles + stutter) plus a low‑frequency
**mosaic/expanded component**, and estimate the **mosaic allele fraction `f`** and
mosaic allele length by maximum likelihood, with a p‑value for `H0: f = 0`. This
is what prancSTR does; it gives a principled, depth‑aware alternative to the
heuristic index and a significance test.

Both views are worth producing — the fragment‑analysis index for continuity with
existing HD literature, the mixture model for statistics and calibration.

---

## 3. Landscape review — existing tools and their logic

| Tool | Purpose | Read evidence used | Diploid genotype? | Somatic/mosaic? | Notes for us |
|---|---|---|---|---|---|
| **ExpansionHunter** (Illumina) | Targeted sizing of *known* loci | spanning + flanking + in‑repeat reads (IRR), sequence‑graph realignment | Yes (ML diploid) | **No** (explicitly cannot resolve mosaicism) | Gold standard for germline sizing; its graph realignment yields per‑read repeat spans we could reuse. Catalog format is a good model. |
| **ExpansionHunter Denovo** | *Discover* novel expansions | anchored IRRs, genome‑wide | No precise size | No | Discovery only; not sizing. Out of scope but good for catalog building. |
| **GangSTR** | Genome‑wide TR genotyping | spanning + mate‑pair insert size + flanking + IRR | Yes | No | Uses insert size to extend range; slow (~10 h/genome). |
| **HipSTR** | Genome‑wide STR genotyping + phasing | spanning reads only (read‑length limited); **learns per‑locus stutter model**; emits **per‑read allele lengths (ALLREADS/MALLREADS)** | Yes | No (but its per‑read output + stutter model feed prancSTR) | The most reuse‑friendly upstream: per‑read lengths + calibrated stutter are precisely what a somatic index needs. |
| **prancSTR** (TRTools, Gymrek lab) | **Mosaic STR detection, single sample, no matched normal** | per‑read lengths from HipSTR VCF + stutter model | — | **Yes**: MLE mosaic fraction `f`, mosaic allele `C`, p‑value `f=0` | Closest existing analog to our goal. Power best at f≈10–20% @ 30–50× PCR‑free; down to ~7%. Works PCR± *if* stutter params are accurate. |
| **STRling / superSTR / STRetch** | Fast screening / discovery | k‑mer / alignment‑free | No | No | Screening, not precise sizing. |
| **REViewer / GraphAlignmentViewer** | Visualize EH pileups | — | — | — | Reuse for QC visualization of tail reads. |
| **Fragment analysis / MiSeq expansion index** (Monckton/Lee/Ciosi) | Somatic index (PCR) | amplicon peaks | — | **Yes** (the index we emulate) | Defines the metric; not sequence‑from‑CRAM. |
| **Mukamel & McCarroll** (*Nature* 2025) | Population somatic expansion from WGS | read‑level CAG lengths at ~15 CAG loci in blood; **base‑quality drop to flag PCR‑slippage reads** | — | **Yes** (population scale) | The proof that this works from plain WGS; source of the key PCR‑artifact filter idea and GWAS‑of‑expansion framing. |

**Two logical families**, then: (a) *genotypers* that produce per‑read length
evidence (EH, GangSTR, HipSTR), and (b) *somatic estimators* that consume that
evidence (prancSTR; the fragment‑analysis index; the McCarroll read‑distribution
approach). **somatiq_re sits in family (b)** and needs a family‑(a) substrate —
either an existing genotyper's per‑read output or our own pysam extraction.

---

## 4. Key methodological building blocks

### 4.1 Read classes (from ExpansionHunter's framework)
- **Spanning reads** — repeat fully inside the read, non‑repeat anchors on *both*
  sides ⇒ **exact per‑read repeat count**. The workhorse for somatic index, because
  somatic increments in blood are usually small (a few units) and stay within
  spanning range for typical germline alleles.
- **Flanking reads** — anchor on one side only ⇒ a **lower bound** on length.
- **In‑repeat reads (IRRs)** — entirely repeat ⇒ evidence of an allele **longer than
  the read**; gives only a coarse count (via IRR abundance). Beyond precise per‑read
  sizing, so used only as a "large‑expansion present" flag, not in the index.

Detectable spanning size ≈ `(read_len − 2·min_anchor) / motif_len`. For 150 bp
reads, CAG (3 bp), 8 bp anchors ⇒ ~44 repeats can be spanned. This bounds where a
per‑read index is meaningful; it fits the small‑somatic‑increment regime well.

### 4.2 Per‑read repeat sizing (our own path)
For each read overlapping the locus:
1. Locate the anchored flanks (align/match the known left/right flank sequence).
2. Count motif copies in the intervening span, tolerant of interruptions and
   sequencing mismatches (motif‑aware scan; a small banded alignment for accuracy).
3. Record: repeat count, read class, MAPQ, per‑base qualities across the repeat,
   soft‑clip/indel flags, strand, whether it is properly paired.

### 4.3 The somatic metrics we output
- **Germline modal allele(s)** — one/two modes from the per‑read histogram.
- **Expansion index** (§2.1) per allele, on the read‑count histogram.
- **Mosaic fraction `f`** and expanded‑allele length via a mixture MLE (§2.2),
  with p‑value.
- **Tail summaries** — fraction of reads > mode+k, mean/skew of the expansion tail,
  max observed length.
- Everything reported **alongside germline length, depth, and QC** (see confounders).

### 4.4 Stutter / error model
- **PCR‑free WGS**: polymerase‑slippage stutter is low but nonzero and asymmetric
  (contraction‑biased). Model it as a per‑locus, per‑motif distribution centered on
  the true allele; the modal ±1 units are largely stutter and must not be counted as
  somatic expansion. Calibrate empirically from homozygous/short alleles and from
  non‑expanding control loci.
- **PCR‑positive WES/amplicon**: stutter is much larger and length‑dependent; borrow
  HipSTR's learned stutter model or fit our own per‑batch. This is the regime where
  "existing methods deal with PCR data" — we support it but flag it as lower‑confidence.

### 4.5 The PCR‑artifact filter (McCarroll trick)
PCR slippage that changes a read's apparent repeat length tends to leave a
**predictable drop in base‑quality** within/after the repeat. Reads showing that
signature are excluded before building the histogram — this is what lets plain WGS
distinguish *biological* somatic expansion from *artefactual* length changes, and is
essential for WES/PCR‑positive data.

---

## 5. Core challenges & confounders (short‑read WGS/WES specifics)

1. **Read‑length ceiling.** Only sub‑read‑length expansions get exact per‑read
   counts. Fine for small somatic increments; not for measuring how far the *rare*
   long‑expanded molecules go (the "armadillo tail" seen in single neurons). Report
   the ceiling explicitly; never silently truncate.
2. **Depth is the resolution limit.** The smallest resolvable mosaic fraction ≈
   1/(spanning‑read depth). 30–50× WGS → a few tens of spanning reads → sensitivity
   floor ~7–20% (matches prancSTR). WES coverage is uneven; some loci will be
   under‑powered — must be depth‑gated and reported.
3. **Depth‑dependence of the index.** More reads sample more tail → naive index
   grows with depth. Must **normalize for depth** (downsample to a fixed count, or
   model‑based) so samples are comparable.
4. **Length‑ & age‑dependence.** Longer inherited alleles are intrinsically more
   unstable; somatic expansion rises with age. The index is only interpretable
   **stratified by / regressed on inherited length and age** — a first‑class output
   caveat, not an afterthought.
5. **Allele assignment (diploid).** Somatic expansion rides mostly on the
   longer/expanded allele. Tail reads must be attributed to the correct parental
   allele (by length separation; easy for well‑separated HD carriers, hard for
   near‑equal alleles).
6. **Batch / chemistry effects.** Stutter and quality profiles vary by prep and
   sequencer. Population studies need per‑batch stutter calibration and batch as a
   covariate.
7. **CRAM specifics.** CRAM needs the exact reference for decode; base qualities may
   be **binned** (e.g. Illumina 8‑bin), which weakens the base‑quality artifact
   filter — detect and adapt.
8. **Reference/catalog correctness.** Off‑target homologous repeats, wrong
   coordinates, or unlisted interruptions corrupt counts. Interruption‑aware motif
   matching + a curated catalog matter.
9. **Scale/throughput.** Hundreds of thousands of CRAMs ⇒ targeted region fetch only
   (never full‑genome scan), streaming, parallel per‑sample/per‑locus, resumable,
   with per‑CRAM cost measured in seconds‑to‑minutes.

---

## 6. Proposed design

### 6.1 Architecture — pluggable front‑ends, one common back‑end (decided)
The somatic‑index/mosaic‑fraction math operates on a **common per‑read length
table** (columns: read_id, allele‑guess, repeat_count, read_class, mapq,
qual_across_repeat, flags). Multiple front‑ends can populate that table, so we are
not locked to one substrate:

- **FE‑1: pysam raw CRAM extraction (primary, self‑contained).** Targeted region
  fetch + our own per‑read sizing. No C++ dependency, fully auditable, full control
  over the QC/artifact filters, easy to scale. Cost: we reimplement + validate
  sizing.
- **FE‑2: parse ExpansionHunter's *realigned* BAM (fast reuse — you already have EH
  results).** EH (with `--analysis-mode`/REViewer output) emits a realigned BAM of
  the reads it used, with graph alignments spanning the repeat. Parsing that gives
  per‑read repeat spans **from EH's validated realignment** — better than naïve
  re‑extraction for reads near/over the read length, and near‑free since you have EH
  outputs. Caveat: only reads EH retained near the locus; must be requested at EH
  run time.
- **FE‑3: HipSTR per‑read output (`ALLREADS`/`MALLREADS`) + its learned stutter
  model.** Feeds a prancSTR‑style mixture estimator directly and supplies a
  calibrated per‑locus stutter model (valuable for the PCR‑WES mode). Run HipSTR only
  where needed.

**Plan:** build FE‑1 as the reference path and **FE‑2 immediately** (reuses your
existing EH results, and cross‑checks FE‑1's sizing on the same reads). FE‑3 is the
stutter‑calibration/PCR‑WES helper. All three emit the identical per‑read table, so
the back‑end and every metric are shared and front‑ends cross‑validate each other.

### 6.2 Pipeline (architecture A)
```
catalog(JSON/BED: locus, motif, ref coords, flanks, inherited-length prior)
        │
   [1] region fetch  ── pysam, CRAM+ref, locus ± window
        │
   [2] read classify ── spanning / flanking / IRR
        │
   [3] per-read size ── motif-aware count + interruption handling
        │
   [4] QC filter    ── MAPQ, proper-pair, base-quality artifact (McCarroll),
        │              soft-clip sanity, strand balance, dup removal
   [5] histogram    ── read-count vs repeat-length, per allele
        │
   [6] germline call ── modal allele(s); stutter deconvolution
        │
   [7] metrics       ── expansion index + mosaic-fraction MLE + tail stats
        │
   [8] emit          ── per-sample × per-locus TSV/Parquet (+ optional per-read dump, QC plot)
```

### 6.3 Interfaces
- **Input:** CRAM/BAM + reference FASTA + locus catalog + sample manifest.
- **Core CLI:** `somatiq call --cram X.cram --ref ref.fa --catalog loci.json --out X.tsv`
- **Runner:** `somatiq run --manifest samples.tsv ...` — parallel (multiprocessing
  locally; array‑job / Nextflow‑ or Snakemake‑friendly on cluster/cloud), resumable
  (skip completed), per‑sample logs, aggregate to Parquet.
- **Output columns (per sample × locus):** modal/inherited length(s), depth,
  spanning/flanking/IRR counts, expansion_index, mosaic_fraction `f`, p‑value,
  tail fraction, QC flags, filter‑pass. Machine‑friendly for downstream GWAS.

### 6.4 Dependencies (keep it lightweight)
`pysam`, `numpy`, `scipy` (MLE/optimization), `pandas`/`pyarrow` (I/O), `click`/`argparse`.
Optional: `matplotlib` for QC pileup/histogram plots. No heavy frameworks in the core.

---

## 7. Decisions (resolved) & v1 loci

**Resolved with the project owner:**
1. **v1 loci:** TCF4, AR, DMPK, ATXN7 (table below).
2. **Data regimes:** both PCR‑free WGS **and** PCR‑positive WES in v1 → the stutter
   model + base‑quality artifact filter are v1 requirements, not later add‑ons.
3. **Architecture:** self‑contained pysam (FE‑1) **plus** consume existing EH output
   (FE‑2, the realigned BAM); run HipSTR (FE‑3) only where a calibrated stutter model
   is needed. (See §6.1.)
4. **Metrics:** both — the CE‑style expansion/instability index (TRACE‑compatible)
   **and** a prancSTR‑style mosaic fraction.
5. **Validation:** NA06075 (known mosaic DMPK control) + reproduce All of Us TCF4
   mosaicism + an internal cohort. (See §8.)
6. **Target:** reproduce the All of Us TCF4 somatic‑mosaicism signal, then apply to
   the internal cohort.

### v1 locus table
| Locus | Disease | Motif | Chrom/ploidy | Typical inherited range | Notes for sizing |
|---|---|---|---|---|---|
| **TCF4** (CTG18.1) | Fuchs endothelial dystrophy | CTG/CAG | chr18, diploid | ~10–40 (common), 40–75+ expanded | **Primary validation target** (AoU high blood mosaicism). Common alleles fit within‑read; large pathogenic alleles exceed 150 bp reads. |
| **AR** | SBMA (Kennedy) | CAG | chrX, **hemizygous in males** | ~9–36 | Sex‑aware ploidy: single germline mode in males (no 2nd allele anchor); diploid in females. |
| **DMPK** | Myotonic dystrophy 1 | CTG | chr19, diploid | 5–34 normal; 35–49 premut; 50+ DM1 | **NA06075 control = 12/56/70 CTG mosaic.** 56–70 CTG ≈ 168–210 bp **> 150 bp read** → tests the read‑length ceiling; expect flanking/IRR, not clean spanning, for the expanded alleles. |
| **ATXN7** | SCA7 | CAG | chr3, diploid | 4–33 normal; 37+ SCA7 | CAG; normal/short‑expanded alleles fit within‑read. |

**Cross‑cutting sizing note:** at 2×150 bp (UKB/AoU/NovaSeq), the spanning ceiling is
~(150−2·anchor)/motif ≈ 40–44 units. Somatic *index* in blood targets small
increments over the modal/common allele, which fits for AR, ATXN7, and common TCF4;
large pathogenic DMPK/TCF4 alleles (and NA06075's 56/70) fall to flanking/IRR and get
only coarse bounds — report the ceiling per read‑length, never silently truncate.

---

## 8. Validation plan
- **Positive control sample — NA06075** (Coriell/NIST DM1 reference, mosaic DMPK
  **12/56/70 CTG**). Expect: germline modal ~12 called cleanly from spanning reads;
  the 56/70 mosaic alleles exceed a 150 bp read → validates the read‑length‑ceiling
  behavior (flanking/IRR flag rather than a false clean call). If long‑read/PCR truth
  on NA06075 is available, compare directly.
- **Primary reproduction target — TCF4 in All of Us.** Reproduce the AoU TCF4 blood
  somatic‑mosaicism signal (common‑allele length mosaicism), then run the internal
  cohort. TCF4 is a high‑mosaicism locus, so it is the strongest signal to confirm the
  pipeline end‑to‑end.
- **CE cross‑check** — where matched fragment‑analysis/MiSeq exists, compare the WGS
  index to a **TRACE** index computed on the same samples (same index definition ⇒
  apples‑to‑apples).
- **Simulation** — spike synthetic mosaic alleles at known `f` (custom or TRTools'
  simuSTR) into real backgrounds; recover `f`/index; establish the sensitivity‑vs‑depth
  curve per locus and read length.
- **Cross‑tool** — germline modal vs ExpansionHunter (FE‑2 reuse); mosaic vs prancSTR
  (FE‑3) on the same CRAMs.
- **WGS↔WES concordance** — for samples with both, compare indices; quantify the
  PCR‑stutter penalty and confirm the base‑quality artifact filter closes the gap.
- **Reproducibility** — replicate CRAMs / technical duplicates; batch/chemistry check.

## 9. Phased roadmap
- **P0** — catalog format + pysam region fetch + per‑read spanning‑read sizing for
  HTT‑CAG; raw histogram out. (Prove the substrate.)
- **P1** — QC/artifact filters + stutter deconvolution + expansion index + depth
  normalization. (Prove the metric.)
- **P2** — mosaic‑fraction mixture MLE + p‑value; sensitivity/depth calibration on
  simulation. (Prove the statistics.)
- **P3** — batch runner, Parquet aggregation, multi‑locus panel, WES/PCR mode.
  (Prove the scale.)
- **P4** — cross‑tool/orthogonal validation + docs. (Prove it's right.)

## 10. Key references
- Dolzhenko et al. 2017 *Genome Res* (ExpansionHunter, PCR‑free long expansions);
  2019 *Bioinformatics* (sequence‑graph EH); 2020 *Genome Biol* (EHdn); 2022 (REViewer).
- Mousavi et al. 2019 *NAR* (GangSTR). Willems et al. 2017 *Nat Methods* (HipSTR).
- Tandon et al. 2024 *Bioinformatics* btae485 (prancSTR / TRTools mosaicism).
- Mukamel & McCarroll 2025 *Nature* (900k biobank somatic expansion); Handsaker/Kashin/Reed
  et al. 2025 *Cell*/*Nat Genet* (single‑cell HTT somatic expansion).
- Ciosi et al. 2019 *eLife* 64674 (somatic expansion index over life); Lee et al.
  2010/2017 (expansion index method); *Nat Med* 2024 (blood somatic CAG ↔ HD biomarkers).
- **TRACE** 2026 bioRxiv (`traceshiny.mgh.harvard.edu`) & TraceTrack 2023 (CE peak‑based
  SI/expansion index software — the metric spec we mirror).
- Trost et al. 2024 *PLOS One* (STR tool comparison).
- Kalman et al. 2013 (NIST/CDC DM1 genomic DNA reference panel; NA06075 = 12/56/70 CTG).
