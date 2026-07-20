# Deploying somatiq_re on a VM without GitHub access

You do not need GitHub to run this — you need the **source archive** or a **container/
conda image**. Pick a path by what the target VM can reach. `somatiq --version` should
report **0.1.0**.

> **Reference build:** GRCh38 only. **hg19 / GRCh37 is not supported or tested in
> v0.1.0** — the shipped `catalogs/somatiq_grch38.json` carries GRCh38 coordinates.
> Running against an hg19 CRAM will mis-locate the loci. An hg19 catalog can be added
> later but will not be validated in this version.

## What every path needs at runtime (on the VM, with your data)
- the CRAM/BAM **and its index** (`.crai`/`.bai`) **co-located** (targeted access;
  the full CRAM is never scanned),
- the **reference FASTA** matching the alignment (GRCh38) for CRAM decode.

---

## Get the source there
Transfer `somatiq_re-src.tar.gz` (git archive of the repo; no top-level folder), then:
```bash
mkdir somatiq_re && tar xzf somatiq_re-src.tar.gz -C somatiq_re && cd somatiq_re
```

## Path A — VM reaches PyPI + conda (just not GitHub)
```bash
conda env create -f environment.yml          # HipSTR, prancSTR(trtools), EH, samtools, pysam
conda activate somatiq_re
pip install --no-deps .                       # the somatiq CLI
somatiq call --cram S.cram --ref GRCh38.fa \
  --catalog catalogs/somatiq_grch38.json --sample S --sex M --out S.tsv
```

## Path B — fully offline VM, via container (recommended)
On a machine **with** internet (the build pulls bioconda):
```bash
docker build -t somatiq_re:0.1.0 .
docker save somatiq_re:0.1.0 | gzip > somatiq_re-image.tar.gz   # move this to the VM
```
On the air-gapped VM:
```bash
docker load < somatiq_re-image.tar.gz
docker run --rm -v /data:/data somatiq_re:0.1.0 call \
  --cram /data/S.cram --ref /data/GRCh38.fa \
  --catalog /app/catalogs/somatiq_grch38.json --sample S --sex M --out /data/S.tsv
```

## Path B′ — fully offline VM, via conda-pack (no Docker)
On a connected machine:
```bash
conda env create -f environment.yml
conda activate somatiq_re && pip install --no-deps .
conda install -c conda-forge conda-pack && conda pack -n somatiq_re -o somatiq_re-env.tar.gz
```
On the offline VM (no conda needed there):
```bash
mkdir -p somatiq_env && tar xzf somatiq_re-env.tar.gz -C somatiq_env
source somatiq_env/bin/activate
conda-unpack
somatiq call --cram S.cram --ref GRCh38.fa --catalog <catalog> --sample S --sex M --out S.tsv
```

## Path C — Python only, no external tools (works anywhere with Python 3.10+)
The core (catalog, parsing, interruption detection, QC, scoring) is **stdlib-only**. If
you already produce HipSTR + prancSTR outputs elsewhere, install core-only and use the
library to parse + QC + score them — no HipSTR/prancSTR/pysam required:
```bash
pip install .                                  # zero third-party deps
somatiq catalog-check --catalog catalogs/somatiq_grch38.json
python - <<'PY'
from somatiq.hipstr import parse_hipstr_vcf
from somatiq.prancstr import parse_prancstr_tab
from somatiq.call import assemble_rows           # merge into scored rows
print("core OK")
PY
```

## Reproducibility / air-gapped mirroring
- Pin exact versions once, on the connected host: `conda env export > env.lock`, and
  record HipSTR/prancSTR/EH `--version` in `docs/tool_flags.md`.
- Verify the install on the VM: `somatiq --version` (→ 0.1.0), `pytest -q` (core tests),
  `somatiq catalog-check --catalog catalogs/somatiq_grch38.json`.
