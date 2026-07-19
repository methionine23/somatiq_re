# Portable P0 image: HipSTR + prancSTR (TRTools) + ExpansionHunter + samtools.
# Built so the pipeline runs in the owner's cloud, next to the CRAMs (no egress here).
FROM mambaorg/micromamba:1.5-jammy

USER root
# Pin exact versions at build time and record them in env.lock (docs/P0_SPEC.md 1).
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        python=3.11 \
        hipstr \
        expansionhunter \
        samtools \
        htslib \
        trtools \
        pysam \
        pandas \
        pyarrow \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /app
COPY . /app
RUN pip install --no-deps -e .

# Smoke test: the CLI and catalog load without the sequencing stack at import time.
RUN somatiq catalog-check --catalog catalogs/somatiq_grch38.json

ENTRYPOINT ["somatiq"]
CMD ["--help"]
