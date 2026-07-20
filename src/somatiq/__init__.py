"""somatiq_re: somatic repeat-expansion index from short-read WGS/WES.

P0 = prancSTR-first pipeline: HipSTR (germline + stutter) -> prancSTR (mosaic f),
with optional annotation from pre-computed ExpansionHunter results.

The core logic (catalog, VCF/MALLREADS/.tab parsing, interruption detection, QC,
score-row assembly) is stdlib-only and unit-tested. Subprocess execution of HipSTR
and prancSTR, and CRAM read access via pysam, are thin wrappers run where the data
lives.
"""

__version__ = "0.1.0"
