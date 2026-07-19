"""somatiq command-line interface (P0)."""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .catalog import load_catalog, validate_catalog
from .hipstr import write_hipstr_regions


def _cmd_catalog_check(args) -> int:
    cat = load_catalog(args.catalog)
    if args.chr_style:
        cat = cat.with_chrom_style(args.chr_style)
    warnings = validate_catalog(cat)
    print(f"catalog OK: {len(cat)} loci, build={cat.genome_build}, chrom_style={cat.chrom_style}")
    for locus in cat.loci:
        r = locus.hipstr_region
        print(
            f"  {locus.locus_id:8s} {locus.chrom}:{locus.ref_start}-{locus.ref_end} "
            f"ref_motif={locus.ref_motif} disease_motif={locus.disease_motif} "
            f"ploidy={locus.ploidy} hipstr[{r.start},{r.end},p{r.period},n{r.ref_copies}]"
        )
    if warnings:
        print("\nwarnings (validate before trusting these loci):", file=sys.stderr)
        for w in warnings:
            print(f"  ! {w}", file=sys.stderr)
    return 0


def _cmd_regions(args) -> int:
    cat = load_catalog(args.catalog)
    if args.chr_style:
        cat = cat.with_chrom_style(args.chr_style)
    write_hipstr_regions(cat, args.out)
    print(f"wrote HipSTR regions: {args.out}")
    return 0


def _cmd_call(args) -> int:
    # Imported lazily so catalog-check/regions work without the CRAM/subprocess stack.
    from .call import RunConfig, run_sample
    from .qc import QCConfig
    from .schema import write_parquet, write_tsv

    cat = load_catalog(args.catalog)
    if args.chr_style:
        cat = cat.with_chrom_style(args.chr_style)
    sample = args.sample or _infer_sample(args.alignment)
    cfg = RunConfig(
        data_mode=args.data_mode,
        def_stutter_model=not args.no_def_stutter_model,
        eh_annotate=args.eh_annotate,
        eh_vcf=args.eh_vcf,
        qc=QCConfig(min_informative_reads=args.min_informative_reads),
        keep_intermediate=args.keep_intermediate,
    )
    rows = run_sample(
        alignment=args.alignment,
        reference=args.ref,
        catalog=cat,
        sample=sample,
        sex=args.sex,
        config=cfg,
        workdir=args.workdir,
    )
    if args.out_format == "parquet":
        write_parquet(rows, args.out)
    else:
        write_tsv(rows, args.out)
    print(f"wrote {len(rows)} rows -> {args.out}")
    return 0


def _infer_sample(path: str) -> str:
    base = path.rsplit("/", 1)[-1]
    for ext in (".cram", ".bam", ".sam"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="somatiq", description="Somatic repeat-expansion index (P0)")
    p.add_argument("--version", action="version", version=f"somatiq {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    cc = sub.add_parser("catalog-check", help="load + validate the locus catalog")
    cc.add_argument("--catalog", required=True)
    cc.add_argument("--chr-style", choices=["chr", "nochr"], default=None)
    cc.set_defaults(func=_cmd_catalog_check)

    rg = sub.add_parser("regions", help="write the HipSTR region file")
    rg.add_argument("--catalog", required=True)
    rg.add_argument("--out", required=True)
    rg.add_argument("--chr-style", choices=["chr", "nochr"], default=None)
    rg.set_defaults(func=_cmd_regions)

    ca = sub.add_parser("call", help="score one sample (HipSTR -> prancSTR)")
    ca.add_argument("--cram", "--bam", dest="alignment", required=True, help="CRAM/BAM path")
    ca.add_argument("--ref", required=True, help="reference FASTA (CRAM decode)")
    ca.add_argument("--catalog", required=True)
    ca.add_argument("--sample", default=None, help="SM tag / VCF sample (default: from filename)")
    ca.add_argument("--sex", default=None, help="M/F (required if any locus is sex_dependent)")
    ca.add_argument("--out", required=True)
    ca.add_argument("--out-format", choices=["tsv", "parquet"], default="tsv")
    ca.add_argument("--data-mode", choices=["wgs"], default="wgs", help="exome mode arrives in P2")
    ca.add_argument("--eh-annotate", action="store_true", default=False,
                    help="attach annotation from a pre-computed EH VCF (default off)")
    ca.add_argument("--eh-vcf", default=None, help="pre-computed ExpansionHunter VCF")
    ca.add_argument("--min-informative-reads", type=int, default=20)
    ca.add_argument("--no-def-stutter-model", action="store_true", default=False,
                    help="disable HipSTR default stutter model (requires a batch)")
    ca.add_argument("--chr-style", choices=["chr", "nochr"], default=None)
    ca.add_argument("--workdir", default=None)
    ca.add_argument("--keep-intermediate", action="store_true", default=False)
    ca.set_defaults(func=_cmd_call)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
