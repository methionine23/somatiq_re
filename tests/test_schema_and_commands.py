from conftest import CATALOG_PATH

from somatiq.call import assemble_rows, match_pranc_to_locus
from somatiq.catalog import load_catalog
from somatiq.hipstr import (
    HipstrCall,
    build_hipstr_command,
    write_hipstr_regions,
)
from somatiq.prancstr import PrancResult, build_prancstr_command
from somatiq.schema import ScoreRow, assemble_score_row, write_tsv


def _dmpk_call(mallreads="0|25;6|20;12|2", alt_copies=22, learned=True):
    return HipstrCall(
        chrom="chr19",
        pos=45770204,
        ref_seq="CAG" * 20,
        alt_seqs=("CAG" * alt_copies,),
        gt=(0, 1),
        gb_bp=(60, alt_copies * 3),
        q=0.99,
        dp=47,
        dstutter=3,
        dflankindel=0,
        mallreads=mallreads,
        allreads=mallreads,
        inframe_up=0.05,
        inframe_down=0.05,
        inframe_pgeom=0.9 if learned else None,
    )


def _dmpk_pranc(f=0.08, c=24.0):
    return PrancResult(
        sample="NA12878", chrom="chr19", pos=45770204, locus="DMPK", motif="CAG",
        allele_a=20.0, allele_b=22.0, mosaic_c=c, f=f, pval=0.001, reads=47,
        mosaic_support=4, read_depth=47, raw={},
    )


def test_assemble_score_row_full():
    cat = load_catalog(CATALOG_PATH)
    dmpk = cat.get("DMPK")
    row = assemble_score_row(
        sample="NA12878", locus=dmpk, sex="M", data_mode="wgs",
        hipstr=_dmpk_call(), pranc=_dmpk_pranc(), read_length=150,
    )
    assert row.prancstr_f == 0.08
    assert row.prancstr_mosaic_allele == 24.0
    assert row.germline_gt_copies == "20|22"
    assert row.informative_reads == 47
    assert row.stutter_model_source == "learned"
    assert row.allele_exceeds_readlen is False
    assert row.interruption_flag is False
    assert row.qc_level1 == "PASS"


def test_ceiling_flag_when_allele_exceeds_readlen():
    cat = load_catalog(CATALOG_PATH)
    dmpk = cat.get("DMPK")
    # a 56-copy allele = 168 bp > 150 bp read (NA06075-like)
    row = assemble_score_row(
        sample="s", locus=dmpk, sex="M", data_mode="wgs",
        hipstr=_dmpk_call(mallreads="0|30;108|8", alt_copies=56), pranc=None,
        read_length=150,
    )
    assert row.allele_exceeds_readlen is True
    assert row.qc_level1 in ("WARN", "FAIL")
    assert "allele_exceeds_readlen" in row.qc_reasons


def test_interruption_annotation_in_row():
    cat = load_catalog(CATALOG_PATH)
    dmpk = cat.get("DMPK")
    hip = _dmpk_call()
    hip = HipstrCall(**{**hip.__dict__, "alt_seqs": ("CAG" * 10 + "CCG" + "CAG" * 9,)})
    row = assemble_score_row(
        sample="s", locus=dmpk, sex="M", data_mode="wgs",
        hipstr=hip, pranc=None, read_length=150,
    )
    assert row.interruption_flag is True
    assert "CCG" in (row.interruption_seq or "")


def test_default_stutter_downgrades_to_warn():
    cat = load_catalog(CATALOG_PATH)
    row = assemble_score_row(
        sample="s", locus=cat.get("DMPK"), sex="M", data_mode="wgs",
        hipstr=_dmpk_call(learned=False), pranc=_dmpk_pranc(), read_length=150,
    )
    assert row.stutter_model_source == "default"
    assert row.qc_level1 == "WARN"


def test_assemble_rows_all_loci_even_when_missing():
    cat = load_catalog(CATALOG_PATH)
    rows = assemble_rows(
        sample="NA12878", sex="M", catalog=cat,
        calls=[_dmpk_call()], pranc_rows=[_dmpk_pranc()],
        read_length=150, data_mode="wgs",
    )
    assert len(rows) == 4  # one row per locus
    by = {r.locus_id: r for r in rows}
    assert by["DMPK"].prancstr_f == 0.08
    # loci with no HipSTR call still get a row (FAIL: no reads)
    assert by["TCF4"].qc_level1 == "FAIL"


def test_match_pranc_by_name():
    cat = load_catalog(CATALOG_PATH)
    assert match_pranc_to_locus([_dmpk_pranc()], cat.get("DMPK")) is not None
    assert match_pranc_to_locus([_dmpk_pranc()], cat.get("TCF4")) is None


def test_write_tsv_roundtrip(tmp_path):
    cat = load_catalog(CATALOG_PATH)
    rows = assemble_rows(
        sample="NA12878", sex="M", catalog=cat,
        calls=[_dmpk_call()], pranc_rows=[_dmpk_pranc()],
        read_length=150, data_mode="wgs",
    )
    out = tmp_path / "out.tsv"
    write_tsv(rows, out)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == ScoreRow.columns()
    assert len(lines) == 5  # header + 4 loci


# ---- command construction ----
def test_build_hipstr_command():
    cmd = build_hipstr_command(["a.cram"], "ref.fa", "r.bed", "o.vcf.gz")
    assert cmd[0] == "HipSTR"
    assert "--def-stutter-model" in cmd
    assert "--str-vcf" in cmd and "o.vcf.gz" in cmd


def test_build_prancstr_command():
    cmd = build_prancstr_command("in.vcf.gz", "pref", samples=["NA12878"])
    assert cmd[0] == "prancSTR"
    assert "--vcftype" in cmd and "hipstr" in cmd
    assert "--readfield" in cmd and "MALLREADS" in cmd
    assert "NA12878" in cmd


def test_write_hipstr_regions(tmp_path):
    cat = load_catalog(CATALOG_PATH)
    p = write_hipstr_regions(cat, tmp_path / "r.bed")
    lines = p.read_text().splitlines()
    assert len(lines) == 4
    cols = lines[0].split("\t")
    assert len(cols) == 6  # chrom start end period ref_copies name
