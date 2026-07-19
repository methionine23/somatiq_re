"""Parser tests: MALLREADS, HipSTR VCF, prancSTR .tab, EH VCF."""
from conftest import CATALOG_PATH

from somatiq.catalog import load_catalog
from somatiq.eh import parse_eh_vcf
from somatiq.hipstr import match_call_to_locus, parse_hipstr_vcf
from somatiq.mallreads import parse_read_field
from somatiq.prancstr import parse_prancstr_lines


# ---- MALLREADS ----
def test_parse_read_field_basic():
    h = parse_read_field("0|25;6|20;12|2", period=3, ref_copies=20)
    assert h.informative_reads == 47
    assert h.by_bpdiff == {0: 25, 6: 20, 12: 2}
    cps = h.by_copies()
    assert cps[20.0] == 25 and cps[22.0] == 20 and cps[24.0] == 2
    assert h.max_copies() == 24.0
    assert h.modal_bpdiff() == 0


def test_parse_read_field_empty():
    for v in (".", "", None):
        h = parse_read_field(v, 3, 20)
        assert h.informative_reads == 0
        assert h.detectable_f_floor() is None


def test_parse_read_field_malformed_tokens_skipped():
    h = parse_read_field("0|10;garbage;6|x;3|5", 3, 20)
    assert h.by_bpdiff == {0: 10, 3: 5}


def test_detectable_floor():
    h = parse_read_field("0|50", 3, 20)
    assert abs(h.detectable_f_floor() - 0.02) < 1e-9


# ---- HipSTR VCF ----
HIPSTR_VCF = "\n".join(
    [
        "##fileformat=VCFv4.1",
        '##INFO=<ID=INFRAME_PGEOM,Number=1,Type=Float,Description="x">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878",
        "\t".join(
            [
                "chr19",
                "45770204",
                "DMPK",
                "CAG" * 20,
                "CAG" * 22,
                ".",
                ".",
                "INFRAME_UP=0.05;INFRAME_DOWN=0.05;INFRAME_PGEOM=0.9;END=45770264",
                "GT:GB:Q:DP:DSTUTTER:DFLANKINDEL:MALLREADS",
                "0|1:60|66:0.99:47:3:0:0|25;6|20;12|2",
            ]
        ),
    ]
)


def test_parse_hipstr_vcf():
    import io
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f:
        f.write(HIPSTR_VCF)
        path = f.name
    calls = parse_hipstr_vcf(path, "NA12878")
    assert len(calls) == 1
    c = calls[0]
    assert c.chrom == "chr19" and c.pos == 45770204
    assert c.gt == (0, 1)
    assert c.q == 0.99 and c.dp == 47 and c.dstutter == 3
    assert c.inframe_pgeom == 0.9 and c.stutter_learned
    assert c.gb_bp == (60, 66)
    assert c.allele_copies(3) == [20.0, 22.0]
    assert len(c.called_allele_seqs()) == 2


def test_match_call_to_locus():
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f:
        f.write(HIPSTR_VCF)
        path = f.name
    calls = parse_hipstr_vcf(path, "NA12878")
    cat = load_catalog(CATALOG_PATH)
    assert match_call_to_locus(calls, cat.get("DMPK")) is not None
    assert match_call_to_locus(calls, cat.get("TCF4")) is None


# ---- prancSTR .tab ----
PRANC_TAB = "\n".join(
    [
        "\t".join(
            [
                "sample", "chrom", "pos", "locus", "motif", "A", "B", "C", "f", "pval",
                "reads", "mosaic_support", "stutter parameter u", "stutter paramter d",
                "stutter paramter rho", "quality factor", "read depth",
            ]
        ),
        "\t".join(
            ["NA12878", "chr19", "45770204", "DMPK", "CAG", "20", "22", "24",
             "0.08", "0.001", "47", "4", "0.05", "0.05", "0.9", "0.99", "47"]
        ),
    ]
)


def test_parse_prancstr_lines():
    rows = parse_prancstr_lines(PRANC_TAB.splitlines())
    assert len(rows) == 1
    r = rows[0]
    assert r.sample == "NA12878" and r.locus == "DMPK"
    assert r.allele_a == 20.0 and r.allele_b == 22.0
    assert r.mosaic_c == 24.0 and r.f == 0.08 and r.pval == 0.001
    assert r.mosaic_support == 4 and r.read_depth == 47


# ---- EH VCF (optional annotation) ----
EH_VCF = "\n".join(
    [
        "##fileformat=VCFv4.1",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878",
        "chr19\t45770204\t.\tC\t<STR22>\t.\tPASS\tREPID=DMPK;VARID=DMPK\tGT:REPCN\t0/1:20/22",
    ]
)


def test_parse_eh_vcf():
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f:
        f.write(EH_VCF)
        path = f.name
    ann = parse_eh_vcf(path, "NA12878", thresholds={"DMPK": 50})
    assert "DMPK" in ann
    assert ann["DMPK"]["gt_short"] == 20 and ann["DMPK"]["gt_long"] == 22
    assert ann["DMPK"]["carrier"] is False
