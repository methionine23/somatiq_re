import pytest
from conftest import CATALOG_PATH

from somatiq.catalog import load_catalog, validate_catalog


def test_shipped_catalog_loads():
    cat = load_catalog(CATALOG_PATH)
    assert len(cat) == 4
    assert {l.locus_id for l in cat} == {"TCF4", "AR", "DMPK", "ATXN7"}
    assert cat.genome_build == "GRCh38"


def test_ref_motif_matches_period():
    cat = load_catalog(CATALOG_PATH)
    for locus in cat:
        assert len(locus.ref_motif) == locus.period


def test_ar_ploidy_sex_dependent():
    cat = load_catalog(CATALOG_PATH)
    ar = cat.get("AR")
    assert ar.resolve_ploidy("M") == 1
    assert ar.resolve_ploidy("male") == 1
    assert ar.resolve_ploidy("F") == 2
    with pytest.raises(ValueError):
        ar.resolve_ploidy(None)
    with pytest.raises(ValueError):
        ar.resolve_ploidy("X")


def test_diploid_ploidy():
    cat = load_catalog(CATALOG_PATH)
    assert cat.get("DMPK").resolve_ploidy(None) == 2


def test_chrom_style_conversion():
    cat = load_catalog(CATALOG_PATH)
    nochr = cat.with_chrom_style("nochr")
    assert nochr.get("DMPK").chrom == "19"
    assert nochr.get("AR").chrom == "X"
    backchr = nochr.with_chrom_style("chr")
    assert backchr.get("DMPK").chrom == "chr19"


def test_validate_returns_warnings_list():
    cat = load_catalog(CATALOG_PATH)
    warns = validate_catalog(cat)
    assert isinstance(warns, list)


def test_duplicate_locus_rejected(tmp_path):
    import json

    bad = {
        "genome_build": "GRCh38",
        "loci": [
            _min_locus("DMPK"),
            _min_locus("DMPK"),
        ],
    }
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(bad))
    with pytest.raises(ValueError):
        load_catalog(p)


def _min_locus(lid):
    return {
        "locus_id": lid,
        "disease": "x",
        "chrom": "chr19",
        "ref_region": {"start": 100, "end": 160, "base": 1},
        "ref_motif": "CAG",
        "disease_motif": "CTG",
        "period": 3,
        "ploidy": "diploid",
        "eh_variant_id": lid,
        "hipstr_region": {"start": 99, "end": 160, "period": 3, "ref_copies": 20, "name": lid},
    }
