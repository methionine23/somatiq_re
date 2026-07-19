import pytest
from conftest import CATALOG_PATH

from somatiq.cli import main


def test_cli_catalog_check(capsys):
    rc = main(["catalog-check", "--catalog", str(CATALOG_PATH)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "4 loci" in out
    assert "DMPK" in out


def test_cli_regions(tmp_path):
    out = tmp_path / "r.bed"
    rc = main(["regions", "--catalog", str(CATALOG_PATH), "--out", str(out)])
    assert rc == 0
    assert len(out.read_text().splitlines()) == 4


def test_cli_regions_nochr(tmp_path):
    out = tmp_path / "r.bed"
    main(["regions", "--catalog", str(CATALOG_PATH), "--out", str(out), "--chr-style", "nochr"])
    assert out.read_text().startswith("18\t") or "\n19\t" in out.read_text()


def test_cli_requires_subcommand():
    with pytest.raises(SystemExit):
        main([])
