from somatiq.interruptions import detect_interruptions, summarize_alleles
from somatiq.qc import QCConfig, level1_qc


def test_pure_tract_no_interruption():
    call = detect_interruptions("CAG" * 20, "CAG")
    assert call.flag is False
    assert call.n_impure_units == 0
    assert call.pure_fraction == 1.0


def test_rotation_treated_as_pure():
    # AGC/GCA are rotations of CAG -> still pure
    call = detect_interruptions("AGC" * 10, "CAG")
    assert call.flag is False


def test_interruption_detected():
    seq = "CAG" * 10 + "CCG" + "CAG" * 9
    call = detect_interruptions(seq, "CAG")
    assert call.flag is True
    assert "CCG" in call.impure_units
    assert call.n_impure_units == 1
    assert 0 < call.pure_fraction < 1


def test_summarize_alleles_any_impure():
    s = summarize_alleles(["CAG" * 20, "CAG" * 5 + "CCG" + "CAG" * 5], "CAG")
    assert s.flag is True
    assert "CCG" in s.impure_units


def test_qc_pass():
    r = level1_qc(
        informative_reads=47, stutter_learned=True, germline_q=0.99,
        allele_exceeds_readlen=False, interruption_flag=False,
    )
    assert r.flag == "PASS"
    assert abs(r.detectable_f_floor - 1 / 47) < 1e-9


def test_qc_fail_low_depth():
    r = level1_qc(
        informative_reads=5, stutter_learned=True, germline_q=0.99,
        allele_exceeds_readlen=False, interruption_flag=False,
        config=QCConfig(min_informative_reads=20),
    )
    assert r.flag == "FAIL"
    assert any("low_depth" in x for x in r.reasons)


def test_qc_warn_default_stutter_and_ceiling():
    r = level1_qc(
        informative_reads=40, stutter_learned=False, germline_q=0.99,
        allele_exceeds_readlen=True, interruption_flag=True,
    )
    assert r.flag == "WARN"
    assert "stutter_default_model" in r.reasons
    assert "allele_exceeds_readlen" in r.reasons
    assert "interruption_present" in r.reasons  # annotation only


def test_qc_interruption_alone_is_not_a_downgrade():
    r = level1_qc(
        informative_reads=40, stutter_learned=True, germline_q=0.99,
        allele_exceeds_readlen=False, interruption_flag=True,
    )
    assert r.flag == "PASS"
    assert "interruption_present" in r.reasons
