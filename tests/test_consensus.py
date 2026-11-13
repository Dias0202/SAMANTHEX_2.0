"""Tests for the consensus merge, focused on low-alignment warnings."""

from samanthex.core.consensus import (
    is_low_alignment_score,
    merge_two_sequences,
)

READ = "ATGCGTACGTTAGCCGATCGATCGTAGCTAGCTAGGCTTACGATCGATCGGATCCATGCAATGCCGTA"


def test_identical_reads_merge_without_warning():
    warnings: list = []
    seq, qual = merge_two_sequences(READ, [30] * len(READ), READ, [30] * len(READ), warnings=warnings)
    assert seq == READ
    assert all(q == 60 for q in qual)
    assert warnings == []


def test_unrelated_reads_produce_visible_warning_but_still_merge():
    other = "TTTTTTTTTTGGGGGGGGGGAAAAAAAAAACCCCCCCCCCTTTTTTTTTTGGGGGGGGGGAAAAAAAAAA"
    warnings: list = []
    seq, _ = merge_two_sequences(
        READ, [30] * len(READ), other, [30] * len(other), warnings=warnings, label="Test merge"
    )
    assert seq  # the existing flow is not blocked
    assert len(warnings) == 1
    assert "LOW ALIGNMENT SCORE" in warnings[0]
    assert warnings[0].startswith("Test merge")


def test_warnings_argument_is_optional():
    seq, _ = merge_two_sequences(READ, [30] * len(READ), READ, [30] * len(READ))
    assert seq == READ


def test_low_score_threshold():
    assert is_low_alignment_score(score=19, min_len=50, match_score=2)
    assert not is_low_alignment_score(score=20, min_len=50, match_score=2)
