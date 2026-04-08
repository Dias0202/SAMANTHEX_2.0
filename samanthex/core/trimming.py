"""Quality trimming for Sanger sequences using Mott's modified algorithm.

Implements the same trimming approach used by phred/phrap and other standard
Sanger sequencing tools. The algorithm finds the longest high-quality region
by maximizing the cumulative quality excess over a threshold.

Reference:
    Ewing B, Green P. Base-calling of automated sequencer traces using
    phred. II. Error probabilities. Genome Research. 1998;8(3):186-194.
"""

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class TrimResult:
    """Result of quality trimming.

    Attributes:
        sequence: Trimmed DNA sequence.
        quality: Trimmed Phred quality scores.
        trim_start: Start position of retained region in original sequence.
        trim_end: End position (exclusive) of retained region.
        original_length: Length of original sequence before trimming.
        bases_trimmed_left: Number of bases trimmed from the 5' end.
        bases_trimmed_right: Number of bases trimmed from the 3' end.
    """
    sequence: str
    quality: List[int]
    trim_start: int
    trim_end: int
    original_length: int

    @property
    def trimmed_length(self) -> int:
        return self.trim_end - self.trim_start

    @property
    def bases_trimmed_left(self) -> int:
        return self.trim_start

    @property
    def bases_trimmed_right(self) -> int:
        return self.original_length - self.trim_end


def mott_trim(
    sequence: str,
    quality: List[int],
    threshold: int = 20,
    min_length: int = 50,
) -> TrimResult:
    """Trim low-quality ends using Mott's modified algorithm.

    Finds the contiguous region that maximizes the sum of
    (quality_i - threshold) across all positions. This is equivalent
    to Kadane's maximum subarray algorithm, which efficiently identifies
    the highest-quality region of the read.

    This is the standard algorithm used by phred for quality trimming.

    Args:
        sequence: DNA sequence string.
        quality: List of Phred quality scores (same length as sequence).
        threshold: Minimum Phred quality score. Bases below this
            contribute negatively to the score.
        min_length: Minimum acceptable length after trimming. If the
            trimmed region is shorter, falls back to window-based trimming.

    Returns:
        TrimResult with the trimmed sequence, quality, and trim metadata.
    """
    n = len(sequence)

    if n == 0 or not quality:
        return TrimResult(
            sequence=sequence,
            quality=quality,
            trim_start=0,
            trim_end=0,
            original_length=0,
        )

    # Mott's algorithm: maximum subarray (Kadane's variant)
    # Score each position as quality_i - threshold
    best_start = 0
    best_end = n
    best_score = float("-inf")

    current_start = 0
    current_score = 0

    for i in range(n):
        current_score += quality[i] - threshold

        if current_score > best_score:
            best_score = current_score
            best_start = current_start
            best_end = i + 1

        if current_score < 0:
            current_score = 0
            current_start = i + 1

    # Validate the result
    if best_score <= 0:
        # No region has net positive quality — try window-based fallback
        best_start, best_end = _window_trim(quality, threshold)

    # Ensure minimum length
    if best_end - best_start < min_length and n >= min_length:
        # If trimmed region is too short, try window-based fallback
        ws, we = _window_trim(quality, threshold)
        if we - ws >= min_length:
            best_start, best_end = ws, we
        else:
            # Keep the entire sequence as last resort
            best_start, best_end = 0, n

    # Final bounds check
    best_start = max(0, best_start)
    best_end = min(n, best_end)
    if best_end <= best_start:
        best_start, best_end = 0, n

    return TrimResult(
        sequence=sequence[best_start:best_end],
        quality=quality[best_start:best_end],
        trim_start=best_start,
        trim_end=best_end,
        original_length=n,
    )


def _window_trim(
    quality: List[int],
    threshold: int,
    window: int = 10,
) -> Tuple[int, int]:
    """Fallback sliding-window trimming.

    Finds the first and last windows where the mean quality exceeds
    the threshold, then returns the region between them.

    Args:
        quality: Phred quality scores.
        threshold: Quality threshold.
        window: Sliding window size.

    Returns:
        Tuple of (start, end) positions.
    """
    n = len(quality)
    if n < window:
        return 0, n

    # Find start: first window where mean quality >= threshold
    start = 0
    for i in range(n - window + 1):
        window_mean = sum(quality[i : i + window]) / window
        if window_mean >= threshold:
            start = i
            break

    # Find end: last window where mean quality >= threshold
    end = n
    for i in range(n - 1, window - 2, -1):
        w_start = max(0, i - window + 1)
        window_mean = sum(quality[w_start : i + 1]) / (i - w_start + 1)
        if window_mean >= threshold:
            end = i + 1
            break

    if end <= start:
        return 0, n

    return start, end
