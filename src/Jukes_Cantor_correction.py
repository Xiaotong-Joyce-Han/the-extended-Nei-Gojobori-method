import math
import sys


def jukes_cantor_correction(p: float, sequence1_name: str, sequence2_name: str, distance_type: str) -> float:
    """
    Apply the Jukes-Cantor correction to an observed proportion of differences.

    Parameters
    ----------
    p : float
        Observed proportion of differences.
        For example, pS = synonymous_differences / synonymous_sites,
        or pN = nonsynonymous_differences / nonsynonymous_sites.
    sequence1_name : str
        Taxon of the first sequence (used for error messages).
    sequence2_name : str
        Taxon of the second sequence (used for error messages).
    distance_type : str
        Type of distance being corrected ("pN" or "pS"), used for error messages.

    Returns
    -------
    float
        Jukes-Cantor corrected evolutionary distance. Returns nan and writes a warning to stderr when p >= 0.75.
    """
    if p < 0:
        raise ValueError(
            f"Error in doing Jukes-Cantor correction for {distance_type} distance between sequences {sequence1_name} and {sequence2_name}.\n"
            f"p distance is expected to >= 0, but got {p}.\n"            
        )
    
    if p >= 0.75:
        print(
            f"Error in doing Jukes-Cantor correction for {distance_type} distance between sequences {sequence1_name} and {sequence2_name}.\n"
            f"Jukes-Cantor correction is undefined for p >= 0.75, but got {p}.\n",
            file=sys.stderr
        )
        return float("nan")

    return -0.75 * math.log(1 - 4 * p / 3)