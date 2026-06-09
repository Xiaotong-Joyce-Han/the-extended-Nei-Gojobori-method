from genetic_code import GENETIC_CODES
from itertools import permutations
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sequences import Sequence


def count_codon_n_s(codon1:str, codon2:str, code_id:int) -> tuple[float, float]:
    """
    Count Nei-Gojobori nonsynonymous (n) and synonymous (s) differences
    between two codons.

    All possible pathways between codon 1 and codon 2 are given equal probability,
    and the average number of differences is returned.

    Pathways that involve stop codons are ignored.

    Returns
    -------
    n_diff : float
        Number of nonsynonymous differences between codon1 and codon2.

    s_diff : float
        Number of synonymous differences between codon1 and codon2.
    """

    code = GENETIC_CODES.get(code_id)

    if code is None:
        raise ValueError(f"NCBI Genetic Code {code_id} does not exist.")

    if codon1 == codon2:
        return 0.0, 0.0


    different_positions = [
    i for i in range(3)
    if codon1[i] != codon2[i]
    ]

    total_n = 0
    total_s = 0
    valid_path_count = 0

    for path in permutations(different_positions):
        current_codon = codon1
        path_n = 0.0
        path_s = 0.0
        valid_path = True

        for position in path:
            next_codon = (
                current_codon[:position]
                + codon2[position]
                + current_codon[position + 1:]
            )

            if code.is_stop(next_codon):
                valid_path = False
                break

            current_aa = code.translate_codon(current_codon)
            next_aa = code.translate_codon(next_codon)

            if current_aa == next_aa:
                path_s += 1.0
            else:
                path_n += 1.0

            current_codon = next_codon

        if valid_path:
            total_n += path_n
            total_s += path_s
            valid_path_count += 1

    if valid_path_count == 0:
        raise ValueError(
            f"No valid mutational path between {codon1} and {codon2} "
            f"under genetic code {code_id}. All paths involve stop codons."
        )

    n_diff = total_n / valid_path_count
    s_diff = total_s / valid_path_count

    return n_diff, s_diff


def count_differences(sequence1: "Sequence", sequence2: "Sequence", code_id: int) -> tuple[list[float], list[float], list[int]]:
    """
    Count Nei-Gojobori nonsynonymous and synonymous differences
    between two sequences of codons.

    Returns
    -------
    n_diffs : list[float]
        List of nonsynonymous differences for each codon pair.

    s_diffs : list[float]
        List of synonymous differences for each codon pair.

    invalid_sites_indices : list[int]
        List of codon indices where there is no valid pathway between the sequence pair.
    """

    n_diffs = []
    s_diffs = []
    invalid_site_indices = []

    codon_index = 0

    for codon1, codon2 in zip(sequence1.codons, sequence2.codons):

        try:

            n_diff, s_diff = count_codon_n_s(codon1, codon2, code_id)
            n_diffs.append(n_diff)
            s_diffs.append(s_diff)

        except ValueError as e:

            invalid_site_indices.append(codon_index)

            original_index = sequence1.original_indices[codon_index]

            print(f"Error processing sequence{sequence1.name} vs sequence{sequence2.name}, codon pair {codon1} vs {codon2}.\n"
                  f"{e}\n"
                  f"They are at codon index {original_index+1} (1-based index).\n"
                  f"For all sequences, skipped this codon site and continued with the next one.\n",
                  file = sys.stderr)
            
            n_diffs.append(-1)
            s_diffs.append(-1)

        codon_index += 1

    return n_diffs, s_diffs, invalid_site_indices

