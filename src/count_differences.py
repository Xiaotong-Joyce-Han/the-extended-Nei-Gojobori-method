from genetic_code import GENETIC_CODES, make_codon_order
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


def make_pair_ns_table(code_id) -> dict[tuple[str,str], tuple[float, float, bool]]:
    """
    For each pair of codons, calculate and memorize their n, s, validity.

    """
    code = GENETIC_CODES.get(code_id)
    codons = make_codon_order(["T", "C", "A", "G"])
    codons = [codon for codon in codons if not code.is_stop(codon)]

    table = {}

    for c1 in codons:
        for c2 in codons:
            try:
                n, s = count_codon_n_s(c1, c2, code_id)
                table[(c1, c2)] = (n, s, True)
            except ValueError:
                table[(c1, c2)] = (0.0, 0.0, False)

    return table


def get_invalid_indices(sequence1: "Sequence", sequence2: "Sequence", table: dict[tuple[str,str], tuple[float, float, bool]], 
                        code_id: int, original_indices: list[int]) -> list[int]:
    """
    First scanning: detect codon positions without valid pathways

    Returns
    -------
    invalid_sites_indices : list[int]
        List of codon indices where there is no valid pathway between this sequence pair.
    """
    invalid_site_indices = []

    for codon_index, (codon1, codon2) in enumerate(zip(sequence1.codons, sequence2.codons)):
        _, _, valid = table[(codon1, codon2)]

        if not valid:

            invalid_site_indices.append(codon_index)
            original_index = original_indices[codon_index]

            print(f"Error processing sequence{sequence1.name} vs sequence{sequence2.name}, codon pair {codon1} vs {codon2}.\n"
                  f"No valid mutational path between {codon1} and {codon2} "
                  f"under genetic code {code_id}. All paths involve stop codons.\n"
                  f"They are at codon index {original_index+1} (1-based index).\n"
                  f"For all sequences, skipped this codon site.\n",
                  file = sys.stderr)
    
    return invalid_site_indices


def count_differences(sequence1: "Sequence", sequence2: "Sequence", 
                      table: dict[tuple[str,str], tuple[float, float, bool]]) -> tuple[float, float]:
    """
    Second scanning: count Nei-Gojobori nonsynonymous and synonymous differences
    between two sequences of codons after invalid sites have been removed.

    Returns
    -------
    n_total : float
        Nonsynonymous differences between this pair of sequence.

    s_total : float
        Synonymous differences between this pair of sequence.
    """

    n_total = 0
    s_total = 0

    for (codon1, codon2) in zip(sequence1.codons, sequence2.codons):
        n_diff, s_diff, valid = table[(codon1, codon2)]

        if not valid:
            raise ValueError("After removing invalid sites without pathways, there should not be such codon pairs.")
        
        n_total += n_diff
        s_total += s_diff

    return n_total, s_total

