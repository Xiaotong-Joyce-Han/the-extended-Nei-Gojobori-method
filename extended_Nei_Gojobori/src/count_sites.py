from genetic_code import GENETIC_CODES
from genetic_code import make_codon_order
from mutation_spectrum import MutationSpectrum
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from sequences import Sequence


def codon_potential_mutations(code_id: int) -> tuple[dict[str, dict[tuple[str, str], int]], dict[str, dict[tuple[str, str], int]]]:
    '''
    Returns
    -------
    nonsynonymous_potential_mutations_dict : dict[str, dict[tuple[str, str], int]]
        Maps each sense codon to a dictionary giving the number of possible
        nonsynonymous point mutations of each mutation type (original nucleotide and target nucleotide).

    synonymous_potential_mutations_dict : dict[str, dict[tuple[str, str], int]]
        Maps each sense codon to a dictionary giving the number of possible
        synonymous point mutations of each mutation type (original nucleotide and target nucleotide).

    For most codons, the sum of possible nonsynonymous point mutations and possible synonymous mutations will be 9.
    However, if the codon can yield stop codons through a point mutation, the sum will be 7 or 8
    '''

    code = GENETIC_CODES[code_id]
    if code is None:
        raise ValueError(f"NCBI Genetic Code {code_id} does not exist.")

    codons = make_codon_order(["T", "C", "A", "G"])
    bases = ["T", "C", "A", "G"]

    mutation_types = [
        (original_base, target_base)
        for original_base in bases
        for target_base in bases
        if original_base != target_base
    ]

    nonsynonymous_potential_mutations_dict = {}
    synonymous_potential_mutations_dict = {}

    for codon in codons:

        if code.is_stop(codon):
            continue

        nonsynonymous_counts = {
            mutation_type: 0 for mutation_type in mutation_types
        }
        synonymous_counts = {
            mutation_type: 0 for mutation_type in mutation_types
        }

        aa = code.translate_codon(codon)

        for position in range(3):
            original_base = codon[position]

            for target_base in bases:
                if target_base == original_base:
                    continue

                mutated_codon = (
                    codon[:position]
                    + target_base
                    + codon[position + 1:]
                )

                if code.is_stop(mutated_codon):
                    continue

                mutated_aa = code.translate_codon(mutated_codon)
                mutation_type = (original_base, target_base)

                if mutated_aa == aa:
                    synonymous_counts[mutation_type] += 1
                else:
                    nonsynonymous_counts[mutation_type] += 1

        nonsynonymous_potential_mutations_dict[codon] = nonsynonymous_counts
        synonymous_potential_mutations_dict[codon] = synonymous_counts

    return (
        nonsynonymous_potential_mutations_dict,
        synonymous_potential_mutations_dict,
    )


def codon_N_S(code_id: int, mutation_spectrum: MutationSpectrum) -> tuple[dict[str, float], dict[str, float]]:
    '''
    Returns
    -------
    codon_N_dict : dict[str, float]
        Maps each sense codon to its N under the mutation spectrum.

    codon_S_dict : dict[str, float]
        Maps each sense codon to its S under the mutation spectrum.
    '''

    nonsynonymous_potential_mutations_dict, synonymous_potential_mutations_dict = codon_potential_mutations(code_id)

    codon_N_dict = dict()
    codon_S_dict = dict()

    for codon, nonsynonymous_counts in nonsynonymous_potential_mutations_dict.items():
        
        N = 0

        for mutation_type, number in nonsynonymous_counts.items():

            if number == 0:
                continue

            original_base = mutation_type[0]
            target_base = mutation_type[1]
            N += mutation_spectrum.get_mutation_rate(original_base, target_base)*number
        
        codon_N_dict[codon] = N
    
    for codon, synonymous_counts in synonymous_potential_mutations_dict.items():

        S = 0

        for mutation_type, number in synonymous_counts.items():

            if number == 0:
                continue

            original_base = mutation_type[0]
            target_base = mutation_type[1]
            S += mutation_spectrum.get_mutation_rate(original_base, target_base)*number
        
        codon_S_dict[codon] = S
    
    return codon_N_dict, codon_S_dict


def sequence_N_S(sequence: "Sequence", codon_N_dict: dict[str, float], codon_S_dict: dict[str, float]) -> None:
    '''
    Parameters
    -------
    sequence : "Sequence"
        An sequence object. The object will update its N, S, N_list, S_list attributes in the function.
    
    codon_N_dict : dict[str, float]
        Maps each sense codon to its N under the mutation spectrum.

    codon_S_dict : dict[str, float]
        Maps each sense codon to its S under the mutation spectrum.
    '''

    sequence.N = 0
    sequence.S = 0
    sequence.N_list = list()
    sequence.S_list = list()

    for codon in sequence.codons:
        sequence.N += codon_N_dict[codon]
        sequence.S += codon_S_dict[codon]
        sequence.N_list.append(codon_N_dict[codon])
        sequence.S_list.append(codon_S_dict[codon])

    expected_sum = sequence.length
    actual_sum = sequence.N + sequence.S
    sequence.N *= expected_sum/actual_sum
    sequence.S *= expected_sum/actual_sum

    for value in sequence.N_list:
        value *= expected_sum/actual_sum
    for value in sequence.S_list:
        value *= expected_sum/actual_sum


