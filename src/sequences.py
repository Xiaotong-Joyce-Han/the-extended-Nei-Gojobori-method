from genetic_code import GENETIC_CODES
import sys
from count_differences import count_differences
from count_sites import sequence_N_S
from Jukes_Cantor_correction import jukes_cantor_correction


class Sequence:
    """
    Attributes
    ----------
    name : str
    sequence : str
    length : int
    codons : list[str]

    original_indices : list[int]
        Maps each current codon position back to its original 0-based codon index
        before pruning ambiguous codons, stop codons, or invalid sites.

    N : float
    S : float

    N_list : list[float]
        Per-codon nonsynonymous sites.

    S_list : list[float]
        Per-codon synonymous sites.
    """

    def __init__(self, name: str, sequence: str):
        self.name = name
        self.sequence = sequence
        self.length = len(sequence)
        self.codons = list()

        for i in range(0, self.length, 3):
            self.codons.append(self.sequence[i: i+3])
        
        # store the 0-based indices 
        self.original_indices = list(range(len(self.codons)))

        self.N = 0
        self.S = 0
        self.N_list = list()
        self.S_list = list()


    def prune(self, indices: list[int]): 
        """
        delete codons at a list of indices from the current sequence
        """       
        for index in sorted(indices, reverse=True):
            self.codons.pop(index)
            self.original_indices.pop(index)

        self.sequence = "".join(self.codons)
        self.length = len(self.sequence)
    

    def contain_ambiguous_base(self) -> list[int]:
        '''
        identify the indices (in the current sequence) of codons with ambiguous base (or anything but "A", "T", "C", "G").
        '''

        indices: list[int] = []

        for i in range(len(self.codons)):
            codon = self.codons[i]

            if any(base not in {"A", "T", "C", "G"} for base in codon):
                indices.append(i)
        
        if indices:
            print(f"Warning: Detected codons containing ambiguous base in {self.name}! Their 1-based codonindices:"
                  f"{[i + 1 for i in indices]}"
                  "\nThose codons are deleted for all sequences", file=sys.stderr)
        
        return indices


    def contain_stop_codon(self, code_id: int) -> list[int]:
        '''
        identify the indices (in the current sequence) of stop codons
        '''

        code = GENETIC_CODES[code_id]
        indices: list[int] = []

        for i in range(len(self.codons)):
            if code.is_stop(self.codons[i]):
                indices.append(i)

        if indices:
            print(f"Warning: Detected stop codons in {self.name}! Their 1-based codon indices:"
                  f"{[i + 1 for i in indices]}"
                  "\nThose codons are deleted for all sequences", file=sys.stderr)
                
        return indices


class Sequences:
    '''
    A group of sequences to be analyzed in a run.
    
    Attributes
    ------
    sequences : list[Sequence]

    pairwise_differences : dict[tuple[str, str], dict[str, list[float] | float]]

        maps each sequence name pair to their pairwise n, s results.

        The inner dictionary contains:

            n_diffs: per-codon n
            s_diffs: per-codon s
            n_total
            s_total

    pairwise_p_distances_and_ratio: dict[tuple[str, str], dict[str, float]] 

        maps each sequence name pair to their pN, pS, pN/pS.

        The inner dictionary contains: "pN", "pS", "pNpS".

    pairwise_JC_distances_and_ratio: dict[tuple[str, str],bdict[str, float]] 

        maps each sequence name pair to their dN, dS, dN/dS.
        
        The inner dictionary contains: "dN", "dS", "dNdS"
    '''


    def __init__(self, sequences: list[Sequence]):

        self.sequences = sequences

        self.pairwise_differences: dict[
            tuple[str, str],
            dict[str, list[float] | float]
        ] = {}

        self.pairwise_p_distances_and_ratio: dict[
            tuple[str, str],
            dict[str, float]
        ] = {}

        self.pairwise_JC_distances_and_ratio: dict[
            tuple[str, str],
            dict[str, float]
        ] = {}


    def __len__(self) -> int:
        '''
        return the number of taxa
        '''
        return len(self.sequences)
    

    def seq_len(self) -> int:
        '''
        return the length of sequences
        '''
        return self.sequences[0].length
    

    def prune_all(self, indices: list[int]):
        '''
        delete codons at specified indices in all of the sequences. Pairwise difference attribute is also updated.
        '''
        
        indices = sorted(indices, reverse=True)

        for seq in self.sequences:
            seq.prune(indices)
        
        if self.pairwise_differences:
            for _, dict_for_pair in self.pairwise_differences.items():
                for _, diffs in dict_for_pair.items():
                    if isinstance(diffs, list):
                        for index in indices:
                            diffs.pop(index)

        if self.seq_len() == 0:
            raise ValueError(
                "All codon sites were removed after filtering ambiguous bases, "
                "stop codons, or sites without valid mutational pathways."
            ) 


    
    def __str__(self):
        '''
        print taxa names in this sequence group.
        '''
        sequence_names = [sequence.name for sequence in self.sequences]
        name_text = ", ".join(sequence_names)
        return (f"Number of taxa: {len(self.sequences)}. They are: {name_text}")

    
    def curate_sequences(self, code_id: int):
        '''
        For all the sequences in the group, identify the indices of stop codons and codons containing ambiguous base.
        Codons are deleted at these indices for all sequences.
        '''

        indices_to_prune: set[int] = set()

        for sequence in self.sequences:
            indices_to_prune.update(sequence.contain_ambiguous_base())
            indices_to_prune.update(sequence.contain_stop_codon(code_id))
        
        indices_to_prune = sorted(indices_to_prune, reverse=True)
        self.prune_all(indices_to_prune)


    def update_total_pairwise_s_n(self):
        '''
        Given the "n_diffs" and "s_diffs" lists in pairwise_differences, update the "n_total" and "s_total"
        '''
        if self.pairwise_differences:
            for _, dict_for_pair in self.pairwise_differences.items():
                dict_for_pair["n_total"] = sum(dict_for_pair["n_diffs"])
                dict_for_pair["s_total"] = sum(dict_for_pair["s_diffs"])


    def count_pairwise_differences(self, code_id: int):
        """
        Count pairwise nonsynonymous and synonymous differences.
        update the n_totals, s_totals, n_diffs, and s_diffs in the pairwise_differences attribute.
        """

        invalid_indices: set[int] = set()

        for i in range(len(self.sequences) - 1):
            for j in range(i + 1, len(self.sequences)):

                seq1 = self.sequences[i]
                seq2 = self.sequences[j]

                n_diffs, s_diffs, invalid_site_indices = count_differences(
                    seq1,
                    seq2,
                    code_id
                )

                invalid_indices.update(invalid_site_indices)

                pair_key = (seq1.name, seq2.name)

                self.pairwise_differences[pair_key] = {
                    "n_diffs": n_diffs,
                    "s_diffs": s_diffs,
                }

        if invalid_indices:
            self.prune_all(sorted(invalid_indices, reverse=True))

        self.update_total_pairwise_s_n()


    def make_pairwise_n_s_matrix_text(self, value_name: str) -> str:
        """
        Return pairwise n_total or s_total as a symmetric matrix in text format.
        The values are rounded to 4 decimal places.

        Parameters
        ----------
        value_name : str
            Either "n_total" or "s_total".
        """

        if value_name not in {"n_total", "s_total"}:
            raise ValueError(
                f"value_name must be 'n_total' or 's_total', but got {value_name!r}."
            )

        names = [seq.name for seq in self.sequences]

        lines: list[str] = []

        # header
        lines.append("\t" + "\t".join(names))

        for name1 in names:
            row = [name1]

            for name2 in names:
                if name1 == name2:
                    value = 0.0
                else:
                    pair_key = (name1, name2)
                    reverse_pair_key = (name2, name1)

                    if pair_key in self.pairwise_differences:
                        value = self.pairwise_differences[pair_key][value_name]
                    elif reverse_pair_key in self.pairwise_differences:
                        value = self.pairwise_differences[reverse_pair_key][value_name]
                    else:
                        value = float("nan")

                row.append(f"{value:.9f}")

            lines.append("\t".join(row))

        return "\n".join(lines)
    

    def count_sequence_N_S(self, codon_N_dict: dict[str, float], codon_S_dict: dict[str, float], codon_Stop_dict: dict[str, float]):
        '''
        count N, S for all the sequences in this group

        Parameters:
        ------
        codon_N_dict : dict[str, float]
            Maps each sense codon to its N under the mutation spectrum.

        codon_S_dict : dict[str, float]
            Maps each sense codon to its S under the mutation spectrum.

        codon_Stop_dict : dict[str, float]
            Maps each sense codon to its Stop under the mutation spectrum.
        '''
        for sequence in self.sequences:
            sequence_N_S(sequence, codon_N_dict, codon_S_dict, codon_Stop_dict)

    
    def make_N_S_text(self) -> str:
        """
        Return total N and S sites for each sequence as text.
        Four decimal places.
        """

        lines: list[str] = []

        header = ["sequence", "N", "S"]
        lines.append("\t".join(header))

        for sequence in self.sequences:
            name = sequence.name
            N = sequence.N
            S = sequence.S

            row = [
                name,
                f"{N:.9f}",
                f"{S:.9f}",
            ]

            lines.append("\t".join(row))

        return "\n".join(lines) 
    

    def calculate_pairwise_pN_pS(self) -> None:
        """
        Calculate pairwise pN and pS.

        pN = n_total / average_N
        pS = s_total / average_S
        pNpS = pN / pS

        Results are stored in self.pairwise_p_distances.
        """

        self.pairwise_p_distances_and_ratio = {}

        name_to_sequence = {
            sequence.name: sequence
            for sequence in self.sequences
        }

        for pair_key, diff_dict in self.pairwise_differences.items():
            
            name1, name2 = pair_key

            seq1 = name_to_sequence[name1]
            seq2 = name_to_sequence[name2]

            n_total = diff_dict["n_total"]
            s_total = diff_dict["s_total"]

            average_N = (seq1.N + seq2.N) / 2
            average_S = (seq1.S + seq2.S) / 2

            if average_N == 0:
                pN = float("nan")
            else:
                pN = n_total / average_N

            if average_S == 0:
                pS = float("nan")
            else:
                pS = s_total / average_S

            if pS == 0:
                pNpS = float("nan")
            else:
                pNpS = pN/pS


            self.pairwise_p_distances_and_ratio[pair_key] = {
                "pN": pN,
                "pS": pS,
                "pNpS": pNpS
            }


    def make_pairwise_pN_pS_text(self, value_name: str) -> str:
        """
        Return pairwise pN or pS or pN/pS as a symmetric matrix in text format.

        Parameters
        ----------
        value_name : str
            Either "pN" or "pS" or "pNpS".
        """

        if value_name not in {"pN", "pS", "pNpS"}:
            raise ValueError(
                f"value_name must be 'pN' or 'pS' pr 'pNpS', but got {value_name!r}."
            )

        names = [seq.name for seq in self.sequences]

        lines: list[str] = []

        # header
        lines.append("\t" + "\t".join(names))

        for name1 in names:
            row = [name1]

            for name2 in names:
                if name1 == name2:
                    value = 0.0
                else:
                    pair_key = (name1, name2)
                    reverse_pair_key = (name2, name1)

                    if pair_key in self.pairwise_p_distances_and_ratio:
                        value = self.pairwise_p_distances_and_ratio[pair_key][value_name]
                    elif reverse_pair_key in self.pairwise_p_distances_and_ratio:
                        value = self.pairwise_p_distances_and_ratio[reverse_pair_key][value_name]
                    else:
                        value = float("nan")

                row.append(f"{value:.9f}")

            lines.append("\t".join(row))

        return "\n".join(lines)

    
    def Jukes_Cantor_correction(self):
        '''
        Conduct Jukes Cantor correction for pN and pS.
        Results are stored in self.pairwise_JC_distances_and_ratio.
        '''

        self.pairwise_JC_distances_and_ratio = {}

        for pair_key, dist_dict in self.pairwise_p_distances_and_ratio.items():

            self.pairwise_JC_distances_and_ratio[pair_key] = {}

            dN = jukes_cantor_correction(dist_dict["pN"], pair_key[0], pair_key[1], "pN")
            self.pairwise_JC_distances_and_ratio[pair_key]["dN"] = dN
            
            dS = jukes_cantor_correction(dist_dict["pS"], pair_key[0], pair_key[1], "pS")
            self.pairwise_JC_distances_and_ratio[pair_key]["dS"] = dS

            if dS == 0:
                dNdS = float("nan")
            else:
                dNdS = dN/dS

            self.pairwise_JC_distances_and_ratio[pair_key]["dNdS"] = dNdS


    def make_pairwise_dN_dS_text(self, value_name: str) -> str:
        """
        Return pairwise dN or dS or dN/dS as a symmetric matrix in text format.

        Parameters
        ----------
        value_name : str
            Either "dN" or "dS" or "dNdS".
        """

        if value_name not in {"dN", "dS", "dNdS"}:
            raise ValueError(
                f"value_name must be 'dN' or 'dS' or 'dNdS', but got {value_name!r}."
            )

        names = [seq.name for seq in self.sequences]

        lines: list[str] = []

        # header
        lines.append("\t" + "\t".join(names))

        for name1 in names:
            row = [name1]

            for name2 in names:
                if name1 == name2:
                    value = 0.0
                else:
                    pair_key = (name1, name2)
                    reverse_pair_key = (name2, name1)

                    if pair_key in self.pairwise_JC_distances_and_ratio:
                        value = self.pairwise_JC_distances_and_ratio[pair_key][value_name]
                    elif reverse_pair_key in self.pairwise_JC_distances_and_ratio:
                        value = self.pairwise_JC_distances_and_ratio[reverse_pair_key][value_name]
                    else:
                        value = float("nan")

                row.append(f"{value:.9f}")

            lines.append("\t".join(row))

        return "\n".join(lines) 