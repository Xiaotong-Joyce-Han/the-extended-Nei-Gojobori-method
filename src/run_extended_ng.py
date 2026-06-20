import argparse
import sys
from contextlib import redirect_stderr
from input_and_output import parse_args
from input_and_output import validate_args
from input_and_output import read_sequence_file
from input_and_output import read_spectrum_file
from input_and_output import append_to_outfile
from genetic_code import GENETIC_CODES
from count_sites import codon_N_S


def main() -> None:
    args = parse_args()

    if args.errorfile is None:
        run_main_logic(args)
    else:
        with open(args.errorfile, "w") as err:
            with redirect_stderr(err):
                run_main_logic(args)


def run_main_logic(args: argparse.Namespace) -> None:

    try:

        validate_args(args)
        
        sequence_file = args.sequence_file
        spectrum_file = args.spectrum_file
        outfile = args.outfile
        code_id = args.genetic_code

        append_to_outfile(outfile, GENETIC_CODES[code_id])
        append_to_outfile(outfile, "")

        seqs = read_sequence_file(sequence_file)
        append_to_outfile(outfile, "Before deleting ambiguous sites and stop codons:")
        append_to_outfile(outfile, seqs)
        append_to_outfile(outfile, f"Length of sequences: {seqs.seq_len()}")
        append_to_outfile(outfile, "")

        seqs.curate_sequences(code_id)
        append_to_outfile(outfile, "After deleting ambiguous sites and stop codons:")
        append_to_outfile(outfile, f"Length of sequences: {seqs.seq_len()}")
        append_to_outfile(outfile, "")

        seqs.count_pairwise_differences(code_id)
        append_to_outfile(outfile, "Further after deleting codons with no valid mutational pathway between any pair of sequences:")
        append_to_outfile(outfile, f"Length of sequences: {seqs.seq_len()}")
        append_to_outfile(outfile, "")

        append_to_outfile(outfile, "Pairwise n between sequences:")
        append_to_outfile(outfile, seqs.make_pairwise_n_s_matrix_text("n_total"))
        append_to_outfile(outfile, "")

        append_to_outfile(outfile, "Pairwise s between sequences:")
        append_to_outfile(outfile, seqs.make_pairwise_n_s_matrix_text("s_total"))
        append_to_outfile(outfile, "")        

        mut_spec = read_spectrum_file(spectrum_file)
        append_to_outfile(outfile, mut_spec)
        append_to_outfile(outfile, "")   

        codon_N_dict, codon_S_dict, codon_Stop_dict = codon_N_S(code_id, mut_spec)
        seqs.count_sequence_N_S(codon_N_dict, codon_S_dict, codon_Stop_dict)
        append_to_outfile(outfile, seqs.make_N_S_text())
        append_to_outfile(outfile, "")   

        seqs.calculate_pairwise_pN_pS()
        append_to_outfile(outfile, "Nonsynonymous p-distances (pN):")
        append_to_outfile(outfile, seqs.make_pairwise_pN_pS_text("pN"))
        append_to_outfile(outfile, "")         

        append_to_outfile(outfile, "Synonymous p-distances (pS):")
        append_to_outfile(outfile, seqs.make_pairwise_pN_pS_text("pS"))
        append_to_outfile(outfile, "")   

        append_to_outfile(outfile, "pN/pS:")
        append_to_outfile(outfile, seqs.make_pairwise_pN_pS_text("pNpS"))
        append_to_outfile(outfile, "") 

        seqs.Jukes_Cantor_correction()
        append_to_outfile(outfile, "Nonsynonymous distances after Jukes-Cantor correction (dN):")
        append_to_outfile(outfile, seqs.make_pairwise_dN_dS_text("dN"))
        append_to_outfile(outfile, "")

        append_to_outfile(outfile, "Synonymous distances after Jukes-Cantor correction (dS):")
        append_to_outfile(outfile, seqs.make_pairwise_dN_dS_text("dS"))
        append_to_outfile(outfile, "")

        append_to_outfile(outfile, "dN/dS:")
        append_to_outfile(outfile, seqs.make_pairwise_dN_dS_text("dNdS"))
        append_to_outfile(outfile, "")

    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()