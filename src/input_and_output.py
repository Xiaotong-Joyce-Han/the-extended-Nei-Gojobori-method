from pathlib import Path
import argparse
from genetic_code import GENETIC_CODES
from sequences import Sequence, Sequences
from mutation_spectrum import MutationSpectrum
import sys
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Nei-Gojobori dN/dS estimation."
    )

    parser.add_argument("--outfile", required=True, help="Path to write output.")
    parser.add_argument("--sequence-file", required=True, help="Path to input coding sequences.")

    parser.add_argument(
        "--spectrum-file",
        default=None,
        help="Path to input mutation spectrum."
    )

    parser.add_argument(
        "--genetic-code",
        type=int,
        default=1,
        help="NCBI genetic code ID. Default: 1."
    )

    parser.add_argument(
        "--errorfile",
        default=None,
        help="Optional path to write warning/error messages. "
             "If not provided, messages are written to stderr."
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    '''
    to check the validate of user-defined arguments.
    '''

    sequence_file = Path(args.sequence_file)
    outfile = Path(args.outfile)
    code_id = args.genetic_code
    spectrum_file = args.spectrum_file

    if not sequence_file.is_file():
        raise ValueError(f"Sequence file does not exist: {sequence_file}")

    if spectrum_file:
        spectrum_file = Path(spectrum_file)
        if not spectrum_file.is_file():
            raise ValueError(f"Spectrum file does not exist: {spectrum_file}")

    if outfile.exists() and outfile.is_dir():
        raise ValueError(f"Output path is a directory, not a file: {outfile}")

    if outfile.parent and not outfile.parent.exists():
        raise ValueError(f"Output directory does not exist: {outfile.parent}")
    
    if outfile.exists() and outfile.stat().st_size > 0:
        print(
            f"Warning: output file '{outfile}' already exists and is not empty. "
            f"New results will be appended to it.",
            file=sys.stderr,
        )
    
    if code_id not in GENETIC_CODES:
        raise ValueError(f"NCBI Genetic Code {code_id} does not exist")
        

def read_sequence_file(sequence_file: str | Path) -> Sequences:
    path = Path(sequence_file)

    with open(path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) < 1:
        raise ValueError(f"Sequence file is empty: {sequence_file}")

    header = lines[0].split()

    if len(header) != 2:
        raise ValueError(
            "The first line of sequence file must contain exactly two numbers: "
            "number of sequences and sequence length."
        )

    try:
        expected_num_sequences = int(header[0])
        expected_sequence_length = int(header[1])
    except ValueError:
        raise ValueError(
            "The first line of sequence file must contain two integers: "
            "number of sequences and sequence length."
        )

    if expected_sequence_length <= 0:
        raise ValueError(
            f"Coding sequence length must be positive, got a length of{expected_sequence_length}."
        )
    
    if expected_sequence_length % 3 != 0:
        raise ValueError(
            f"Coding sequence length must be a multiple of 3, "
            f"but got {expected_sequence_length}."
        )

    sequence_lines = lines[1:]

    if len(sequence_lines) % 2 != 0:
        raise ValueError(
            "Sequence file format error: "
            "each sequence must have one name line and one sequence line."
        )

    actual_num_sequences = len(sequence_lines) // 2

    if actual_num_sequences != expected_num_sequences:
        raise ValueError(
            f"Number of sequences does not match the header: "
            f"expected {expected_num_sequences}, got {actual_num_sequences}."
        )
    
    if actual_num_sequences <= 1:
        raise ValueError(
            f"Number of sequences is not valid: "
            f"expected at least 2, got {actual_num_sequences}"
        )

    sequences: list[Sequence] = []
    seen_names: set[str] = set()

    for i in range(0, len(sequence_lines), 2):
        name = sequence_lines[i]
        sequence = sequence_lines[i + 1].upper()

        if len(sequence) != expected_sequence_length:
            raise ValueError(
                f"Sequence length mismatch for {name}: "
                f"expected {expected_sequence_length}, got {len(sequence)}."
            )
        
        if name in seen_names:
            raise ValueError(f"Duplicated sequence name: {name}")
        
        seen_names.add(name) 

        sequences.append(Sequence(name, sequence))

    return Sequences(sequences)


SIX_RATE_KEYS = [
    "ATtoTA",
    "ATtoCG",
    "ATtoGC",
    "CGtoTA",
    "CGtoAT",
    "CGtoGC",
]


DEFAULT_SPECTRUM_RATES = [1.0 / 3.0] * 6


def read_spectrum_file(spectrum_file: str | Path | None) -> MutationSpectrum:
    """
    Read mutation spectrum from a file and return a MutationSpectrum object.

    Accepted formats
    ----------------
    Format 1: ordered 6-parameter format in arbitrary unit

        0.33 (per AT site per generation * some constant)
        0.33 (per AT site per generation * some constant)
        0.34 (per AT site per generation * some constant)
        0.20 (per CG site per generation * some constant)
        0.50 (per CG site per generation * some constant)
        0.30 (per CG site per generation * some constant)

    Only numbers should be present in the mutation spectrum file. The units here are only for illustration.

    The order is:

        ATtoTA
        ATtoCG
        ATtoGC
        CGtoTA
        CGtoAT
        CGtoGC

    Format 2: labeled 6-parameter format in arbitrary unit

        ATtoTA 0.33
        CGtoGC 0.30
        ATtoCG 0.33
        CGtoAT 0.50
        ATtoGC 0.34
        CGtoTA 0.20

    The order does not matter.

    Format 3: ordered 8-parameter format

        1
        1
        1
        1
        1
        1
        0.8
        1.2

    The order is:

        ATtoTA relative rate among AT-origin mutations 
        ATtoCG relative rate among AT-origin mutations
        ATtoGC relative rate among AT-origin mutations
            (The above three parameters should be normalized to the same scale)
        CGtoTA relative rate among CG-origin mutations
        CGtoAT relative rate among CG-origin mutations
        CGtoGC relative rate among CG-origin mutations
            (The above three parameters should be normalized to the same scale)
        total relative AT-origin mutation rate
        total relative CG-origin mutation rate
            (The above two parameters should be normalized to the same scale)

    Format 4: labeled 8-parameter format

        ATtoTA 1
        ATtoCG 1
        ATtoGC 1
        CGtoTA 0.4
        CGtoAT 0.5
        CGtoGC 0.3
        ATrate 0.8
        CGrate 1.2
    
    The order doesn't matter.

    Notes
    -----
    If spectrum_file is None, the default spectrum is used:

        ATtoTA = ATtoCG = ATtoGC = CGtoTA = CGtoAT = CGtoGC = 1/3

    Returns
    -------
    MutationSpectrum
        A MutationSpectrum object constructed from 6 final mutation-type rates.
    """

    if spectrum_file is None:
        return MutationSpectrum(DEFAULT_SPECTRUM_RATES)

    path = Path(spectrum_file)

    lines: list[str] = []

    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.split("#", maxsplit=1)[0].strip()

            if line:
                lines.append(line)

    if not lines:
        raise ValueError(f"Spectrum file '{path}' is empty.")

    # Decide whether this is an ordered file or a labeled file.
    first_tokens = lines[0].replace("=", " ").replace(":", " ").split()

    if len(first_tokens) == 1:
        is_labeled = False
    elif len(first_tokens) == 2:
        is_labeled = True
    else:
        raise ValueError(
            f"Spectrum file '{path}' has an invalid line: {lines[0]!r}."
        )

    # Case 1 or 3: ordered format.
    if not is_labeled:

        values: list[float] = []

        for line in lines:
            tokens = line.replace("=", " ").replace(":", " ").split()

            if len(tokens) != 1:
                raise ValueError(
                    f"Spectrum file '{path}' mixes ordered and labeled lines: {line!r}."
                )

            try:
                value = float(tokens[0])
            except ValueError:
                raise ValueError(
                    f"Spectrum file '{path}' contains a non-numeric value: {line!r}."
                )
            
            values.append(value)


        if len(values) == 6:
            rates = values

        elif len(values) == 8:
            at_relative_rates = values[0:3]
            cg_relative_rates = values[3:6]
            at_total_rate = values[6]
            cg_total_rate = values[7]

            at_sum = sum(at_relative_rates)
            cg_sum = sum(cg_relative_rates)            

            if at_sum <= 0:
                raise ValueError(
                    f"Spectrum file '{path}' has invalid AT relative rates. "
                    f"Their sum must be positive."
                )

            if cg_sum <= 0:
                raise ValueError(
                    f"Spectrum file '{path}' has invalid CG relative rates. "
                    f"Their sum must be positive."
                )

            rates = [
                at_total_rate * at_relative_rates[0] / at_sum,
                at_total_rate * at_relative_rates[1] / at_sum,
                at_total_rate * at_relative_rates[2] / at_sum,
                cg_total_rate * cg_relative_rates[0] / cg_sum,
                cg_total_rate * cg_relative_rates[1] / cg_sum,
                cg_total_rate * cg_relative_rates[2] / cg_sum,
            ]

        else:
            raise ValueError(
                f"Spectrum file '{path}' contains {len(values)} values. "
                f"Expected either 6 or 8 values."
            )

    # Case 2 or 4: labeled format.
    else:
        label_to_value: dict[str, float] = {}

        for line in lines:
            tokens = line.replace("=", " ").replace(":", " ").split()

            if len(tokens) != 2:
                raise ValueError(
                    f"Spectrum file '{path}' has an invalid labeled line: {line!r}."
                )

            label = tokens[0]
            value_text = tokens[1]

            try:
                value = float(value_text)
            except ValueError:
                raise ValueError(
                    f"Spectrum file '{path}' contains a non-numeric value: {line!r}."
                )

            if label in label_to_value:
                raise ValueError(
                    f"Spectrum file '{path}' contains duplicated label: {label}."
                )

            label_to_value[label] = value

        six_labels = set(SIX_RATE_KEYS)
        eight_labels = set(SIX_RATE_KEYS + ["ATrate", "CGrate"])
        given_labels = set(label_to_value)

        if given_labels == six_labels:
            rates = [label_to_value[name] for name in SIX_RATE_KEYS]

        elif given_labels == eight_labels:
            at_relative_rates = [
                label_to_value["ATtoTA"],
                label_to_value["ATtoCG"],
                label_to_value["ATtoGC"],
            ]

            cg_relative_rates = [
                label_to_value["CGtoTA"],
                label_to_value["CGtoAT"],
                label_to_value["CGtoGC"],
            ]

            at_total_rate = label_to_value["ATrate"]
            cg_total_rate = label_to_value["CGrate"]

            at_sum = sum(at_relative_rates)
            cg_sum = sum(cg_relative_rates)

            if at_sum <= 0:
                raise ValueError(
                    f"Spectrum file '{path}' has invalid AT relative rates. "
                    f"Their sum must be positive."
                )

            if cg_sum <= 0:
                raise ValueError(
                    f"Spectrum file '{path}' has invalid CG relative rates. "
                    f"Their sum must be positive."
                )

            rates = [
                at_total_rate * at_relative_rates[0] / at_sum,
                at_total_rate * at_relative_rates[1] / at_sum,
                at_total_rate * at_relative_rates[2] / at_sum,
                cg_total_rate * cg_relative_rates[0] / cg_sum,
                cg_total_rate * cg_relative_rates[1] / cg_sum,
                cg_total_rate * cg_relative_rates[2] / cg_sum,
            ]

        else:
            raise ValueError(
                f"Spectrum file '{path}' has invalid labels.\n"
                f"For 6-parameter format, expected labels: {sorted(six_labels)}\n"
                f"For 8-parameter format, expected labels: {sorted(eight_labels)}\n"
                f"Observed labels: {sorted(given_labels)}"
            )

    return MutationSpectrum(rates)



def append_to_outfile(outfile: str | Path, content: Any) -> None:
    """
    Append text to outfile. Create outfile if it does not exist.
    """
    outfile = Path(outfile)

    with outfile.open("a", encoding="utf-8") as f:
        f.write(str(content) + "\n")