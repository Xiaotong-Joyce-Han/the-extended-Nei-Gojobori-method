from pathlib import Path
import re


class GeneticCode:
    def __init__(self, name, code_id, codon_to_aa):
        self.name = name
        self.id = code_id
        self.codon_to_aa = codon_to_aa

    def translate_codon(self, codon) -> str:
        codon = codon.upper()
        return self.codon_to_aa[codon]

    def is_stop(self, codon):
        '''
        detects whether a codon is stop codon under the current genetic.
        This function won't detect invalid codons (containing anything but "A", "T", "C", "G") and should be used together with ambiguous base detecting
        '''
        if codon not in self.codon_to_aa:
            return False
        return self.translate_codon(codon) == "*"
    
    def __str__(self):
        codon_mapping = "\n".join(
            f"  {codon} -> {aa}"
            for codon, aa in self.codon_to_aa.items()
        )

        return (
            f"GeneticCode {self.id}: {self.name}\n"
            f"Codon to AA mapping:\n"
            f"{codon_mapping}"
        )
    

def make_codon_order(bases: list[str]) -> list[str]:

    """
    Parameters
    ----------
    bases : list[str]
        List of bases, e.g. ["T", "C", "A", "G"]

    Returns
    -------
    codons : list[str]
        List of codons in the order specified by the NCBI genetic code
    """

    codons = [
    first + second + third
    for first in bases
    for second in bases
    for third in bases
    ]

    return codons


def parse_one_code_block(block):

    names = re.findall(r'name\s+"(.*?)"', block, flags=re.DOTALL)
    code_id = re.search(r"id\s+(\d+)", block).group(1)
    aa_string = re.search(r'ncbieaa\s+"([A-Z*]+)"', block).group(1)

    name = " ".join(names[0].split())
    code_id = int(code_id)

    codon_order = make_codon_order(["T", "C", "A", "G"])
    codon_to_aa = dict(zip(codon_order, aa_string))

    return GeneticCode(name, code_id, codon_to_aa)


def load_genetic_codes(path: str) -> dict[int, GeneticCode]:

    text = Path(path).read_text()

    blocks = re.findall(
        r"\{\s*(name.*?ncbieaa.*?--\s*Base3\s+[TCAG]+.*?)\s*\}",
        text,
        flags=re.DOTALL,
    )

    genetic_codes = {}

    for block in blocks:
        code = parse_one_code_block(block)
        genetic_codes[code.id] = code

    return genetic_codes


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
FILE_PATH = DATA_DIR / "genetic_code.prt"

GENETIC_CODES = load_genetic_codes(FILE_PATH)
