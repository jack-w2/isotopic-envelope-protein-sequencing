import csv
from sequence import Seq
from molecular_formula import MolecularFormula


def get_aa_one_letter_codes(aa_file):
    """Get amino acids one letter codes from given file."""
    with open(aa_file, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        return [row[1] for row in reader]


aa_one_letter_codes = get_aa_one_letter_codes('amino_acids.csv')
print(aa_one_letter_codes)

for aa in aa_one_letter_codes:
    seq = Seq(aa, 'full')
    formula = seq.convert_to_molecular_formula()
    for aa2 in aa_one_letter_codes:
        for aa3 in aa_one_letter_codes:
            seq_to_test = Seq(aa2, 'full') + Seq(aa3, 'full')
            formula_to_test = seq_to_test.convert_to_molecular_formula() - MolecularFormula({'H': 2, 'O': 1})
            if formula == formula_to_test:
                print('hurra', aa2, '+', aa3, '=', aa)
            for aa4 in aa_one_letter_codes:
                seq_to_test = Seq(aa2, 'full') + Seq(aa3, 'full') + Seq(aa4, 'full')
                formula_to_test = seq_to_test.convert_to_molecular_formula() - MolecularFormula({'H': 2, 'O': 1}) * 2
                if formula == formula_to_test:
                    print('hurra', aa2, '+', aa3, '+', aa4, '=', aa)
                for aa5 in aa_one_letter_codes:
                    seq_to_test = Seq(aa2, 'full') + Seq(aa3, 'full') + Seq(aa4, 'full') + Seq(aa5, 'full')
                    formula_to_test = seq_to_test.convert_to_molecular_formula() - MolecularFormula({'H': 2, 'O': 1}) * 3
                    if formula == formula_to_test:
                        print('hurra', aa2, '+', aa3, '+', aa4, '+', aa5, '=', aa)
                    for aa6 in aa_one_letter_codes:
                        seq_to_test = Seq(aa2, 'full') + Seq(aa3, 'full') + Seq(aa4, 'full') + Seq(aa5, 'full') + Seq(aa6, 'full')
                        formula_to_test = seq_to_test.convert_to_molecular_formula() - MolecularFormula({'H': 2, 'O': 1}) * 4
                        if formula == formula_to_test:
                            print('hurra', aa2, '+', aa3, '+', aa4, '+', aa5, '+', aa6, '=', aa)
                        for aa7 in aa_one_letter_codes:
                            seq_to_test = Seq(aa2, 'full') + Seq(aa3, 'full') + Seq(aa4, 'full') + Seq(aa5, 'full') + Seq(aa6, 'full') + Seq(aa7, 'full')
                            formula_to_test = seq_to_test.convert_to_molecular_formula() - MolecularFormula({'H': 2, 'O': 1}) * 5
                            if formula == formula_to_test:
                                print('hurra', aa2, '+', aa3, '+', aa4, '+', aa5, '+', aa6, '+', aa7, '=', aa)

