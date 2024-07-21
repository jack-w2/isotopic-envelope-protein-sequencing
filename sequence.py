import csv
from molecular_formula import MolecularFormula


class Seq:
    def __init__(self, seq, type):
        self.seq = seq
        self.type = type

    def __str__(self):
        return f'Seq({self.seq}, {self.type})'

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        if isinstance(other, str):
            return Seq(self.seq+other, self.type)
        assert self.type == other.type
        return Seq(self.seq + other.seq, type=self.type)

    def convert_to_molecular_formula(self):
        def convert_formula_str_to_obj(formula_str):
            elements_dict = {}
            for i, char in enumerate(formula_str):
                if char.isalpha():
                    chemical_element = char
                    if i + 1 >= len(formula_str) or formula_str[i + 1].isalpha():
                        amount = 1
                    else:
                        if i + 2 < len(formula_str) and formula_str[i + 2].isdigit():
                            amount = int(formula_str[i + 1:i + 3])
                        else:
                            amount = int(formula_str[i + 1])
                    elements_dict[chemical_element] = amount
            return MolecularFormula(elements_dict)

        with open('amino_acids.csv', 'r') as f:
            reader = csv.reader(f, delimiter=';')
            formulas_dict = {row[1]: convert_formula_str_to_obj(row[2]) for row in reader}
        result_formula = MolecularFormula({})
        for char in self.seq:
            result_formula += formulas_dict[char]
        return result_formula
