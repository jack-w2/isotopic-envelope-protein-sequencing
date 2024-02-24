from masserstein import Spectrum
from collections import Counter


class MolecularFormula:
    def __init__(self, elements_dict):
        cnt = Counter(elements_dict)
        self.c = cnt['C']
        self.h = cnt['H']
        self.n = cnt['N']
        self.o = cnt['O']

    def __add__(self, other):
        # consider using Counter() + Counter()
        updated_elements_dict = {
            'C': self.c + other.c,
            'H': self.h + other.h,
            'N': self.n + other.n,
            'O': self.o + other.o
        }
        return MolecularFormula(updated_elements_dict)

    def __str__(self):
        s = ''
        for element, amount in self.__dict__.items():
            if amount > 0:
                s += f'{element.upper()}{amount}'
        return s

    def generate_spectrum(self):
        Spectrum(str(self)).plot()


# tests
formula1 = MolecularFormula({'C': 2, 'H': 6, 'O': 1})
print(formula1)
