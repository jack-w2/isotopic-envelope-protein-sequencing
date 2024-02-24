from masserstein import Spectrum


class MolecularFormula:
    def __init__(self, elements_dict):
        self.c = elements_dict.get('C', 0)
        self.h = elements_dict.get('H', 0)
        self.n = elements_dict.get('N', 0)
        self.o = elements_dict.get('O', 0)
        self.s = elements_dict.get('S', 0)

    def __add__(self, other):
        updated_elements_dict = {
            'C': self.c + other.c,
            'H': self.h + other.h,
            'N': self.n + other.n,
            'O': self.o + other.o,
            'S': self.s + other.s
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
