from masserstein import Spectrum


class MolecularFormula:
    def __init__(self, c, h, n, o):
        self.c = c
        self.h = h
        self.n = n
        self.o = o

    def __add__(self, other):
        return MolecularFormula(self.c + other.c, self.h + other.h, self.n + other.n, self.o + other.o)

    def generate_spectrum(self):
        Spectrum(f'C{self.c}H{self.h}N{self.n}O{self.o}').plot()
