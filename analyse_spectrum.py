import csv
from icecream import ic
from masserstein import estimate_proportions, Spectrum
from create_spectrum import create_spectrum


def create_aa_spectrum_objects(aa_file):
    with open(aa_file, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        aa_dict = {row[0]: row[1] for row in reader}
    spectras = [Spectrum(formula, adduct='H', label=label) for label, formula in aa_dict.items()]
    for spectrum in spectras:
        spectrum.normalize()
    return spectras


def analyse_spectrum(exp_spectrum, aa_file):
    exp_spectrum.normalize()
    return estimate_proportions(exp_spectrum, create_aa_spectrum_objects(aa_file))


# tests
exp_spectrum = create_spectrum()
estimations = analyse_spectrum(exp_spectrum, 'amino_acids.csv')
ic(estimations['proportions'])
ic(estimations['noise'][:11])
