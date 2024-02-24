from IsoSpecPy import IsoDistribution, IsoTotalProb
from masserstein import Spectrum
from analyse_spectrum import analyse_spectrum
from icecream import ic
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt
import csv
import itertools


def generate_prefixes_and_suffixes(seq):
    """Split given seq to create all possible suffixes and prefixes."""
    seqs = itertools.chain.from_iterable([(seq[:i], 'pref'), (seq[i:], 'suf')] for i in range(1, len(seq)))
    return list(seqs)


def scale_intensities(intensities):
    """Scale given intensities to make sure they sum to 1."""
    intensities_sum = sum(intensities)
    return [i/intensities_sum for i in intensities]


def load_config_file(config_file_path):
    """Load given toml config file. Returns a dictionary with parameter names as keys."""
    with open(config_file_path, 'rb') as tf:
        config = tomllib.load(tf)
    return config


def add_noise(masserstein_spectrum, nb_of_noise_peaks=100, noise_fraction=0.1, sd=0.01):
    plt.figure()
    plt.title('raw spectrum')
    masserstein_spectrum.plot()
    masserstein_spectrum.normalize()    # maybe not necessary
    plt.figure()
    plt.title('normalized spectrum')
    masserstein_spectrum.plot()
    masserstein_spectrum.add_chemical_noise(nb_of_noise_peaks, noise_fraction)
    plt.figure()
    plt.title('spectrum with chemical noise')
    masserstein_spectrum.plot()
    masserstein_spectrum.gaussian_smoothing()
    plt.figure()
    plt.title('spectrum with gaussian smoothing')
    masserstein_spectrum.plot()
    masserstein_spectrum.add_gaussian_noise(sd)
    plt.figure()
    plt.title('spectrum with gaussian noise')
    masserstein_spectrum.plot()


def scoring_function(seq, aa_index):
    """Scoring function idea to use analyse_spectrum"""
    spectrum_to_test = create_raw_spectrum_from_fasta(seq)
    proportions = analyse_spectrum(spectrum_to_test)['proportions']
    return proportions[aa_index]


def find_next_best_letter(seq, aa_file='amino_acids.csv'):
    """For given seq find next best amino acid."""
    def get_aa_one_letter_codes(aa_file):
        """Get amino acids one letter codes from given file."""
        with open(aa_file, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            return [row[1] for row in reader]
    aa_one_letter_codes = get_aa_one_letter_codes(aa_file)
    best_seq = ''
    best_score = 0
    for aa_index, aa in enumerate(aa_one_letter_codes):
        test_seq = seq + aa
        test_score = scoring_function(test_seq, aa_index)
        if test_score > best_score:
            best_score = test_score
            best_seq = test_seq
    return best_seq


def create_spectrum():
    config = load_config_file('config.toml')

    seq = config['fasta']
    seqs = generate_prefixes_and_suffixes(seq)

    envelopes = [IsoTotalProb(0.999, fasta=seqs[s][0]) for s in range(len(seqs))]
    intensities = []
    for i in config['probs']:
        intensities.extend([i * 0.5] * 2)

    intensities = scale_intensities(intensities)

    spectre = IsoDistribution.LinearCombination(envelopes, intensities)
    masses_and_intensities = list(zip(spectre.masses, spectre.probs))

    print(masses_and_intensities)

    # spectre.plot()
    masserstein_spectrum = Spectrum(confs=masses_and_intensities)
    add_noise(masserstein_spectrum, config['noise']['nb_of_noise_peaks'], config['noise']['noise_fraction'], config['noise']['gaussian_noise_sd'])
    return masserstein_spectrum


def create_raw_spectrum_from_fasta(seq):
    """Create raw spectrum for given fasta string (without noise and other stuff)."""
    spectre = IsoTotalProb(0.999, fasta=seq)
    masses_and_intensities = list(zip(spectre.masses, spectre.probs))
    return Spectrum(confs=masses_and_intensities)

