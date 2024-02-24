from IsoSpecPy import IsoDistribution, IsoTotalProb
from masserstein import Spectrum
from analyse_spectrum import analyse_spectrum
from collections import Counter
from sequence import Seq
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
    seqs = itertools.chain.from_iterable([Seq(seq[:i], 'pref'), Seq(seq[i:], 'suf')] for i in range(1, len(seq)))
    ret = list(seqs)
    print("AAAAA", ret)
    return ret


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
    return
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


def scoring_function(seqs_to_test, experimental_spectrum, aa_one_letter_codes):
    """Scoring function idea to use analyse_spectrum"""
    spectra_to_test = [create_raw_spectrum_from_fasta(seq) for seq in seqs_to_test]
    proportions = analyse_spectrum(experimental_spectrum, spectra_to_test)['proportions']

    def select_n_best(n):
        counter = Counter(dict(zip(aa_one_letter_codes, proportions)))
        return counter.most_common(n)
    return select_n_best(3)


def find_next_best_letter(seq_to_test, experimental_spectrum, aa_file='amino_acids.csv'):
    """For given seq find next best amino acid."""
    def get_aa_one_letter_codes(aa_file):
        """Get amino acids one letter codes from given file."""
        with open(aa_file, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            return [row[1] for row in reader]
    aa_one_letter_codes = get_aa_one_letter_codes(aa_file)
    seqs_to_test = [Seq(seq_to_test.seq + aa, seq_to_test.type) for aa in aa_one_letter_codes]
    best_letters = scoring_function(seqs_to_test, experimental_spectrum, aa_one_letter_codes)
    return best_letters


def create_spectrum():
    config = load_config_file('config.toml')

    seq = config['fasta']
    seqs = generate_prefixes_and_suffixes(seq)

    envelopes = [IsoTotalProb(0.999, fasta=seqs[s].seq, formula='OH' if seqs[s].type == 'pref' else 'H') for s in range(len(seqs))]
    intensities = []
    for i in config['probs']:
        intensities.extend([i * 0.5] * 2)

    intensities = scale_intensities(intensities)

    spectre = IsoDistribution.LinearCombination(envelopes, intensities)
    masses_and_intensities = list(zip(spectre.masses, spectre.probs))

    #print(masses_and_intensities)

    # spectre.plot()
    masserstein_spectrum = Spectrum(confs=masses_and_intensities, label='experimental')
    add_noise(masserstein_spectrum, config['noise']['nb_of_noise_peaks'], config['noise']['noise_fraction'], config['noise']['gaussian_noise_sd'])
    return masserstein_spectrum


def create_raw_spectrum_from_fasta(seq):
    """Create raw spectrum for given fasta string (without noise and other stuff)."""
    spectre = IsoTotalProb(0.999, fasta=seq.seq, formula='OH' if seq.type == 'pref' else 'H')
    masses_and_intensities = list(zip(spectre.masses, spectre.probs))
    masserstein_spectrum = Spectrum(confs=masses_and_intensities, label='theoretical')
    masserstein_spectrum.normalize()
    return masserstein_spectrum


# tests
exp_spectrum = create_spectrum()
exp_spectrum.normalize(target_value=100000.0)
#Spectrum.plot_all([exp_spectrum, create_raw_spectrum_from_fasta(Seq('MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGG', 'pref'))], cmap=['green', 'blue'])
simulated_seq = Seq('MALW', 'pref')
for i in range(110):
    best_letters = find_next_best_letter(simulated_seq, exp_spectrum)
    print(best_letters)
    # plt.figure()
    # plt.title('result')
    # exp_spectrum.plot(show=False, color='red')
    # create_raw_spectrum_from_fasta(simulated_seq).plot(show=False, color='blue')
    plt.close()
    print("Norm:", sum([x[1] for x in exp_spectrum.confs]))
    exp_spectrum.normalize(1000.0)
    next_steps = [simulated_seq + letter[0] for letter in best_letters]
    print(next_steps)
    next_steps = list(map(create_raw_spectrum_from_fasta, next_steps))
    Spectrum.plot_all([exp_spectrum] + next_steps, cmap=['black', 'blue', "red", "yellow"])
    simulated_seq.seq += best_letters[0][0]
    # plt.legend()
    # plt.show()
    print(simulated_seq)
