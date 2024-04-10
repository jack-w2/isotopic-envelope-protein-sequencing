from IsoSpecPy import IsoDistribution, IsoTotalProb
from masserstein import Spectrum
from analyse_spectrum import analyse_spectrum
from collections import Counter
from sequence import Seq
from icecream import ic
from datetime import datetime
from pathlib import Path
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
    # print("AAAAA", ret)
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


def scoring_function(seqs_to_test, experimental_spectrum, aa_one_letter_codes, parameters_set):
    """Scoring function idea to use analyse_spectrum"""
    spectra_to_test = [create_raw_spectrum_from_fasta(seq) for seq in seqs_to_test]
    proportions = analyse_spectrum(experimental_spectrum, spectra_to_test, mtd=parameters_set[0], mdc=parameters_set[1], mmd=parameters_set[2], mtd_th=parameters_set[3])['proportions']

    def select_n_best(n):
        counter = Counter(dict(zip(aa_one_letter_codes, proportions)))
        return counter.most_common(n)
    return select_n_best(3)


def find_next_best_letter(seq_to_test, experimental_spectrum, parameters_set, aa_file='amino_acids.csv'):
    """For given seq find next best amino acid."""
    def get_aa_one_letter_codes(aa_file):
        """Get amino acids one letter codes from given file."""
        with open(aa_file, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            return [row[1] for row in reader]
    aa_one_letter_codes = get_aa_one_letter_codes(aa_file)
    seqs_to_test = [Seq(seq_to_test.seq + aa, seq_to_test.type) for aa in aa_one_letter_codes]
    best_letters = scoring_function(seqs_to_test, experimental_spectrum, aa_one_letter_codes, parameters_set)
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


def get_replacements_dict(replacements_file='amino_acids_replacements.csv'):
    replacements_dict = {}
    with open(replacements_file, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            replacements_dict[row[0]] = row[1:]
    return replacements_dict


# tests
model_seq = load_config_file('config.toml')['fasta']


# Functions for testing
def check_if_matches_model_seq(seq, model_seq):
    """Check if last letter from guessed seq is the same in the model_seq."""
    return seq[-1] == model_seq[len(seq) - 1]


def give_helping_hand(seq, model_seq):
    """Return one next correct letter from model_seq."""
    return model_seq[len(seq) - 1]


mtd = [0.1, 0.01, 0.001, 0.0001]
mdc = [1e-7, 1e-8, 1e-9, 1e-10]
mmd = [-1, 0.2, 0.4, 0.6]
mtd_th = [None, 0.1, 0.4, 0.6]
parameters_matrix = itertools.product(mtd, mdc, mmd, mtd_th)

log_file_name = f'log_file_{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}.txt'
for parameters_set in parameters_matrix:
    lines_to_file = []
    exp_spectrum = create_spectrum()
    exp_spectrum.normalize(target_value=100000.0)
    simulated_seq = Seq('MALW', 'pref')
    parameters_info = f'{simulated_seq}, parameters: {parameters_set}'
    lines_to_file.append(parameters_info)
    print(parameters_info)
    for i in range(110):
        best_letters = find_next_best_letter(simulated_seq, exp_spectrum, parameters_set)
        print(best_letters)
        # plt.close()
        print("Norm:", sum([x[1] for x in exp_spectrum.confs]))
        exp_spectrum.normalize(1000.0)
        next_steps = [simulated_seq + letter[0] for letter in best_letters]
        print(next_steps)
        # next_steps = list(map(create_raw_spectrum_from_fasta, next_steps))
        # Spectrum.plot_all([exp_spectrum] + next_steps, cmap=['black', 'blue', "red", "yellow"])
        replacements_dict = get_replacements_dict()
        if best_letters[0][0] in replacements_dict.keys():
            for letter in [l[0] for l in best_letters]:
                if letter in replacements_dict[best_letters[0][0]]:
                    # select letter as next letter in simulated_seq
                    # break if letter found
                    simulated_seq.seq += letter
                    break
        else:
            simulated_seq.seq += best_letters[0][0]
        print(simulated_seq)
        if not check_if_matches_model_seq(simulated_seq.seq, model_seq):
            lines_to_file.append(best_letters)
            lines_to_file.append(simulated_seq)
            break

    Path('tests').mkdir(exist_ok=True)
    with open(f'tests/{log_file_name}', 'a') as log_file:
        log_file.writelines([f'{str(line)}\n' for line in lines_to_file])
        log_file.write('\n\n')
