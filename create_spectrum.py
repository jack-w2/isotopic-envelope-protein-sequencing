from IsoSpecPy import IsoDistribution, IsoTotalProb
from masserstein import Spectrum
from analyse_spectrum import analyse_spectrum
from sequence import Seq
from priority_queue import Queue, QueueItem
from icecream import ic
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt
import csv
import itertools

log_lines = []


def generate_prefixes_and_suffixes(seq):
    """Split given seq to create all possible suffixes and prefixes."""
    seqs = itertools.chain.from_iterable([Seq(seq[:i], 'pref'), Seq(seq[i:], 'suf')] for i in range(1, len(seq)))
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
    """Add chemical and electronic noise to the given spectrum. Works in place."""
    # plt.figure()
    # plt.title('raw spectrum')
    # masserstein_spectrum.plot()
    masserstein_spectrum.normalize()    # maybe not necessary
    # plt.figure()
    # plt.title('normalized spectrum')
    # masserstein_spectrum.plot()
    masserstein_spectrum.add_chemical_noise(nb_of_noise_peaks, noise_fraction)
    # plt.figure()
    # plt.title('spectrum with chemical noise')
    # masserstein_spectrum.plot()
    masserstein_spectrum.gaussian_smoothing(sd=0.3)
    # plt.figure()
    # plt.title('spectrum with gaussian smoothing')
    # masserstein_spectrum.plot()
    masserstein_spectrum.add_gaussian_noise(sd)
    # plt.figure()
    # plt.title('spectrum with electronic noise')
    # masserstein_spectrum.plot()


def scoring_function(seqs_to_test, experimental_spectrum, aa_one_letter_codes, parameters_set):
    """Scoring function for the analyse_spectrum function."""
    spectra_to_test_hydrogen_removed = [create_raw_spectrum_from_fasta(seq, remove_one_hydrogen=True) for seq in seqs_to_test]
    proportions_hydrogen_removed = analyse_spectrum(experimental_spectrum, spectra_to_test_hydrogen_removed, mtd=parameters_set[0], mdc=parameters_set[1], mmd=parameters_set[2], mtd_th=parameters_set[3])['proportions']
    spectra_to_test_standard = [create_raw_spectrum_from_fasta(seq) for seq in seqs_to_test]
    proportions_standard = analyse_spectrum(experimental_spectrum, spectra_to_test_standard, mtd=parameters_set[0], mdc=parameters_set[1], mmd=parameters_set[2], mtd_th=parameters_set[3])['proportions']
    spectra_to_test_extra_hydrogen = [create_raw_spectrum_from_fasta(seq, add_extra_hydrogen=True) for seq in seqs_to_test]
    proportions_extra_hydrogen = analyse_spectrum(experimental_spectrum, spectra_to_test_extra_hydrogen, mtd=parameters_set[0], mdc=parameters_set[1], mmd=parameters_set[2], mtd_th=parameters_set[3])['proportions']

    score = []
    for proportion in range(len(proportions_standard)):
        numerator = proportions_standard[proportion]
        denominator = max(proportions_extra_hydrogen[proportion], proportions_hydrogen_removed[proportion])
        log_lines[-1].append(aa_one_letter_codes[proportion])
        log_lines[-1].append(numerator)
        log_lines[-1].append(denominator)
        if denominator == 0.0:
            denominator = 1
        score.append(numerator / denominator**parameters_set[4])
        log_lines[-1].append(numerator / denominator**parameters_set[4])

    best_letters_all = dict(zip(aa_one_letter_codes, score))
    max_proportion = max(best_letters_all.values())
    best_letters = [k for k,v in best_letters_all.items() if v > 0.5 * max_proportion]
    return best_letters, best_letters_all


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


def create_spectrum(config):
    """Create experimental spectrum."""
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

    # spectre.plot()
    masserstein_spectrum = Spectrum(confs=masses_and_intensities, label='experimental')
    add_noise(masserstein_spectrum, config['noise']['nb_of_chemical_noise_peaks'], config['noise']['chemical_noise_amount'], config['noise']['electronic_noise_sd'])
    return masserstein_spectrum


def create_raw_spectrum_from_fasta(seq, add_extra_hydrogen=False, remove_one_hydrogen=False):
    """Create raw spectrum for given fasta string (without noise and other stuff)."""
    if add_extra_hydrogen:
        formula_pref = 'H2O'
        formula_suf = 'H2'
    elif remove_one_hydrogen:
        formula_pref = 'O'
        formula_suf = ''
    else:
        formula_pref = 'OH'
        formula_suf = 'H'
    spectre = IsoTotalProb(0.999, fasta=seq.seq, formula=formula_pref if seq.type == 'pref' else formula_suf)
    masses_and_intensities = list(zip(spectre.masses, spectre.probs))
    masserstein_spectrum = Spectrum(confs=masses_and_intensities, label='theoretical')
    masserstein_spectrum.normalize()
    return masserstein_spectrum


def get_replacements_dict(replacements_file='amino_acids_replacements.csv'):
    """Prepare the amino acids replacements dictionary (some amino acids can "sum" to a different amino acid)."""
    replacements_dict = {}
    with open(replacements_file, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            replacements_dict[row[0]] = row[1:]
    return replacements_dict


def calculate_priority(simulated_seq, config, best_letters_all, cost_so_far):
    """
    Calculate priority for priority queue item. Similar to A* algorithm function.
    Uses alpha factor defined in the config file to balance the proportion of each component.
    """
    alpha = config['heuristic_factor']
    proportion = best_letters_all[simulated_seq.seq[-1]] + cost_so_far
    formula_string = str(Seq(config['fasta'], 'full').convert_to_molecular_formula())
    simulated_seq_formula_string = str(simulated_seq.convert_to_molecular_formula())
    heuristic = int(formula_string[1:formula_string.index('H')]) - int(simulated_seq_formula_string[1:simulated_seq_formula_string.index('H')])
    priority = alpha * proportion + (1 - alpha) * heuristic
    log_lines[-1].append(priority)
    log_lines[-1].append(alpha)
    log_lines[-1].append(proportion)
    log_lines[-1].append(heuristic)
    return priority, simulated_seq_formula_string, proportion


def main():
    config = load_config_file('config.toml')
    exp_spectrum = create_spectrum(config)
    exp_spectrum.normalize(target_value=100000.0)
    simulated_seq = Seq('MALW', 'pref')
    cost_so_far = 0
    parameters_set = [0.1, 1e-06, 0.4, 0.4, 0.3]
    parameters_info = f'{simulated_seq}, parameters: {parameters_set}'
    print(parameters_info)
    log_lines.append(parameters_info)
    q = Queue()
    for i in range(110):
        #print('queue:', q)
        if q:
            considered_state = q.dequeue()
            simulated_seq = considered_state.seq
            cost_so_far = considered_state.cost_so_far
        simulated_seq = Seq(config['fasta'][:i], 'pref')
        print('simulated_seq:', simulated_seq)
        print('expected letter:', config['fasta'][i])
        log_lines.append([simulated_seq])
        best_letters, best_letters_all = find_next_best_letter(simulated_seq, exp_spectrum, parameters_set)
        print('best_letters:', best_letters)
        print()
        log_lines[-1].append(best_letters)
        if config['fasta'][i] not in best_letters:
            break
        exp_spectrum.normalize(1000.0)
        # next_steps = [simulated_seq + letter[0] for letter in best_letters]
        # print(next_steps)
        # next_steps = list(map(create_raw_spectrum_from_fasta, next_steps))
        # Spectrum.plot_all([exp_spectrum] + next_steps, cmap=['black', 'blue', "red", "yellow"])
        replacements_dict = get_replacements_dict()
        for letter in best_letters:
            possible_simulated_seq = Seq('', simulated_seq.type)
            if letter in replacements_dict.keys():
                set_to_test = set(best_letters) & set(replacements_dict[letter])
                if len(set_to_test) > 0:
                    found_good_replacement = False
                    for letter2 in best_letters:
                        letter_proportion = best_letters_all[letter2]
                        replacements_for_best_letter = replacements_dict[letter]
                        replacement1_proportion = best_letters_all[replacements_for_best_letter[0]]
                        replacement2_proportion = best_letters_all[replacements_for_best_letter[1]]
                        conditions = [
                            letter2 in replacements_dict[letter],
                            letter_proportion > 0.2 * replacement1_proportion
                            or
                            letter_proportion > 0.2 * replacement2_proportion
                        ]
                        if all(conditions):
                            # select letter2 as next letter in simulated_seq
                            # break if letter found
                            possible_simulated_seq.seq = simulated_seq.seq + letter2
                            found_good_replacement = True
                            break
                    if not found_good_replacement:
                        possible_simulated_seq.seq = simulated_seq.seq + letter
                else:
                    possible_simulated_seq.seq = simulated_seq.seq + letter
            else:
                possible_simulated_seq.seq = simulated_seq.seq + letter
            priority, simulated_seq_formula_string, cost_so_far = calculate_priority(possible_simulated_seq, config, best_letters_all, cost_so_far)
            q.enqueue(QueueItem(cost_so_far + priority, possible_simulated_seq, simulated_seq_formula_string, cost_so_far))


if __name__ == "__main__":
    main()
    with open('parameters.csv', 'w') as lf:
        writer = csv.writer(lf, dialect='excel')
        writer.writerows(log_lines)
