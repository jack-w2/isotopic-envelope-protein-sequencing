"""
Supplementary ablation experiment for Chapter 4 (Section 4.6) of the thesis.

Reproduces the search of create_spectrum.py with three knobs exposed that are
hardcoded in the main pipeline: the priority weight alpha, whether the isobaric
substitution mechanism (Section 3.6) is active, and the iteration budget.

Unlike the main autonomous corpus (tests/logs, budget T=110, historical protocol),
this experiment uses a smaller, self-contained budget (T=30) chosen to keep the
runtime of a small parameter sweep tractable; it is not intended to reproduce the
main corpus exactly, only to let alpha and the isobaric mechanism be compared
against each other under matched conditions.

Usage:
    python3 ablation_experiment.py -f <fasta> -a <alpha> --no-isobaric -l <logfile> -T <budget>
"""
import argparse
import copy
import csv

import numpy.random

from create_spectrum import (
    load_config_file, create_spectrum, find_next_best_letter,
    calculate_priority, get_replacements_dict, log_lines
)
from priority_queue import Queue, QueueItem
from sequence import Seq


def run(fasta, alpha, isobaric_enabled, budget, logfile, seed=None):
    log_lines.clear()
    config = load_config_file('config.toml')
    config = copy.deepcopy(config)
    config['fasta'] = fasta
    config['heuristic_factor'] = alpha

    if seed is not None:
        numpy.random.seed(seed)
    exp_spectrum = create_spectrum(config)
    exp_spectrum.normalize(target_value=100000.0)
    simulated_seq = Seq('', 'pref')
    cost_so_far = 0
    parameters_set = [0.1, 1e-06, 0.4, 0.4, 0.3]
    parameters_set.extend([
        config['factors_from_linear_model']['numerator_factor'],
        config['factors_from_linear_model']['denominator_factor'],
        config['factors_from_linear_model']['score_factor'],
    ])

    q = Queue()
    replacements_dict = get_replacements_dict()

    for i in range(budget):
        if q:
            considered_state = q.dequeue()
            simulated_seq = considered_state.seq
            cost_so_far = considered_state.cost_so_far
        if i >= len(fasta):
            break
        expected_letter = fasta[i]
        best_letters, best_letters_all = find_next_best_letter(simulated_seq, exp_spectrum, parameters_set)
        if expected_letter not in best_letters:
            break
        exp_spectrum.normalize(1000.0)

        for letter in best_letters:
            possible_simulated_seq = Seq('', simulated_seq.type)
            if isobaric_enabled and letter in replacements_dict.keys():
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
                            or letter_proportion > 0.2 * replacement2_proportion
                        ]
                        if all(conditions):
                            possible_simulated_seq.seq = simulated_seq.seq + letter2
                            found_good_replacement = True
                            break
                    if not found_good_replacement:
                        possible_simulated_seq.seq = simulated_seq.seq + letter
                else:
                    possible_simulated_seq.seq = simulated_seq.seq + letter
            else:
                possible_simulated_seq.seq = simulated_seq.seq + letter
            priority, simulated_seq_formula_string, cost_so_far_new = calculate_priority(
                possible_simulated_seq, config, best_letters_all, cost_so_far)
            q.enqueue(QueueItem(cost_so_far_new + priority, possible_simulated_seq, simulated_seq_formula_string, cost_so_far_new))
            if possible_simulated_seq.seq == fasta[:i + 1]:
                log_lines[possible_simulated_seq.seq].append(True)

    headers = ['seq', 'numerator', 'denominator', 'score', 'priority', 'alpha', 'proportion', 'heuristic', 'is_correct_option']
    with open(logfile, 'w') as lf:
        writer = csv.writer(lf, dialect='excel')
        writer.writerow(headers)
        for key, val in log_lines.items():
            values_count = len(val)
            if values_count == 7:
                writer.writerow([key, *val, False])
            elif values_count == 8:
                writer.writerow([key, *val])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-a', '--alpha', type=float, default=0.6)
    parser.add_argument('--no-isobaric', action='store_true')
    parser.add_argument('-T', '--budget', type=int, default=30)
    parser.add_argument('-l', '--logfile', type=str, required=True)
    parser.add_argument('-s', '--seed', type=int, default=None)
    args = parser.parse_args()
    run(args.fasta, args.alpha, not args.no_isobaric, args.budget, args.logfile, args.seed)


if __name__ == '__main__':
    main()
