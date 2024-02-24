from IsoSpecPy import IsoDistribution, IsoTotalProb
from masserstein import Spectrum
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt


def generate_prefixes_and_suffixes(seq):
    """Split given seq to create all possible suffixes and prefixes."""
    return [(seq[:i], seq[i:]) for i in range(1, len(seq))]


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


def main():
    config = load_config_file('config.toml')

    seq = config['fasta']
    seqs = generate_prefixes_and_suffixes(seq)

    envelopes = [IsoTotalProb(0.999, fasta=j) for i in seqs for j in i]
    intensities = []
    for i in config['probs']:
        intensities.extend([i * 0.5] * 2)

    intensities = scale_intensities(intensities)

    spectre = IsoDistribution.LinearCombination(envelopes, intensities)
    masses_and_intensities = [(m, p) for (m, p) in zip(spectre.masses, spectre.probs)]

    print(masses_and_intensities)

    # spectre.plot()
    masserstein_spectrum = Spectrum(confs=masses_and_intensities)
    add_noise(masserstein_spectrum)


if __name__ == '__main__':
    main()
