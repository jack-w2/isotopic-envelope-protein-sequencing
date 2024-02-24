from IsoSpecPy import IsoDistribution, IsoTotalProb
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


def generate_prefixes_and_suffixes(seq):
    return [(seq[:i], seq[i:]) for i in range(1, len(seq))]


def scale_intensities(intensities):
    intensities_sum = sum(intensities)
    return [i/intensities_sum for i in intensities]


def load_config_file(config_file_path):
    with open(config_file_path, 'rb') as tf:
        config = tomllib.load(tf)
    return config


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
    print([(m, p) for (m, p) in zip(spectre.masses, spectre.probs)])

    spectre.plot()


if __name__ == '__main__':
    main()
