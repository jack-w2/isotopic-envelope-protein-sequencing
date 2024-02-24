from IsoSpecPy import IsoDistribution, IsoTotalProb
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


def generate_prefixes_and_suffixes(seq):
    return [(seq[:i], seq[i:]) for i in range(1, len(seq))]


with open('config.toml', 'rb') as tf:
    config = tomllib.load(tf)


with open('insulin.fasta', 'r') as fasta_file:
    seq = ''.join(fasta_file.read().splitlines()[1:])

seqs = generate_prefixes_and_suffixes(seq)

envelopes = [IsoTotalProb(0.999, fasta=j) for i in seqs for j in i]
intensities = []
for i in config['probs']:
    intensities.extend([i * 0.5] * 2)


spectre = IsoDistribution.LinearCombination(envelopes, intensities)
print([(m, p) for (m, p) in zip(spectre.masses, spectre.probs)])

spectre.plot()
