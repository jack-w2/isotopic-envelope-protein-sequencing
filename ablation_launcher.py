"""Launcher for the Section 4.6 ablation experiment: runs a small protein subset
under alpha in {0.4, 0.8} and under isobaric-resolution disabled (alpha=0.6),
all with the reduced budget T=30 of ablation_experiment.py."""
import multiprocessing
import subprocess
from pathlib import Path


def load_fasta():
    def parse_fasta(path):
        seqs, header, seq = [], None, []
        with open(path) as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if header is not None:
                        seqs.append((header, ''.join(seq)))
                    header = line[1:]; seq = []
                else:
                    seq.append(line)
            if header is not None:
                seqs.append((header, ''.join(seq)))
        return seqs
    return {h.split('|')[1]: s.replace('I', 'L') for h, s in parse_fasta('tests/proteins.fasta')}

PROTEIN_IDS = ['A6NNB3', 'O00422', 'O14508', 'O15263', 'P01308', 'O95298', 'P15090', 'P10620', 'O43181', 'P07305']

CONFIGS = [
    ('alpha04', ['-a', '0.4']),
    ('alpha06_seeded', ['-a', '0.6']),
    ('alpha08', ['-a', '0.8']),
    ('noisobaric', ['-a', '0.6', '--no-isobaric']),
]


def run_job(args):
    protein_id, seq, config_name, extra_args, seed = args
    outdir = Path(f'tests/ablation/{config_name}')
    outdir.mkdir(parents=True, exist_ok=True)
    logfile = outdir / f'{protein_id}.csv'
    cmd = ['python3', 'ablation_experiment.py', '-f', seq, '-l', str(logfile), '-s', str(seed)] + extra_args
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    return protein_id, config_name, result.returncode, result.stderr[-500:] if result.returncode != 0 else ''


def main():
    fasta = load_fasta()
    jobs = []
    for idx, pid in enumerate(PROTEIN_IDS):
        seq = fasta[pid]
        seed = 1000 + idx  # same seed shared across all configs for this protein
        for config_name, extra_args in CONFIGS:
            jobs.append((pid, seq, config_name, extra_args, seed))
    print(f'Running {len(jobs)} jobs...')
    with multiprocessing.Pool(8) as pool:
        for pid, config_name, rc, err in pool.imap_unordered(run_job, jobs):
            status = 'OK' if rc == 0 else f'FAIL({rc}): {err}'
            print(f'{pid} / {config_name}: {status}', flush=True)


if __name__ == '__main__':
    main()
