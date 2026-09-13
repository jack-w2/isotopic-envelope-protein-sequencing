import argparse
import multiprocessing
import subprocess
from pathlib import Path


def parse_fasta(path):
    seqs, header, seq = [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    seqs.append((header, ''.join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            seqs.append((header, ''.join(seq)))
    return [(header, sequence.replace('I', 'L')) for header, sequence in seqs]


def launcher(job):
    protein, logs_directory, budget, teacher_force = job
    accession = protein[0].split('|')[1]
    logfile = Path(logs_directory) / f'{accession}.csv'
    Path(logs_directory).mkdir(parents=True, exist_ok=True)
    cmd = [
        'python', 'create_spectrum.py',
        '-f', protein[1],
        '-l', str(logfile),
        '-T', str(budget),
    ]
    if teacher_force:
        cmd.append('--teacher-force')
    print(logfile)
    subprocess.run(cmd, check=False)


def main():
    parser = argparse.ArgumentParser(description='Run the pipeline on every protein in a FASTA file.')
    parser.add_argument('--fasta-file', default='tests/proteins.fasta', help='Multi-record FASTA of target proteins.')
    parser.add_argument('--outdir', default='tests/logs_run', help='Directory for per-protein CSV logs.')
    parser.add_argument('-T', '--budget', type=int, default=110, help='Iteration budget passed to create_spectrum.py.')
    parser.add_argument('--teacher-force', action='store_true', help='Calibration regime: force the true prefix at every step.')
    args = parser.parse_args()

    proteins = parse_fasta(args.fasta_file)
    jobs = [(protein, args.outdir, args.budget, args.teacher_force) for protein in proteins]
    with multiprocessing.Pool() as pool:
        pool.map(launcher, jobs)


if __name__ == '__main__':
    main()
