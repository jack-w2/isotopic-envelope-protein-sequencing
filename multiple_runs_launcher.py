import os
import multiprocessing
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser


def prepare_proteins_list(file_path):
    with open(file_path, 'r') as handle:
        return [(protein[0], protein[1].replace('I', 'L')) for protein in SimpleFastaParser(handle)]

def launcher(protein):
    logs_directory = 'logs_25_params_from_linear_model'
    logfile_name = f'{protein[0].split("|")[1]}.csv'
    print(logfile_name)
    Path(f'tests/{logs_directory}').mkdir(exist_ok=True)
    os.system(f'python create_spectrum.py -f {protein[1]} -l tests/{logs_directory}/{logfile_name}')

def main():
    proteins = prepare_proteins_list('tests/proteins.fasta')
    with multiprocessing.Pool() as pool:
        pool.map(launcher, proteins)


if __name__ == "__main__":
    main()
