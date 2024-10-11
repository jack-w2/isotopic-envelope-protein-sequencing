import csv
import heapq
from pathlib import Path
from itertools import groupby
from prettytable import PrettyTable

logs = list(Path('tests/logs').glob('*'))

summary = PrettyTable(field_names=['protein', 'algorithm_performance', 'random_performance'])

for log_file in logs:
    algorithm_performance = 0
    random_performance = 0
    with open(log_file, 'r') as f:
        reader = csv.reader(f, dialect='excel')
        headers = next(reader, None)
        log = {row[0]: row[1:] for row in reader}
    headers_dict = {header: number for number, header in enumerate(headers[1:])}
    seqs_groups = groupby(log.keys(), key=len)
    for seq_len, seqs in seqs_groups:
        queue = []
        for seq in list(seqs):
            heapq.heappush(queue, (log[seq][headers_dict['priority']], seq))
            if log[seq][-1]:
                random_performance += 1 / seq_len
        if log[queue[0][1]][-1]:
            algorithm_performance += 1
    summary.add_row([log_file.stem, algorithm_performance, random_performance])
print(summary)
