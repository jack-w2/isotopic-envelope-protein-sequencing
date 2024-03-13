with open('log_file_13-03-2024-00-29-57.txt', 'r') as log_file:
    logs = [line for line in log_file.read().splitlines() if line]

print(logs)

highest_len = 0
best_seq = ''
best_letters = ''
best_parameters = ''
for index, line in enumerate(logs):
    if line.startswith('Seq('):
        seq_start = len('Seq(')
        comma_index = line.index(',')
        seq_to_test = line[seq_start:comma_index]
        if len(seq_to_test) > highest_len:
            best_seq = seq_to_test
            best_letters = logs[index - 1]
            best_parameters = logs[index - 2]

print('best_seq:', best_seq)
print('best_letters:', best_letters)
print('best_parameters:', best_parameters)
