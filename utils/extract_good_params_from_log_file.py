from pathlib import Path
log_file_path = Path('../tests/good_log_file.txt')
model_seq = 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGLVEQCCTSLCSLYQLENYCN'
with open(log_file_path, 'r') as log_file:
    log_text = log_file.read()
log_lines = log_text.splitlines()
good_params = []
for index, line in enumerate(log_lines):
    if line.startswith('Seq(') and line[-2].isalpha():
        seq = line[len('Seq('):line.index(',')]
        if seq == model_seq or seq[:-1] == model_seq[:-1]:
            params_line = log_lines[index - 2]
            params = params_line[params_line.index('parameters: ') + len('parameters: ') + 1:-1].replace(' ', '')
            print(params)
            good_params.append(params)

with open(log_file_path.parent / f'good_params_{log_file_path.name}', 'w') as good_params_file:
    good_params_file.write('\n'.join(good_params))
