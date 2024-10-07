import csv
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from sklearn.model_selection import train_test_split

logs = list(Path('tests/logs').glob('*'))

seqs = []
x = []
y = []
for log in logs:
    with open(log, 'r') as f:
        reader = csv.reader(f, dialect='excel')
        headers = next(reader, None)[1:-1]
        for row in reader:
            seqs.append(row[0])
            x.append(map(float,row[1:-1]))
            y.append(row[-1])

seqs = pd.Series(seqs)
x = pd.DataFrame(x, columns=headers)
y = pd.Series(1 if i == 'True' else 0 for i in y)

x['numerator^2'] = x['numerator'] ** 2
x['numerator^3'] = x['numerator'] ** 3
x['denominator^2'] = x['denominator'] ** 2
x['denominator^3'] = x['denominator'] ** 3
x.drop(columns=['alpha', 'heuristic', 'priority'], inplace=True)
cols = x.columns.to_list()
cols = ['numerator', 'numerator^2', 'numerator^3', 'denominator', 'denominator^2', 'denominator^3', *cols[2:-4]]
x = x[cols]

pd.set_option('display.max_columns', None)
print(seqs)
print(x)
print(y)

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2)

x_train = pd.DataFrame(x_train)
y_train = pd.Series(y_train)

print(x_train)

sm.add_constant(x_train)
mod = sm.OLS(y_train, x_train)
res = mod.fit()
print(res.summary())

