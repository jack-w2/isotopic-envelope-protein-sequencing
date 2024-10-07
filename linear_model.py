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
        next(reader, None)
        for row in reader:
            seqs.append(row[0])
            x.append(map(float,row[1:-1]))
            y.append(row[-1])

seqs = pd.Series(seqs)
x = pd.DataFrame(x)
y = pd.Series(1 if i == 'True' else 0 for i in y)

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

