import csv
from pathlib import Path
from sklearn import linear_model
from sklearn.model_selection import train_test_split

logs = list(Path('tests/logs').glob('*'))

x = []
y = []
for log in logs:
    with open(log, 'r') as f:
        reader = csv.reader(f, dialect='excel')
        next(reader, None)
        for row in reader:
            x.append(row[1:-1])
            y.append(row[-1])
print(x)
y = [1 if i == 'True' else 0 for i in y]
print(y)

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2)

reg = linear_model.LinearRegression()
import numpy as np
reg.fit(np.array(x_train, dtype=np.float64), y_train)

print(reg.coef_)
