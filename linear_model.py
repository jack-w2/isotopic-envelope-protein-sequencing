import csv
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_X_y, check_is_fitted, check_array

logs = list(Path('tests/logs').glob('*'))

x = []
y = []
for log in logs:
    with open(log, 'r') as f:
        reader = csv.reader(f, dialect='excel')
        next(reader, None)
        for row in reader:
            x.append(map(float,row[1:-1]))
            y.append(row[-1])
print(x)
y = [1 if i == 'True' else 0 for i in y]
print(y)

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2)

x_train = pd.DataFrame(x_train)
y_train = pd.Series(y_train)


class MyLinearRegression(BaseEstimator, RegressorMixin):
    def __init__(self, fit_intercept=True):

        self.fit_intercept = fit_intercept


    """
    Parameters
    ------------
    column_names: list
            It is an optional value, such that this class knows 
            what is the name of the feature to associate to 
            each column of X. This is useful if you use the method
            summary(), so that it can show the feature name for each
            coefficient
    """
    def fit(self, X, y, column_names=None):

        if column_names is None:
            column_names = []
        if self.fit_intercept:
            X = sm.add_constant(X)

        # Check that X and y have correct shape
        X, y = check_X_y(X, y)


        self.X_ = X
        self.y_ = y

        if len(column_names) != 0:
            X = pd.DataFrame(X)
            column_names.insert(0,'intercept')
            print('X ', X)
            X.columns = column_names

        self.model_ = sm.OLS(y, X)
        self.results_ = self.model_.fit()
        return self


    def predict(self, X):
        # Check is fit had been called
        check_is_fitted(self, 'model_')

        # Input validation
        X = check_array(X)

        if self.fit_intercept:
            X = sm.add_constant(X)
        return self.results_.predict(X)


    def get_params(self, deep = False):
        return {'fit_intercept':self.fit_intercept}


    def summary(self):
        print(self.results_.summary())


model = MyLinearRegression()
model.fit(x_train, y_train)
model.summary()
model.predict(x_train)

scores = cross_val_score(MyLinearRegression(), x_train, y_train, cv=10, scoring='neg_mean_squared_error', n_jobs=-1)
print(scores)
print(scores.mean())
