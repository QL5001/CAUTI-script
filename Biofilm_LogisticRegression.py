#Logistic regression classification between 13 high and 16 low biofilm E. coli isolates using sparse partial least squares discriminant analysis (sPLSDA) Component-1 values
#Same method applied to the logistic regression classfication between two B2 sub-clades using sparse principal component analysis (sPCA) PC1 values
#Begin

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

TrainData = pd.read_csv('Biofilm_SPLSDA-PC1.csv', header=0)  #Open the file for logistic regression analysis including the 13 low and 16 high biofilm isolates' sPLSDA Component-1 values
StrainNum = TrainData['Strain']
BiofilmLevel = TrainData['Biofilm']

X_train = TrainData['PC1']
X_train_ = X_train.to_numpy()
X_train_ = X_train_.reshape(-1, 1)
sc = StandardScaler(with_mean=True, with_std=True)
X_train_std = sc.fit_transform(X_train_)

biofilm_le = LabelEncoder()
y_train = biofilm_le.fit_transform(TrainData['Biofilm'])
y_train_ = y_train

lr = LogisticRegression(penalty='none', random_state=1)
scores = []
kfold = list(StratifiedKFold(n_splits=5, random_state=1, shuffle=True).split(X_train_std, y_train_))  #Five-fold cross validations
for k, (train, test) in enumerate(kfold):
    lr.fit(X_train_std[train], y_train_[train])
    score = lr.score(X_train_std[test], y_train_[test])
    scores.append(score)

scores_df = pd.DataFrame(scores)  #Prediction accuracy
scores_df.to_csv('Biofilm_LogisticRegression_Accuracy.csv')

#Exported csv files are used for further data analyses and drawing figures in Excel or Prism
#End
