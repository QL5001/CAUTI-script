#Receiver operating characteristic (ROC) curve of Logistic regression classification between 13 high and 16 low biofilm E. coli isolates using sparse partial least squares discriminant analysis (sPLSDA) Component-1 values
#Same method applied to the Receiver operating characteristic (ROC) curve of logistic regression classfication between two B2 sub-clades using sparse principal component analysis (sPCA) PC1 values
#Begin

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc

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
kfold = list(StratifiedKFold(n_splits=5, random_state=1, shuffle=True).split(X_train_std, y_train_))  #Five-fold cross validations

lr_fpr = []
lr_tpr = []
lr_auc = []
lr_thresholds = []

for k, (train, test) in enumerate(kfold):
    probas = lr.fit(X_train_std[train], y_train_[train]).predict_proba(X_train_std[test])
    fpr, tpr, thresholds = roc_curve(y_train_[test], probas[:, 0], pos_label=0, drop_intermediate=False)
    lr_fpr.append(fpr)  #False positive rate
    lr_tpr.append(tpr)  #True positive rate
    roc_auc = auc(fpr, tpr)  #Area under the ROC curve
    lr_auc.append(roc_auc)
    lr_thresholds.append(thresholds)  #Threshold in logistic regression

lr_fpr_df = pd.DataFrame(lr_fpr)
lr_fpr_df.to_csv('Biofilm_LogisticRegression_ROC_fpr.csv')

lr_tpr_df = pd.DataFrame(lr_tpr)
lr_tpr_df.to_csv('Biofilm_LogisticRegression_ROC_tpr.csv')

lr_auc_df = pd.DataFrame(lr_auc)
lr_auc_df.to_csv('Biofilm_LogisticRegression_AUC.csv')

lr_thresholds_df = pd.DataFrame(lr_thresholds)
lr_thresholds_df.to_csv('Biofilm_LogisticRegression_Thresholds.csv')

#Exported csv files are used for further data analyses and drawing figures in Excel or Prism
#End
