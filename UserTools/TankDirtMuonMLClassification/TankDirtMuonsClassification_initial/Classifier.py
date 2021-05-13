import sys
import glob
import numpy as np
import pandas as pd
import tensorflow as tf
import tempfile
import random
import csv
import matplotlib
#matplotlib.use('Agg')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from array import array
from sklearn import datasets
from sklearn import metrics
from sklearn import model_selection
from sklearn import preprocessing
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier

from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from sklearn.utils import shuffle
from sklearn import svm, datasets
from sklearn.model_selection import train_test_split 
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels

#Create dataframe of both tank and dirt muon data
Tdf = pd.read_csv("tank_muons_data_1000evts.csv") 
Ddf = pd.read_csv("dirt_muons_data_1000evts.csv")

T_rows = len(Tdf.index) 
D_rows = len(Ddf.index)
Total_rows = T_rows + D_rows 

Tank_1 = [1]*T_rows
Dirt_0 = [0]*D_rows
Tdf.insert(1, "T or D", Tank_1, True)
Ddf.insert(1, "T or D", Dirt_0, True)

frames = [Tdf, Ddf]
df = pd.concat(frames, ignore_index = True) 
df = df.sample(frac=1)
rows = int(len(df.index))
print(df.head())

#Isolate input and output variables for algorithm 
#data = df . values
Input_data0 = df.loc[:,"pmt_charge_0":"pmt_charge_132"]
Input_data1 = df.loc[:,"T_0":"T_9999"]
Input_data2 = df.loc[:,"numberofPEs"]
Input_data_frames = [ Input_data0 , Input_data1 , Input_data2 ]
Input_data = pd.concat(Input_data_frames, axis = 1)
Output_data = df.loc[:,"T or D"]
print("Defined input and output data")
print(Output_data)

#Train and Test Splitter
trainX = np.array(Input_data[:1400])
testX = np.array(Input_data[1400:])

dataY0 = np.array(Output_data[:1400])
print (dataY0)
print(len(dataY0))
Y = dataY0
print("before Y.shape ",Y.shape)
Y=Y.reshape(-1)
print("after: ",Y.shape)
print("Y: ",Y[0])

Ytest = np.array(Output_data[1400:])
Ytest = Ytest.reshape(-1)
print(Ytest)
print("Ytest length: " + str(len(Ytest)))

# Scale data (training set) to 0 mean and unit standard deviation.
scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(trainX)
print("X.shape: ", X.shape)
print("types: ",type(X), " ",type(Y))
Xtest = scaler.fit_transform(testX)
print("testX.shape: ", testX.shape," Ytest.shape ",Ytest.shape)

# encode class values as integers
encoder = LabelEncoder()
encoder.fit(Output_data[:1400])
encoded_Y = encoder.transform(Output_data[:1400])
print("encoded ..")

def create_model():
    # create model
    model = Sequential()
    model.add(Dense(1000, input_dim=10134, activation='relu'))
    model.add(Dense(60, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))   
    # Compile model
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    return model

estimator = KerasClassifier(build_fn=create_model, epochs=100, batch_size=5, verbose=0)
estimator.fit(X,encoded_Y,verbose=0)
print("Model defined..")

'''
# evaluate baseline model using a pipeline with: i)stratified cross validation and ii)a standardized dataset
estimators = []
estimators.append(('standardize', StandardScaler()))
estimators.append(('mlp', KerasClassifier(build_fn=create_model, epochs=100, batch_size=5, verbose=0)))
pipeline = Pipeline(estimators)
kfold = StratifiedKFold(n_splits=10, shuffle=True)
results = cross_val_score(pipeline, X, encoded_Y, cv=kfold)
print("Standardized: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))
'''
'''
#evaluate this model using stratified cross validation in the scikit-learn
#estimator = KerasClassifier(build_fn=create_model, epochs=100, batch_size=5, verbose=0)
kfold = StratifiedKFold(n_splits=10, shuffle=True)
results = cross_val_score(estimator, X, encoded_Y, cv=kfold)
print("Baseline: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))
#pred_results = cross_val_score(estimator, Xtest, Ytest, cv=kfold)
#print("Baseline: %.2f%% (%.2f%%)" % (pred_results.mean()*100, pred_results.std()*100))
'''
# checkpoint
filepath="weights_best.hdf5"
checkpoint = ModelCheckpoint(filepath, monitor='val_accuracy', verbose=1, save_best_only=True, save_weights_only=True, mode='auto')
callbacks_list = [checkpoint]
# Fit the model
history = estimator.fit(X, encoded_Y, validation_split=0.33, epochs=50, batch_size=5, callbacks=callbacks_list, verbose=0)
#I can add: validation_data=(testX,testy) instead of validation_split

# summarize history for accuracy:
f, ax2 = plt.subplots(1,1)
ax2.plot(history.history['accuracy'])
ax2.plot(history.history['val_accuracy'])
ax2.set_title('model accuracy')
ax2.set_ylabel('Performance')
ax2.set_xlabel('Epochs')
#ax2.set_xlim(0.,10.)
ax2.legend(['train', 'test'], loc='upper left')
plt.savefig("keras_train_testAcc.pdf")

# summarize history for loss
f, ax2 = plt.subplots(1,1)
ax2.plot(history.history['loss'], label='train')
ax2.plot(history.history['val_loss'], label='test')
ax2.set_title('model loss')
ax2.set_ylabel('Performance')
ax2.set_xlabel('Epochs')
#ax2.set_xlim(0.,10.)
#ax2.legend(['train', 'test'], loc='upper left')
plt.savefig("keras_train_testLoss.pdf")

#ROC curve:
from sklearn.metrics import roc_curve

Ypred = estimator.predict(Xtest).ravel()
Ypred =np.vstack(Ypred)
Ytest=np.vstack(Ytest)
#store predictions to csv:
df1 = pd.DataFrame(data=Ytest,columns=['Yvalues'])
df2 = pd.DataFrame(data=Ypred,columns=['Prediction'])
print("check df head: ",df1.head()," ", df2.head())
df = pd.concat([df1,df2], axis=1)
print(df.head())
df.to_csv("predictionsKeras.csv")

Y_probs = estimator.predict_proba(Xtest)
#print("Ytest: ",Ytest[:10]," and predicts ", Ypred[:10])
#print(type(Xtest)," ",type(Ytest)," Ytest.shape ",Ytest.shape," Ypred.shape ",Ypred.shape)
fpr_keras, tpr_keras, thresholds_keras = roc_curve(Ytest, Y_probs[:, 1], pos_label =1) #Ypred)
#AUC:
from sklearn.metrics import auc
auc_keras = auc(fpr_keras, tpr_keras)

f1, ax1 = plt.subplots(1,1)
ax1.plot([0, 1], [0, 1], 'k--')
ax1.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
ax1.set_xlabel('False positive rate')
ax1.set_ylabel('True positive rate')
ax1.set_title('ROC curve')
ax1.legend(loc='best')
plt.savefig("ROC_curve.pdf")
#plt.show()

print(metrics.classification_report(Ytest, Ypred))


