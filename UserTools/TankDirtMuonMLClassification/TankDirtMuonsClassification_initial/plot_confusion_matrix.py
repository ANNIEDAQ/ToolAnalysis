import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels

data = pd.read_csv("predictionsKeras.csv")
print(data.head())
class_names=["Background","Signal"] #= data0["11"].values
class_types=["Background","Signal"]

positives = data.loc[(data['Yvalues']==1)]
print("P: ", data.loc[(data['Yvalues']==1)].shape)
print("Pred P (TP): ", positives.loc[(data['Prediction']==1)].shape)
print("False N: ",positives.loc[(data['Prediction']==0)].shape)
negatives= data.loc[(data['Yvalues']==0)]
print("N: ", data.loc[(data['Yvalues']==0)].shape)
print("Pred N (TN): ", negatives.loc[(data['Prediction']==0)].shape)
print("False P: ",negatives.loc[(data['Prediction']==1)].shape)

#convert strings to numbers: 
#data1 = data0.replace("muon", 0)
#data = data1.replace("electron", 1)
#print(data.head())

def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    print(unique_labels(y_true,y_pred))
    #classes2 = class_types[unique_labels(y_true, y_pred)]
    #print(classes2)
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=class_types, yticklabels=class_types,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax


np.set_printoptions(precision=2)

#print("data[11]:")
#print(data["11"])
#print("data[Prediction]:")
#print(data["Prediction"])
#print(class_names[0])
#print(class_names[1])
#print(class_names[2])
#print(class_names[3])
#print(class_names[4])
#print(class_names[5])
# Plot non-normalized confusion matrix
plot_confusion_matrix(data["Yvalues"], data["Prediction"], classes=class_names,
                      title='Confusion matrix without normalization')

# Plot normalized confusion matrix
plot_confusion_matrix(data["Yvalues"], data["Prediction"], classes=class_names, normalize=True,
                      title=' Normalized confusion matrix')

plt.show()
