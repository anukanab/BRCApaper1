#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib as mlt
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
from sklearn import svm
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier


# Read the PSI matrix

# In[2]:


PSI = pd.read_csv("/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana-New-Analysis/SF/TrueSplice/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-75p.txt", sep = "\t")
print PSI.head()
print PSI.shape


# In[3]:


PSIDf = PSI.drop(['Symbol', 'Description', 'Examined-Junction', 'Background-Major-Junction', 'AltExons', 'ProteinPredictions', 'dPSI', 'ClusterID', 'Coordinates', 'EventAnnotation'], axis = 1)


# In[4]:


print PSIDf.shape
print PSIDf.head()


# Read and process the training dataset

# In[5]:


metaDataAnalysis = pd.read_csv("/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana-New-Analysis/SF/MutationAll/Events-dPSI_0.1_adjp/PSI.SF3B1_Mut_vs_Others.txt", sep='\t')
metaDataAnalysisDf = metaDataAnalysis.drop(['Inclusion-Junction', 'Event-Direction', 'ClusterID', "UpdatedClusterID", 'AltExons', 'EventAnnotation', 'Coordinates', 'ProteinPredictions', 'dPSI', 'rawp', 'adjp'], axis = 1)
print metaDataAnalysisDf.head()


# In[6]:


#Filter for only the common features
dfMerged = pd.merge(PSIDf, metaDataAnalysisDf, how = "inner", on = "UID")
print dfMerged.shape
print PSIDf.shape
print metaDataAnalysisDf.shape
FeatureList = dfMerged['UID'].tolist()
metaDataAnalysisDf_F = metaDataAnalysisDf[metaDataAnalysisDf['UID'].isin(FeatureList)] #filter training dataset
print metaDataAnalysisDf_F.shape

PSIDf_F = PSIDf[PSIDf['UID'].isin(FeatureList)] #filter test dataset
print PSIDf_F.shape


# In[15]:


print metaDataAnalysisDf_F.head()


# In[7]:


training = metaDataAnalysisDf_F.transpose()
training.reset_index(inplace=True)
training.columns=training.iloc[0]
trainingDf = training[1:]
trainingDf.iloc[0, trainingDf.columns.get_loc('UID')] = 1
trainingDf.iloc[1, trainingDf.columns.get_loc('UID')] = 0

print trainingDf.head()


# In[16]:


print trainingDf.shape


# In[19]:


X_train = trainingDf.drop("UID", axis=1)
y_train = trainingDf['UID']
y_train=y_train.astype('int')
X_train = X_train.values
y_train = y_train.values
print y_train


# In[20]:


print X_train
print X_train.shape


# Process the test dataset. Replace all samples with 0 (Others/no SF3B1 mutation) or 1 (with SF3B1 mutation).

# In[9]:


# read the mutation file
Mut = pd.read_csv("/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana-New-Analysis/SF/TrueSplice/All-Combined-Mutation-cBIO-BROAD-Illumina-Removed-Duplicates.txt", sep = "\t")
print Mut.head()
MutList = Mut.loc[Mut['Gene'] == 'SF3B1', 'Sample'].tolist()
print MutList


# In[10]:


###transpose the PSIDf_F
PSIDf_FT = PSIDf_F.transpose()
PSIDf_FT.reset_index(inplace=True)
PSIDf_FT.columns=PSIDf_FT.iloc[0]
PSIDf_FT = PSIDf_FT[1:]
print PSIDf_FT.head()
print PSIDf_FT.shape


# In[11]:


###Replace UID/ samples with 0 or 1
PSIDf_FT.loc[~PSIDf_FT['UID'].isin(MutList), "UID"] = 0
PSIDf_FT.loc[PSIDf_FT['UID'].isin(MutList), "UID"] = 1
df = PSIDf_FT.fillna(PSIDf_FT.mean())
print df.head(10)


# In[23]:


print df.shape


# In[24]:


X_test = df.drop("UID", axis=1)
y_test = df['UID']
y_test=y_test.astype('int')
X_test = X_test.values
y_test = y_test.values
print X_test.shape
print y_test.shape


# In[30]:


import pickle
with open('/Users/bha1ov/Documents/BRCAPaper1/X_train.p','wb') as f:
    pickle.dump(X_train,f)
with open('/Users/bha1ov/Documents/BRCAPaper1/y_train.p','wb') as f:
    pickle.dump(y_train,f)
with open('/Users/bha1ov/Documents/BRCAPaper1/y_test.p','wb') as f:
    pickle.dump(y_test,f)
with open('/Users/bha1ov/Documents/BRCAPaper1/X_test.p','wb') as f:
    pickle.dump(X_test,f)


# In[34]:


svc = LinearSVC(loss="hinge", C=0.001)
svc.fit(X_train, y_train)
y_dec = svc.decision_function(X_test)
y_pred = svc.predict(X_test)
print y_pred
print y_dec


# Altanalyze code
# regr = LinearSVC()
#     regr.fit(Xobs,X[:,0])
#     q=regr.predict(Y)

# from sklearn.svm import SVC
# svclassifier = SVC(kernel='linear')
# svclassifier.fit(X_train, y_train)
# y_pred = svclassifier.predict(X_test)

# In[35]:


from sklearn.metrics import classification_report, confusion_matrix
print(confusion_matrix(y_test,y_pred))
print(classification_report(y_test,y_pred))


# In[42]:


svc1 = LinearSVC(penalty = 'l1', dual=False, C=0.000001)
svc1.fit(X_train, y_train)
y_dec1 = svc1.decision_function(X_test)
y_pred1 = svc1.predict(X_test)
print y_pred1
print y_dec1
print(confusion_matrix(y_test,y_pred1))
print(classification_report(y_test,y_pred1))


# In[43]:


# using point biserial correlation
data = np.concatenate([X_test,y_test.reshape(-1,1)],axis=1)
import scipy.stats as sc
corr = []
for j in range(X_test.shape[1]):
    corr.append(sc.pointbiserialr(y_test,X_test[:,j])[0])
corr = np.array(corr).reshape(-1,1)
print corr


# In[44]:


# cosian similarity
from sklearn.metrics.pairwise import cosine_similarity

def cosine(sample,centroid1,centroid2):
    to1 = cosine_similarity(sample,centroid1)
    to2 = cosine_similarity(sample,centroid2)
    print(to1,to2)
    if to1 < to2:
        return 1
    else:
        return 0

result = []
for i in range(X_test.shape[0]):
    label = cosine(X_test[i,:].reshape(1,-1),X_train[1].reshape(1,-1),X_train[0].reshape(1,-1))
    result.append(label)
result = np.array(result).reshape(-1,1)

print result


# In[47]:


from sklearn.cluster import AgglomerativeClustering
model = AgglomerativeClustering(n_clusters = 2).fit(X_test)
a = model.labels_
from sklearn.metrics import accuracy_score, precision_score, recall_score
print accuracy_score(y_test, a)
print confusion_matrix(y_test, a)


# In[49]:


import seaborn as sns
data = pd.DataFrame(data = X_test.T, index = np.arange(X_test.shape[1]), columns = np.arange(X_test.shape[0]))
col_colorbar = pd.Series(y_test, name = 'mutation').map({0: "orange", 1:"blue"})
sns.clustermap(data=data, method='ward', col_colors=col_colorbar, cmap='viridis')


# D1 = numpy.ma.corrcoef(x) # ma indicates missing value + correlation
# 
# coefr=numpy.ma.corrcoef(values1,values2)
