# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 21:26:54 2021

@author: Maria Ines
"""
import sklearn
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import make_scorer
from sklearn.svm import SVC
from sklearn import metrics
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_predict
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import LabelEncoder
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
import scipy
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.0"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.1.0\bin"
os.environ["R_USER"]= r"C:\Users\Maria Ines\Anaconda3\Lib\site-packages\rpy2"
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri # install any dependency package if you get error like "module not found"
pandas2ri.activate()
base=rpackages.importr("base")
utils= rpackages.importr("utils")
utils.chooseCRANmirror(ind=1)
utils.install_packages("readxl")
readxl = rpackages.importr("readxl")
utils.install_packages("openxlsx")
openxlsx = rpackages.importr("openxlsx")
# import rpy2.robjects.lib.readxl as readxl
#The data should be already prepared and normalized!

class Read_data:
    '''
    The omics data should have the features in the rows and samples as columns,
    the following parameters should be provided:
        omics_files: a dictionary with the name (type of omic, ex. metadata) and file (the path to the dataset: for now only format .txt, .xlsx and .csv are available)
        number_omics: the number of omics used (at least two)
        skip_row: the number of lines the function should skip to start reading the xlsx, also a dictionary.
        feature_name_col_index: a dictionary that contains the name of the omic and the index of the column that contains the features names.
    
    '''
    
    def __init__(self, omics_files, feature_name_col_index,number_omics=2, skip_row=1, header=0, dec=","): 
        
        self.number_omics=number_omics
        self.omics_files= omics_files
        self.skip_row=skip_row
        self.feature_name_col_index= feature_name_col_index
        self.header= header
        self.dec=dec

    
    def read_data(self):
        data={}
        for name,file in self.omics_files.items():
            # print(name, ":")
            if ".xlsx" in file:
                # robjects.r(''' read <- function(file,skip){
                #      file=read.xlsx(file, startRow= skip, rowNames=TRUE)
                #      return (file)}''')
                # read= robjects.r["read"]
                # df_r= read(file, self.skip_row) 
                
                
                feature_col= self.feature_name_col_index[name]
                df=pd.read_excel(file, header=self.header, index_col=feature_col, skiprows=self.skip_row)
                df_r= robjects.conversion.py2rpy(df)
                # print(df.head(6))
                print(utils.head(df_r))
                data[name]=df_r
            if ".txt" in file:
                feature_col= self.feature_name_col_index[name]
                df= pd.read_table(file, header=self.header, sep="\t", index_col=feature_col)
                # df.to_csv('Trasncriptomics_CS1.csv', index=True)
                df_r= robjects.conversion.py2rpy(df)
                # print(df.head(6))
                print(utils.head(df_r))

                data[name]=df_r
            if ".csv" in file:
                feature_col= self.feature_name_col_index[name]
                df= pd.read_csv(file,header=0, index_col=feature_col, skiprows=self.skip_row, decimal=self.dec)
                df_r= robjects.conversion.py2rpy(df)
                # print(df.head(6))
                print(utils.head(df_r))

                data[name]=df_r
        # print(utils.head(data["Metadata"]))
        return data 



# transcriptomics_r=robjects.conversion.py2rpy(transcriptomics_py)



























        
    # def read_files(self): #We obtain 2 dataframes for each omics (one in python and another in r)
    #     # metadata= pd.read_excel(self.metadata, index_col=0, header=1) #Metadata 
        
    #     for name, file in self.omics_files.items():
    #         data_r={}
    #         robjects.r(''' read <- function(file,skip){
    #          # file=read_excel(file,skip= skip)
    #          # return (file)}''')
             
    #         # robjects.r('''rownames  <-  function(df,col_index){
    #         #     # row.names(df) = df[,col_index]
    #         #     # return (df)}''')
                
    #         read= robjects.r["read"]
    #         # rownames= robjects.r["rownames"]
                
    #         if ".xlsx" in file:
    #             # print(name, ".xlsx file")
    #           # #    robjects.r(''' read <- function(file,skip){
    #           #  file=read_excel(file,skip= skip)
    #           #  return (file)}''')
    #             # read= robjects.r["read"]
    #             df_r= read(file, self.skip_row) 
               
    #             print(name, ":")
    #             print("Dimension of", name, ": ", base.dim(df_r))
    #             print(utils.head(df_r))
    #         data_r[name]=df_r
                
                
    #             # if name.lower() == "metadata":
                   
                    
    #             #     robjects.r('''rownames  <-  function(df,col_index){
    #             # row.names(df) = df[,col_index]
    #             # # return (df)}''')   
    #             #     rownames= robjects.r["rownames"]
    #             #   df_r=rownames(df_r,1)
    #             #   print(utils.head(base.row_names(df_r)))
    #             #   data_r[name]=df_r
    #             #   continue
    #             # data_r[name]=df_r
                
                 
    #         elif ".txt" in file:
    #         # else:
                
    #             df_r= utils.read_table(file,header = True, sep = "\t", dec = ".", row_names=1)
    #             # transcriptomics_py= pd.read_csv(self.omics_files.get(omics), sep="\t")
    #             print(name, ":")
    #             print("Dimension of", name, ": ", base.dim(df_r))
    #             # print((transcriptomics_py.head(6)))
    #             print(utils.head(df_r))
    #             # transcriptomics_r=robjects.conversion.py2rpy(transcriptomics_py)
                
    #             # print("Dimension of Transcriptomics Py: ", transcriptomics_py.shape)
    #             # data["transcriptomics"]=transcriptomics_py
    #             data_r[name]=df_r
    #             continue
    #         return data_r
                
                
                
                
                
                
                
            # data_r={}
            # robjects.r(''' read <- function(metadata){
            #      metadata=read_excel(metadata,skip=1)
            #      return (metadata)}''')
                                       
            # read= robjects.r["read"]
            # metadata_r=read(self.metadata)
            # metadata_py=robjects.conversion.rpy2py(metadata_r)
            # print("Dimension of Metadata R: ", base.dim(metadata_r))
            # print("Dimension of Metadata Py: ",metadata_py.shape)
            # print("Metadata: ")
            # print(utils.head(metadata_r))
            # names= metadata_py.loc[:,"Sample ID"]
            # robjects.r('''rownames  <-  function(df,col_index){
            #     row.names(df) = df[,col_index]
            #     return (df)}''')
            # rownames= robjects.r["rownames"]
            # metadata_r=rownames(metadata_r,1)
            # print(utils.head(base.row_names(metadata_r)))
            # print((metadata_py.head(6)))
            # data["metadata"]=metadata_py
            # data_r["metadata"]=metadata_r
        #     if omics == "Transcriptomics":
        #         transcriptomics_r= utils.read_table(self.omics_files.get(omics),header = True, sep = "\t", dec = ".", row_names=1)
        #         # transcriptomics_py= pd.read_csv(self.omics_files.get(omics), sep="\t")
        #         print("Transcriptomics: ")
        #         # print((transcriptomics_py.head(6)))
        #         print(utils.head(transcriptomics_r))
        #         # transcriptomics_r=robjects.conversion.py2rpy(transcriptomics_py)
        #         print("Dimension of Transcriptomics R: ", base.dim(transcriptomics_r))
        #         # print("Dimension of Transcriptomics Py: ", transcriptomics_py.shape)
        #         # data["transcriptomics"]=transcriptomics_py
        #         data_r["transcriptomics"]=transcriptomics_r
        #     if omics == "Metabolomics":
        #         metabolomics_r=read(self.omics_files.get(omics))
        #         # metabolomics_py=robjects.conversion.rpy2py(metabolomics_r)
        #         print("Metabolomics: ")
        #         print(utils.head(metadata_r))
        #         # print((metabolomics_py.head(6))) 
        #         print("Dimension of Metabolomics R: ", base.dim(metabolomics_r))
        #         # print("Dimension of Metabolomics Py: ",metabolomics_py.shape)
        #         # data["metabolomics"]=metabolomics_py
        #         data_r["metabolomics"]=metabolomics_r
        # return data_r
        
        


    