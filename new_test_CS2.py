# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 13:31:49 2021

@author: Maria Ines
"""

from Read_data import Read_data
from Preprocessing import Preprocessing
from Analysis import Analysis
from Annotation import Annotation
from Machine_Learning import Machine_Learning
from Multiomics_integration import Multiomics_integration
from Unsupervised_Learning import Unsupervised_Learning

class test:
    
    '''Read Data'''
    
    '''Case Study 2- Arabidopsis thaliana '''
    
    omics_files= {"Metadata":"C:/Users/Maria Ines/Desktop/Case Studies/Case_Study2/metadata_GSE65046.xlsx",
                  "Transcriptomics": "C:/Users/Maria Ines/Desktop/Case Studies/Case_Study2/transcriptomics_GSE65046.csv",
                  "Fluxomics":"C:/Users/Maria Ines/Desktop/Case Studies/Case_Study2/fluxomics.xlsx"}
    
    feature_name_col_index={"Metadata":0,"Transcriptomics":0, "Fluxomics":0}
    
    skip_rows={"Metadata":0, "Transcriptomics": 0, "Fluxomics": 0}
    
    scaled={"Transcriptomics": True, "Fluxomics": False}
    
    read_data= Read_data(omics_files,feature_name_col_index,number_omics=2, skip_row=0, header=0, dec=".")
    
    data=read_data.read_data()
    
    '''Preprocessing:  1) Missing Values
                       2) Filter
                       3) Normalization 
                       4) Data Discretization 
    
    '''
  
    
    preprocessing= Preprocessing(data, variance_threshold = 0.3)
    preprocessing.na_delete()
    data_filter= preprocessing.filtering()
    data_norm=preprocessing.normalization(scaled)
    
    # # ''' Exploratory Analysis '''
    items={"barplot": ["treatment"]}
    analysis= Analysis(data_norm,items)
    # exploratory= analysis.exploratory()
    # # # summary=analysis.summary_data()
    # pheatmap =analysis.pheatmap(which_data= "ALL", columns=["treatment"]) #which_data - name of omics data or ALL 
    # pca=analysis.pca(which_data="ALL", variable= {"shape": None,"color":"treatment"}) #variable= {"shape":"Variety","color":"Vintage"})#, pcas=["PC1","PC2"]) #which data - same as pheatmap function
    # # # variable = only one variable or a dictionary of different variables and what you want, eg. "Variety" or {"shape":"Variety" ,"color":"Vintage"} 
    # # #                                                         # pcas don't matter, this function can only do PC1 vs PC2                                                            

    annotation= Annotation(default=False, annotation_file="C:/Users/Maria Ines/Desktop/Case Studies/Case_Study2/Annotation-Transcriptomics.xlsx")
    annotation=annotation.read_annotation(header=0, feature_col=0 , skip_rows=0) 
    # # # print(annotation) 
    top= analysis.differential_expression(annotation,which_data="Transcriptomics", y_pred="treatment", filter_results= 407, default="CS2") #annotation: the directory for the annotation file.
    # # #y_pred= name of variable we want to predict from metadata.
    #filter_results= 212, for example, tells that the data will be filtered in 212 transcripts that are more diferentially expressed.
    
    # ''' Machine Learning '''
    
    # ''' Individual Omic Analysis '''
    
    Machine_Learning = Machine_Learning(top, data_norm, y_pred="treatment")
    train_test_datasets = Machine_Learning.train_models(which_data="ALL")
    # svm_res= Machine_Learning.SVM(train_test_datasets, which_data="ALL")
    
    # rf_res=Machine_Learning.RF(train_test_datasets, which_data="ALL")
    # nn_res=Machine_Learning.ANN(train_test_datasets, which_data="ALL")
       
    
    ''' Multiomics Integration '''
    
    ''' omics = (omics1, omics2)'''
    res_concat, res_multi, omic_1, omic_2, Y, dataset_concat= Machine_Learning.save_data(train_test_datasets, omics=("Transcriptomics","Fluxomics")) #In multiomics integration there must be at least 2 diferent omics, therefore this should be a tuple
    # print(res_concat)   # res_multi= list(train_omic1,test_omic1, train_omic2, test_omic2, Y_train, Y_test)
    
    ''' a) Concatenation - Based Integration '''
    
    MI= Multiomics_integration(res_concat, res_multi, omic_1, omic_2, Y, train_test_datasets, y_pred="treatment", omics=("Transcriptomics","Fluxomics"))                 
    # diablo_res= MI.DIABLO()
    # smspl_res=MI.SMSPL()
    # stack_gen_res= MI.Stack_Generalisation()
    # lasso_res=MI.Lasso_regression()
    # svm_res=MI.SVM()
    # nn_res= MI.NNs()
    # rfs_res= MI.RF()
    
    
    ''' b) Transformation - Based Integration '''
    snftool_res= MI.SNFtool()
    # CAN_res= MI.CAN_TBI()
    # RVM_Ada_res= MI.RVM_Ada_TBI ()
    
    
    ''' c) Model - Based Integration '''
    # ensemble_class_res= MI.Ensemble_classifier(option= "op1", voting="ALL") #voting= "Soft", "Hard", "All"
    
    # Option  =   op1  <- Ensemble Classifier with 2 models SVM
                # op2  <- Ensemble Classifier with 2 different models (Decision Tree and Guassian Naive Bayes)
                # op3  <- Ensemble Classifier with 2 Neural Networks
                # op4  <- Ensemble Classifier with ensemble of 2 Naive Bayes (combination of models using voting classifier with NB as recommended for a MBI approach)
                # ALL  <- Executes all the above
                
                
        
    ''' Unsupervised Learning - Multiomics Integration'''
    UnL= Unsupervised_Learning(train_test_datasets, omics=("Transcriptomics", "Metabolomics"), y_pred="Berry")
    
    ''' Concatenation-Based'''
    # MFA_res= UnL.MFA(type_=["c","c","n"]) #multiple factor analysis (MFA)
    
    ''' Transformation - Based '''
    # NEMO_res= UnL.NEMO( num_clusters= 2, num_neighbors=2, k=20) #NEighborhood based Multi-Omics clustering)
    
    ''' Model - Based '''
    # BCC_res=UnL.BCC(K=2, a=1,b=1,IndivAlpha = "FALSE",Concentration=1,NumDraws = 1000)
    
    
    