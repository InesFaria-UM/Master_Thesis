# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 13:53:55 2021

@author: Maria Ines
"""
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
import numpy as np
from rpy2.robjects import pandas2ri # install any dependency package if you get error like "module not found"
pandas2ri.activate()
base=rpackages.importr("base")
utils= rpackages.importr("utils")
utils.chooseCRANmirror(ind=1)

utils.install_packages("dplyr")
dplyr = rpackages.importr("dplyr")

utils.install_packages("simfinapi")
simfinapi=rpackages.importr("simfinapi")


# utils.remove_packages("Rcpp")
# utils.install_packages("C:/Users/Maria Ines/Desktop/Artigos/Tese/Rcpp_1.0.7.zip", type = "source")
# Rcpp=rpackages.importr("Rcpp")
utils.update_packages(ask= False)
utils.install_packages("writexl")
writexl=rpackages.importr("writexl")

utils.install_packages("caret")
caret=rpackages.importr("caret")

class Machine_Learning:
    
    def __init__(self, top,data_norm_final, y_pred):
        '''top is a dictionary that contains the results from the differential expression for a specific omic, if wanted'''
        self.top=top
        self.data_norm_final=data_norm_final
        self.y_pred=y_pred
        self.res_concat=None
        self.res_multi=None
        self.y= None
        
    
    def train_models(self, which_data="ALL"):
        robjects.r('''create_dataset <- function(data,y,y_pred, top_mv){
                   
                      dataset= cbind(t(data[top_mv,]), y) #ou top_mv_50

                      final=length(top_mv)+1
                      # print(head(dataset))
                      colnames(dataset)[final]= y_pred
                      print(head(dataset[,y_pred]))
                      print(dim(dataset))
                      ind=sample(2, nrow(dataset), replace = TRUE, prob=c(0.7,0.3))
                      trainData= dataset[ind==1,]
                      testData=  dataset[ind==2,]
                      print("Train Data")
                      print(dim(trainData))
                      print(table(trainData[,y_pred]))
                      print("Test Data")
                      print( dim(testData))
                      print(table(testData[,y_pred])) 
                      dados=list(dataset, trainData, testData)
                      return (dados) }''')
            
        robjects.r('''create_dataset_2 <- function(data,y,y_pred){
               
                  dataset= cbind(t(data), y) #ou top_mv_50

                  final=length(dataset)
                  # print(head(dataset))
                  colnames(dataset)[final]= y_pred
                  print(head(dataset[,y_pred]))
                  print(dim(dataset))
                  ind=sample(2, nrow(dataset), replace = TRUE, prob=c(0.7,0.3))
                  trainData= dataset[ind==1,]
                  testData=  dataset[ind==2,]
                  print("Train Data")
                  print(dim(trainData))
                  print(table(trainData[,y_pred]))
                  print("Test Data")
                  print( dim(testData))
                  print(table(testData[,y_pred])) 
                  dados=list(dataset, trainData, testData)
                  return (dados)}''')
            
            
            
            
       
        # print(utils.head(self.data_norm_final["Metabolomics"]))            
        # data=robjects.DataFrame(self.data_norm_final["Metabolomics"])
        # print(utils.head(data))
        train_test_datasets={}
        if which_data == "ALL":
            for name,data in self.data_norm_final.items():
                # print(name)
                if name.lower() =="metadata":
                    print("Executing y_ALL")
                    data=robjects.DataFrame(data)
                    y_2= data.rx(self.y_pred)
                    # print(utils.str(y_2)) 
                    for name,data in self.data_norm_final.items():
                          if name.lower()!= "metadata":
                             create_dataset=robjects.r["create_dataset"]
                             create_dataset_2=robjects.r["create_dataset_2"]
                             if self.top[name] is not None:
                                dados=create_dataset(data,y_2,self.y_pred, self.top[name])
                                print(base.dim(dados[0]))
                                # print(utils.head(dados[0]))
                                dados_py= robjects.conversion.rpy2py(dados[0])
                                dados_py.to_excel(name+".xlsx", index=True)
                                # dados_=np.array(dados) 
                                # print(dados_)
                                train_test_datasets[name]=dados
                                
                                
                             if self.top[name] is None:
                                dados=create_dataset_2(data,y_2,self.y_pred)
                                print(base.dim(dados[0]))
                                dados_py= robjects.conversion.rpy2py(dados[0])
                                dados_py.to_excel(name+".xlsx", index=True)
                                # dados_=np.array(dados) 
                                # print(dados_)
                                train_test_datasets[name]=dados
                      # dataset= base.cbind(t(data[self.dif_exp_top,]),y)
                                
                              
                                
                              
                                
        if which_data != "ALL":
            for name,data in self.data_norm_final.items():
                # print(name)
                if name.lower() =="metadata":
                    print("Executing y_ not ALL")
                    data=robjects.DataFrame(data)
                    y_2= data.rx(self.y_pred)
                     
                    for name,data in self.data_norm_final.items():
                          if name.lower()== which_data.lower():
                             create_dataset=robjects.r["create_dataset"]
                             create_dataset_2=robjects.r["create_dataset_2"]
                             top_res=self.top[name]
                             print(top_res)
                             if top_res is not None:
                                dados=create_dataset(data,y_2,self.y_pred, self.top[name])
                                dados_py= robjects.conversion.rpy2py(dados[0])
                                dados_py.to_excel(name+".xlsx", index=True)
                                # dados_=np.array(dados) 
                                # print(dados_)
                                train_test_datasets[name]=dados
                                
                             if self.top[name] is None:
                                dados=create_dataset_2(data,y_2,self.y_pred)
                                dados_py= robjects.conversion.rpy2py(dados[0])
                                dados_py.to_excel(name+".xlsx", index=True)
                                # dados_=np.array(dados) 
                                # print(dados_)
                                train_test_datasets[name]=dados
                      # dataset= base.cbind(t(data[self.dif_exp_top,]),y)
        
        return train_test_datasets
    
    def SVM (self,train_test_datasets, which_data):

         robjects.r('''svm_res <- function(dataset, trainData, testData, y_pred){
                    install.packages("e1071")
                    library(e1071)
                    install.packages("ROCR")
                    library(ROCR)
                    print("Number of NAs:")
                    print(sum(sapply(sapply(dataset,is.na),sum))) 
        
                    factors <-factor(trainData[,y_pred])
                    
                    modelsvm= svm(factors~., trainData)
                    svm_pred=predict(modelsvm, testData)
                    # print(svm_pred)
                    print(table(svm_pred, testData[,y_pred]))
                    a=table(svm_pred, testData[,y_pred])

                    library(caret) 
                    a=confusionMatrix(a)
                    recall=a$byClass["Recall"]
                    precision=a$byClass["Precision"]
                    # rmse=function (obs,pred) sqrt(mean((obs-pred)^2))
                    # mad=function(obs,pred) mean(abs(obs-pred))
        
                    print("PECC:")
                    print(sum(svm_pred==testData[,y_pred])/length(testData[,y_pred])) #(TP + TN) / (TP + TN + FP + FN)
                    b=sum(svm_pred==testData[,y_pred])/length(testData[,y_pred])

                    # print("Recall:")
                    # c=rmse(as.numeric(factor(testData[,y_pred])), as.numeric(factor(svm_pred)))
                    
                    res= list(a,b,recall,precision)
                    return (res)}''') 
             
             
             
             
             
         robjects.r('''svm_caret <- function(dataset, trainData, testData, y_pred){
                    install.packages("e1071")
                    library(e1071)
                    library(pROC)
                    library(caret)
                    
                    install.packages("MLeval")
                    library(MLeval)
                    
                    set.seed(123451)
                    print("Number of NAs:")
                    print(sum(sapply(sapply(dataset,is.na),sum))) 
                    print(sum(sapply(sapply(trainData,is.na),sum))) 
                    print(sum(sapply(sapply(testData,is.na),sum))) 
                    
                    factors <-factor(trainData[,"treatment"])
                    # print(factors)
                    testData[,y_pred]=factor(testData[,y_pred])
                    
                    index= which(colnames(trainData) == y_pred)
                    print(index)
                    # print(head(testData))
                    # print(testData[,y_pred])
                    trainData$treatment=as.factor(trainData$treatment)
                    
                    # Set up Repeated k-fold Cross Validation
                    train_control <- trainControl(method="repeatedcv", number=10, repeats=3, 
                     classProbs=TRUE,
                     savePredictions = TRUE)
                   
                    
                   # Fit the model 
                    # model.svm <- train(trainData[,!(names(trainData) %in% y_pred)],factors, method = "svmLinear", trControl = train_control,  metric="Accuracy",
                    # tuneLength=15, preProcess = c("center","scale"),tuneGrid = expand.grid(C = seq(0, 2, length = 20)))
                    
                    # model.svm <- train(trainData[,-index],trainData[,index], method = "svmLinear", trControl = train_control,  metric="Accuracy",
                    # tuneLength=15, preProcess = c("center","scale"),tuneGrid = expand.grid(C = seq(0, 2, length = 20)))
                    
                    
                    
                    model.svm <- train(treatment~.,data=trainData, method = "svmLinear",trControl = train_control,  metric="Accuracy",
                                        tuneLength=15,tuneGrid = expand.grid(C = seq(0, 2, length = 4)))
                    print(model.svm)
                    svm_pred=predict(model.svm, testData)
                    
                    
                    
                                        
                                        
                    ## run MLeval
                    
                    res <- evalm(model.svm)
                    
                    ## get ROC

                    roc=res$roc

                    # coefs <- model.svm$finalModel@coef[[1]]
                    # mat <- model.svm$finalModel@xmatrix[[1]]
                    # print(coefs %*% mat)
                    # fs= coefs %*% mat
                    
                    varImp= varImp(model.svm, scale = FALSE)
                    
                    
                    # png(file="plot_svm_VarImp.png", width=512, height=512)
                    print(plot(varImp))
                    # dev.off()
                    
                    # png(file="svm_roc_curve.png", width=512, height=512)
                    # print(roc(response=testData[,y_pred], predictor = svm_pred,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE))
                    
                    print("Here!")
                    print(svm_pred)
                    # print(table(svm_pred, testData[,y_pred]))
                    # a=table(svm_pred, testData[,y_pred])

                    library(caret) 
                    a=confusionMatrix(svm_pred,testData[,y_pred])
                    print(a)
                    recall=a$byClass["Recall"]
                    precision=a$byClass["Precision"]
                    # rmse=function (obs,pred) sqrt(mean((obs-pred)^2))
                    # mad=function(obs,pred) mean(abs(obs-pred))
        
                    print("PECC:")
                    print(sum(svm_pred==testData[,y_pred])/length(testData[,y_pred])) #(TP + TN) / (TP + TN + FP + FN)
                    b=sum(svm_pred==testData[,y_pred])/length(testData[,y_pred])

                    # print("Recall:")
                    # c=rmse(as.numeric(factor(testData[,y_pred])), as.numeric(factor(svm_pred)))
                    
                    res= list(a,b,recall,precision,roc, varImp)
                    return (res)}''')
            

             

             
             
             
             
         if which_data=="ALL":
            for name, data in self.data_norm_final.items():
                if name.lower() != "metadata":
                   print(name)
                   lista=train_test_datasets[name] 
                   file= open("SVM_"+ name+"_.txt", "w")
                   dataset=lista[0]
                   print("Dimension Dataset:")
                   print(base.dim(dataset))
                   trainData=lista[1]
                   print("Dimension of TrainData:")
                   print(base.dim(trainData))
                   testData=lista[2]  
                   print("Dimension of TestData:")
                   print(base.dim(testData))
                   svm_res=robjects.r["svm_caret"]
                   res_svm= svm_res(dataset, trainData, testData, self.y_pred)
                   file.write(f'Dimension Dataset: {base.dim(dataset)} \n')
                   file.write(f'Dimension Train Data: {base.dim(trainData)} \n')
                   file.write(f'Dimension Test Data: {base.dim(testData)} \n')
                   file.write(f'Confusion Matrix: \n {res_svm[0]} \n')
                   file.write(f'PECC: {res_svm[1]} \n')
                   file.write(f'Recall: {res_svm[2]} \n')
                   file.write(f'Precision: {res_svm[3]} \n')
                   file.write(f'Relevant features: \n {res_svm[5]} \n')
                   file.write(f'See figure svm_roc_{name}.png for the roc curve.')
                   file.close()
                   grdevices = importr('grDevices')
    
                   grdevices.png(file="svm_roc_"+name+".png", width=512, height=512)
                   print(res_svm[4])
                   grdevices.dev_off()
                   # file.write(f'RMSE: {res_svm[2]} \n')
                   
                
         if which_data != "ALL": 
            for name, data in self.data_norm_final.items():
                if name.lower() == which_data.lower():
                   print(name)
                   lista=train_test_datasets[name]        
                   dataset=lista[0]
                   print("Dimension Dataset:")
                   print(base.dim(dataset))
                   trainData=lista[1]
                   print("Dimension of TrainData:")
                   print(base.dim(trainData))
                   testData=lista[2]  
                   print("Dimension of TestData:")
                   print(base.dim(testData))
                   svm_res=robjects.r["svm_caret"]
                   res_svm= svm_res(dataset, trainData, testData, self.y_pred)
                   file= open("SVM_"+ name+"_.txt", "w")
                   file.write(f'Dimension Dataset: {base.dim(dataset)} \n')
                   file.write(f'Dimension Train Data: {base.dim(trainData)} \n')
                   file.write(f'Dimension Test Data: {base.dim(testData)} \n')
                   file.write(f'Confusion Matrix: {res_svm[0]} \n')
                   file.write(f'PECC: {res_svm[1]} \n')
                   file.write(f'Recall: {res_svm[2]} \n')
                   file.write(f'Precision: {res_svm[3]} \n')
                   file.write(f'Relevant features: \n {res_svm[5]} \n')
                   file.close()
                   grdevices = importr('grDevices')
    
                   grdevices.png(file="svm_roc_"+name+".png", width=512, height=512)
                   print(res_svm[4])
                   grdevices.dev_off()
                   
                   grdevices.png(file="svm_varImp_"+name+".png", width=512, height=512)
                   robjects.r(f'plot({res_svm[5]}')
                   grdevices.dev_off()
         return          
                   
                   
                   
                   
                   
                   
         
    def RF (self,train_test_datasets, which_data):

        robjects.r('''rf_res <- function(dataset, trainData, testData, y_pred){
                   install.packages("randomForest")
                   library(randomForest)
                   
                   # install.packages("caret")
                   # library(caret)
                   
                   set.seed(12345)
                   micro.rf=randomForest(factor(trainData[,y_pred])~.,trainData,importance=TRUE)
                   pred.rf=predict(micro.rf,testData)
                   pecc=function(obs,pred) sum(obs==pred)/length(pred)
                   print(pecc(pred.rf,testData[,y_pred]))
                   b=pecc(pred.rf,testData[,y_pred])
                   print(table(pred.rf,testData[,y_pred]))
                   a=table(pred.rf,testData[,y_pred])
                   
                   library(caret) 
                   a=confusionMatrix(a)
                   recall=a$byClass["Recall"]
                   precision=a$byClass["Precision"]
                   
                   
                   
                   
                   # print(round(importance(micro.rf),2))
                   c=round(importance(micro.rf),2)
                   model_rf=train(trainData[,1:10], trainData[,11], method="rf")
                   c=model_rf$results
                   c=varImp(model_rf)
                   res=list(a,b,recall, precision, c)
                   return (res)}''') 
            
            
        robjects.r(''' rf_caret <- function(dataset, trainData, testData, y_pred){
           
                install.packages("MLeval")
                    library(MLeval)
                   
                library(caret)
                set.seed(123451)
                trainData$treatment=as.factor(trainData$treatment)
                # factors <-as.factor(trainData[,y_pred])
                
                res.ctrl=trainControl(method="repeatedcv",number=10, repeats=3, search="grid", classProbs=TRUE,savePredictions = TRUE)
                # tunegrid <- expand.grid(.mtry=c(sqrt(ncol(trainData[,!(names(trainData) %in% y_pred)]))))
                
                tunegrid <- expand.grid(.mtry=3)
                model.rf=train(trainData[,!(names(trainData) %in% y_pred)],trainData[,y_pred], method="rf",tuneLength  = 20, trControl=res.ctrl, tuneGrid=tunegrid)
               
                
               ## run MLeval
                    
                    res <- evalm(model.rf)
                    
                    ## get ROC

                    roc=res$roc

                varImp= varImp(model.rf, scale = FALSE)
               
                png(file="plot_rf_VarImp.png", width=512, height=512)
                print(plot(varImp))
                dev.off()
                # rf_gridsearch <- train(trainData[,!(names(trainData) %in% y_pred)],trainData[,y_pred], method="rf", metric="Accuracy", tuneGrid=tunegrid, trControl=res.ctrl, ntree=80)
                # print(rf_gridsearch)
                # png(file="gridsearch_rf.png")            
                
                # print((plot(rf_gridsearch))
                # dev.off()
                
                
                # control <- trainControl(method="repeatedcv", number=3, repeats=5, search="grid")
                # tunegrid <- expand.grid(.mtry=c(sqrt(ncol(trainData[,!(names(trainData) %in% y_pred)]))))
                # modellist <- list()
                # for (ntree in c(20,40,60,80,100)) {
                #  	set.seed(123451)
                #  	fit <- train(trainData[,!(names(trainData) %in% y_pred)],trainData[,y_pred], method="rf", metric="Accuracy", tuneGrid=tunegrid, trControl=control, ntree=ntree)
                #  	key <- toString(ntree)
                #  	modellist[[key]] <- fit
                # }
                # # compare results
                # results <- resamples(modellist)
                # print(summary(results))
                # png(file="dotplot_rf_ntree.png")
                # print(dotplot(results))
                #   dev.off()               
                
                
                print(testData[,y_pred])
                
                
                
                pred.rf= predict(model.rf, testData)
                a=confusionMatrix(pred.rf, factor(testData[,y_pred]))
                print(a)
                recall=a$byClass["Recall"]
                precision=a$byClass["Precision"]
                pecc=function(obs,pred) sum(obs==pred)/length(pred)
                PECC= pecc(pred.rf,factor(testData[,y_pred]))
                res_rf=list(a,PECC,recall,precision, pred.rf,roc,varImp)
            
                
                return(res_rf)}''')
                   
                   
        
        if which_data=="ALL":
           for name, data in self.data_norm_final.items():
               if name.lower() != "metadata":
                  print(name)
                  lista=train_test_datasets[name] 
            
                  
                  dataset=lista[0]
                  print("Dimension Dataset:")
                  print(base.dim(dataset))
                  trainData=lista[1]
                  print("Dimension of TrainData:")
                  print(base.dim(trainData))
                  testData=lista[2]  
                  print("Dimension of TestData:")
                  print(base.dim(testData))
                  rf_res=robjects.r["rf_caret"]
                  res_rf= rf_res(dataset, trainData, testData, self.y_pred)
                  file= open("RF_"+ name+"_.txt", "w")
                  file.write(f'Dimension Dataset: {base.dim(dataset)} \n')
                  file.write(f'Dimension Train Data: {base.dim(trainData)} \n')
                  file.write(f'Dimension Test Data: {base.dim(testData)} \n')
                  file.write(f'Confusion Matrix: {res_rf[0]} \n')
                  file.write(f'PECC: {res_rf[1]} \n')
                  file.write(f'Recall: {res_rf[2]} \n')
                  file.write(f'Precision: {res_rf[3]} \n')
                  file.write(f'Importance: {res_rf[4]} \n')
                  file.write(f'Relevant features: \n {res_rf[6]} \n')
                  file.write(f'See figure rf_roc_{name}.png for the roc curve.')
                  file.close()
                  grdevices = importr('grDevices')
    
                  grdevices.png(file="rf_roc_"+name+".png", width=512, height=512)
                  print(res_rf[5])
                  grdevices.dev_off()
                   
        if which_data != "ALL": 
            for name, data in self.data_norm_final.items():
                if name.lower() == which_data.lower():
                   print(name)
                   lista=train_test_datasets[name]        
                   dataset=lista[0]
                   print("Dimension Dataset:")
                   print(base.dim(dataset))
                   trainData=lista[1]
                   print("Dimension of TrainData:")
                   print(base.dim(trainData))
                   testData=lista[2]  
                   print("Dimension of TestData:")
                   print(base.dim(testData))
                   rf_res=robjects.r["rf_caret"]
                   res_rf= rf_res(dataset, trainData, testData, self.y_pred)
                   file= open("RF_"+ name+"_.txt", "w")
                   file.write(f'Dimension Dataset: {base.dim(dataset)} \n')
                   file.write(f'Dimension Train Data: {base.dim(trainData)} \n')
                   file.write(f'Dimension Test Data: {base.dim(testData)} \n')
                   file.write(f'Confusion Matrix: {res_rf[0]} \n')
                   file.write(f'PECC: {res_rf[1]} \n')
                   file.write(f'Recall: {res_rf[2]} \n')
                   file.write(f'Precision: {res_rf[3]} \n')
                   file.write(f'Importance: {res_rf[4]} \n')
                   file.write(f'Relevant features: \n {res_rf[6]} \n')
                   file.write(f'See figure rf_roc_{name}.png for the roc curve.')
                   file.close()
                   grdevices = importr('grDevices')
    
                   grdevices.png(file="rf_roc_"+name+".png", width=512, height=512)
                   print(res_rf[5])
                   grdevices.dev_off()
            return
         
            
    def ANN (self,train_test_datasets, which_data):

        robjects.r('''nns_res <- function(dataset, trainData, testData, y_pred){
                   install.packages("nnet")
                   library(nnet)
                   set.seed(123451)
                   
                   factors <-factor(trainData[,y_pred])
                   nn=nnet(factors~.,trainData,size=3,MaxNWts=18000)
                    
                   nn_prev=predict(nn,testData,type="class")
                   print(nn_prev)
                   c= nn_prev
                   print(table(nn_prev,testData[,y_pred]))
                   a=table(nn_prev,testData[,y_pred])
                   
                   library(caret) 
                   a=confusionMatrix(table(nn_prev,testData[,y_pred]))
                   recall=a$byClass["Recall"]
                   precision=a$byClass["Precision"]
                   
                   
                   print(sum(nn_prev==testData[,y_pred])/length(testData[,y_pred]))
                   b= sum(nn_prev==testData[,y_pred])/length(testData[,y_pred])
                    # pecc=function(obs,pred) sum(obs==pred)/length(pred)
                    # print(pecc(nn_prev,testData[,y_pred]))
                   res=list(a,b,recall, precision, c)
                   return (res)}''')  
            
            
        robjects.r('''ann_caret_res <- function(dataset, trainData, testData, y_pred){
            library(caret)
            set.seed(123451)
            install.packages("MLeval")
                    library(MLeval)
            
            factors <-as.factor(trainData[,y_pred])
            
            res.ctrl=trainControl(method="repeatedcv",number=10,repeats=3, classProbs = TRUE,summaryFunction = twoClassSummary,
                     savePredictions = TRUE)
            netGrid <-  expand.grid(size = seq(1,100,10),
                        decay = c(0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7))

            model.ann=train(trainData[,!(names(trainData) %in% y_pred)],trainData[,y_pred], method="nnet",metric="Accuracy", tuneLength=20, trControl=res.ctrl,tuneGrid = netGrid)
            
            ## run MLeval
                    
                    res <- evalm(model.ann)
                    
                    ## get ROC

                    roc=res$roc

            varImp= print(varImp(model.ann, scale = FALSE))
            pred.ann= predict(model.ann, testData)
            a=confusionMatrix(pred.ann, factor(testData[,y_pred]))
            print(a)
            recall=a$byClass["Recall"]
            precision=a$byClass["Precision"]
            pecc=function(obs,pred) sum(obs==pred)/length(pred)
            PECC= pecc(pred.ann,testData[,y_pred])
            res_ann=list(a,PECC,recall,precision, pred.ann,roc,varImp)
            
            
            return(res_ann)}''')
        
        if which_data=="ALL":
           for name, data in self.data_norm_final.items():
               if name.lower() != "metadata":
                  print(name)
                  lista=train_test_datasets[name] 
            
                  
                  dataset=lista[0]
                  print("Dimension Dataset:")
                  print(base.dim(dataset))
                  trainData=lista[1]
                  print("Dimension of TrainData:")
                  print(base.dim(trainData))
                  testData=lista[2]  
                  print("Dimension of TestData:")
                  print(base.dim(testData))
                  nns_res=robjects.r["ann_caret_res"]
                  res_nns= nns_res(dataset, trainData, testData, self.y_pred)
                  file= open("NNs_"+ name+"_.txt", "w")
                  file.write(f'Dimension Dataset: {base.dim(dataset)} \n')
                  file.write(f'Dimension Train Data: {base.dim(trainData)} \n')
                  file.write(f'Dimension Test Data: {base.dim(testData)} \n')
                  file.write(f'Confusion Matrix: {res_nns[0]} \n')
                  file.write(f'PECC: {res_nns[1]} \n')
                  file.write(f'Recall: {res_nns[2]} \n')
                  file.write(f'Precision: {res_nns[3]} \n')
                  file.write(f'Predicted: {res_nns[4]} \n')
                  file.write(f'Relevant features: \n {res_nns[6]} \n')
                  file.write(f'See figure ann_roc_{name}.png for the roc curve.')
                  file.close()
                  grdevices = importr('grDevices')
    
                  grdevices.png(file="ann_roc_"+name+".png", width=512, height=512)
                  print(res_nns[5])
                  grdevices.dev_off()
        if which_data != "ALL": 
            for name, data in self.data_norm_final.items():
                if name.lower() == which_data.lower():
                   print(name)
                   lista=train_test_datasets[name]        
                   dataset=lista[0]
                   print("Dimension Dataset:")
                   print(base.dim(dataset))
                   trainData=lista[1]
                   print("Dimension of TrainData:")
                   print(base.dim(trainData))
                   testData=lista[2]  
                   print("Dimension of TestData:")
                   print(base.dim(testData))
                   nns_res=robjects.r["ann_caret_res"]
                   res_nns= nns_res(dataset, trainData, testData, self.y_pred)
                   file= open("NNs_"+ name+"_.txt", "w")
                   file.write(f'Dimension Dataset: {base.dim(dataset)} \n')
                   file.write(f'Dimension Train Data: {base.dim(trainData)} \n')
                   file.write(f'Dimension Test Data: {base.dim(testData)} \n')
                   file.write(f'Confusion Matrix: {res_nns[0]} \n')
                   file.write(f'PECC: {res_nns[1]} \n')
                   file.write(f'Recall: {res_nns[2]} \n')
                   file.write(f'Precision: {res_nns[3]} \n')
                   file.write(f'Predicted: {res_nns[4]} \n')
                   file.write(f'Relevant features: \n {res_nns[6]} \n')
                   file.write(f'See figure ann_roc_{name}.png for the roc curve.')
                   file.close()
                   grdevices = importr('grDevices')
    
                   grdevices.png(file="ann_roc_"+name+".png", width=512, height=512)
                   print(res_nns[5])
                   grdevices.dev_off()
            return
        
 
    
 
    def save_data(self, train_test_datasets, omics):
        robjects.r('''train_test_datasets_concat <- function(omic1, omic2){
             install.packages("caret")
             library(caret)
             
             dataset_concat= cbind(omic1[,1:(length(omic1)-1)],omic2[,1:length(omics2)])
            
            
            } ''')
                   
                   

        
        robjects.r('''dataset_concat <- function(omic1, omic2, y_pred){
                    bol=all(omic1[,length(omic1)] == omic2[,length(omic2)]) 
                    if (bol == TRUE) {
                            install.packages("caret")
                            library(caret)
                            omic1_new=omic1
                            omic1_new[,y_pred]=NULL
                            dataset_concat= cbind(omic1[,1:length(omic1_new)],omic2[,1:(length(omic2))])
                            row.names(dataset_concat)=row.names(omic1)
                            print("Dimension of Dataset Concatenated:")
                            print(dim(dataset_concat)) 
                            dataset=as.data.frame(dataset_concat)
                            
                            write_xlsx(dataset,"dataset_concat.xlsx") 
                            print("The dataset concatenated was saved in your directory!")
                            
                            # omic_1=dataset[,1:(length(omic1)-1)]
                            
                            # start_2 = length(dataset)- length(omic1)
                            # omic_2=dataset_concat[,(start_2+1):(length(dataset)-1)]
                            # Y= dataset[,length(dataset)]
                            # print(head(Y))
                            
                            
                            # inTrain=createDataPartition(y=Y,p=0.7,list=F)
                            # trainDa_multi=dataset[inTrain,]
                            # testDa_multi=dataset[-inTrain,]
                            
                            # train_test_concat= list(trainDa_multi, testDa_multi)
                            
                            #  #Concatenated Datasets
                            # print(dim(trainDa_multi))
                            # print(dim(testDa_multi))
                            
                            
                            
                    } else {
                        print("Y is not the same")
                    }
                
                    return (dataset)}''')
                        
                        
        robjects.r('''train_test_concat <- function(y_pred, dataset_concat){
                        install.packages("caret")            
                        library(caret)
                        Y= dataset_concat[,y_pred]
                        inTrain=createDataPartition(y=Y,p=0.7,list=F)
                        trainDa_concat=dataset_concat[inTrain,]
                        testDa_concat=dataset_concat[-inTrain,]
                        # print("Done!")
                        train_test_concat= list(trainDa_concat, testDa_concat)
                        
               
                    return (train_test_concat)}
                   
                    ''')
                        
        robjects.r('''train_test_multiomics <- function(len_concat,omic1,omic2,trainDa_multi, testDa_multi, y_pred){
            
                        #Datasets for other types of integration
                        # trainDa_multi= train_test_concat[0]
                        # testDa_multi=train_test_concat[1]
                        # print(head(train_test_concat))
                        # print(head(trainDa_multi))
                        # print(dim(testDa_multi))                        
                        train_omic1=trainDa_multi[,1:(length(omic1))]
                        # print(head(train_omic1))
                        start_2 = len_concat- length(omic1)
                        
                        train_omic2=trainDa_multi[,(start_2):(len_concat-1)]
                        # print(head(train_omic2))
                        Y_train= trainDa_multi[,y_pred]
                        
                        test_omic1=testDa_multi[,1:(length(omic1))]
                        test_omic2=testDa_multi[,(start_2):(len_concat-1)]
                        Y_test= testDa_multi[,y_pred]
                       
                        train_test_omics= list(train_omic1,test_omic1, train_omic2, test_omic2, Y_train, Y_test)
                        return (train_test_omics)
                                                        
                            
                            
            
            }''')
            

               
        
        
        omics_=[]
        names=[]
        for name, lista in train_test_datasets.items():
            if name in omics:
                #confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        if len(omics_) ==2:
            omics1= omics_[0]
            omics2= omics_[1]
            dataset_concat=robjects.r["dataset_concat"]
            dataset_concat_= dataset_concat(omics1,omics2, self.y_pred)
            # print(utils.head(dataset_concat))
            dataset_concat_py= robjects.conversion.rpy2py(dataset_concat_)
            dataset_concat_py.to_excel("dataset_concat.xlsx", index=True)   
        
            print("Dimension of "+ names[0] + ":")
            dim1=base.dim(omics1)
            print(dim1)
            range_omic1= np.arange(1, dim1[1], 1)
            dataset_concat=robjects.DataFrame(dataset_concat_)
            omic_1=dataset_concat.rx(True, range_omic1)    
            # print(utils.head(omic_1)) #Perfect!
            
            print("Dimension of "+ names[1]+ ":")
            dim2=base.dim(omics2)
            print(dim2)
            dim_concat=base.dim(dataset_concat)
            len_concat=dim_concat[1]
            final_2= dim_concat[1]
            range_omic2= np.arange(dim1[1],final_2, 1)
            omic_2=dataset_concat.rx(True,range_omic2)
            # print(utils.head(omic_2))
            
            Y= dataset_concat.rx(True, self.y_pred)
            # print(Y)
            tt=robjects.r["train_test_concat"]
            res_concat=tt(self.y_pred, dataset_concat_)
            trainDa_multi= res_concat[0]
            testDa_multi=res_concat[1]
            # print(utils.head(trainDa_multi))
            # print(utils.head(testDa_multi))
            # Y_2= robjects.FactorVector(Y)
            # # print(utils.head(trainDa_concat))
            # ''' Test/Train Concat '''
            # inTrain=caret.createDataPartition(y=Y_2,p=0.7)
            # trainDa_concat=dataset_concat[inTrain,]
            # testDa_concat=dataset_concat[-inTrain,]
            # print(utils.head(trainDa_concat))
            # train_test_concat=robjects.r["train_test_concat"]
            # res_concat= train_test_concat(self.y_pred, dataset_concat)
        
            train_test_multiomics=robjects.r["train_test_multiomics"]
        
            res_multi = train_test_multiomics(len_concat,omic_1,omic_2,trainDa_multi, testDa_multi, self.y_pred)
            self.res_concat = res_concat
            self.res_multi=res_multi
        return res_concat, res_multi, omic_1, omic_2, Y, dataset_concat
        