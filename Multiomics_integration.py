# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 22:52:17 2021

@author: Maria Ines
"""
import warnings 
import sklearn
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import make_scorer
from sklearn.svm import SVC
from sklearn import metrics
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc, plot_roc_curve
from sklearn.model_selection import KFold
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
from lime import lime_text
from lime.lime_text import LimeTextExplainer
from tabulate import tabulate

import scipy
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

utils.install_packages("graphics")
graphics= rpackages.importr("graphics")

utils.install_packages("grDevices")
grdevices = importr('grDevices')
utils.install_packages("dplyr")
dplyr = rpackages.importr("dplyr")

# utils.install_packages("simfinapi")
# simfinapi=rpackages.importr("simfinapi")


# utils.remove_packages("Rcpp")
# utils.install_packages("C:/Users/Maria Ines/Desktop/Artigos/Tese/Rcpp_1.0.7.zip", type = "source")
# Rcpp=rpackages.importr("Rcpp")
utils.update_packages(ask= False)
utils.install_packages("writexl")
writexl=rpackages.importr("writexl")

utils.install_packages("caret")
caret=rpackages.importr("caret")

robjects.r('install.packages("BiocManager", '
    'repos="http://cran.r-project.org")')
robjects.r('BiocManager::install("mixOmics")')
mixOmics= rpackages.importr("mixOmics")
# robjects.r('library(mixOmics)')

utils.install_packages("simEd")
simEd=rpackages.importr("simEd")

class Multiomics_integration:
    
    def __init__(self, res_concat, res_multi, omic_1, omic_2, Y, train_test_datasets, y_pred, omics=("Transcriptomics","Metabolomics")):
        self.res_concat= res_concat
        self.res_multi=res_multi
        self.y_pred=y_pred
        self.omic_1=omic_1
        self.omic_2=omic_2
        self.Y=Y
    #list(train_omic1,test_omic1, train_omic2, test_omic2, Y_train, Y_test)
        self.train_omic1 = res_multi[0]
        self.test_omic1=res_multi[1]
        self.train_omic2=res_multi[2]
        self.test_omic2=res_multi[3]
        self.Y_train=res_multi[4]
        self.Y_test=res_multi[5]
        self.train_test_datasets= train_test_datasets
        self.omics= omics


    '''Multivariate-Based Integration'''
    
    def DIABLO (self):
        
        robjects.r(''' evaluate.DIABLO.performance <- function(confusion.mat, true_lable, predict_label){
          install.packages("ROCR")
          library(ROCR)
          
          tp <- confusion.mat[2,2]
          tn <- confusion.mat[1,1]
          fp <- confusion.mat[2,1]
          fn <- confusion.mat[1,2]
          
          Accuracy <- (tp + tn)/(tp + tn + fp + fn)
          Sensitivity <- tp/(tp + fn)
          Specificity <- tn/(tn + fp)
          Recall = tp/(tp + fp)
          
          predict_label <- as.factor(predict_label)
          true_lable <- as.numeric(true_lable)
          predict_label <- as.numeric(predict_label)
          pred <- prediction(predict_label, true_lable)
          perf <- performance(pred, measure = "tpr", x.measure = "fpr")
          auc <- performance(pred, measure = "auc")
          AUC <- auc@y.values[[1]]
          
          perf <- c(Accuracy, Sensitivity, Specificity, Recall, AUC)
          return(perf)
          
        }''')
        
        robjects.r(''' exec_diablo <- function(train_omic1,train_omic2, Y_train,y_pred, trainDa_multi){
            
            
                data = list(omics1 = train_omic1, omics2= train_omic2)
    
                print(lapply(data, dim))
            
                #outcome
                trainDa_multi=as.data.frame(trainDa_multi)
                Y= factor(trainDa_multi[,y_pred])
                print(Y)
                print(summary(Y))
                
                # Y_diablo = factor(Y_train[,y_pred])
                # print(summary(Y_diablo))
                # Y_diablo= Y_train
            
            pls.res=pls(X=train_omic1,Y=train_omic2,ncomp=3)
            a=cor(pls.res$variates$X,pls.res$variates$Y) 
            print(diag(a))
            
            pca.res=pca(X=train_omic1,ncomp=3) 
            print(diag(a))
            print(tune.pca(train_omic1, ncomp = 10, center = TRUE, scale = FALSE))
            
            spls.res=spls(X=train_omic1,Y=train_omic2,ncomp=3)
            b=cor(spls.res$variates$X,pls.res$variates$Y) 
            print(diag(b))
            
            # Parameter choice
            # Design
            design = matrix(0.84, ncol = length(data), nrow = length(data), 
                            dimnames = list(names(data), names(data)))
            diag(design) = 0
            
            print(design) 
            
            # Tuning the number of components
            
            sgccda.res = block.splsda(X = data, Y =Y, ncomp = 5, 
                                       design = design)
            
            cor(sgccda.res$variates$omics1, sgccda.res$variates$omics2)
            
            set.seed(123) # for reproducibility, only when the `cpus' argument is not used
            # this code takes a couple of min to run
            perf.diablo = perf(sgccda.res, validation = 'loo', folds = 10, nrepeat = 10) #Leave One Out
            
            
            print("Done!")
            png(file="perf_diablo.png", width=512, height=512)
            print(plot(perf.diablo))
            dev.off()
            
            
            perf.diablo$choice.ncomp$WeightedVote
            ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
            print(ncomp)
            
            
            #Tuning keepX
            
            #test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
             #                  proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
            
            #tune.TCGA = tune.block.splsda(X = data, Y = Y_diablo, ncomp = ncomp, 
             #                             test.keepX = test.keepX, design = design,
              #                            validation = 'Mfold', folds = 10, nrepeat = 1,
               #                           cpus = 2, dist = "centroids.dist")
            
            
            #list.keepX = tune.TCGA$choice.keepX
            #list.keepX
            
            print("Start sgccda.res")
            sgccda.res = block.splsda(X = data, Y = Y, ncomp= 2, #está assim para testar (tem de ser ncomp=ncomp)
                                      design = design)
            #sgccda.res   # list the different functions of interest related to that object
            
            print(sgccda.res$design)
            
            print("Done sgccda.res!")
            
            DIABLO_transcriptomics=selectVar(sgccda.res, block = "omics1", comp = 1)$omics1$name
            DIABLO_metabolomics=selectVar(sgccda.res, block ="omics2", comp = 1)$omics2$name
            
            write.table(DIABLO_transcriptomics,"DIABLO_transcriptomics.txt", sep=",",col.names = FALSE, row.names = FALSE)
            write.table(DIABLO_metabolomics,"DIABLO_metabolomics.txt", sep=",",col.names = FALSE, row.names = FALSE)
            
            # png(file="select_var.png", width=512, height=512)
            
            # # selectVar(sgccda.res, block = 'omics1', comp = 1)$omics1$name 
            # selectVar(sgccda.res)$transcriptomics$name 
            # dev.off()
            
            # png(file="plot_diablo.png", width=512, height=512)
            # plotDiablo(sgccda.res, ncomp =1)
            # dev.off()
    
            png(file="plotIndiv.png", width=512, height=512)
            plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO', ellipse = TRUE)
            dev.off()
            
            png(file="plotArrow.png", width=512, height=512)
            plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
            dev.off()
            
            png(file="plotVar.png", width=512, height=512)
            plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, pch = c(16, 17), cex = c(2,2), col = c('darkorchid', 'lightgreen')) 
            dev.off()
            
            object=circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, # não posso correr isto
            color.blocks= c('darkorchid', 'brown1'),
            color.cor = c("chocolate3","grey20"), size.labels = 1.5)
            png(file="circosPlot.png", width=1000, height=1000)
            circosPlot(sgccda.res, cutoff =0.85, line = TRUE, # não posso correr isto
            color.blocks= c('darkorchid', 'brown1'),
            color.cor = c("chocolate3","grey20"), size.variables = 0.6)
            dev.off()
            
            print("Attributes")
            print(head(object))
            
            # write.table(object,"DIABLO_circosplotobject.txt", sep=",",col.names = TRUE, row.names = TRUE)
            write.csv(object,"DIABLO_circosplotobject.csv",col.names = TRUE, row.names = TRUE)
            png(file="network.png", width=1000, height=1000)
            network(sgccda.res, blocks = c(1,2), color.node = c('darkorchid', 'lightgreen'), cutoff = 0.5)
            dev.off()

            png(file="plotLoadings.png",width=600, height=600)
            plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median', ndisplay=20, size.name=1, size.legend=1)
            
            dev.off()
            
            return (sgccda.res)}  ''')
            
            
        robjects.r(''' predict_diablo <- function(sgccda.res, y_pred, test_omic1, test_omic2, testDa_multi){
            
            set.seed(123)# for reproducibility, only when the `cpus' argument is not used
            
            
            perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, 
            dist = 'centroids.dist')
            
            print(perf.diablo)  # lists the different outputs
            
            # Performance with Majority vote
            perf.diablo$MajorityVote.error.rate
            
            
            perf.diablo$WeightedVote.error.rate
            
            auc.splsda = auroc(sgccda.res, roc.block = "omics1", roc.comp = 2)
            print(auc.splsda)
            
            png(file="roc_diablo_final.png")
            print(auc.splsda)
            dev.off()
            
            auc.splsda2 = auroc(sgccda.res, roc.block = "omics2", roc.comp = 2)
            print(auc.splsda)
            
            png(file="roc_diablo_final_omics2.png")
            print(auc.splsda2)
            dev.off()
            
            
            # prepare test set data: here one block (proteins) is missing
            data.test = list(omics1 = test_omic1, 
                                  omics2 = test_omic2)
            
            testDa_multi=as.data.frame(testDa_multi)
            Y_test= factor(testDa_multi[,y_pred])
            
            print(dim(Y_test))
            print(summary(Y_test))
            
            predict.diablo = predict(sgccda.res, newdata = data.test)
            # the warning message will inform us that one block is missing
            #predict.diablo # list the different outputs
            
            confusion.mat = get.confusion_matrix(truth = Y_test, 
                                 predicted = predict.diablo$WeightedVote$centroids.dist[,2])
            
            print(confusion.mat)
            
            perf=evaluate.DIABLO.performance(confusion.mat, Y_test, predict.diablo$WeightedVote$centroids.dist[,2])
            
            print(perf)
            # table=as.data.frame(confusion.mat)
            # # table_2=table(table)
            # colnames(table)=rownames(table)
            # rownames(table)=colnames(table)
            # print(table)
            # # # a=table(predict.diablo$WeightedVote$centroids.dist[,2], Y_test)
            # # # print(a)
            # # library(caret) 
            # # a=confusionMatrix(table(table))
            # # recall=a$byClass["Recall"]
            # # precision=a$byClass["Precision"]
            # pecc_2= (table[2,2]+table[1,1])/(table[2,2]+table[1,1]+table[1,2]+table[2,1])
            # precision= table[1,1]/(table[1,1]+table[2,1])
            # recall= table[1,1]/(table[1,1]+ table[1,2])
            
            # print(pecc_2)
            # print(precision)
            # print(recall)
            
            
            
            pecc=function(obs,pred) sum(obs==pred)/length(pred)
                  
            print(get.BER(confusion.mat))
            
            AUC= 1-get.BER(confusion.mat)
            print(pecc(predict.diablo$WeightedVote$centroids.dist[,2],Y_test))
         
            diablo_res= list(auc.splsda,confusion.mat,get.BER(confusion.mat),pecc(predict.diablo$WeightedVote$centroids.dist[,2],Y_test), perf,AUC)   
            
            return (diablo_res)}''')
       
        trainDa_multi= self.res_concat[0] 
        trainDa_multi_=robjects.DataFrame(trainDa_multi)
        print(utils.head(trainDa_multi_))
        Y_train=trainDa_multi_.rx(self.y_pred)
        # print(Y_train)
        Y_train_=robjects.FactorVector(self.Y_train)
        train_omic1_=robjects.DataFrame(self.train_omic1)
        train_omic2_=robjects.DataFrame(self.train_omic2)
        
        test_omic1_=robjects.DataFrame(self.test_omic1)
        test_omic2_=robjects.DataFrame(self.test_omic2)
        
        testDa_multi= self.res_concat[1] 
        testDa_multi_=robjects.DataFrame(testDa_multi)
        Y_test= testDa_multi_.rx(self.y_pred)
        # print(Y_train)
        Y_test_=robjects.FactorVector(self.Y_test)
        
        data=base.list(omic1= train_omic1_, omic2= train_omic2_)
        
        print(base.lapply(data,base.dim))
        
        #outcome
        trainDa_multi= base.as_data_frame(trainDa_multi_)
        # Y_i= factor(trainDa_multi$Berry)
        print(Y_train_)
        print(base.summary(Y_train_))
        
        
        # Y_train_= robjects.DataFrame(self.Y_train)
        exec_diablo= robjects.r["exec_diablo"]
        sgccda_res= exec_diablo(train_omic1_,train_omic2_, Y_train,self.y_pred,trainDa_multi)
        predict_diablo=robjects.r["predict_diablo"]
        res_pred= predict_diablo(sgccda_res,self.y_pred, test_omic1_, test_omic2_, testDa_multi)
        
        auc_splsda= res_pred[0]
        confusion_mat= res_pred[1]
        BER= res_pred[2]
        PECC= res_pred[3]
        # perf=res_pred[4]
        Accuracy=res_pred[4][0]
        Sensitivity=res_pred[4][1]
        Specificity=res_pred[4][2]
        Recall=res_pred[4][3]
        AUC=res_pred[5]
        # precision=res_pred[4]
        # recall=res_pred[5]

        
        file= open("Diablo_results.txt", "w")
        file.write(f'Model (different outputs): {auc_splsda} \n')
        file.write(f'Dimension Test Data: {base.dim(testDa_multi)} \n')
        file.write(f'Summary Y_test:\n {base.summary(Y_test)} \n')
        file.write(f'Confusion Matrix: \n {confusion_mat} \n')
        file.write(f'BER: {BER} \n')
        file.write(f'PECC: {PECC} \n')
        # file.write(f'Accuracy, Sensitivity, Specificity, Recall, AUC: {perf} \n')
        file.write((f'Accuracy: {Accuracy} \n'))
        file.write((f'Sensitivity: {Sensitivity} \n'))
        file.write((f'Specificity: {Specificity} \n'))
        file.write((f'Recall: {Recall} \n'))
        file.write((f'AUC: {AUC} \n'))
        
        # file.write(f'Precision: {precision} \n')
        # file.write(f'Recall: {recall} \n')

        file.close()
    
        print("Done DIABLO")
        
    
    def SMSPL(self): 
        
        robjects.r(''' packages <- function( ) {
            install.packages("Matrix")
            library(Matrix)
            
            install.packages("tseries")
            library(tseries)
            
            install.packages("glmnet")
            library(glmnet)
            
            library(caret)
            
            install.packages("ncvreg")
            library(ncvreg)
            
            install.packages("pROC")
            library(pROC)
            
            install.packages("ROCR")
            library(ROCR)
            
            
            library(ggplot2)
            
            install.packages("reticulate")
            library(reticulate)
            return (print("Packages installed!"))}''')
            
            
            
            
        robjects.r(''' evaluate.ensemble <- function(predictlabels,truelabels){
  
          prelable_final = matrix(0, nrow = length(truelabels), ncol = 1)
          for(i in 1:length(truelabels)){
            freq_table <- as.data.frame(table(predictlabels[i,]))
            vote_index <- which.max(freq_table$Freq)
            prelable_final[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
          }
          return(prelable_final)
        } ''')
        
        robjects.r(''' evaluate.DIABLO.performance <- function(confusion.mat, true_lable, predict_label){
          tp <- confusion.mat[2,2]
          tn <- confusion.mat[1,1]
          fp <- confusion.mat[2,1]
          fn <- confusion.mat[1,2]
          
          Accuracy <- (tp + tn)/(tp + tn + fp + fn)
          Sensitivity <- tp/(tp + fn)
          Specificity <- tn/(tn + fp)
          Recall = tp/(tp + fp)
          
          predict_label <- as.factor(predict_label)
          true_lable <- as.numeric(true_lable)
          predict_label <- as.numeric(predict_label)
          pred <- prediction(predict_label, true_lable)
          perf <- performance(pred, measure = "tpr", x.measure = "fpr")
          auc <- performance(pred, measure = "auc")
          AUC <- auc@y.values[[1]]
          
          perf <- c(Accuracy, Sensitivity, Specificity, Recall, AUC)
          return(perf)
          
        }''')
        
        
        robjects.r('''mvselfpace.rank.multiclass <- function(dev_prob, true_label, lambda, 
                                               gamma, v_iter, View_id, View_num, 
                                               sample_select) {
          
          # initialize the loss function
          loss = matrix(0, nrow = length(true_label), ncol = View_num)
          label = matrix(1, nrow = length(true_label), ncol = 1)
          
          #calculate the loss
          for(m in 1: View_num){
            if(m != View_id){
              loss[,m] = (dev_prob[,m] - label)^2    #squared error
            }else{
              next;
            }
          }
          
          # Update v(View_num-j)
          for(m in 1:View_num){
            if(m != View_id){
              for(i in 1:length(true_label)){
                if(loss[i,m] < lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m])){
                  v_iter[i,m] = 1
                }else{
                  v_iter[i,m] = 0
                }
              }
            }
          }
          
          # Update vj
          loss[,View_id] = (dev_prob[,View_id] - label)^2
          for(i in 1:length(true_label)){
            if(loss[i,View_id] < lambda[View_id] + gamma * (sum(v_iter[i,])-v_iter[i,View_id])){
              v_iter[i,View_id] = 1
            }else{
              v_iter[i,View_id] = 0
            }
          }
          
          ## sort sample
          class.idx <- list()
          sample.thr <- list()
          selectedidx <- list()
          selectedsample <- list()
          V_iter = matrix(0, nrow = length(true_label), ncol = View_num)
          for(i in 1:length(unique(true_label))){
            sample.thr[[i]] <- matrix(0, nrow = 1, ncol = View_num)
          }
          
          for(i in 1:View_num){
            for(j in 1:length(unique(true_label))){
              
              class.idx <- which(true_label==j)
              
              if(length(which(v_iter[class.idx,i]==1))!=0){
                sample.thr[[j]][,i] <- sort(loss[class.idx,i][which(v_iter[class.idx,i]==1)])[min(length(class.idx), 
                                                                                                  sample_select[[j]],
                                                                                                  length(which(v_iter[class.idx,i]==1)))]
              }
              if(length(unique(loss[class.idx,i]))!=1){
                selectedidx[[j]] <- intersect(class.idx[which(v_iter[class.idx,i] == 1)],   ## v_iter = 1 && loss is small
                                              class.idx[which(loss[class.idx,i] <= sample.thr[[j]][,i])])
              }
            }
            selectedsample[[i]] <- unlist(selectedidx)
            V_iter[selectedsample[[i]],i] = 1
          }
        
          cat("The ",View_id, "-th modality select ",length(selectedsample[[i]])," samples.\n", sep = "")
          
          #return the result
          return(V_iter)
          
        }''')
        
        robjects.r('''calculate.final.label <- function(pred_label, true_label){
          final_label <- matrix(0, nrow = length(true_label), ncol = 1)
          for(i in 1:length(true_label)){
            freq_table <- as.data.frame(table(pred_label[i,]))
            vote_index <- which.max(freq_table$Freq)
            final_label[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
          }
          return(final_label)
          
        }''')
        
        robjects.r('''selfpaced.muliticlass <- function(X_train, Y_train, X_test, Y_test, lambda, uplambda, Iter_num) {
          
          #the list storing the result for each iteration
          valpredmatrix = list()
          valprobmatrix = list()
          evlpredmatrix = list()
          evlprobmatrix = list()
          coefmatrix = list()
          coefnummatrix = list()
          coefnamematrix = list()
          valmaps <- replicate(Iter_num,0)
          evlmaps <- replicate(Iter_num,0)
          label_train <- matrix(1, nrow = length(Y_train), ncol = 1)
          label_test <- matrix(1, nrow = length(Y_test), ncol = 1)
          
          #the starting values
          cvfit<-cv.glmnet(X_train,Y_train,alpha=1,family="multinomial",type.multinomial="grouped")
          val.pred <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
          val.prob <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
          
          v.idx = selfpace.rank.multiple(dev_decval = val.prob, 
                                         dev_labels = Y_train, 
                                         lambda = lambda)
          this.training.vidx = v.idx
          iters.vidx = list()	
          
          for(iter in 1:Iter_num) {
            if(length(this.training.vidx) == length(Y_train)){break}
            cat("Starting the ",iter,"-th iteration.\t", sep = "")
            # iters.vidx[[iter]] = this.training.vidx
            
            # glmnet (Lasso, Elastic net & L2)
            cvfit<-cv.glmnet(data.matrix(X_train[this.training.vidx,]),
                             Y_train[this.training.vidx],
                             alpha=1,
                             family="multinomial",
                             type.multinomial="grouped",
                             nfolds = 5, 
                             lambda = seq(0.06,0.11,by=0.01))
            
            valprediction <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
            valprobiction <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
            
            tsprediction <- as.numeric(predict(cvfit,X_test,type="class",s="lambda.min"))
            tsprobiction <- apply(predict(cvfit,X_test,type="response",s="lambda.min"),1,max)
            
            coefprediction <- predict(cvfit,type="coefficients",s="lambda.min")  # coef
            coef.idx <- which(coefprediction$`1`[-1]!=0)
            coef.name <- rownames(coefprediction$`1`)[coef.idx]
            coef.number <- length(coef.idx)
            
            #evaluate the training and test error
            val_loss <- sum((valprobiction - label_train)^2)
            evl_loss <- sum((tsprobiction - label_test)^2)
            
            #self-paced learning
            selectedidx = selfpace.rank.multiple(dev_decval = valprobiction, 
                                                 dev_labels = Y_train, 
                                                 lambda = lambda)
            this.training.vidx = selectedidx
            cat("Select ", length(selectedidx), " samples.\t", sep = "")
            
            #change the parameter accoding to the step size
            for(i in 1:length(lambda)){
              lambda[[i]] = lambda[[i]] + uplambda[[i]]
            }
            
            
            #store the prediction for this iteration
            valpredmatrix[[iter]] = valprediction
            evlpredmatrix[[iter]] = tsprediction
            coefmatrix[[iter]]= coefprediction
            coefnummatrix[[iter]] = coef.number
            coefnamematrix[[iter]] = coef.name
            valmaps[iter] <- val_loss
            evlmaps[iter] <- evl_loss
            
            cat("Finish the ",iter,"-th iteration.\n", sep = "")
          }
          
          results <- list("valpredmatrix" = valpredmatrix, 
                          "evlpredmatrix" = evlpredmatrix, 
                          "valmaps" = valmaps,
                          "evlmaps" = evlmaps,
                          "Coef" = coefmatrix, 
                          "NumbernonzeroCoef" = coefnummatrix,
                          "Coef.name" = coefnamematrix)
          return(results)
        }''')
        
        
        robjects.r('''selfpace.rank.multiple <- function(dev_decval, dev_labels, lambda) {
          
          #calculate the loss
          label = matrix(1, nrow = length(dev_labels), ncol = 1)
          loss = (dev_decval-label)^2	
          
          class.idx <- list()
          sample.thr <- list()
          selectedposidx <- list()
          for(i in 1:length(unique(dev_labels))){
            class.idx[[i]] <- which(dev_labels==i)
            sample.thr[[i]] <- sort(loss[class.idx[[i]]])[min(length(class.idx[[i]]), lambda[[i]])]
            
            if(length(unique(loss[class.idx[[i]]]))!=1){
              selectedposidx[[i]] <- class.idx[[i]][which(loss[class.idx[[i]]] <= sample.thr[[i]])]
            }else{
              selectedposidx[[i]] <- sample(class.idx[[i]], 
                                            size = min(lambda[[i]], length(class.idx[[i]])), 
                                            replace = FALSE)
            }
          }
          
          selectedidx = unlist(selectedposidx)
          
          return(selectedidx)
        }''')
        
        
        
        robjects.r('''selfpaced.EN.multiclass <- function(X_train, Y_train, lambda, uplambda, Iter_num) {
          
          #the list storing the result for each iteration
          cvfitmatrix = list()
          valmaps <- replicate(Iter_num,0)
          label_train <- matrix(1, nrow = length(Y_train), ncol = 1)
          
          #the starting values
          cvfit<-cv.glmnet(X_train,Y_train,alpha=1,family="multinomial",type.multinomial="grouped")
          valpred <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
          valprob <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
          
          this.training.vidx = selfpace.rank1.multiclass(dev_decval = valprob, 
                                                         dev_labels = Y_train, 
                                                         lambda = lambda)
          
          for(iter in 1:Iter_num){
            
            if(length(this.training.vidx) == length(Y_train)){break}
            cat("Starting the ",iter,"-th iteration.\t", sep = "")
            
            # glmnet (Lasso, Elastic net & L2)
            cvfit<-cv.glmnet(data.matrix(X_train[this.training.vidx,]),
                             Y_train[this.training.vidx],
                             alpha=1,
                             family="multinomial",
                             type.multinomial="grouped",
                             nfolds = 5, 
                             lambda = seq(0.06,0.11,by=0.01))
            
            valprediction <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
            valprobiction <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
            
            #self-paced learning
            selectedidx = selfpace.rank1.multiclass(dev_decval = valprobiction, dev_labels = Y_train, lambda)
            this.training.vidx = selectedidx
            cat("Select ", length(selectedidx), " samples.\t", sep = "")
            
            #change the parameter accoding to the step size
            for(i in 1:length(lambda)){
              lambda[[i]] = lambda[[i]] + uplambda[[i]]
            }
            
            #store the prediction for this iteration
            cvfitmatrix[[iter]] <- cvfit
            
            #evaluate the training and test error
            val_loss <- sum((valprobiction - label_train)^2)
            valmaps[iter] <- val_loss
            
            cat("Finish the ",iter,"-th iteration.\n", sep = "")
          }
          
          ## best results ##
          best.iter <- which(valmaps == min(valmaps[1:length(cvfitmatrix)]))
          best.cvfit <- cvfitmatrix[[best.iter]]
          
          return(best.cvfit)
        }''')
        
        
        robjects.r('''selfpace.rank1.multiclass <- function(dev_decval, dev_labels, lambda){
          #calculate the loss
          label = matrix(1, nrow = length(dev_labels), ncol = 1)
          loss = (dev_decval-label)^2	
          
          class.idx <- list()
          sample.thr <- list()
          selectedposidx <- list()
          for(i in 1:length(unique(dev_labels))){
            class.idx[[i]] <- which(dev_labels==i)
            sample.thr[[i]] <- sort(loss[class.idx[[i]]])[min(length(class.idx[[i]]), lambda[[i]])]
            
            if(length(unique(loss[class.idx[[i]]]))!=1){
              selectedposidx[[i]] <- class.idx[[i]][which(loss[class.idx[[i]]] <= sample.thr[[i]])]
            }else{
              selectedposidx[[i]] <- sample(class.idx[[i]], 
                                            size = min(lambda[[i]], length(class.idx[[i]])), 
                                            replace = FALSE)
            }
          }
          
          selectedidx = unlist(selectedposidx)
          
          return(selectedidx)
        }''')
        
        
        robjects.r(''' evaluate.ensemble.sPLSDA <- function(predictlabels,truelabels){
          prelable_final = matrix(0, nrow = length(truelabels), ncol = 1)
          
          Prelabel_matrix = matrix(0, nrow = length(truelabels), ncol = dim(predictlabels)[2])
          colnames(Prelabel_matrix) <- colnames(predictlabels)
          
          for(i in 1:length(truelabels)){
            for(j in 1:dim(predictlabels)[2]){
              Prelabel_matrix[i,j] <- as.integer(as.character(predictlabels[i,j]))
            }
          }
          
          for(i in 1:length(truelabels)){
            freq_table <- as.data.frame(table(Prelabel_matrix[i,]))
            vote_index <- which.max(freq_table$Freq)
            prelable_final[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
          }
          return(prelable_final)
        }''')
        
        robjects.r('''SMSPL.multiclass.rank <- function(dev_prob, dev_labels, v_iter, lambda, gamma, View_num, num_sample) {
          
          # Initialize the loss function
          loss = matrix(0, nrow = length(dev_labels), ncol = View_num)
          label = matrix(1, nrow = length(dev_labels), ncol = 1)
          v_iter_new = matrix(0, nrow = length(dev_labels), ncol = View_num)
          
          # Calculate the loss
          for(m in 1:View_num){
            loss[,m] = (dev_prob[,m] - label)^2    #squared error
          }
          
          # Update the weight of samples
          for(m in 1:View_num){
            for(i in 1:length(dev_labels)){
              if(loss[i,m] < lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m]) - gamma){
                v_iter_new[i,m] <- 1
              }else if(loss[i,m] > lambda[m] + gamma * (sum(v_iter[i,]) - v_iter[i,m])){
                v_iter_new[i,m] <- 0
              }else{
                v_iter_new[i,m] <- (lambda[m] - loss[i,m])/gamma + sum(v_iter[i,]) - v_iter[i,m]
              }
            }
          }
          
          ## sort sample
          class.idx <- list()
          sample.thr <- list()
          selectedidx <- list()
          selectedsample <- list()
          V_iter = matrix(0, nrow = length(dev_labels), ncol = View_num)
          for(i in 1:length(unique(dev_labels))){
            sample.thr[[i]] <- matrix(0, nrow = 1, ncol = View_num)
          }
          
          for(i in 1:View_num){
            for(j in 1:length(unique(dev_labels))){
              
              class.idx <- which(dev_labels==j)
              
              if(length(which(v_iter_new[class.idx,i] != 0)) != 0){
                sample.thr[[j]][,i] <- sort(loss[class.idx,i][which(v_iter_new[class.idx,i]!=0)])[min(length(class.idx), 
                                                                                                  num_sample[[j]],
                                                                                                  length(which(v_iter_new[class.idx,i]!=0)))]
              }
              
              if((length(unique(loss[class.idx,i]))!=1)  &&  (length(class.idx[which(loss[class.idx,i] <= sample.thr[[j]][,i])]) > 1)){
                selectedidx[[j]] <- intersect(class.idx[which(v_iter_new[class.idx,i] == 1)],   ## v_iter = 1 && loss is small
                                              class.idx[which(loss[class.idx,i] <= sample.thr[[j]][,i])])
              }
            }
            selectedsample[[i]] = unlist(selectedidx)
            V_iter[selectedsample[[i]],i] <-v_iter_new[selectedsample[[i]],i]
          }
        
          cat("Selected ",length(selectedsample[[i]])," samples.\n", sep = "")
          #return the result
          return(V_iter)
          
        } ''')
           
                  
        packages= robjects.r["packages"]
        install_packages= packages()
    
    

        robjects.r(''' run_smspl <- function(train_omic1, train_omic2, test_omic1, test_omic2, trainDa_multi_, testDa_multi_, y_pred){
            library(caret)
            library(glmnet)
            
            ##--------------
            ## 1. Load data
            ##--------------
            # source("Function_performance.R")
            # load("/1-data/1-datatrainTestDatasetsNormalized.RDATA")
            
            
            print(dim(train_omic1))
            print(dim(train_omic2))
            
            print(dim(test_omic1))
            print(dim(test_omic2))
            
            ## Preparing data
            data.train <- list("omic1" = as.matrix(train_omic1),
                               "omic2" = as.matrix(train_omic2))
            colnames(data.train$omic1) <- paste("omic1", colnames(data.train$omic1), sep = "_")
            colnames(data.train$omic2) <- paste("omic2", colnames(data.train$omic2), sep = "_")
            
            train_group <- factor(trainDa_multi_[, y_pred])
            print(train_group)
            
            data.test <- list("omic1" = as.matrix(test_omic1),
                              "omic2" = as.matrix(test_omic2))
            colnames(data.test$omic1) <- paste("omic1", colnames(data.test$omic1), sep = "_")
            colnames(data.test$omic2) <- paste("omic2", colnames(data.test$omic2), sep = "_")
            test_group <- factor(testDa_multi_[,y_pred])
            
            y_train <-as.numeric(train_group)
            y_test <-as.numeric(test_group)
            y_train
            
            
            ##--------------
            # 2. SMSPL
            ##--------------
            ##-----------------
            ## Step 1 : Initialization parameters
            ##-----------------
            View_num = 2
            iter_num = 50
            gamma <- 0.56                        # adjusting parameters
            lambda <- c(0.95, 0.76)        # adjusting parameters
            uplambda <- 0.02                     # adjusting parameters
            
            ## setting selected sample number in each iteration
            sample.select <- list()
            sample.add <- list()
            Times <- 50         #approximate iteration times
            
            for(i in 1: length(unique(y_train))){
              sample.select[[i]] <- 10
              sample.add[[i]] <- ceiling(length(which(y_train==i))/Times)
            }
            
            valpredmatrix = list()
            valprobmatrix = list()
            evlpredmatrix = list()
            evlprobmatrix = list()
            coefmatrix = list()
            nonzerocoefmatrix = list()
            coef_idx = list()
            coef_coef = list()
            coef_value = list()
            coef_name = list()
            iter.weight = list()
            selectedidx <- list()
            sample.weight <- list()
            
            loss = matrix(0, nrow = length(y_train), ncol = View_num)
            v_iter = matrix(0, nrow = length(y_train), ncol = View_num)
            
            for(iter in 1:iter_num) {
              valpredmatrix[[iter]] = matrix(0, nrow = length(y_train), ncol = View_num)
              valprobmatrix[[iter]] = matrix(0, nrow = length(y_train), ncol = View_num)
              evlpredmatrix[[iter]] = matrix(0, nrow = length(y_test), ncol = View_num)
              evlprobmatrix[[iter]] = matrix(0, nrow = length(y_test), ncol = View_num)
              coefmatrix[[iter]] =  list()
              nonzerocoefmatrix[[iter]] = matrix(0, nrow = 1, ncol = View_num)
              iter.weight[[iter]] = list()
            }
            
            val_labels <- matrix(1, nrow = length(y_train), ncol = 3)
            evl_labels <- matrix(1, nrow = length(y_test), ncol = 3)
            valmaps <- replicate(iter_num,0)
            evlmaps <- replicate(iter_num,0)
            
            
            ##------------------------------------
            ## Step 2.1: Initialization classifier
            ##------------------------------------
            # using all samples to initialize classifiers
            for(i in 1:View_num){
              cvfit <- cv.glmnet(x = data.train[[i]],
                                 y = y_train,
                                 alpha = 1,
                                 family = "multinomial",
                                 type.multinomial="grouped")
              valpredmatrix[[1]][,i] <- as.numeric(predict(cvfit, data.train[[i]], type = "class", s = "lambda.min"))
              valprobmatrix[[1]][,i] <- apply(predict(cvfit,data.train[[i]], type = "response", s = "lambda.min"), 1, max)
              evlpredmatrix[[i]][,i] <- as.numeric(predict(cvfit, data.test[[i]], type = "class", s = "lambda.min"))
              evlprobmatrix[[1]][,i] <- apply(predict(cvfit,data.test[[i]], type = "response", s = "lambda.min"), 1, max)
            }
            
            v_iter = SMSPL.multiclass.rank(dev_prob = valprobmatrix[[1]], 
                                           dev_labels = y_train, 
                                           v_iter = v_iter,
                                           lambda = lambda,
                                           gamma = gamma,
                                           View_num = View_num,
                                           num_sample = sample.select)
            
            for(i in 1:View_num){
              selectedidx[[i]] <- which(v_iter[,i]!=0)              # selected samples index
              sample.weight[[i]] <- v_iter[which(v_iter[,i]!=0),i]  # the weight of selected samples
              iter.weight[[1]][[i]] <- sample.weight[[i]]
            }
            val_loss <- sum((as.vector(valprobmatrix[[1]]) - val_labels)^2)
            evl_loss <- sum((as.vector(evlprobmatrix[[1]]) - evl_labels)^2)
            valmaps[1] <- val_loss
            evlmaps[1] <- evl_loss
            
            # ##--------------------------
            # ## Step 2.2: Optimization
            # ##--------------------------
            for(iter in 1:iter_num){
            
              cat("Starting the ",iter, "-th iterations.\n", sep = "")
              ## Training Ensemble classifier ##
              for(i in 1:View_num){
                cvfit <- cv.glmnet(x = data.train[[i]][selectedidx[[i]],],
                                   y = y_train[selectedidx[[i]]],
                                   weights = sample.weight[[i]],
                                   alpha = 1,
                                   family = "multinomial",
                                   type.multinomial = "grouped")
                
                valpredmatrix[[iter]][,i] <- as.numeric(predict(cvfit, data.train[[i]], type = "class", s = "lambda.min"))
                valprobmatrix[[iter]][,i] <- apply(predict(cvfit,data.train[[i]], type = "response", s = "lambda.min"), 1, max)
                evlpredmatrix[[iter]][,i] <- as.numeric(predict(cvfit, data.test[[i]], type = "class", s = "lambda.min"))
                evlprobmatrix[[iter]][,i] <- apply(predict(cvfit,data.test[[i]], type = "response", s = "lambda.min"), 1, max)
                coefmatrix[[iter]][[i]] <- predict(cvfit, type = "coefficients", s = "lambda.min")
                nonzerocoefmatrix[[iter]][[i]] <- length(which(coefmatrix[[iter]][[i]]$`1`!=0)-1)
              }
              
              ## evaluate the training and test error
              val_loss <- sum((as.vector(valprobmatrix[[iter]]) - val_labels)^2)
              evl_loss <- sum((as.vector(evlprobmatrix[[iter]]) - evl_labels)^2)
              valmaps[iter] <- val_loss
              evlmaps[iter] <- evl_loss
              
              ## if all the sample used to tranining the classifiers, then stop the iteration.
              if(length(unlist(selectedidx)) == (length(y_train)*View_num)){break}
              
              ## update lambda and valpredmatrix for next iteriation
              lambda <- uplambda + lambda
              for(i in 1:length(sample.select)){
                sample.select[[i]] <- sample.select[[i]] + sample.add[[i]]
              }
              
              ## select samples
              v_iter = SMSPL.multiclass.rank(dev_prob = valprobmatrix[[iter]], 
                                             dev_labels = y_train, 
                                             v_iter = v_iter,
                                             lambda = lambda,
                                             gamma = gamma,
                                             View_num = View_num,
                                             num_sample = sample.select)
              
              for(i in 1:View_num){
                selectedidx[[i]] <- which(v_iter[,i]!=0)              # selected samples index
                sample.weight[[i]] <- v_iter[which(v_iter[,i]!=0),i]  # the weight of selected samples
                iter.weight[[iter]][[i]] <- sample.weight[[i]]
              }
            
            }
            
            ##----------------------------------------------------
            # Step 3: Find the run with the best valudation map
            ##----------------------------------------------------
            ## best results ##
            best.iter <- which(valmaps == min(valmaps[1:length(which(valmaps!=0))]))
            best_valperf <- valpredmatrix[[best.iter]]
            best_evlperf <- evlpredmatrix[[best.iter]]
            best_coef <- coefmatrix[[best.iter]]
            best_numcoef <- nonzerocoefmatrix[[best.iter]]
            
            ## record label
            final_val_label <- calculate.final.label(pred_label = best_valperf, true_label = y_train)
            final_evl_label <- calculate.final.label(pred_label = best_evlperf, true_label = y_test)
            
            ## record selected features
            
            for(i in 1:View_num){
              coef_idx[[i]] <- which(best_coef[[i]]$`1`!=0)[-1]
              coef_coef[[i]] <- best_coef[[i]]$`1`[coef_idx[[i]]]
              coef_name[[i]] <- rownames(best_coef[[i]]$`1`)[coef_idx[[i]]]
              
            }
            coef.omic1 <- cbind(coef_name[[1]], coef_coef[[1]])
            coef.omic2 <- cbind(coef_name[[2]], coef_coef[[2]])
            
            
            ## calculate features ##
            num_feature <- matrix(0, nrow = 1, ncol = 2)
            for(i in 1:View_num){
              num_feature[,i] <- length(coef_name[[i]])
            }
            
            ##-------------------------------------------
            ## Step 4: Evaluate the best performance
            ##-------------------------------------------
            conMatrix.train <- confusionMatrix(as.factor(y_train),
                                               as.factor(final_val_label),
                                               mode = "prec_recall")
            conMatrix.test <- confusionMatrix(as.factor(y_test),
                                              as.factor(final_evl_label),
                                              mode = "prec_recall")
            
            
            
            png(file="SMSPL_ROC_curves.png", width=512, height=512)
            par(mfcol = c(1,2))
            roc(response = y_train, predictor =final_val_label ,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,print.auc=TRUE)
            title(main = "SMSPL_TRAIN")
            
            roc(response = y_test, predictor =final_evl_label ,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,print.auc=TRUE)
            title(main = "SMSPL_TEST")
            dev.off()
            
            
            
            ## record results
            perf.SMSPL <- list("coef.omic1" = coef.omic1,
                               "coef.omic2" = coef.omic2,
                               "feature.num" = num_feature,
                               "Perf.Train" = conMatrix.train,
                               "Perf.Test" = conMatrix.test)
                    
                    
                    
                   return (perf.SMSPL) } ''')    
            
        run_smspl=robjects.r["run_smspl"]
        
        trainDa_multi= self.res_concat[0] 
        trainDa_multi_=robjects.DataFrame(trainDa_multi)
        train_omic1_=robjects.DataFrame(self.train_omic1)
        train_omic2_=robjects.DataFrame(self.train_omic2)
        Y_train=trainDa_multi_.rx(self.y_pred)
        
        test_omic1_=robjects.DataFrame(self.test_omic1)
        test_omic2_=robjects.DataFrame(self.test_omic2)
        testDa_multi= self.res_concat[1] 
        testDa_multi_=robjects.DataFrame(testDa_multi)
        Y_test= testDa_multi_.rx(self.y_pred)
        # print(Y_train)
        # Y_test_=robjects.FactorVector(self.Y_test)
        perf_SMSPL= run_smspl(train_omic1_, train_omic2_, test_omic1_, test_omic2_, trainDa_multi_, testDa_multi, self.y_pred)
        print(perf_SMSPL)
        
        file= open("SMSPL_results.txt", "w")
        file.write(f'SMSPL results: {perf_SMSPL} \n')
        file.close()
        
        return print("SMSPL excecuted successfully!")
    
    
    
    def Stack_Generalisation(self):
        
        ''' Stack Generalisation with only glms (Generalised Linear Model) '''
        robjects.r('''stacked <- function(trainDa_multi, testDa_multi, y_pred){
                   
                   library(h2o)
                   h2o.init(nthreads = -1)
                   h2o.removeAll()  

                    # Import a sample binary outcome train/test set into H2O
                    train <- as.h2o(trainDa_multi)
                    test <- as.h2o(testDa_multi)
            
                    y <- y_pred
                    x <- setdiff(names(train), y)

                    # For binary classification, response should be a factor
                    train[, y] <- as.factor(train[, y])
                    test[, y] <- as.factor(test[, y])

                    # Number of CV folds (to generate level-one data for stacking)
                    family <- "binomial"
                    nfolds <- 15

                    # There are a few ways to assemble a list of models to stack toegether:
                    # 1. Train individual models and put them in a list
                    # 2. Train a grid of models
                    # 3. Train several grids of models
                    # Note: All base models must have the same cross-validation folds and
                    # the cross-validated predicted values must be kept.


                    # 1. Generate a 2-model ensemble (GLM + GLM)

                    # Train & Cross-validate a GBM
                    glm1 <- h2o.glm(  
              x = x, y = y,
              family = family,
              training_frame = train,
              nfolds = nfolds,
              fold_assignment = "Modulo",
              keep_cross_validation_predictions = TRUE,
              seed=1
              )
            
            
            glm2 <- h2o.glm(  
              x = x, y = y,
              family = family,
              training_frame = train,
              nfolds = nfolds,
              fold_assignment = "Modulo",
              keep_cross_validation_predictions = TRUE,
              seed=1
              )


            # Train a stacked ensemble using the GBM and RF above
            ensemble <- h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = train,
                                base_models = list(glm1, glm2))

            # Eval ensemble performance on a test set
            perf <- h2o.performance(ensemble, newdata = test)

            # Compare to base learner performance on the test set
            perf_gbm_test <- h2o.performance(glm1, newdata = test)
            perf_rf_test <- h2o.performance(glm2, newdata = test)
            baselearner_best_auc_test <- max(h2o.auc(perf_gbm_test), h2o.auc(perf_rf_test))
            ensemble_auc_test <- h2o.auc(perf)
            print(sprintf("Best Base-learner Test AUC:  %s", baselearner_best_auc_test))
            print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
            # [1] "Best Base-learner Test AUC:  0.76979821502548"
            # [1] "Ensemble Test AUC:  0.773501212640419"
            
            # Generate predictions on a test set (if neccessary)
            pred <- h2o.predict(ensemble, newdata = test)
            
            
            # find importance
            my_varimp1 <- h2o.varimp(glm1)
            my_varimp2 <- h2o.varimp(glm2)
            # print(my_varimp1)
            # print(my_varimp2)
            write.table(my_varimp1,"Stack_varImp1.txt", sep=",",col.names = TRUE, row.names = FALSE)
            write.table(my_varimp2,"Stack_varImp2.txt", sep=",",col.names = TRUE, row.names = FALSE)
            list=list(perf,pred, my_varimp1, my_varimp2)
            #
                   
                   
                   
             return (list) }''')
                
                
                
        robjects.r('''run_stack <- function(trainDa_multi, testDa_multi, y_pred){
            
            install.packages("devtools")
            library(devtools)
            install_github("h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package")
            library(h2oEnsemble)
            h2o.init(nthreads = -1)
            h2o.removeAll()
            # import data
            
            # convert to h2o objects
            train <- as.h2o(trainDa_multi)
            test <- as.h2o(testDa_multi)
            
            
            
            y <- y_pred
            x <- setdiff(names(train), y)
            family <- "binomial"
            nfolds <- 5  
            
            train[, y] <- as.factor(train[, y])  
            test[, y] <- as.factor(test[, y])
            
            glm1 <- h2o.glm(  
              x = x, y = y,
              family = family,
              training_frame = train,
              nfolds = nfolds,
              fold_assignment = "Modulo",
              keep_cross_validation_predictions = TRUE
              )
            
            
            glm2 <- h2o.glm(  
              x = x, y = y,
              family = family,
              training_frame = train,
              nfolds = nfolds,
              fold_assignment = "Modulo",
              keep_cross_validation_predictions = TRUE
              )
            
            
            models <- list(glm1, glm2)
            metalearner <- "h2o.glm.wrapper"
            
            
            
            # stack existing models
            stack <- h2o.stack(
              models = models,
              response_frame = train[,y],
              metalearner = metalearner, 
              seed = 1,
              keep_levelone_data = TRUE
              )
            
            
            # Compute test set performance:
            perf <- h2o.ensemble_performance(stack, newdata = test)
            # h2o.accuracy(perf)
            # h2o.auc(perf)
            
            # png(file="Stack_Generalisation_ROC.png", width=512, height=512)
            # plot(perf, type = "roc")
            # title(main = "Stack_Generalisation_ROC.png")
            # dev.off()
            
            # a=confusionMatrix(table(perf, test[, y]))
            # print(a)
            # recall=a$byClass["Recall"]
            # precision=a$byClass["Precision"]
            
            # print(a)
            # print(recall)            
            # print(precision)
            
            # png(file="SNFtool_ROC.png", width=512, height=512)
            # roc(response = develop[-trainSample], predictor =newLabel[-c(1:n)] ,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE ,print.auc=TRUE)
            # title(main = "SNFtool")
            # dev.off()
            
            pred <- h2o.predict(stack, newdata = test)
            
            
            return (perf)}''')
        run_stack= robjects.r["stacked"]
        
        trainDa_multi= self.res_concat[0] 
        trainDa_multi_=robjects.DataFrame(trainDa_multi)
        
        testDa_multi= self.res_concat[1] 
        testDa_multi_=robjects.DataFrame(testDa_multi)
        
        
        lista= run_stack(trainDa_multi, testDa_multi, self.y_pred)
        perf= lista[0]
        pred=lista[1]
        varImp1=lista[2]
        varImp2=lista[3]
        print(varImp1)
        print(varImp2)
        file= open("Stack_results_CB.txt", "w")
        file.write(f'Stack Generalisation results:\n {perf} \n')
        file.write(f'Stack Generalisation Test Prediction:\n {pred} \n')
        file.write(f'Variable Importance Model1: Ckeck "Stack_varImp1.txt" e "Stack_varImp2.txt" to see full results \n {varImp1} \n')
        file.write(f'Variable Importance Model2: \n {varImp2} \n')
        file.close()
        
        return print("Stack Generalisation excecuted successfully!")
    
    def Lasso_regression(self):
        
        
        
        
        
        
        robjects.r(''' run_lasso  <- function(trainDa_multi, testDa_multi, y_pred){
            install.packages("glmnet")
            library(glmnet)
            eval_results <- function(true, predicted, df) {
              SSE <- sum((predicted - true)^2)
              SST <- sum((true - mean(true))^2)
              R_square <- 1 - SSE / SST
              RMSE = sqrt(SSE/nrow(df))
            
               # Model performance metrics
            data.frame(
              RMSE = RMSE,
              Rsquare = R_square
            )
              
            }
              
            Y_train= trainDa_multi[,y_pred]
            Y_test= testDa_multi[,y_pred]
            lambdas <- 10^seq(2, -3, by = -.1)
            
            # Setting alpha = 1 implements lasso regression
            library(glmnet)
            
            
            lasso_reg <- cv.glmnet(as.matrix(trainDa_multi[,1:(length(trainDa_multi)-1)]), as.numeric(factor(Y_train)), alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5, family="binomial")
            
            # Best 
            lambda_best <- lasso_reg$lambda.min 
            lambda_best
            
            
            lasso_model <- glmnet(as.matrix(trainDa_multi[,1:(length(trainDa_multi)-1)]),  as.numeric(factor(Y_train)), alpha = 1, lambda = lambda_best, standardize = TRUE, family="binomial")
            
            assess= assess.glmnet(lasso_model,newx= as.matrix(testDa_multi[,1:(length(testDa_multi)-1)]), newy=testDa_multi[,y_pred])
            print(assess)
            cnf <- confusion.glmnet(lasso_model, newx= as.matrix(testDa_multi[,1:(length(testDa_multi)-1)]), newy=testDa_multi[,y_pred])
            print(cnf)
            png(file="Lasso_ROC.png")
            plot(roc.glmnet(lasso_model,newx= as.matrix(testDa_multi[,1:(length(testDa_multi)-1)]), newy=testDa_multi[,y_pred] ), type="l")
            dev.off()
                   
            co= coef(lasso_model, s = "lambda.1se")
            print(co)
            # write.table(co,"LASSO_varImp.txt", sep=",",col.names = TRUE, row.names = FALSE)
            inds<-which(co!=0)
            variables<-row.names(co)[inds]
            variables<-variables[!(variables %in% '(Intercept)')];
            print(variables)
            
            predictions_train <- predict(lasso_model, s = lambda_best, newx = as.matrix(trainDa_multi[,1:(length(trainDa_multi)-1)]))
            print(eval_results(as.numeric(factor(Y_train)), predictions_train, trainDa_multi))
            
            predictions_test <- predict(lasso_model, s = lambda_best, newx = as.matrix(testDa_multi[,1:(length(testDa_multi)-1)]))
            
            # print(table(predictions_test,Y_test))
            # print("Accuracy:")
            # print(mean(lasso.pred==testy))
            # # recall=a$byClass["Recall"]
            # # precision=a$byClass["Precision"]
            
            print(eval_results(as.numeric(factor(Y_test)), predictions_test, testDa_multi))
            
            lasso_res= list(eval_results(as.numeric(factor(Y_train)), predictions_train, trainDa_multi),eval_results(as.numeric(factor(Y_test)), predictions_test, testDa_multi), assess, cnf, variables)

            
                   
                   
            return (lasso_res)}''')
        run_lasso= robjects.r["run_lasso"]
        trainDa_multi= self.res_concat[0] 
        trainDa_multi_=robjects.DataFrame(trainDa_multi)
        
        testDa_multi= self.res_concat[1] 
        testDa_multi_=robjects.DataFrame(testDa_multi)
        
        lasso_res= run_lasso(trainDa_multi_, testDa_multi_, self.y_pred)
        
        Train_res= lasso_res[0]
        Test_res= lasso_res[1]
        
        Assess= lasso_res[2]
        cnf=lasso_res[3]
        variables=lasso_res[4]
        file= open("Lasso_regression.txt", "w")
        file.write('Lasso Regression results: \n')
        file.write(f'Train Results: \n {Train_res} \n')
        file.write(f'Test Results: \n {Test_res} \n')
        file.write(f'Performance asssessment: \n {Assess} \n')
        file.write(f'Confusion Matrix: \n {cnf} \n')
        file.write(f'Most important features: \n {variables} \n')
        file.write("See Lasso_ROC.png for the ROC curve.")
        file.close()
        
        return print("Lasso Regresion excecuted successfully!")
    
    
    
    def SVM (self):
        
        omics_=[]
        names=[]
        for name, lista in self.train_test_datasets.items():
            if name in self.omics:
                # confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        print(omics_)
        if len(omics_) ==2:
            
            omics1= omics_[0]
            omics2= omics_[1]
               
            print("Dimension of Dataset Concatenated:")
            dataset_concat= pd.read_excel("dataset_concat.xlsx", index_col=0) 
            print(dataset_concat.shape)
            
            
            
            dim1=base.dim(omics1)
            # print(dim1)
            # range_omic1= np.arange(1, dim1[1], 1)
            
            
            
            # dim2=base.dim(omics2)
            # print(dim2)
            dim_concat=base.dim(dataset_concat)
            # len_concat=dim_concat[1]
            final_2= dim_concat[1]
            # range_omic2= np.arange(dim1[1],final_2, 1)
            # omic_2=dataset_concat.rx(True,range_omic2)
            # print(utils.head(omic_2))
       
        omic1=self.omics[0]
        # omic2=omics[1]
        omic1=dataset_concat.iloc[:,np.arange(0, dim1[1]-1, 1)]
        print("Dimension of "+ names[0] + ":")
        print(omic1.shape)
        
        print("Dimension of "+ names[1]+ ":")
        omics2=dataset_concat.iloc[:,np.arange(dim1[1]-1,final_2-1, 1)]
        print(omics2.shape)
        # Y_svm= dataset_concat.loc[:, self.y_pred]

        X= dataset_concat.drop(self.y_pred, axis=1)
        print(X.shape)

        y= dataset_concat[self.y_pred]
        print(y.shape)
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)
        
        print(X_train.shape, y_train.shape)
        
        
        print(X_test.shape, y_test.shape)
        
        
        from sklearn.model_selection import GridSearchCV
        
        from matplotlib import pyplot as plt
        from sklearn import svm

        # def f_importances(coef, names):
        #     imp = coef
        #     imp,names = zip(*sorted(zip(imp,names)))
        #     plt.barh(range(len(names)), imp, align='center')
        #     plt.yticks(range(len(names)), names)
        #     plt.show()



        features_names = np.array(dataset_concat.columns)
        print(features_names)



        print("Grid Search: \n") 
        parameters = {'kernel':['linear'], 'C':[1, 3, 10, 100], 'gamma':[0.01, 0.001]}
        svm_model_d = SVC( )
        
        opt_model_d = GridSearchCV(svm_model_d, parameters, cv = 2)
        opt_model_d.fit(X, y)
        print (opt_model_d.best_estimator_)
        # print(opt_model_d.best_estimator_.feature_importances_)
        importance = opt_model_d.best_estimator_.coef_[0]
        # # summarize feature importance
        # importance= importance.sort
        # file=open('SVM_variImp_GridSearch.txt','w')
        # for i,v in enumerate(importance):
        #     print(f'Feature: {features_names[i]}, Score: %.5f' % (v))
        #     file.writelines(f'Feature: {features_names[i]}, Score: %.5f \n' % (v))

        # file.close()
        headers = ["name", "score"]
        values = sorted(zip(X.columns, importance), key=lambda x: x[1] * -1)
        print(tabulate(values, headers, tablefmt="plain"))
        importance= tabulate(values, headers, tablefmt="plain")
        
        
        scores_gs = cross_val_score(opt_model_d, X, y, cv = 2)
        print("CV average score: %.2f" % scores_gs.mean())
        
       # precision, recall and F1
        from sklearn.preprocessing import LabelBinarizer
        
        lb = LabelBinarizer()
        y = np.array([number[0] for number in lb.fit_transform(y)])
        y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
        # recall = cross_val_score(svm_model_d, X_train, y_train, cv=5, scoring='recall')
        # print('Recall', np.mean(recall), recall)
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', np.mean(precision), precision)
       
        # print(scores["precision"])
        # lb = preprocessing.LabelBinarizer()
        # lb.fit(y)
        # recall = cross_val_score(svm_model_d, X, y, cv=5, scoring='recall')
        # print('Recall', recall.mean())
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', precision.mean())
        
        
        # y_pred = cross_val_predict(opt_model_d, X, y, cv=5)
        # print(confusion_matrix(y, y_pred))
        # print(classification_report(y, y_pred))
        # print("Accuracy: ",metrics.accuracy_score(y, y_pred))
        # print("Recall: ", metrics.recall_score(y, y_pred))
        # print("Precision: ", metrics.precision_score(y, y_pred))
        # false_positive_rate, true_positive_rate, thresholds = roc_curve(y, y_pred)
        # roc_auc = auc(false_positive_rate, true_positive_rate) 
        # print("ROC_AUC:", roc_auc)  
        
        
        
        y_pred = cross_val_predict(opt_model_d, X_test, y_test, cv=2)
        rs_cm= confusion_matrix(y_test, y_pred)
        print(rs_cm)
        cl_rs=classification_report(y_test, y_pred)
        print(cl_rs)
        ac_rs=metrics.accuracy_score(y_test, y_pred)
        print("Accuracy: ", ac_rs)
        re_rs= metrics.recall_score(y_test, y_pred)
        print("Recall: ", metrics.recall_score(y_test, y_pred))
        prec_rs=metrics.precision_score(y_test, y_pred)
        print("Precision: ", prec_rs)
        false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
        roc_auc = auc(false_positive_rate, true_positive_rate) 
        print("ROC_AUC:", roc_auc)  
        
        
        
        
        
        plt.figure()
        plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Grid Search')
        plt.legend(loc="lower right")
        plt.savefig('ROC_SVM_gridsearch.png')
        plt.show()
        print("Random Search: \n")
        
        def report(results, n_top=3):
            for i in range(1, n_top + 1):
                candidates = np.flatnonzero(results['rank_test_score'] == i)
                for candidate in candidates:
                    print("Model with rank: {0}".format(i))
                    print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                          results['mean_test_score'][candidate],
                          results['std_test_score'][candidate]))
                    print("Parameters: {0}".format(results['params'][candidate]))
                    print("")
        
        

        svm_model = SVC( )
        
        parameters = {'kernel':['linear'], 'C':[1, 3, 10, 100], 'gamma':[0.01, 0.001]}
        
        rand_search = RandomizedSearchCV(svm_model, param_distributions=parameters, cv = 2)
        
        rand_search.fit(X, y)
        
        
        print (rand_search.best_estimator_)
        report(rand_search.cv_results_)
        # print(rand_search.best_estimator_.feature_importances_)
        # get importance
        # importance = rand_search.best_estimator_.coef_[0]
        # summarize feature importance
        # for i,v in enumerate(importance):
        #     print('Feature: %0d, Score: %.5f' % (i,v))
        
        importance2 = rand_search.best_estimator_.coef_[0]
        # # summarize feature importance
        # importance= importance.sort
        # file=open('SVM_variImp_GridSearch.txt','w')
        # for i,v in enumerate(importance):
        #     print(f'Feature: {features_names[i]}, Score: %.5f' % (v))
        #     file.writelines(f'Feature: {features_names[i]}, Score: %.5f \n' % (v))

        # file.close()
        headers = ["name", "score"]
        values = sorted(zip(X.columns, importance2), key=lambda x: x[1] * -1)
        print(tabulate(values, headers, tablefmt="plain"))
        importance2= tabulate(values, headers, tablefmt="plain")
        
        
        scores_rs = cross_val_score(rand_search, X, y, cv = 2)
        print("CV average score: %.2f" % scores_rs.mean())
        
        
        from sklearn.preprocessing import LabelBinarizer
        
        lb = LabelBinarizer()
        y = np.array([number[0] for number in lb.fit_transform(y)])
        
        y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
        # recall = cross_val_score(svm_model_d, X_train, y_train, cv=5, scoring='recall')
        # print('Recall', np.mean(recall), recall)
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', np.mean(precision), precision)
       
        # print(scores["precision"])
        # lb = preprocessing.LabelBinarizer()
        # lb.fit(y)
        # recall = cross_val_score(svm_model_d, X, y, cv=5, scoring='recall')
        # print('Recall', recall.mean())
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', precision.mean())
        
        
        y_pred_rs = cross_val_predict(rand_search, X_test, y_test, cv=2)
        rs_cm= confusion_matrix(y_test, y_pred_rs)
        print(rs_cm)
        cl_rs=classification_report(y_test, y_pred_rs)
        print(cl_rs)
        ac_rs=metrics.accuracy_score(y_test, y_pred_rs)
        print("Accuracy: ", ac_rs)
        re_rs= metrics.recall_score(y_test, y_pred_rs)
        print("Recall: ", metrics.recall_score(y_test, y_pred_rs))
        prec_rs=metrics.precision_score(y_test, y_pred_rs)
        print("Precision: ", prec_rs)
        false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_rs)
        roc_auc_rs = auc(false_positive_rate, true_positive_rate) 
        print("ROC_AUC:", roc_auc_rs)  
        
        plt.figure()
        plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_rs)
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Random Search')
        plt.legend(loc="lower right")
        plt.savefig('ROC_SVM_randomsearch.png')
        plt.show()
        
      

        
        file= open("SVM_CBI.txt", "w")
        file.write('SVM (Concatenation_Based Integration results: \n')
        file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
        file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
        file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
        file.write("\n")
        file.write("Grid Search: \n")
        file.write(f"Best Estimator: {opt_model_d.best_estimator_}\n")
        file.write(f'Accuracy: \n {metrics.accuracy_score(y_test, y_pred)} \n')
        file.write(f'Precision: \n {metrics.precision_score(y_test, y_pred)} \n')
        file.write(f'Recall: \n { metrics.recall_score(y_test, y_pred)} \n')
        file.write(f'Confusion Matrix: \n {confusion_matrix(y_test, y_pred)} \n')
        file.write(f'Classification report: \n {classification_report(y_test, y_pred)} \n')
        file.write(f'ROC_AUC: \n {roc_auc} \n')
        file.write(f'Feature Importance: \n {importance} \n')
        file.write("\n")        
        file.write("Random Search: \n")
        file.write(f"Best Estimator: {rand_search.best_estimator_}\n")
        file.write(f'Accuracy: \n {ac_rs} \n')
        file.write(f'Precision: \n {prec_rs} \n')
        file.write(f'Recall: \n { re_rs} \n')
        file.write(f'Confusion Matrix: \n {rs_cm} \n')
        file.write(f'Classification report: \n {cl_rs} \n')
        file.write(f'ROC_AUC: \n {roc_auc_rs} \n')
        file.write(f'Feature Importance: \n {importance2} \n')
        file.write("See ROC_SVM_gridsearch.png and ROC_SVM_randomsearch for the ROC curves.")
        
        
        file.close()
        
        # return print("SVM excecuted successfully!")
    
    

    
    def NNs(self):
        omics_=[]
        names=[]
        for name, lista in self.train_test_datasets.items():
            if name in self.omics:
                # confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        if len(omics_) ==2:
            omics1= omics_[0]
            omics2= omics_[1]
               
            print("Dimension of Dataset Concatenated:")
            dataset_concat= pd.read_excel("dataset_concat.xlsx", index_col=0) 
            print(dataset_concat.shape)
            
            
            
            dim1=base.dim(omics1)
            # print(dim1)
            range_omic1= np.arange(1, dim1[1], 1)
            
            
            
            dim2=base.dim(omics2)
            # print(dim2)
            dim_concat=base.dim(dataset_concat)
            len_concat=dim_concat[1]
            final_2= dim_concat[1]
            range_omic2= np.arange(dim1[1],final_2, 1)
            # omic_2=dataset_concat.rx(True,range_omic2)
            # print(utils.head(omic_2))
       
        # omic1=omics[0]
        # omic2=omics[1]
        omics1=dataset_concat.iloc[:,np.arange(0, dim1[1]-1, 1)]
        print("Dimension of "+ names[0] + ":")
        print(omics1.shape)
        
        print("Dimension of "+ names[1]+ ":")
        omics2=dataset_concat.iloc[:,np.arange(dim1[1]-1,final_2-1, 1)]
        print(omics2.shape)
        # Y_nns= dataset_concat.loc[:, self.y_pred]

        X= dataset_concat.drop(self.y_pred, axis=1)
        print(X.shape)

        y= dataset_concat[self.y_pred]
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)
        # mlp = MLPClassifier(hidden_layer_sizes=(8,8,8), activation='relu', solver='adam', max_iter=500)
        # mlp.fit(X_train,y_train)
        # y_pred = mlp.predict(X_test)

        # print(confusion_matrix(y_test,y_pred))
        # print(classification_report(y_test,y_pred))
        
        # # from sklearn.neural_network import MLPClassifier
        # # from sklearn.model_selection import KFold
        # # from sklearn.model_selection import cross_val_score
        
        # mlp = MLPClassifier(hidden_layer_sizes=(8,8,8), activation='relu', solver='adam', max_iter=1000)
        # kf = KFold(n_splits=5)
        
        # result = cross_val_score(mlp , X, y, cv = kf)
        # print("Avg accuracy: {}".format(result.mean()))
        
        # file= open("NNs_CBI.txt", "w")
        # file.write('NNs (Concatenation_Based Integration results: \n')
        # file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
        # file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
        # file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
        # file.write(f'Confusion Matrix: \n {confusion_matrix(y_test, y_pred)} \n')
        # file.write(f'Classification report: \n {classification_report(y_test, y_pred)} \n')
        # file.write(f'With Cross Validation to prevent overfitting: \n {"Avg accuracy: {}".format(result.mean())} \n')
        # file.close()
        
        # return print("NNs excecuted successfully!")
       #  from sklearn.model_selection import GridSearchCV

       #  print("Grid Search: \n") 
       #  parameters = {"solver":('lbfgs','sgd','adam') ,"activation": ("logistic","relu"),"alpha":(0.0001, 0.001, 0.01) , "hidden_layer_sizes": ((15,),(25,),(50,),(75,),(100,))}
       #  mlp_model_d = MLPClassifier(max_iter= 900)
        
       #  opt_model_d = GridSearchCV(mlp_model_d, parameters, cv = 5)
       #  opt_model_d.fit(X, y)
       #  print (opt_model_d.best_estimator_)
       #  scores_gs = cross_val_score(opt_model_d, X, y, cv = 5)
       #  print("CV average score: %.2f" % scores_gs.mean())
        
       # # precision, recall and F1
       #  from sklearn.preprocessing import LabelBinarizer
        
       #  lb = LabelBinarizer()
       #  y = np.array([number[0] for number in lb.fit_transform(y)])
        
       #  # recall = cross_val_score(svm_model_d, X_train, y_train, cv=5, scoring='recall')
       #  # print('Recall', np.mean(recall), recall)
       #  # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
       #  # print('Precision', np.mean(precision), precision)
       
       #  # print(scores["precision"])
       #  # lb = preprocessing.LabelBinarizer()
       #  # lb.fit(y)
       #  # recall = cross_val_score(svm_model_d, X, y, cv=5, scoring='recall')
       #  # print('Recall', recall.mean())
       #  # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
       #  # print('Precision', precision.mean())
        
        
       #  y_pred = cross_val_predict(opt_model_d, X, y, cv=5)
       #  print(confusion_matrix(y, y_pred))
       #  print(classification_report(y, y_pred))
       #  print("Accuracy: ",metrics.accuracy_score(y, y_pred))
       #  print("Recall: ", metrics.recall_score(y, y_pred))
       #  print("Precision: ", metrics.precision_score(y, y_pred))
       #  false_positive_rate, true_positive_rate, thresholds = roc_curve(y, y_pred)
       #  roc_auc = auc(false_positive_rate, true_positive_rate) 
       #  print("ROC_AUC:", roc_auc)  
        
       #  plt.figure()
       #  plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
       #  plt.plot([0, 1], [0, 1], 'k--')
       #  plt.xlim([0.0, 1.0])
       #  plt.ylim([0.0, 1.05])
       #  plt.xlabel('False Positive Rate')
       #  plt.ylabel('True Positive Rate')
       #  plt.title('ROC Grid Search')
       #  plt.legend(loc="lower right")
       #  plt.show()
       #  plt.savefig('ROC_ANN_gridsearch.png')
        print("Random Search: \n")
        
        def report(results, n_top=3):
            for i in range(1, n_top + 1):
                candidates = np.flatnonzero(results['rank_test_score'] == i)
                for candidate in candidates:
                    print("Model with rank: {0}".format(i))
                    print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                          results['mean_test_score'][candidate],
                          results['std_test_score'][candidate]))
                    print("Parameters: {0}".format(results['params'][candidate]))
                    print("")
        
    
        # print(help(MLP()))
        mlp_model = MLPClassifier(max_iter= 900)
            
        parameters = {"solver":('lbfgs','sgd','adam') ,"activation": ("identity","logistic","tanh","relu"),"alpha":(0.0001, 0.001, 0.01) , "hidden_layer_sizes": ((15,),(25,),(50,),(75,),(100,))}

        rand_search = RandomizedSearchCV(mlp_model, param_distributions=parameters, cv = 2)
        
        rand_search.fit(X, y)
        # print (dir(rand_search))
        # print(rand_search.feature_importances_)
    
        # pred = mlp_model2.predict(X)
        # print (rand_search.best_estimator_.feature_importances_)
        report(rand_search.cv_results_)

        # features_names = np.array(dataset_concat.columns)
        # importance = rand_search.best_estimator_.coefs_[1]
        # print(importance)
        # # summarize feature importance
        # file=open('ANN_variImp_GridSearch.txt','w')
        # for i,v in enumerate(importance):
        #     print(f'Feature: {features_names[i]}, Score: %.5f' % (v))
        #     file.writelines(f'Feature: {features_names[i]}, Score: %.5f \n' % (v))

        # file.close()
        
        scores_rs = cross_val_score(rand_search, X, y, cv = 2)
        print("CV average score: %.2f" % scores_rs.mean())
        
        
        from sklearn.preprocessing import LabelBinarizer
        
        lb = LabelBinarizer()
        y = np.array([number[0] for number in lb.fit_transform(y)])
        y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
        
        # print(y)
        # print(y_test)
        # recall = cross_val_score(svm_model_d, X_train, y_train, cv=5, scoring='recall')
        # print('Recall', np.mean(recall), recall)
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', np.mean(precision), precision)
       
        # print(scores["precision"])
        # lb = preprocessing.LabelBinarizer()
        # lb.fit(y)
        # recall = cross_val_score(svm_model_d, X, y, cv=5, scoring='recall')
        # print('Recall', recall.mean())
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', precision.mean())
        
        
        # y_pred_rs = cross_val_predict(rand_search, X, y, cv=5)
        # rs_cm= confusion_matrix(y, y_pred_rs)
        # print(rs_cm)
        # cl_rs=classification_report(y, y_pred_rs)
        # print(cl_rs)
        # ac_rs=metrics.accuracy_score(y, y_pred_rs)
        # print("Accuracy: ", ac_rs)
        # re_rs= metrics.recall_score(y, y_pred_rs)
        # print("Recall: ", metrics.recall_score(y, y_pred_rs))
        # prec_rs=metrics.precision_score(y, y_pred_rs)
        # print("Precision: ", prec_rs)
        # false_positive_rate, true_positive_rate, thresholds = roc_curve(y, y_pred_rs)
        # roc_auc_rs = auc(false_positive_rate, true_positive_rate) 
        # print("ROC_AUC:", roc_auc_rs)  
        
        
        y_pred_rs = cross_val_predict(rand_search, X_test, y_test, cv=2)
        rs_cm= confusion_matrix(y_test, y_pred_rs)
        print(rs_cm)
        cl_rs=classification_report(y_test, y_pred_rs)
        print(cl_rs)
        ac_rs=metrics.accuracy_score(y_test, y_pred_rs)
        print("Accuracy: ", ac_rs)
        re_rs= metrics.recall_score(y_test, y_pred_rs)
        print("Recall: ", metrics.recall_score(y_test, y_pred_rs))
        prec_rs=metrics.precision_score(y_test, y_pred_rs)
        print("Precision: ", prec_rs)
        false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_rs)
        roc_auc_rs = auc(false_positive_rate, true_positive_rate) 
        print("ROC_AUC:", roc_auc_rs)  
        
        import lime 
        from lime import lime_tabular
        lime_explainer = lime_tabular.LimeTabularExplainer(
            training_data=np.array(X_train),
            feature_names=X_train.columns,
            class_names=["PreV", "PostV"],
            mode='classification')
        warnings.simplefilter(action='ignore', category=FutureWarning)
        i= 1
        idx=4
        # test_1=X_test.iloc[1,:]
        # print(test_1)
        
        
        exp = lime_explainer.explain_instance(X_test.iloc[i,:].astype(int).values, rand_search.predict_proba, num_features=20)
        exp.show_in_notebook()
        exp.save_to_file('ANN_LIME_'+ str(i)+'.html')
        # exp.show_in_notebook(show_table=True)
        varImp1= exp.as_map()
        fig = exp.as_pyplot_figure()
        fig.savefig('lime_report_ANN_'+str(i)+'.jpg')

        exp = lime_explainer.explain_instance(X_test.iloc[idx,:].astype(int).values, rand_search.predict_proba, num_features=20)
        exp.show_in_notebook()
        exp.save_to_file('ANN_LIME_'+ str(idx)+'.html')
        # exp.show_in_notebook(show_table=True)
        varImp2= exp.as_map()
        fig = exp.as_pyplot_figure()
        fig.savefig('lime_report_ANN_'+str(idx)+'.jpg')
    
         
        # def get_feature_importance(j, n):
        #     s = accuracy_score(y_test, y_pred_rs) # baseline score
        #     total = 0.0
        #     for i in range(n):
        #       perm = np.random.permutation(range(X_test.shape[0]))
        #       X_test_ = X_test.copy()
        #       X_test_.iloc[:, j] = X_test.iloc[perm, j]
        #       y_pred_ = rand_search.best_estimator_.predict(X_test_)
        #       s_ij = accuracy_score(y_test, y_pred_)
        #       total += s_ij
        #     return s - total / n
     
    
        # # Feature importances
        # f = []
        # for j in range(X_test.shape[1]):
        #   f_j = get_feature_importance(j, 100)
        #   f.append(f_j)
            
        # plt.figure(figsize=(10, 5))
        # plt.bar(range(X_test.shape[1]), f, color="r", alpha=0.7)
        # plt.xticks(ticks=range(X_test.shape[1]))
        # plt.xlabel("Feature")
        # plt.ylabel("Importance")
        # plt.title("Feature importances ANN (CBI)")
        # plt.show()
    
 
    
        plt.figure()
        plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_rs)
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Random Search')
        plt.legend(loc="lower right")
        
        plt.savefig('ROC_ANN_randomsearch.png')
        plt.show()
    
        
        file= open("ANN_CBI.txt", "w")
        file.write('ANN (Concatenation_Based Integration results: \n')
        file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
        file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
        file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
        file.write("\n")
        # file.write("Grid Search: \n")
        # file.write(f"Best Estimator: {opt_model_d.best_estimator_}\n")
        # file.write(f'Accuracy: \n {metrics.accuracy_score(y, y_pred)} \n')
        # file.write(f'Precision: \n {metrics.precision_score(y, y_pred)} \n')
        # file.write(f'Recall: \n { metrics.recall_score(y, y_pred)} \n')
        # file.write(f'Confusion Matrix: \n {confusion_matrix(y, y_pred)} \n')
        # file.write(f'Classification report: \n {classification_report(y, y_pred)} \n')
        # file.write(f'ROC_AUC: \n {roc_auc} \n')
        # file.write("\n")        
        file.write("Random Search: \n")
        file.write(f"Best Estimator: {rand_search.best_estimator_}\n")
        file.write(f'Accuracy: \n {ac_rs} \n')
        file.write(f'Precision: \n {prec_rs} \n')
        file.write(f'Recall: \n { re_rs} \n')
        file.write(f'Confusion Matrix: \n {rs_cm} \n')
        file.write(f'Classification report: \n {cl_rs} \n')
        # file.write(f'Feature Importance for two samples: \n {varImp1} \n')
        file.write(f'ROC_AUC: \n {roc_auc_rs} \n')
        file.write("See ROC_ANN_gridsearch.png and ROC_ANN_randomsearch for the ROC curves and 'lime_report_ANN.png' and ANN_LIME_.html for the feature importance plot. ")
        
        
        file.close()


    
    
    
    
    
    
    def RF(self):
        omics_=[]
        names=[]
        for name, lista in self.train_test_datasets.items():
            if name in self.omics:
                # confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        if len(omics_) ==2:
            omics1= omics_[0]
            omics2= omics_[1]
               
            print("Dimension of Dataset Concatenated:")
            dataset_concat= pd.read_excel("dataset_concat.xlsx", index_col=0) 
            print(dataset_concat.shape)
            
            
            
            dim1=base.dim(omics1)
            # print(dim1)
            range_omic1= np.arange(1, dim1[1], 1)
            
            
            
            dim2=base.dim(omics2)
            # print(dim2)
            dim_concat=base.dim(dataset_concat)
            len_concat=dim_concat[1]
            final_2= dim_concat[1]
            range_omic2= np.arange(dim1[1],final_2, 1)
            # omic_2=dataset_concat.rx(True,range_omic2)
            # print(utils.head(omic_2))
       
        # omic1=omics[0]
        # omic2=omics[1]
        omics1=dataset_concat.iloc[:,np.arange(0, dim1[1]-1, 1)]
        print("Dimension of "+ names[0] + ":")
        print(omics1.shape)
        
        print("Dimension of "+ names[1]+ ":")
        omics2=dataset_concat.iloc[:,np.arange(dim1[1]-1,final_2-1, 1)]
        print(omics2.shape)
        # Y_nns= dataset_concat.loc[:, self.y_pred]

        X= dataset_concat.drop(self.y_pred, axis=1)
        print(X.shape)

        y= dataset_concat[self.y_pred]
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)
        
        # rf = RandomForestClassifier(n_estimators=100, oob_score=True, random_state=123456)
        # rf.fit(X_train, y_train)
        
        
        # predicted = rf.predict(X_test)
        # accuracy = accuracy_score(y_test, predicted)
        # print(f'Out-of-bag score estimate: {rf.oob_score_:.3}')
        # print(f'Mean accuracy score: {accuracy:.3}')
         
        # kf = KFold(n_splits=5)
        
        # result = cross_val_score(rf , X, y, cv = kf)
        # print("Avg accuracy: {}".format(result.mean()))
        
        # file= open("RFs_CBI.txt", "w")
        # file.write('NNs (Concatenation_Based Integration results: \n')
        # file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
        # file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
        # file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
        # file.write(f'Out-of-bag score estimate: {rf.oob_score_:.3} \n')
        # file.write(f'Mean accuracy score: {accuracy:.3} \n')
        # file.write(f'With Cross Validation to prevent overfitting: \n {"Avg accuracy: {}".format(result.mean())} \n')
        # file.close()
        
        # return print("RFs excecuted successfully!")
        
        
       #  from sklearn.model_selection import GridSearchCV

       #  print("Grid Search: \n")
       #  parameters={"max_depth": [2, 3, None], "max_features": [2,4,6], "min_samples_split": [2,4,6],
       #        "min_samples_leaf": [2,4,6], "bootstrap": [True, False], "criterion": ["gini", "entropy"]}
       #  rf_model_d = RandomForestClassifier( )
        
       #  opt_model_d = GridSearchCV(rf_model_d, parameters, cv = 5)
       #  opt_model_d.fit(X, y)
       #  print (opt_model_d.best_estimator_)
       #  scores_gs = cross_val_score(opt_model_d, X, y, cv = 5)
       #  print("CV average score: %.2f" % scores_gs.mean())
        
       # # precision, recall and F1
       #  from sklearn.preprocessing import LabelBinarizer
        
       #  lb = LabelBinarizer()
       #  y = np.array([number[0] for number in lb.fit_transform(y)])
        
       #  # recall = cross_val_score(svm_model_d, X_train, y_train, cv=5, scoring='recall')
       #  # print('Recall', np.mean(recall), recall)
       #  # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
       #  # print('Precision', np.mean(precision), precision)
       
       #  # print(scores["precision"])
       #  # lb = preprocessing.LabelBinarizer()
       #  # lb.fit(y)
       #  # recall = cross_val_score(svm_model_d, X, y, cv=5, scoring='recall')
       #  # print('Recall', recall.mean())
       #  # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
       #  # print('Precision', precision.mean())
        
        
       #  y_pred = cross_val_predict(opt_model_d, X, y, cv=5)
       #  print(confusion_matrix(y, y_pred))
       #  print(classification_report(y, y_pred))
       #  print("Accuracy: ",metrics.accuracy_score(y, y_pred))
       #  print("Recall: ", metrics.recall_score(y, y_pred))
       #  print("Precision: ", metrics.precision_score(y, y_pred))
       #  false_positive_rate, true_positive_rate, thresholds = roc_curve(y, y_pred)
       #  roc_auc = auc(false_positive_rate, true_positive_rate) 
       #  print("ROC_AUC:", roc_auc)  
        
       #  plt.figure()
       #  plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
       #  plt.plot([0, 1], [0, 1], 'k--')
       #  plt.xlim([0.0, 1.0])
       #  plt.ylim([0.0, 1.05])
       #  plt.xlabel('False Positive Rate')
       #  plt.ylabel('True Positive Rate')
       #  plt.title('ROC Grid Search')
       #  plt.legend(loc="lower right")
       #  plt.show()
       #  plt.savefig('ROC_RF_gridsearch.png')
        
        print("Random Search: \n")
        
        def report(results, n_top=3):
            for i in range(1, n_top + 1):
                candidates = np.flatnonzero(results['rank_test_score'] == i)
                for candidate in candidates:
                    print("Model with rank: {0}".format(i))
                    print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                          results['mean_test_score'][candidate],
                          results['std_test_score'][candidate]))
                    print("Parameters: {0}".format(results['params'][candidate]))
                    print("")
        
        

        rf_model = RandomForestClassifier( )
        
        parameters={"max_depth": [2, 3, None], "max_features": [2,4,6], "min_samples_split": [2,4,6],
              "min_samples_leaf": [2,4,6], "bootstrap": [True, False], "criterion": ["gini", "entropy"]}
       
        rand_search = RandomizedSearchCV(rf_model, param_distributions=parameters, cv = 2)
        
        rand_search.fit(X, y)
        
        print (rand_search.best_estimator_)
        report(rand_search.cv_results_)
        print(rand_search.best_estimator_.feature_importances_)
        varImp=rand_search.best_estimator_.feature_importances_
        
        headers = ["name", "score"]
        values = sorted(zip(X_train.columns, varImp), key=lambda x: x[1] * -1)
        print(tabulate(values, headers, tablefmt="plain"))
        importance= tabulate(values, headers, tablefmt="plain")
        
        # sorted_idx = varImp.argsort()
        # plt.barh(dataset_concat.columns[sorted_idx], rand_search.best_estimator_.feature_importances_[sorted_idx])
        # plt.xlabel("Random Forest Feature Importance")
        # plt.savefig('RF_varImp.png')
        
        # indices = np.argsort(varImp)

        # # customized number 
        # num_features = 10 
        
        # plt.figure(figsize=(10,100))
        # plt.title('Feature Importances')
        
        # # only plot the customized number of features
        # plt.barh(range(num_features), varImp[indices[-num_features:]], color='b', align='center')
        # plt.yticks(range(num_features), [dataset_concat.columns[i] for i in indices[-num_features:]])
        # plt.xlabel('Relative Importance')
        # plt.show()
                
        
        
        scores_rs = cross_val_score(rand_search, X, y, cv = 2)
        print("CV average score: %.2f" % scores_rs.mean())
        
        
        from sklearn.preprocessing import LabelBinarizer
        
        lb = LabelBinarizer()
        y = np.array([number[0] for number in lb.fit_transform(y)])
        y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
        # recall = cross_val_score(svm_model_d, X_train, y_train, cv=5, scoring='recall')
        # print('Recall', np.mean(recall), recall)
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', np.mean(precision), precision)
       
        # print(scores["precision"])
        # lb = preprocessing.LabelBinarizer()
        # lb.fit(y)
        # recall = cross_val_score(svm_model_d, X, y, cv=5, scoring='recall')
        # print('Recall', recall.mean())
        # precision = cross_val_score(svm_model_d, X, y, cv=5, scoring='precision')
        # print('Precision', precision.mean())
        
        
        # y_pred_rs = cross_val_predict(rand_search, X, y, cv=5)
        # rs_cm= confusion_matrix(y, y_pred_rs)
        # print(rs_cm)
        # cl_rs=classification_report(y, y_pred_rs)
        # print(cl_rs)
        # ac_rs=metrics.accuracy_score(y, y_pred_rs)
        # print("Accuracy: ", ac_rs)
        # re_rs= metrics.recall_score(y, y_pred_rs)
        # print("Recall: ", metrics.recall_score(y, y_pred_rs))
        # prec_rs=metrics.precision_score(y, y_pred_rs)
        # print("Precision: ", prec_rs)
        # false_positive_rate, true_positive_rate, thresholds = roc_curve(y, y_pred_rs)
        # roc_auc_rs = auc(false_positive_rate, true_positive_rate) 
        # print("ROC_AUC:", roc_auc_rs)  
        
        
        
        y_pred_rs = cross_val_predict(rand_search, X_test, y_test, cv=2)
        rs_cm= confusion_matrix(y_test, y_pred_rs)
        print(rs_cm)
        cl_rs=classification_report(y_test, y_pred_rs)
        print(cl_rs)
        ac_rs=metrics.accuracy_score(y_test, y_pred_rs)
        print("Accuracy: ", ac_rs)
        re_rs= metrics.recall_score(y_test, y_pred_rs)
        print("Recall: ", metrics.recall_score(y_test, y_pred_rs))
        prec_rs=metrics.precision_score(y_test, y_pred_rs)
        print("Precision: ", prec_rs)
        false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_rs)
        roc_auc_rs = auc(false_positive_rate, true_positive_rate) 
        print("ROC_AUC:", roc_auc_rs)  
        
        
        
        
        
        
        plt.figure()
        plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_rs)
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Random Search')
        plt.legend(loc="lower right")
        plt.savefig('ROC_RF_randomsearch.png')
        plt.show()
        
      

        
        file= open("RF_CBI.txt", "w")
        file.write('SVM (Concatenation_Based Integration results: \n')
        file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
        file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
        file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
        file.write("\n")
        # file.write("Grid Search: \n")
        # file.write(f"Best Estimator: {opt_model_d.best_estimator_}\n")
        # file.write(f'Accuracy: \n {metrics.accuracy_score(y, y_pred)} \n')
        # file.write(f'Precision: \n {metrics.precision_score(y, y_pred)} \n')
        # file.write(f'Recall: \n { metrics.recall_score(y, y_pred)} \n')
        # file.write(f'Confusion Matrix: \n {confusion_matrix(y, y_pred)} \n')
        # file.write(f'Classification report: \n {classification_report(y, y_pred)} \n')
        # file.write(f'ROC_AUC: \n {roc_auc} \n')
        # file.write("\n")        
        file.write("Random Search: \n")
        file.write(f"Best Estimator: {rand_search.best_estimator_}\n")
        file.write(f'Accuracy: \n {ac_rs} \n')
        file.write(f'Precision: \n {prec_rs} \n')
        file.write(f'Recall: \n { re_rs} \n')
        file.write(f'Confusion Matrix: \n {rs_cm} \n')
        file.write(f'Classification report: \n {cl_rs} \n')
        file.write(f'ROC_AUC: \n {roc_auc_rs} \n')
        file.write(f'Feature Importance: \n {importance} \n')
        file.write("See ROC_RF_gridsearch.png and ROC_SVM_randomsearch for the ROC curves.")
        
        
        file.close()
        
        
    ''' Transformation - Based Integration (TBI)'''    
        
    def SNFtool (self):
        
        robjects.r(''' run_packages <- function() {
            require(devtools)
            install_version("heatmap.plus", version = "1.3", repos = "http://cran.us.r-project.org")
            library(heatmap.plus)
            install_version("SNFtool", version = "2.3.0", repos = "http://cran.us.r-project.org")
            library(SNFtool) 
            install.packages("pROC")
            library(pROC)
            return (print("Package SNFtool installed!"))}''')
            
        #Function Dist2   
        robjects.r('''dist2 <- function(X,C) {
                # Calculates the squared euclidean distance between two matrices where rows 
                #   represent a single data point or patient.
                #
                # Args:
                #   X: Matrix with each row representing a single data point (or patient)
                #   C: Matrix with each row representing a single data point (or patient)
                #
                # Returns:
                #   res: A NxM matrix where nrow(X) == N and nrow(C) == M. Element [n,m] 
                #       is the squared euclidean distance between rows N[n,] and C[m,].
            
                ndata <- nrow(X)
                ncentres <- nrow(C)
                
                sumsqX <- rowSums(X^2)
                sumsqC <- rowSums(C^2)
                  
                XC <- 2 * (X %*% t(C))
                res <- matrix(rep(sumsqX, times=ncentres), ndata, ncentres) + 
                    t(matrix(rep(sumsqC, times=ndata), ncentres, ndata)) - XC
                res[res < 0] <- 0
            
                return(res)
            }  ''')   
            
                
                
        robjects.r('''run_snftool  <- function(K, alpha, T, omic1, omic2, Y){
            
            ## First, set all the parameters:
            K = K		# number of neighbors, usually (10~30)
            alpha = alpha  	# hyperparameter, usually (0.3~0.8)
            T = T	# Number of Iterations, usually (10~20)
            
            ## Data1 is of size n x d_1, 
            ## where n is the number of patients, d_1 is the number of genes, 
            ## Data2 is of size n x d_2, 
            ## where n is the number of patients, d_2 is the number of methylation
            omic1_= cbind(omic1, Y)
            omic2_= cbind(omic2,Y)
            
            print(dim(omic1))
            print(dim(omic2))
            
            ## Here, the simulation data (SNFdata) has two data types. They are complementary to each other. 
            ## And two data types have the same number of points. 
            ## The first half data belongs to the first cluster; the rest belongs to the second cluster.
            
            factors <- factor(Y)
            print(factors)
            # Convert the factor to numbers:
            factors=as.numeric(factors)
            develop=as.matrix(factors)
            
            
            truelabel = c(develop,develop) ## the ground truth of the simulated data
            
            ## Calculate distance matrices
            ## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
            
            ## If the data are all continuous values, we recommend the users to perform 
            ## standard normalization before using SNF, 
            ## though it is optional depending on the data the users want to use.  
            # Data1 = standardNormalization(Data1);
            # Data2 = standardNormalization(Data2);
            
            
            
            ## Calculate the pair-wise distance; 
            ## If the data is continuous, we recommend to use the function "dist2" as follows 
            
            Dist1 = (dist2(as.matrix(omic1),as.matrix(omic1)))^(1/2)
            Dist2 = (dist2(as.matrix(omic2),as.matrix(omic2)))^(1/2)
            
            ## next, construct similarity graphs
            W1 = affinityMatrix(Dist1, K, alpha)
            W2 = affinityMatrix(Dist2, K, alpha)
            
            ## These similarity graphs have complementary information about clusters.
            
            png(file="displayCl1.png", width=512, height=512)
            displayClusters(W1,develop)
            dev.off()
            
            png(file="displayCl2.png", width=512, height=512)
            displayClusters(W2,develop)
            dev.off()
            
            
            ## next, we fuse all the graphs
            ## then the overall matrix can be computed by similarity network fusion(SNF):
            W = SNF(list(W1,W2), K, T)
            
            NMI_scores <- rankFeaturesByNMI(list(omic1, omic2), W)

            print(NMI_scores)
            
            ## With this unified graph W of size n x n, 
            ## you can do either spectral clustering or Kernel NMF. 
            ## If you need help with further clustering, please let us know. 
            
            ## You can display clusters in the data by the following function
            ## where C is the number of clusters.
            C = 2 								# number of clusters
            group = spectralClustering(W,C) 	# the final subtypes information
            png(file="displayClusters.png", width=512, height=512)
            displayClusters(W, group)
            dev.off()
            NMI = calNMI(group,develop)
            
            ## You can get cluster labels for each data point by spectral clustering
            labels = spectralClustering(W, C)
            
            # dev.off()
            # png(file="plots.png", width=512, height=512)
            # plot(omic1, col=labels, main='Data type 1')
            # plot(omic2, col=labels, main='Data type 2')
            # }
            
            ################################################################################
            # Provide an example of predicting the new labels with label propagation
            
            # Load views into list "dataL" and the cluster assignment into vector "label"
             #data(Digits)
            
            
            dataL=list(omic1,omic2)
            # Create the training and test data
            n = floor(0.8*length(develop)) # number of training cases
            trainSample = sample.int(length(develop), n)
            train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
            test = lapply(dataL, function(x) x[-trainSample, ]) # Test the rest of the data set
            groups = develop[trainSample]
            
            # Set the other
            # K = 20
            # alpha = 0.5
            # t = 20
            method = TRUE
            
            # Apply the prediction function to the data
            newLabel = groupPredict(train,test,groups,K,alpha,T,method)
            print(newLabel)
            # Compare the prediction accuracy
            accuracy = sum(develop[-trainSample] == newLabel[-c(1:n)])/(length(develop) - n)
            print("Accuracy:")
            print( accuracy)
            
            a=confusionMatrix(table(newLabel[-c(1:n)], develop[-trainSample]))
            recall=a$byClass["Recall"]
            precision=a$byClass["Precision"]
            
            print(a)
            print(recall)            
            print(precision)
            
            png(file="SNFtool_ROC.png", width=512, height=512)
            roc(response = develop[-trainSample], predictor =newLabel[-c(1:n)] ,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
            title(main = "SNFtool")
            dev.off()
            
            SNFtool_res= list(develop[-trainSample], newLabel[-c(1:n)], accuracy, a, recall, precision)
            return (SNFtool_res)}''')
                
                
                
                
        # robjects.r('''run_snftool_CS2  <- function(K, alpha, T, omic1, omic2, Y){
        #     install.packages("pROC")
        #     library(pROC)
        #     ## First, set all the parameters:
        #     K = K		# number of neighbors, usually (10~30)
        #     alpha = alpha  	# hyperparameter, usually (0.3~0.8)
        #     T = T	# Number of Iterations, usually (10~20)
            
        #     ## Data1 is of size n x d_1, 
        #     ## where n is the number of patients, d_1 is the number of genes, 
        #     ## Data2 is of size n x d_2, 
        #     ## where n is the number of patients, d_2 is the number of methylation
        #     omic1_= cbind(omic1, Y)
        #     omic2_= cbind(omic2,Y)
            
        #     print(dim(omic1))
        #     print(dim(omic2))
            
        #     ## Here, the simulation data (SNFdata) has two data types. They are complementary to each other. 
        #     ## And two data types have the same number of points. 
        #     ## The first half data belongs to the first cluster; the rest belongs to the second cluster.
            
        #     factors <- factor(Y)
        #     print(factors)
        #     # Convert the factor to numbers:
        #     factors=as.numeric(factors)
        #     develop=as.matrix(factors)
            
            
        #     truelabel = c(develop,develop) ## the ground truth of the simulated data
            
        #     ## Calculate distance matrices
        #     ## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
            
        #     ## If the data are all continuous values, we recommend the users to perform 
        #     ## standard normalization before using SNF, 
        #     ## though it is optional depending on the data the users want to use.  
        #     # Data1 = standardNormalization(Data1);
        #     # Data2 = standardNormalization(Data2);
            
            
            
        #     ## Calculate the pair-wise distance; 
        #     ## If the data is continuous, we recommend to use the function "dist2" as follows 
            
        #     Dist1 = (dist2(as.matrix(omic1),as.matrix(omic1)))^(1/2)
        #     Dist2 = (dist2(as.matrix(omic2),as.matrix(omic2)))^(1/2)
            
        #     ## next, construct similarity graphs
        #     W1 = affinityMatrix(Dist1, K, alpha)
        #     W2 = affinityMatrix(Dist2, K, alpha)
            
        #     ## These similarity graphs have complementary information about clusters.
            
        #     png(file="displayCl1.png", width=512, height=512)
        #     displayClusters(W1,develop)
        #     dev.off()
            
        #     png(file="displayCl2.png", width=512, height=512)
        #     displayClusters(W2,develop)
        #     dev.off()
            
            
        #     ## next, we fuse all the graphs
        #     ## then the overall matrix can be computed by similarity network fusion(SNF):
        #     W = SNF(list(W1,W2), K, T)
            
            
        #     NMI_scores <- rankFeaturesByNMI(list(omic1, omic2), W)

        #     print(NMI_scores)
            
        #     ## With this unified graph W of size n x n, 
        #     ## you can do either spectral clustering or Kernel NMF. 
        #     ## If you need help with further clustering, please let us know. 
            
        #     ## You can display clusters in the data by the following function
        #     ## where C is the number of clusters.
        #     C = 2 								# number of clusters
        #     group = spectralClustering(W,C) 	# the final subtypes information
        #     png(file="displayClusters.png", width=512, height=512)
        #     displayClusters(W, group)
        #     dev.off()
        #     NMI = calNMI(group,develop)
            
        #     ## You can get cluster labels for each data point by spectral clustering
        #     labels = spectralClustering(W, C)
            
        #     # dev.off()
        #     # png(file="plots.png", width=512, height=512)
        #     # plot(omic1, col=labels, main='Data type 1')
        #     # plot(omic2, col=labels, main='Data type 2')
        #     # }
            
        #     ################################################################################
        #     # Provide an example of predicting the new labels with label propagation
            
        #     # Load views into list "dataL" and the cluster assignment into vector "label"
        #      #data(Digits)
            
            
        #     dataL=list(omic1,omic2)
        #     # Create the training and test data
        #     n = floor(1*length(develop)) # number of training cases
        #     trainSample = sample.int(length(develop), n)
        #     train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
        #     test = lapply(dataL, function(x) x[trainSample, ]) # Test the rest of the data set
        #     groups = develop[trainSample]
            
        #     # Set the other
        #     # K = 20
        #     # alpha = 0.5
        #     # t = 20
        #     method = TRUE
            
        #     # Apply the prediction function to the data
        #     newLabel = groupPredict(train,test,groups,K,alpha,T,method)
        #     print(newLabel)
        #     # Compare the prediction accuracy
        #     accuracy = sum(develop[trainSample] == newLabel[c(1:n)])/(length(develop) - n)
        #     print("Accuracy:")
        #     print( accuracy)
            
        #     a=confusionMatrix(table(newLabel[c(1:n)], develop[trainSample]))
        #     recall=a$byClass["Recall"]
        #     precision=a$byClass["Precision"]
            
        #     print(a)
        #     print(recall)            
        #     print(precision)
            
        #     png(file="SNFtool_ROC.png", width=512, height=512)
        #     roc(response = develop[trainSample], predictor =newLabel[c(1:n)] ,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
        #     title(main = "SNFtool")
        #     dev.off()
            
        #     SNFtool_res= list(develop[trainSample], newLabel[c(1:n)], accuracy, a, recall, precision)
        #     return (SNFtool_res)}''')
        run_packages= robjects.r["run_packages"]
        run_packages()
        
        run_snftool= robjects.r["run_snftool"]
        
        print(utils.head(self.omic_1))
        print(utils.head(self.omic_2))
        print(self.Y)
        Y= robjects.FactorVector(self.Y)
        K = 20		# number of neighbors, usually (10~30)
        alpha = 0.5  	# hyperparameter, usually (0.3~0.8)
        T = 20 	# Number of Iterations, usually (10~20)
        
        
        snftool_res= run_snftool(K, alpha, T, self.omic_1, self.omic_2,Y)
        trueLabel= snftool_res[0]
        newLabel= snftool_res[1]
        accuracy=snftool_res[2]
        cm=snftool_res[3]
        recall=snftool_res[4]
        precision=snftool_res[5]
        file= open("SNFtool_TBI.txt", "w")
        file.write('SNFtool results: \n')
        file.write(f'Dimension of dataset {self.omics[0]}: {base.dim(self.omic_1)} \n')
        file.write(f'Dimension of dataset {self.omics[1]}: {base.dim(self.omic_2)} \n')
        file.write(f'True Label:\n {trueLabel} \n')
        file.write(f'New Label:\n {newLabel} \n')
        file.write(f'Accuracy: {accuracy} \n')
        file.write(f'Confusion Matrix: {cm} \n')
        file.write(f'Recall: {recall} \n')
        file.write(f'Precision: {precision} \n')
        file.close()
        
        return print("SNFtool excecuted successfully!")
    
    
    
    def CAN_TBI(self):
        robjects.r(''' predict.CAN <-function(object,...,score=FALSE){
          LMatrix<-list(...)
          index<-length(LMatrix)
          if(index!=length(object@Weights)){
            stop("The number of associate networks is inconsistent")
          }
        
          for(i in 1: index){
            M<-LMatrix[[i]]
            if(unique(diag(M)!=0)){
              diag(LMatrix[[i]])<-0
              warning(paste0("The diagonal elements of ",i,"th Weight Matrix should be all Zero"))
            }
            if(ncol(M)!=nrow(M)){
              Warn<-paste0("The dimension of the ",i,"th Weight Matrix doesn't match")
              stop(Warn)
            }
          }
        
          if(object@Numeric){
            Status<-c(object@Label,rep('unknown',(ncol(M)-length(object@Label))))
            LEVEL<-as.numeric(object@Level)
          }else{
            Status<-c(as.character(object@Label),rep('unknown',(ncol(M)-length(object@Label))))
            LEVEL<-object@Level
          }
        
          Status[Status!=LEVEL[1]&Status!=LEVEL[2]]<-0
          Status[Status==LEVEL[1]]<--1
          Status[Status==LEVEL[2]]<-1
          Status<-as.numeric(Status)
        
          L<-matrix(0,nrow=length(Status),ncol=length(Status))
          for(i in 1:index){
            L<-L+LMatrix[[i]]*object@Weights[i]
          }
          E<-diag(nrow(L))
          L<-E+diag(rowSums(L))-L
          Score<-solve(L)%*%Status
        
          Pos<-median(Score[Status==1])
          Neg<-median(Score[Status==-1])
          Result<-ifelse(abs(Score-Pos)<=abs(Score-Neg),1,-1)
        
          Prediction<-Result[Status==0]
          Prediction<-ifelse(Prediction==-1,LEVEL[1],LEVEL[2])
          if(!score){
            return(Prediction)
          }else{
            H<-data.frame(Prediction=Prediction,Score=Score[Status==0])
          }
        
        }''')
        
        
        robjects.r('''predict.crvm <-function(object,newdata,prob=FALSE){
          P<-length(object@Level)-1
          if(is.matrix(newdata)&P>1){
            stop("The test data should be a list of kernel matrices")
          }
        
          if(is.list(newdata)&length(newdata)!=P){
            stop("The number of kernel matrices in newdata should equla to number of classer -1")
          }
        
          if(is.matrix(newdata)){
            N<-nrow(newdata)
          } else{
            N<-nrow(newdata[[1]])
          }
        
          KernMat<-newdata
          if(object@bias){
            for(p in 1:P){
              KernMat[[p]]<-cbind(rep(1,nrow(KernMat[[p]])),KernMat[[p]])
            }
          }
        
          probabilities<-multinomial(KernMat,object@Weights,N,P)
          probabilities<-cbind(probabilities,(1-rowSums(probabilities)))
        
          Prediction<-apply(probabilities,1,which.max)
          if(is.numeric(object@Prediction)){
            for(i in 1:(P+1)){
              Prediction[Prediction==i]<-as.numeric(object@Level)[i]
            }
          }else{
            for(i in 1:(P+1)){
              Prediction[Prediction==i]<-object@Level[i]
            }
          }
        
          if(prob){
            colnames(probabilities)<-object@Level
            return(probabilities)
          }else{
            return(Prediction)
          }
        }''')
        
        robjects.r(''' predict.Ada <-function(object,newdata,prob=FALSE){
          P<-length(object@Level)-1
          if(is.matrix(newdata)&P>1){
            stop("The test data should be a list of kernel matrices")
          }
        
          if(is.list(newdata)&length(newdata)!=P){
            stop("The number of kernel matrices in newdata should equla to number of classer -1")
          }
        
          if(is.matrix(newdata)){
            N<-nrow(newdata)
          } else{
            N<-nrow(newdata[[1]])
          }
        
          KernMat<-newdata
          if(object@bias){
            for(p in 1:P){
              KernMat[[p]]<-cbind(rep(1,nrow(KernMat[[p]])),KernMat[[p]])
            }
          }
        
          probabilities<-multinomial(KernMat,object@Weights,N,P)
          probabilities<-cbind(probabilities,(1-rowSums(probabilities)))
        
          Prediction<-apply(probabilities,1,which.max)
          if(is.numeric(object@Prediction)){
            for(i in 1:(P+1)){
              Prediction[Prediction==i]<-as.numeric(object@Level)[i]
            }
          }else{
            for(i in 1:(P+1)){
              Prediction[Prediction==i]<-object@Level[i]
            }
          }
        
          if(prob){
            colnames(probabilities)<-object@Level
            return(probabilities)
          }else{
            return(Prediction)
          }
        }''')
        
        
        
        robjects.r('''predict.MultR <-function(object,...,prob=FALSE){
        
          tmpdata<-list(...)
          Index<-length(tmpdata)
        
          newdata<-vector("list",Index)
          for(i in 1:Index){
            tmp<-vector("list",1)
            tmp[[1]]<-tmpdata[[i]]
            newdata[[i]]<-tmp
          }
          N<-length(object@Model)
        
          if(N!=Index){
            stop("The number of data sources doesn't match")
          }
        
          P<-length(object@Level)-1
          dimension<-c()
          for(i in 1:Index){
            if(length(newdata[[i]])!=P){
              stop(paste0("The number of kernel matrices in ",i,"th newdata should equla to ",P))
            }
            for(j in 1:P){
              dimension<-c(dimension,nrow(newdata[[i]][[j]]))
            }
          }
        
          if(length(unique(dimension))!=1){
            stop("The dimension of the test data is not consistent")
          }
        
          D<-unique(dimension)
        
          if(object@Para@Bias){
            for(i in 1:Index){
              for(p in 1:P){
                newdata[[i]][[p]]<-cbind(rep(1,nrow(KernMat[[p]])),KernMat[[p]])
              }
            }
          }
        
          for(i in 1:Index){
            if(i==1){
              probabilities<-multinomial(newdata[[i]],object@Model[[i]]@Weights,D,P)
              probabilities<-cbind(probabilities,(1-rowSums(probabilities)))
            }else{
              tmp<-multinomial(newdata[[i]],object@Model[[i]]@Weights,D,P)
              tmp<-cbind(tmp,(1-rowSums(tmp)))
              probabilities<-probabilities+tmp
            }
          }
          probabilities<-probabilities/Index
        
          Prediction<-apply(probabilities,1,which.max)
        
          if(is.numeric(object@Prediction)){
            for(i in 1:(P+1)){
              Prediction[Prediction==i]<-as.numeric(object@Level)[i]
            }
          }else{
            for(i in 1:(P+1)){
              Prediction[Prediction==i]<-object@Level[i]
            }
          }
        
          if(prob){
            colnames(probabilities)<-object@Level
            return(probabilities)
          }else{
            return(Prediction)
          }
        }''')
              
        robjects.r(''' packages <- function() {
            Required    <- c("caret","signal","MASS","kernlab","pROC")
            New.Package <- ! Required %in% installed.packages()[,"Package"]
            if ( sum(New.Package)!=0){
                install.packages(Required[New.Package])}

            install.packages('MDIntegration_1.0.tar.gz', repos = NULL, type = "source") 
            library(MDIntegration)
        
            return(print("MDIntegration package installed successfully!"))}
                     ''')
                     
                     
        packages= robjects.r["packages"]
        packages()
        
        ''' Composite Association Network'''
        
        robjects.r(''' run_CAN <- function(omic1, omic2, Y){
            library(caret)
            # library("MDIntegration")
            
            install.packages('lime')
            library(lime)
            
            W_omic1= cor(t(omic1),method="pearson")
            W_omic2=cor(t(omic2),method="pearson")
            # Y= as.matrix(Y)
            Label= factor(Y)
            print(Label)
            set.seed(1024)
            Train       <- createDataPartition(Label, p = 0.75, list = FALSE)
            Test        <- setdiff(1:length(Label),Train)
            
            ### Training Data 
            Train_omic1 <- W_omic1[Train,Train]
            Train_omic2 <- W_omic2[Train,Train]
            Train_Label <- Label[Train]
            
            #Training the Model
            Model       <- CANetwork(Train_omic1,Train_omic2,Status = Train_Label)
            
            #Prediction with CANetwork Model
            #The data should be ordered before get the prediction. Otherwise,you may get the wrong results.
            
            Index       <- c(Train,Test)
            Prediction  <- predict(Model,W_omic1[Index,Index],W_omic2[Index,Index],score = FALSE)
            
            #Get the Contingency Table
            Cont_table=table(Label[Test],Prediction)
            
            library(caret)
            a=confusionMatrix(table(Label[Test],Prediction))
            print(a)
            recall=a$byClass["Recall"]
            precision=a$byClass["Precision"]
            
            # print(table(Label[Test],Prediction))
            pecc=function(obs,pred) sum(obs==pred)/length(pred)
            
                     
            #Plot AUC Curve
            library(pROC)
            Prediction.Score  <- predict(Model,W_omic1[Index,Index],W_omic2[Index,Index],score = TRUE)
            # print(Prediction.Score$Score)
            print(pecc(Prediction,Label[Test]))
            
            png(file="CAN_res_TBI.png", width=512, height=512)
            roc(response = Label[Test], predictor = Prediction.Score$Score,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
            title(main = "Composite Association Network")
            dev.off()
            
            res_CAN= list(a, pecc(Prediction,Label[Test]), precision, recall)
            
            return(res_CAN)}''')
            
            
        robjects.r(''' run_CAN_CS2 <- function(omic1, omic2, Y){
            library(caret)
            # library("MDIntegration")
            W_omic1= cor(t(omic1),method="pearson")
            W_omic2=cor(t(omic2),method="pearson")
            # Y= as.matrix(Y)
            Label= factor(Y)
            print(Label)
            set.seed(1024)
            Train       <- 1:length(Label)
            Test        <- 1:length(Label)
            
            ### Training Data 
            Train_omic1 <- W_omic1[Train,Train]
            Train_omic2 <- W_omic2[Train,Train]
            Train_Label <- Label[Train]
            
            #Training the Model
            Model       <- CANetwork(Train_omic1,Train_omic2,Status = Train_Label)
            
            #Prediction with CANetwork Model
            #The data should be ordered before get the prediction. Otherwise,you may get the wrong results.
            
            Index       <- c(Train,Test)
            Prediction  <- predict(Model,W_omic1[Index,Index],W_omic2[Index,Index],score = FALSE)
            
            #Get the Contingency Table
            Cont_table=table(Label[Test],Prediction)
            
            library(caret)
            a=confusionMatrix(table(Label[Test],Prediction))
            print(a)
            recall=a$byClass["Recall"]
            precision=a$byClass["Precision"]
            
            # print(table(Label[Test],Prediction))
            pecc=function(obs,pred) sum(obs==pred)/length(pred)
            
            #Plot AUC Curve
            library(pROC)
            Prediction.Score  <- predict(Model,W_omic1[Index,Index],W_omic2[Index,Index],score = TRUE)
            # print(Prediction.Score$Score)
            print(pecc(Prediction,Label[Test]))
            
            png(file="CAN_res_TBI.png", width=512, height=512)
            roc(response = Label[Test], predictor = Prediction.Score$Score,plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
            title(main = "Composite Association Network")
            dev.off()
            
            res_CAN= list(a, pecc(Prediction,Label[Test]), precision, recall)
            
            return(res_CAN)}''')
     
        print(utils.head(self.omic_1))
        print(utils.head(self.omic_2))
        print(self.Y)
        Y= robjects.FactorVector(self.Y)
        run_CAN= robjects.r["run_CAN"]
        CAN_res= run_CAN(self.omic_1, self.omic_2, Y)
        Cont_table=CAN_res[0]
        PECC=CAN_res[1]
        precision=CAN_res[2]
        recall=CAN_res[3]
        
        file= open("CAN_TBI.txt", "w")
        file.write('CAN results: \n')
        file.write(f'Dimension of dataset {self.omics[0]}: {base.dim(self.omic_1)} \n')
        file.write(f'Dimension of dataset {self.omics[1]}: {base.dim(self.omic_2)} \n')
        file.write(f'Confusion Matrix:\n {Cont_table} \n')
        file.write(f'PECC:\n {PECC} \n')
        file.write(f'Precision:\n {precision} \n')
        file.write(f'Recall:\n {recall} \n')
        file.write(' See image CAN_res_TBI.png for the AUC curve plot. \n')
        
        file.close()
        
        return(print("CAN successfully executed!"))
    
    def RVM_Ada_TBI(self ):
        robjects.r('''setGeneric("predict.crvm",function(object,newdata,prob=FALSE){
                  P<-length(object@Level)-1
                  if(is.matrix(newdata)&P>1){
                    stop("The test data should be a list of kernel matrices")
                  }
                
                  if(is.list(newdata)&length(newdata)!=P){
                    stop("The number of kernel matrices in newdata should equla to number of classer -1")
                  }
                
                  if(is.matrix(newdata)){
                    N<-nrow(newdata)
                  } else{
                    N<-nrow(newdata[[1]])
                  }
                
                  KernMat<-newdata
                  if(object@bias){
                    for(p in 1:P){
                      KernMat[[p]]<-cbind(rep(1,nrow(KernMat[[p]])),KernMat[[p]])
                    }
                  }
                
                  probabilities<-multinomial(KernMat,object@Weights,N,P)
                  probabilities<-cbind(probabilities,(1-rowSums(probabilities)))
                
                  Prediction<-apply(probabilities,1,which.max)
                  if(is.numeric(object@Prediction)){
                    for(i in 1:(P+1)){
                      Prediction[Prediction==i]<-as.numeric(object@Level)[i]
                    }
                  }else{
                    for(i in 1:(P+1)){
                      Prediction[Prediction==i]<-object@Level[i]
                    }
                  }
                
                  if(prob){
                    colnames(probabilities)<-object@Level
                    return(probabilities)
                  }else{
                    return(Prediction)
                  }
                })''')
        
        robjects.r(''' setGeneric("predict.Ada",function(object,newdata,prob=FALSE){
                  P<-length(object@Level)-1
                  if(is.matrix(newdata)&P>1){
                    stop("The test data should be a list of kernel matrices")
                  }
                
                  if(is.list(newdata)&length(newdata)!=P){
                    stop("The number of kernel matrices in newdata should equla to number of classer -1")
                  }
                
                  if(is.matrix(newdata)){
                    N<-nrow(newdata)
                  } else{
                    N<-nrow(newdata[[1]])
                  }
                
                  KernMat<-newdata
                  if(object@bias){
                    for(p in 1:P){
                      KernMat[[p]]<-cbind(rep(1,nrow(KernMat[[p]])),KernMat[[p]])
                    }
                  }
                
                  probabilities<-multinomial(KernMat,object@Weights,N,P)
                  probabilities<-cbind(probabilities,(1-rowSums(probabilities)))
                
                  Prediction<-apply(probabilities,1,which.max)
                  if(is.numeric(object@Prediction)){
                    for(i in 1:(P+1)){
                      Prediction[Prediction==i]<-as.numeric(object@Level)[i]
                    }
                  }else{
                    for(i in 1:(P+1)){
                      Prediction[Prediction==i]<-object@Level[i]
                    }
                  }
                
                  if(prob){
                    colnames(probabilities)<-object@Level
                    return(probabilities)
                  }else{
                    return(Prediction)
                  }
                })
                ''')
            
        
        
        robjects.r('''setGeneric("predict.MultR",function(object,...,prob=FALSE){

              tmpdata<-list(...)
              Index<-length(tmpdata)
            
              newdata<-vector("list",Index)
              for(i in 1:Index){
                tmp<-vector("list",1)
                tmp[[1]]<-tmpdata[[i]]
                newdata[[i]]<-tmp
              }
              N<-length(object@Model)
            
              if(N!=Index){
                stop("The number of data sources doesn't match")
              }
            
              P<-length(object@Level)-1
              dimension<-c()
              for(i in 1:Index){
                if(length(newdata[[i]])!=P){
                  stop(paste0("The number of kernel matrices in ",i,"th newdata should equla to ",P))
                }
                for(j in 1:P){
                  dimension<-c(dimension,nrow(newdata[[i]][[j]]))
                }
              }
            
              if(length(unique(dimension))!=1){
                stop("The dimension of the test data is not consistent")
              }
            
              D<-unique(dimension)
            
              if(object@Para@Bias){
                for(i in 1:Index){
                  for(p in 1:P){
                    newdata[[i]][[p]]<-cbind(rep(1,nrow(KernMat[[p]])),KernMat[[p]])
                  }
                }
              }
            
              for(i in 1:Index){
                if(i==1){
                  probabilities<-multinomial(newdata[[i]],object@Model[[i]]@Weights,D,P)
                  probabilities<-cbind(probabilities,(1-rowSums(probabilities)))
                }else{
                  tmp<-multinomial(newdata[[i]],object@Model[[i]]@Weights,D,P)
                  tmp<-cbind(tmp,(1-rowSums(tmp)))
                  probabilities<-probabilities+tmp
                }
              }
              probabilities<-probabilities/Index
            
              Prediction<-apply(probabilities,1,which.max)
            
              if(is.numeric(object@Prediction)){
                for(i in 1:(P+1)){
                  Prediction[Prediction==i]<-as.numeric(object@Level)[i]
                }
              }else{
                for(i in 1:(P+1)){
                  Prediction[Prediction==i]<-object@Level[i]
                }
              }
            
              if(prob){
                colnames(probabilities)<-object@Level
                return(probabilities)
              }else{
                return(Prediction)
              }
            })
            ''')
                          
        robjects.r(''' packages <- function() {
            Required    <- c("caret","signal","MASS","kernlab","pROC")
            New.Package <- ! Required %in% installed.packages()[,"Package"]
            if ( sum(New.Package)!=0){
                install.packages(Required[New.Package])}

            install.packages('MDIntegration_1.0.tar.gz', repos = NULL, type = "source") 
            library(MDIntegration)
        
            return(print("MDIntegration package installed successfully!"))}
                     ''')

        
        robjects.r(''' run_RVM <- function(omic1, omic2, Y){
            
        
            install.packages("kernlab")
            library(kernlab)
            library(signal)
            library(MASS)
            library(MDIntegration)
            library(caret)
            library(pROC)
            # Y= as.matrix(dataset_concat[,425])
            Label= factor(Y)
            print(Label)
            
            rbf <- rbfdot(sigma = 0.05)
            rbf
            
            ## calculate kernel matrix
            
            
            Kern_omic1= kernelMatrix(rbf, as.matrix(omic1)) 
            Kern_omic2= kernelMatrix(rbf, as.matrix(omic2))
            print(dim(Kern_omic1))
            print(dim(Kern_omic2))
            set.seed(1024)
            Train          <- createDataPartition(Label, p = 0.75, list = FALSE)
            Test           <- setdiff(1:length(Label),Train)
            
            ### Training Data 
            Train_omic1    <- Kern_omic1[Train,Train]
            Train_omic2   <- Kern_omic2[Train,Train]
            Train_Label    <- Label[Train]
            
            
            #Training an RVM Model
            Para1           <- ParaRVM(Bias = FALSE, Boost = FALSE, MaxIts = 20)
            RVM.Model       <- RVMInt(Train_omic1,Train_omic2,classMat = Train_Label,Para = Para1)
            
            # #Training a Boosted-RVM Model
            Para2           <- ParaRVM(Bias = FALSE, Boost = TRUE, MaxIts = 10, BoostIts = 10,resample_size = 15) 
            AdaRVM.Model    <- RVMInt(Train_omic1,Train_omic2,classMat = Train_Label,Para = Para2)
            
            
            #Get the relevance vector index of RVM and Boosted-RVM model
            #The relevance index returned by function RVMIndex follows the order of input data sources.
            
            ### Indices for RVM Model
            RVM.omic1         <- Train[RVMIndex(RVM.Model)[[1]]]
            RVM.omic2         <- Train[RVMIndex(RVM.Model)[[2]]]
            
            
            # # ### Indices for Boosted-RVM Model
            AdaRVM.omic1     <- Train[RVMIndex(AdaRVM.Model)[[1]]]
            AdaRVM.omic2     <- Train[RVMIndex(AdaRVM.Model)[[2]]]
            
            #Predict with RVM model and Boosted-RVM model
            ### Prediction with RVM model 
            RVM.Test.omic1   <- Kern_omic1[Test,RVM.omic1]
            RVM.Test.omic2   <- Kern_omic2[Test,RVM.omic2]
            RVM.Pred         <- predict(RVM.Model,RVM.Test.omic1,RVM.Test.omic2)
            
            # ### Prediction with RVM model 
            Ada.Test.omic1  <- Kern_omic1[Test,AdaRVM.omic1]
            Ada.Test.omic2  <- Kern_omic2[Test,AdaRVM.omic2]
            Ada.Pred        <- predict(AdaRVM.Model,Ada.Test.omic1,Ada.Test.omic2)
            
            #Get the Contingency Table
            RVM_conf_matrix= table(RVM.Pred,Label[Test])
            # print(table(RVM.Pred,Label[Test]))
            
            library(caret)
            a=confusionMatrix(table(Label[Test],RVM.Pred))
            print(a)
            recall_rvm=a$byClass["Recall"]
            prec_rvm=a$byClass["Precision"]
            
            
            Ada_conf_matrix=table(Ada.Pred,Label[Test])
            # print(table(Ada.Pred,Label[Test]))
            
            
            library(caret)
            b=confusionMatrix(table(Label[Test],Ada.Pred))
            print(b)
            recall_ada=a$byClass["Recall"]
            prec_ada=a$byClass["Precision"]
            
            #Plot AUC Curve
            RVM.Pred         <- predict(RVM.Model, RVM.Test.omic1, RVM.Test.omic2, prob = TRUE)
            
            Ada.Pred        <- predict(AdaRVM.Model, Ada.Test.omic1, Ada.Test.omic2, prob =TRUE)
            
            pecc=function(obs,pred) sum(obs==pred)/length(pred)
            # RVM_PECC= pecc(RVM.Pred,Label[Test])
            # ADA_PECC= pecc(Ada.Pred[,2],Label[Test])
            # print("RVM PECC:")
            # print(RVM_PECC)
            
            # print("ADA PECC:")
            # print(ADA_PECC)
            
            png(file="AUC_RVM_Ada_TBI.png", width=512, height=512)
            par(mfcol = c(1,2))
            
            roc(response = Label[Test], predictor = RVM.Pred[,2],plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
            
            title(main = "Relevance Vector Machine")
            
            roc(response = Label[Test], predictor = Ada.Pred[,2],plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
            title(main = "Boosted Relevance Vector Machine")
            dev.off()
            
            
            RVM_res=list(RVM.Model, a, prec_rvm, recall_rvm, AdaRVM.Model, b, prec_ada, recall_ada) 
            return(RVM_res)
            }''')
            
            
        robjects.r(''' run_RVM_CS2 <- function(omic1, omic2, Y){
        
    
        install.packages("kernlab")
        library(kernlab)
        library(signal)
        library(MASS)
        library(MDIntegration)
        library(caret)
        library(pROC)
        # Y= as.matrix(dataset_concat[,425])
        Label= factor(Y)
        print(Label)
        
        rbf <- rbfdot(sigma = 0.05)
        rbf
        
        ## calculate kernel matrix
        
        
        Kern_omic1= kernelMatrix(rbf, as.matrix(omic1)) 
        Kern_omic2= kernelMatrix(rbf, as.matrix(omic2))
        print(dim(Kern_omic1))
        print(dim(Kern_omic2))
        set.seed(1024)
        Train          <- 1:length(Label)
        Test           <- 1:length(Label)
        
        ### Training Data 
        Train_omic1    <- Kern_omic1
        Train_omic2   <- Kern_omic2
        Train_Label    <- Label
        
        
        #Training an RVM Model
        Para1           <- ParaRVM(Bias = FALSE, Boost = FALSE, MaxIts = 20)
        RVM.Model       <- RVMInt(Train_omic1,Train_omic2,classMat = Train_Label,Para = Para1)
        
        print("Done RVM.Model")
        # #Training a Boosted-RVM Model
        Para2           <- ParaRVM(Bias = FALSE, Boost = TRUE, MaxIts = 5, BoostIts = 5, resample_size = 20) 
        AdaRVM.Model    <- RVMInt(Train_omic1,Train_omic2,classMat = Train_Label,Para = Para2)
        
        print("Done AdaRVM.Model!")
        #Get the relevance vector index of RVM and Boosted-RVM model
        #The relevance index returned by function RVMIndex follows the order of input data sources.
        
        ### Indices for RVM Model
        RVM.omic1         <- Train[RVMIndex(RVM.Model)[[1]]]
        RVM.omic2         <- Train[RVMIndex(RVM.Model)[[2]]]
        
        
        # # ### Indices for Boosted-RVM Model
        AdaRVM.omic1     <- Train[RVMIndex(AdaRVM.Model)[[1]]]
        AdaRVM.omic2     <- Train[RVMIndex(AdaRVM.Model)[[2]]]
        
        #Predict with RVM model and Boosted-RVM model
        ### Prediction with RVM model 
        RVM.Test.omic1   <- Kern_omic1[,RVM.omic1]
        RVM.Test.omic2   <- Kern_omic2[,RVM.omic2]
        RVM.Pred         <- predict(RVM.Model,RVM.Test.omic1,RVM.Test.omic2)
        
        # ### Prediction with RVM model 
        Ada.Test.omic1  <- Kern_omic1[,AdaRVM.omic1]
        Ada.Test.omic2  <- Kern_omic2[,AdaRVM.omic2]
        Ada.Pred        <- predict(AdaRVM.Model,Ada.Test.omic1,Ada.Test.omic2)
        
        #Get the Contingency Table
        RVM_conf_matrix= table(RVM.Pred,Label)
        # print(table(RVM.Pred,Label[Test]))
        
        library(caret)
        a=confusionMatrix(table(Label,RVM.Pred))
        print(a)
        recall_rvm=a$byClass["Recall"]
        prec_rvm=a$byClass["Precision"]
        
        
        Ada_conf_matrix=table(Ada.Pred,Label)
        # print(table(Ada.Pred,Label))
        
        
        library(caret)
        b=confusionMatrix(table(Label,Ada.Pred))
        print(b)
        recall_ada=a$byClass["Recall"]
        prec_ada=a$byClass["Precision"]
        
        #Plot AUC Curve
        RVM.Pred         <- predict(RVM.Model, RVM.Test.omic1, RVM.Test.omic2, prob = TRUE)
        
        Ada.Pred        <- predict(AdaRVM.Model, Ada.Test.omic1, Ada.Test.omic2, prob =TRUE)
        
        pecc=function(obs,pred) sum(obs==pred)/length(pred)
        # RVM_PECC= pecc(RVM.Pred,Label[Test])
        # ADA_PECC= pecc(Ada.Pred[,2],Label[Test])
        # print("RVM PECC:")
        # print(RVM_PECC)
        
        # print("ADA PECC:")
        # print(ADA_PECC)
        
        png(file="AUC_RVM_Ada_TBI.png", width=512, height=512)
        par(mfcol = c(1,2))
        
        roc(response = Label[Test], predictor = RVM.Pred[,2],plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
        
        title(main = "Relevance Vector Machine")
        
        roc(response = Label[Test], predictor = Ada.Pred[,2],plot=TRUE,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE)
        title(main = "Boosted Relevance Vector Machine")
        dev.off()
        
        
        RVM_res=list(RVM.Model, a, prec_rvm, recall_rvm, AdaRVM.Model, b, prec_ada, recall_ada) 
        return(RVM_res)
        }''')
            
        print(utils.head(self.omic_1))
        print(utils.head(self.omic_2))
        print(self.Y)
        Y= robjects.FactorVector(self.Y)
        packages=robjects.r["packages"]
        packages()
        run_RVM= robjects.r["run_RVM"]
        RVM_res= run_RVM(self.omic_1, self.omic_2, Y)
        Model_RVM=RVM_res[0]
        cm_rvm=RVM_res[1]
        precision_rvm=RVM_res[2]
        recall_rvm=RVM_res[3]
        
        Model_Ada=RVM_res[3]
        cm_Ada=RVM_res[4]
        precision_ada=RVM_res[5]
        recall_ada=RVM_res[6]
        file= open("RVM_Ada_TBI.txt", "w")
        
        file.write(f'Dimension of dataset {self.omics[0]}: {base.dim(self.omic_1)} \n')
        file.write(f'Dimension of dataset {self.omics[1]}: {base.dim(self.omic_2)} \n')
        file.write('RVM results: \n')
        file.write(f'RVM Model:\n {Model_RVM} \n')
        file.write(f'Confusion Matrix:\n {cm_rvm} \n')
        file.write(f'Precision:\n {precision_rvm} \n')
        file.write(f'Recall:\n {recall_rvm} \n')
        
        file.write("\n")
        file.write('Ada results: \n')
        file.write(f'Ada Model:\n {Model_Ada} \n')
        file.write(f'Confusion Matrix:\n {cm_Ada} \n')
        file.write(f'Precision:\n {precision_ada} \n')
        file.write(f'Recall:\n {recall_ada} \n')
  
        file.write(' See image AUC_RVM_Ada_TBI_res.png for the AUC curve plot. \n')
        
        file.close()
        
        return(print(" RVM and AdaBoost executed successfully!"))
    
    
    
    '''Model - Based Integration '''
    
    def Ensemble_classifier(self, option, voting):
        '''
        Executes Ensemble Classifier using both Soft and Hard Voting, depending 
        on the chosen option.

        Parameters
        ----------
        train_test_datasets : TYPE
            DESCRIPTION.
        omics : TYPE
            DESCRIPTION.
        option : STR
            Option:
                op1  <- Soft Voting with 2 models (SVM)
            
                op2  <- Soft Voting with 2 models (Decision Tree and Guassian Naive Bayes)
            
                op3  <- Soft Voting with 2 Neural Networks
                        
                op4  <- Hard Voting with ensemble of  2 Naive Bayes (combination of models using voting classifier with NB as recommended for a MBI approach)
                        
                ALL  <- Executes all the above.

        Returns
        -------
        .txt files with the results obtain from the executed option.
        

        '''
        
        def fit_multiple_estimators(classifiers, X_list, y, X_test_list, opt,clas,sample_weights = None):
            
            # Convert the labels `y` using LabelEncoder, because the predict method is using index-based pointers
            # which will be converted back to original data later.
            le_ = LabelEncoder()
            le_.fit(y)
            transformed_y = le_.transform(y)
            print(transformed_y)
            
            
           # Fit all estimators with their respective feature arrays
            # estimators_ = [clf.fit(X, y) if sample_weights is None else clf.fit(X, y, sample_weights) for clf, X in zip([clf for _, clf in classifiers], X_list)]
        
            estimators_=[]
            importances={}
            for clf, X, X_test,clas in zip([clf for _, clf in classifiers], X_list, X_test_list, clas):
                if sample_weights is None:
                    print(clas)
                    
                    clf=clf.fit(X,y)
                    estimators_.append(clf.fit(X,y))
                    if "coef_" in dir(clf):
                        print(clas+": Has coef_")
                        varImp= clf.coef_[0]
                        headers = ["name", "score"]
                        values = sorted(zip(X.columns, varImp), key=lambda x: x[1] * -1)
                        print(tabulate(values, headers, tablefmt="plain"))
                        importance= tabulate(values, headers, tablefmt="plain")
                        importances[clas]=importance
                    elif "feature_importances_" in dir(clf):
                        print(clas+": Has feature_importances_")
                        varImp= clf.feature_importances_
                        headers = ["name", "score"]
                        values = sorted(zip(X.columns, varImp), key=lambda x: x[1] * -1)
                        print(tabulate(values, headers, tablefmt="plain"))
                        importance= tabulate(values, headers, tablefmt="plain")
                        importances[clas]=importance
                                            
                    elif "coef_" and "feature_importances_"  not in dir(clf):
                        print(clas+": Doesn't have coef_ nor feature_importances_")
                        import lime 
                        from lime import lime_tabular
                        lime_explainer = lime_tabular.LimeTabularExplainer(
                            training_data=np.array(X),
                            feature_names=X.columns,
                            class_names=["Control", "Drought"],
                            mode='classification')
                        warnings.simplefilter(action='ignore', category=FutureWarning)
                        
                        idx=1
                        # test_1=X_test.iloc[1,:]
                        # print(test_1)
                        
                      
                        exp = lime_explainer.explain_instance(X_test.iloc[idx,:].astype(int).values, clf.predict_proba, num_features=20)
                        exp.show_in_notebook()
                        exp.save_to_file(clas+'_LIME_'+ str(idx)+ opt +'.html')
                        # exp.show_in_notebook(show_table=True)
                       
                        from sklearn.inspection import permutation_importance
                        imps = permutation_importance(clf, X, y)
                        print(imps.importances_mean)
                        headers = ["name", "score"]
                        values = sorted(zip(X.columns, imps.importances_mean), key=lambda x: x[1] * -1)
                        print(tabulate(values, headers, tablefmt="plain"))
                        importance= tabulate(values, headers, tablefmt="plain")
                        importances[clas]=importance
                       
                else:
                    estimators_.append(clf.fit(X,y,sample_weights))
                
                
                
                
                
        
        
        
        
            return estimators_, le_,importances


        def predict_from_multiple_estimator(estimators, label_encoder, X_list, weights = None):
        
            # Predict 'soft' voting with probabilities
        
            pred1 = np.asarray([clf.predict_proba(X) for clf, X in zip(estimators, X_list)])
            # print(pred1)
            pred2 = np.average(pred1, axis=0, weights=weights)
            # print(pred2)
            pred = np.argmax(pred2, axis=1)
            # print(pred)
            return label_encoder.inverse_transform(pred)


        def predict_from_multiple_estimator_hard(estimators, label_encoder, X_list, weights=None):
            # Predict 'Hard' voting with probabilities
             
              pred2 = np.asarray([clf.predict(X) for clf, X in zip(estimators, X_list)])
              # print(pred2)
              lista=[]
              for i in pred2:
                le_ = LabelEncoder()
                le_.fit(i)
                transformed_y = le_.transform(i)
                # print(transformed_y)
                lista.append(transformed_y)
              # return print(np.asarray(lista))
              pred3= np.asarray(lista)
              # print("Lista:",lista)
              # pred1= np.bincount(np.asarray(lista), weights=weights)
              pred= np.apply_along_axis(lambda x: np.argmax(np.bincount(x, weights=weights)), axis=0, arr=pred3.astype('int'))
              # print(pred)
        
            # Convert integer predictions to original labels:
              return label_encoder.inverse_transform(pred)
                
        def option1(dataset_concat,Y,dim1):
            ''' Case Study 1'''
            opt="option1"
            X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
             
            
        #Divide the X into different feature datas:
            print(X_train.iloc[:,0:(dim1[1]-1)].head(6))
            print(X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1].head(6))
            X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1] #transcriptomics and metabolomics
            X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1])-1]
            
            X_train_list = [X_train1, X_train2]
            X_test_list = [X_test1, X_test2]
            
            
        #     '''Case Study 2'''
            
        #     # print(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)])
        #     # X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
            
        #     # Y_1=dataset_concat.loc[:,"treatment"].astype('category')
        #     X_train, X_test, y_train, y_test= dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.loc[:,"treatment"] ,dataset_concat.loc[:,"treatment"]
            
        # #Divide the X into different feature datas:
        
        #     X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)] #transcriptomics and metabolomics
        #     X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)]
            
        #     X_train_list = [X_train1, X_train2]
        #     X_test_list = [X_test1, X_test2]
            
            print(X_train1.shape)
            print(X_train2.shape)
            print(X_test1.shape)
            print(X_test2.shape)
            
           
            #Get list of classifiers:
            
            # from sklearn.neighbors import KNeighborsClassifier
            # from sklearn.svm import SVC
            
            # Make sure the number of estimators here are equal to number of different feature datas
            # classifiers = [('knn',  KNeighborsClassifier(3)),
            #     ('svc', SVC(kernel="linear", C=0.025, probability=True))]
            
            classifiers = [('svc', SVC(kernel="linear", C=0.025, probability=True)),
                ('svc', SVC(kernel="linear", C=0.025, probability=True))]
            #Fit the classifiers with the data:
            clas=["scv","svc1"]
            fitted_estimators, label_encoder,importances = fit_multiple_estimators(classifiers, X_train_list,y_train,X_test_list,opt,clas) 
            
                # for j,v in enumerate(varImp):
                #     print(f'Feature: {j}, Score: %.5f' % (j,v))
                #     file.writelines(f'Feature: {features_names[i]}, Score: %.5f \n' % (v))

                
                
                
            # print(dir(fitted_estimators[0]))
            
            #Predict using the test data:
            if voting.lower() == "soft":
                y_pred = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                acc_soft= accuracy_score(y_test, y_pred)
                cm_soft= confusion_matrix(y_test, y_pred)
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                roc_auc = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc)  
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt1_soft.png')
                plt.show()
                
                
                res= [acc_soft, cm_soft, roc_auc, cl_rs]
                
                
                
            if voting.lower()== "hard":
                y_pred = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                acc_hard= accuracy_score(y_test, y_pred)
                cm_hard= confusion_matrix(y_test, y_pred)
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                
                # from sklearn.preprocessing import LabelBinarizer
                # lb = LabelBinarizer()
                # # y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                # y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                # false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                # roc_auc = auc(false_positive_rate, true_positive_rate) 
                
                # print("ROC_AUC:", roc_auc)  
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt1_hard.png')
                # plt.show()
                
                
                
                res= [acc_hard,cm_hard, roc_auc,cl_rs]
            
            #Get accuracy of predictions:
            if voting.lower()=="all":
                y_pred_soft = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                y_pred_hard = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                
                acc_soft=  accuracy_score(y_test, y_pred_soft) 
                acc_hard=  accuracy_score(y_test, y_pred_hard)
                cm_soft= confusion_matrix(y_test, y_pred_soft)
                cm_hard= confusion_matrix(y_test, y_pred_hard)
                cl_rs_soft=classification_report(y_test, y_pred_soft, zero_division=0)
                cl_rs_hard=classification_report(y_test, y_pred_hard, zero_division=0)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_soft = np.array([number[0] for number in lb.fit_transform(y_pred_soft)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                
            
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_soft)
                roc_auc_soft = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_soft)  
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_soft)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC ')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt1_all_soft.png')
                plt.show()
                
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_hard = np.array([number[0] for number in lb.fit_transform(y_pred_hard)])
                
               
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_hard)
                roc_auc_hard = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_hard)  
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_hard)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt1_all_hard.png')
                # plt.show()
                
                
                
                res= [acc_soft, acc_hard, cm_soft,cm_hard, roc_auc_soft, roc_auc_hard, cl_rs_soft, cl_rs_hard]
            print("Soft Accuracy: ", acc_soft)
            return res,importances
                
        def option2(dataset_concat,Y,dim1):
            X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
            
            ''' Case Study 1'''
            opt="option2"
            X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
             
            
        #Divide the X into different feature datas:
            print(X_train.iloc[:,0:(dim1[1]-1)].head(6))
            print(X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1].head(6))
            X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1] #transcriptomics and metabolomics
            X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1])-1]
            
            X_train_list = [X_train1, X_train2]
            X_test_list = [X_test1, X_test2]
            
            
            '''Case Study 2'''
            
        #     # print(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)])
        #     # X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
            
        #     # Y_1=dataset_concat.loc[:,"treatment"].astype('category')
        #     X_train, X_test, y_train, y_test= dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.loc[:,"treatment"] ,dataset_concat.loc[:,"treatment"]
            
        # #Divide the X into different feature datas:
        
        #     X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)] #transcriptomics and metabolomics
        #     X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)]
            
        #     X_train_list = [X_train1, X_train2]
        #     X_test_list = [X_test1, X_test2]
            
            # print(X_train1.shape)
            # print(X_train2.shape)
            # print(X_test1.shape)
            # print(X_test2.shape)
         
            classifiers = [
               ("dtree", DecisionTreeClassifier(random_state=0, max_depth=2)),
               ("gnb", GaussianNB())]

            #Fit the classifiers with the data:
            clas=["dtree","gnb"]
            fitted_estimators, label_encoder,importances = fit_multiple_estimators(classifiers, X_train_list,y_train,X_test_list,opt,clas) 
            
            
            #Predict using the test data:
            
            if voting.lower() == "soft":
                y_pred = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                acc_soft= accuracy_score(y_test, y_pred)
                cm_soft= confusion_matrix(y_test, y_pred)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                roc_auc = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc)  
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Grid Search')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt2_soft.png')
                plt.show()
                
                res= [acc_soft, cm_soft, roc_auc, cl_rs]
                
            if voting.lower()== "hard":
                y_pred = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                acc_hard= accuracy_score(y_test, y_pred)
                cm_hard= confusion_matrix(y_test, y_pred)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                roc_auc = auc(false_positive_rate, true_positive_rate) 
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                print("ROC_AUC:", roc_auc)  
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt2_hard.png')
                # plt.show()
                
                
                res= [acc_hard,cm_hard,cl_rs]
            
            #Get accuracy of predictions:
            if voting.lower()=="all":
                y_pred_soft = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                y_pred_hard = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                
                acc_soft=  accuracy_score(y_test, y_pred_soft) 
                acc_hard=  accuracy_score(y_test, y_pred_hard)
                cm_soft= confusion_matrix(y_test, y_pred_soft)
                cm_hard= confusion_matrix(y_test, y_pred_hard)
                cl_rs_soft=classification_report(y_test, y_pred_soft, zero_division=0)
                cl_rs_hard =classification_report(y_test, y_pred_hard, zero_division=0)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_soft = np.array([number[0] for number in lb.fit_transform(y_pred_soft)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_soft)
                roc_auc_soft = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_soft)  
                
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_soft)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Grid Search')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt2_all_soft.png')
                plt.show()
                
                
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_hard = np.array([number[0] for number in lb.fit_transform(y_pred_hard)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_hard)
                roc_auc_hard = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_hard)  
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_hard)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt2_all_hard.png')
                # plt.show()
                
                
            
                res= [acc_soft, acc_hard, cm_soft,cm_hard, roc_auc_soft, roc_auc_hard, cl_rs_soft, cl_rs_hard]
            print("Soft Accuracy: ", acc_soft)

            return res, importances          
                        
        def option3(dataset_concat,Y,dim1):
            ''' Case Study 1'''
            opt="option3"
            X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
             
            
        #Divide the X into different feature datas:
            print(X_train.iloc[:,0:(dim1[1]-1)].head(6))
            print(X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1].head(6))
            X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1] #transcriptomics and metabolomics
            X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1])-1]
            
            X_train_list = [X_train1, X_train2]
            X_test_list = [X_test1, X_test2]
            
            
            '''Case Study 2'''
            
        #     # print(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)])
        #     # X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
            
        #     # Y_1=dataset_concat.loc[:,"treatment"].astype('category')
        #     X_train, X_test, y_train, y_test= dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.loc[:,"treatment"] ,dataset_concat.loc[:,"treatment"]
            
        # #Divide the X into different feature datas:
        
        #     X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)] #transcriptomics and metabolomics
        #     X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)]
            
        #     X_train_list = [X_train1, X_train2]
        #     X_test_list = [X_test1, X_test2]
            
            # print(X_train1.shape)
            # print(X_train2.shape)
            # print(X_test1.shape)
            # print(X_test2.shape)
            
            classifiers = [ ("mlp1",MLPClassifier(solver='lbfgs', hidden_layer_sizes=30,
                        random_state=0)), ("mlp",MLPClassifier(solver='lbfgs', max_iter= 15000,hidden_layer_sizes=30,
                        random_state=0))]
            #Fit the classifiers with the data:
            clas=["mlp1","mlp2"]
            fitted_estimators, label_encoder,importances = fit_multiple_estimators(classifiers, X_train_list,y_train,X_test_list,opt,clas) 
            
            
            if voting.lower() == "soft":
                y_pred = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                acc_soft= accuracy_score(y_test, y_pred)
                cm_soft= confusion_matrix(y_test, y_pred)
                
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                roc_auc = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc) 
                
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Grid Search')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt3_soft.png')
                plt.show()
                
                
                
                res= [acc_soft, cm_soft, roc_auc, cl_rs]
            if voting.lower()== "hard":
                y_pred = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                acc_hard= accuracy_score(y_test, y_pred)
                cm_hard= confusion_matrix(y_test, y_pred)
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                roc_auc = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc)  
                
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt3_hard.png')
                # plt.show()
                
                res= [acc_hard,cm_hard, roc_auc, cl_rs]
                
            
            #Get accuracy of predictions:
            if voting.lower()=="all":
                y_pred_soft = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                y_pred_hard = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                
                acc_soft=  accuracy_score(y_test, y_pred_soft) 
                acc_hard=  accuracy_score(y_test, y_pred_hard)
                cm_soft= confusion_matrix(y_test, y_pred_soft)
                cm_hard= confusion_matrix(y_test, y_pred_hard)
                
                cl_rs_soft=classification_report(y_test, y_pred_soft, zero_division=0)
                cl_rs_hard=classification_report(y_test, y_pred_hard, zero_division=0)
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_soft = np.array([number[0] for number in lb.fit_transform(y_pred_soft)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_soft)
                roc_auc_soft = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_soft)  
                
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_soft)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Grid Search')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt3_all_soft.png')
                plt.show()
                
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_hard = np.array([number[0] for number in lb.fit_transform(y_pred_hard)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_hard)
                roc_auc_hard = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_hard)  
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_hard)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt3_all_hard.png')
                # plt.show()
                
                
                
                
                res= [acc_soft, acc_hard, cm_soft,cm_hard, roc_auc_soft, roc_auc_hard,cl_rs_soft, cl_rs_hard]
            print("Soft Accuracy: ", acc_soft)

            return res, importances    
        
        def option4(dataset_concat,Y,dim1):
            
        #     ''' Case Study 1'''
            opt="option4"
            X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
             
            
        #Divide the X into different feature datas:
            print(X_train.iloc[:,0:(dim1[1]-1)].head(6))
            print(X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1].head(6))
            X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1])-1] #transcriptomics and metabolomics
            X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1])-1]
            
            X_train_list = [X_train1, X_train2]
            X_test_list = [X_test1, X_test2]
            
            ''' Case Study 2 '''
            
            # print(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)])
          # X_train, X_test, y_train, y_test = train_test_split(dataset_concat.iloc[:,0:(dataset_concat.shape[1]-1)], Y, test_size=0.20)
            
          # Y_1=dataset_concat.loc[:,"treatment"].astype('category')
        #    X_train, X_test, y_train, y_test= dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.iloc[:,0:(dataset_concat.shape[1]-2)], dataset_concat.loc[:,"treatment"] ,dataset_concat.loc[:,"treatment"]
            
        # #Divide the X into different feature datas:
        
        #    X_train1, X_train2 = X_train.iloc[:,0:(dim1[1]-1)], X_train.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)] #transcriptomics and metabolomics
        #    X_test1, X_test2 = X_test.iloc[:,0:(dim1[1]-1)], X_test.iloc[:,dim1[1]:(dataset_concat.shape[1]-2)]
            
        #    X_train_list = [X_train1, X_train2]
        #    X_test_list = [X_test1, X_test2]
            
            print(X_train1.shape)
            print(X_train2.shape)
            print(X_test1.shape)
            print(X_test2.shape)
            
            classifiers = [ ("gnb_micro",GaussianNB()), ("gnb_prot", GaussianNB())]
            
            clas=["gnb_micro","gnb_prot"]
            #Fit the classifiers with the data:
            
            fitted_estimators, label_encoder,importances = fit_multiple_estimators(classifiers, X_train_list,y_train,X_test_list,opt, clas) 
            
            if voting.lower() == "soft":
                y_pred = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                acc_soft= accuracy_score(y_test, y_pred)
                cm_soft= confusion_matrix(y_test, y_pred)
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                roc_auc = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc) 
                
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Grid Search')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt4_soft.png')
                plt.show()
                
                
                
                res= [acc_soft, cm_soft, roc_auc, cl_rs]
            if voting.lower()== "hard":
                y_pred = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                acc_hard= accuracy_score(y_test, y_pred)
                cm_hard= confusion_matrix(y_test, y_pred)
                cl_rs=classification_report(y_test, y_pred, zero_division=0)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred = np.array([number[0] for number in lb.fit_transform(y_pred)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
                roc_auc = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc)  
                
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt4_hard.png')
                # plt.show()
                
                res= [acc_hard,cm_hard, roc_auc, cl_rs]
            
            #Get accuracy of predictions:
            if voting.lower()=="all":
                y_pred_soft = predict_from_multiple_estimator(fitted_estimators, label_encoder, X_test_list)
                y_pred_hard = predict_from_multiple_estimator_hard(fitted_estimators, label_encoder, X_test_list)
                
                acc_soft=  accuracy_score(y_test, y_pred_soft) 
                acc_hard=  accuracy_score(y_test, y_pred_hard)
                cm_soft= confusion_matrix(y_test, y_pred_soft)
                cm_hard= confusion_matrix(y_test, y_pred_hard)
                
                cl_rs_soft=classification_report(y_test, y_pred_soft, zero_division=0)
                cl_rs_hard=classification_report(y_test, y_pred_hard, zero_division=0)
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_soft = np.array([number[0] for number in lb.fit_transform(y_pred_soft)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_soft)
                roc_auc_soft = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_soft)  
                
                
                
                plt.figure()
                plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_soft)
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Grid Search')
                plt.legend(loc="lower right")
                plt.savefig('ROC_opt4_all_soft.png')
                plt.show()
                
                
                from sklearn.preprocessing import LabelBinarizer
                lb = LabelBinarizer()
                y_pred_hard = np.array([number[0] for number in lb.fit_transform(y_pred_hard)])
                y_test= np.array([number[0] for number in lb.fit_transform(y_test)])
                false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred_hard)
                roc_auc_hard = auc(false_positive_rate, true_positive_rate) 
                print("ROC_AUC:", roc_auc_hard)  
                
                # plt.figure()
                # plt.plot(false_positive_rate, true_positive_rate, label='ROC curve (area = %0.2f)' % roc_auc_hard)
                # plt.plot([0, 1], [0, 1], 'k--')
                # plt.xlim([0.0, 1.0])
                # plt.ylim([0.0, 1.05])
                # plt.xlabel('False Positive Rate')
                # plt.ylabel('True Positive Rate')
                # plt.title('ROC Grid Search')
                # plt.legend(loc="lower right")
                # plt.savefig('ROC_opt4_all_hard.png')
                # plt.show()
                
                
                res= [acc_soft, acc_hard, cm_soft,cm_hard, roc_auc_soft, roc_auc_hard, cl_rs_soft, cl_rs_hard]
            
            print("Soft Accuracy: ", acc_soft)
            return res, importances           
            
            
            
            
            
        omics_=[]
        names=[]
        for name, lista in self.train_test_datasets.items():
            if name in self.omics:
                # confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        if len(omics_) ==2:
            omics1= omics_[0]
            omics2= omics_[1]
               
            print("Dimension of Dataset Concatenated:")
            dataset_concat= pd.read_excel("dataset_concat.xlsx", index_col=0) 
            print(dataset_concat.shape)
            
            # print(dataset_concat.head(6))
            
            dim1=base.dim(omics1)
            # print(dim1)
            range_omic1= np.arange(1, dim1[1], 1)
            
            range_omic1_py=np.arange(0,dim1[1],1)
            
            dim2=base.dim(omics2)
            # print(dim2)
            dim_concat=base.dim(dataset_concat)
            len_concat=dim_concat[1]
            final_2= dim_concat[1]
            range_omic2= np.arange(dim1[1],final_2, 1)
            # omic_2=dataset_concat.rx(True,range_omic2)
            # print(utils.head(omic_2))
       
        # omic1=omics[0]
        # omic2=omics[1]
        omics1=dataset_concat.iloc[:,np.arange(0, dim1[1]-1, 1)]
        print("Dimension of "+ names[0] + ":")
        print(omics1.shape)
        
        print("Dimension of "+ names[1]+ ":")
        range_omic2_py=np.arange(dim1[1]-1,final_2-1, 1)
        omics2=dataset_concat.iloc[:,np.arange(dim1[1]-1,final_2-1, 1)]
        print(omics2.shape)
        
        Y= dataset_concat.loc[:, self.y_pred]
        # print(Y)
        #The logic is taken from VotingClassifier source.
        #Now test the above methods. First get some data:
        
        #Split the data into train and test:
        
        # from sklearn.model_selection import train_test_split
       
            
        if option == "op1":
            res, importances= option1(dataset_concat,Y,dim1)
            if voting.lower()=="soft":
                acc_soft=res[0]
                cm_soft=res[1]
                roc_auc=res[2]
                cl_rs=res[3]
                
                name= "EnsembleC_MBI_soft"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix: \n {cm_soft} \n')
                file.write(f'Accuracy: \n {acc_soft} \n')
                file.write(f'ROC_AUC: \n {roc_auc} \n')
                for name, data in importances.items():
                    print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                file.write(f'Classification Report: \n {cl_rs} \n')
                file.close()
                               
                
            if voting.lower()=="hard":
                acc_hard=res[0]
                # cm_hard=[1]
                roc_auc=res[2]
                cl_rs=res[3]
                name= "EnsembleC_MBI_hard"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                # file.write(f'Confusion Matrix: \n {cm_hard} \n')
                file.write(f'Accuracy: \n {acc_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                # file.write(f'ROC_AUC: \n {roc_auc} \n')
                file.write(f'Classification Report: \n {cl_rs} \n ')
                file.close()
                
            if voting.lower()== "all":
                acc_soft=res[0]
                acc_hard=res[1]
                cm_soft=res[2]
                roc_auc_soft=res[4]
                roc_auc_hard=res[5]
                cl_rs_soft=res[6]
                cl_rs_hard=res[7]
                # cm_hard=[3]
                
                name= "EnsembleC_MBI_all"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix Soft: \n {cm_soft} \n')
                file.write(f'Accuracy Soft: \n {acc_soft} \n')
                file.write(f'ROC_AUC Soft: \n {roc_auc_soft} \n')
                # file.write(f'Confusion Matrix Hard: \n {cm_hard} \n')
                file.write(f'Accuracy Hard: \n {acc_hard} \n')
                # file.write(f'ROC_AUC Hard: \n {roc_auc_hard} \n')
                file.write(f'Classification Report Soft: \n {cl_rs_soft} \n')
                file.write(f'Classification Report Hard: \n {cl_rs_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                # file.write(f'Variable Importance: \n {importances} \n')
                file.close()
                
        if option== "op2":
            res,importances= option2(dataset_concat,Y,dim1)
            if voting.lower()=="soft":
                acc_soft=res[0]
                cm_soft=res[1]
                roc_auc=res[2]
                cl_rs=res[3]
                name= "EnsembleC_MBI_soft"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix: \n {cm_soft} \n')
                file.write(f'Accuracy: \n {acc_soft} \n')
                file.write(f'ROC_AUC: \n {roc_auc} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                file.write(f'Classification Report: \n {cl_rs} \n ')
                file.close()
                               
                
            if voting.lower()=="hard":
                acc_hard=res[0]
                cm_hard=[1]
                roc_auc=res[2]
                cl_rs=res[3]
                name= "EnsembleC_MBI_hard"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix: \n {cm_hard} \n')
                file.write(f'Accuracy: \n {acc_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                # file.write(f'ROC_AUC: \n {roc_auc} \n')
                file.write(f'Classification Report: \n {cl_rs} \n ')
                file.close()
                
            if voting.lower()== "all":
                acc_soft=res[0]
                acc_hard=res[1]
                cm_soft=res[2]
                cm_hard=[3]
                roc_auc_soft=res[4]
                roc_auc_hard=res[5]
                cl_rs_soft=res[6]
                cl_rs_hard=res[7]
                
                name= "EnsembleC_MBI_all"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix Soft: \n {cm_soft} \n')
                file.write(f'Accuracy Soft: \n {acc_soft} \n')
                
                file.write(f'Confusion Matrix Hard: \n {cm_hard} \n')
                file.write(f'Accuracy Hard: \n {acc_hard} \n')
                file.write(f'Classification Report Soft: \n {cl_rs_soft} \n')
                file.write(f'Classification Report Hard: \n {cl_rs_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                file.close()
                
        if option =="op3":
            res,importances= option3(dataset_concat,Y,dim1)
            if voting.lower()=="soft":
                acc_soft=res[0]
                cm_soft=res[1]
                roc_auc_soft=res[2]
                
                cl_rs_hard=res[3]
                
                name= "EnsembleC_MBI_soft"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix: \n {cm_soft} \n')
                file.write(f'Accuracy: \n {acc_soft} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                file.write(f'Classification Report: \n {cl_rs} ')
                file.close()
                               
                
            if voting.lower()=="hard":
                acc_hard=res[0]
                cm_hard=[1]
                roc_auc_hard=res[2]
                cl_rs_hard=res[3]
                name= "EnsembleC_MBI_hard"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix: \n {cm_hard} \n')
                file.write(f'Accuracy: \n {acc_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                file.write(f'Classification Report: \n {cl_rs} ')
                file.close()
                
            if voting.lower()== "all":
                acc_soft=res[0]
                acc_hard=res[1]
                cm_soft=res[2]
                cm_hard=[3]
                roc_auc_soft=res[4]
                roc_auc_hard=res[5]
                cl_rs_soft=res[6]
                cl_rs_hard=res[7]
                name= "EnsembleC_MBI_all"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix Soft: \n {cm_soft} \n')
                file.write(f'Accuracy Soft: \n {acc_soft} \n')
                file.write(f'ROC AUC Soft: \n {roc_auc_soft}')
                file.write(f'Confusion Matrix Hard: \n {cm_hard} \n')
                file.write(f'Accuracy Hard: \n {acc_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                # file.write(f'ROC AUC Hard: \n {roc_auc_hard} \n')
                file.write(f'Classification Report Soft: \n {cl_rs_soft}  \n')
                file.write(f'Classification Report Hard: \n {cl_rs_hard} ')
                file.close()
  
        if option =="op4":

            res,importances= option4(dataset_concat,Y,dim1)
            if voting.lower()=="soft":
                acc_soft=res[0]
                cm_soft=res[1]
                roc_auc=res[2]
                cl_rs=res[3]
    
                name= "EnsembleC_MBI_soft"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix: \n {cm_soft} \n')
                file.write(f'Accuracy: \n {acc_soft} \n')
                file.write(f'Roc Auc: \n {roc_auc} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                file.write(f'Classification Report: \n {cl_rs} ')
                file.close()
                               
                
            if voting.lower()=="hard":
                acc_hard=res[0]
                # cm_hard=[1]
                roc_auc=res[2]
                cl_rs=res[3]
                name= "EnsembleC_MBI_hard"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                # file.write(f'Confusion Matrix: \n {cm_hard} \n')
                # file.write(f' Roc Auc: \n {roc_auc} \n')
                file.write(f'Accuracy: \n {acc_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                file.write(f'Classification Report: \n {cl_rs} \n')
                file.close()
                
            if voting.lower()== "all":
                acc_soft=res[0]
                acc_hard=res[1]
                cm_soft=res[2]
                # cm_hard=[3]
                roc_auc_soft=res[4]
                roc_auc_hard=res[5]
                cl_rs_soft=res[6]
                cl_rs_hard=res[7]
                
                name= "EnsembleC_MBI_all"+ option + ".txt"
                file= open(name,"w")
                file.write('Ensemble Classifier (Model_Based Integration results: \n')
                file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                file.write(f'Confusion Matrix Soft: \n {cm_soft} \n')
                file.write(f'Accuracy Soft: \n {acc_soft} \n')
                file.write(f' Roc Auc Soft: \n {roc_auc_soft} \n')
                # file.write(f'Confusion Matrix Hard: \n {cm_hard} \n')
                file.write(f'Accuracy Hard: \n {acc_hard} \n')
                for name, data in importances.items():
                    # print(len(importances))
                    file.write(f'Variable Importance: \n {data} \n')
                # file.write(f' Roc Auc Hard: \n {roc_auc_hard} \n')
                file.write(f'Classification Report Soft: \n {cl_rs_soft} \n')
                file.write(f'Classification Report Hard: \n {cl_rs_hard} \n ')
                file.close()
                
        if option == "ALL":
            res1, importances1= option1(dataset_concat,Y,dim1) 
            res2, importances2= option2(dataset_concat,Y,dim1)
            res3, importances3= option3(dataset_concat,Y,dim1)
            res4, importances4= option4(dataset_concat,Y,dim1)
            
            res_all=[res1,res2,res3,res4]
            importances_all=[importances1,importances2,importances3,importances4]
            str_res=["op1","op2","op3","op4"]
            for i in range(len(res_all)):
                
                if voting.lower()=="soft":
                    acc_soft=res_all[i][0]
                    cm_soft=res_all[i][1]
                    roc_auc=res_all[i][2]
                    cl_rs=res_all[i][3]
                    
                    name= "EnsembleC_MBI_soft"+ option +"_"+ str_res[i] + ".txt"
                    file= open(name,"w")
                    file.write('Ensemble Classifier (Model_Based Integration results: \n')
                    file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                    file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                    file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                    file.write(f'Confusion Matrix: \n {cm_soft} \n')
                    file.write(f'Accuracy: \n {acc_soft} \n')
                    file.write(f'Roc Auc : \n {roc_auc} \n')
                    file.write(f'Classification Report: \n {cl_rs} ')
                    
                    for name, data in importances_all[i].items():
                        # print(len(importances))
                        file.write(f'Variable Importance: \n {data} \n')
                    file.close()
                                   
                    
                if voting.lower()=="hard":
                    acc_hard=res_all[i][0]
                    # cm_hard=[1]
                    roc_auc=res_all[i][2]
                    cl_rs=res_all[i][3]

                    name= "EnsembleC_MBI_hard"+ option + "_"+ str_res[i] + ".txt"
                    file= open(name,"w")
                    file.write('Ensemble Classifier (Model_Based Integration results: \n')
                    file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                    file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                    file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                    # file.write(f'Confusion Matrix: \n {cm_hard} \n')
                    file.write(f'Accuracy: \n {acc_hard} \n')
                    # file.write(f'Roc Auc : \n {roc_auc} \n')
                    file.write(f'Classification Report: \n {cl_rs} \n ')
                    for name, data in importances_all[i].items():
                        # print(len(importances))
                        file.write(f'Variable Importance: \n {data} \n')
                    file.close()
                    
                if voting.lower()== "all":
                    acc_soft=res_all[i][0]
                    acc_hard=res_all[i][1]
                    cm_soft=res_all[i][2]
                    roc_auc_soft=res_all[i][4]
                    roc_auc_hard=res_all[i][5]
                    cl_rs_soft=res_all[i][6]
                    cl_rs_hard=res_all[i][7]
                    # print(acc_soft)
                    # print(acc_hard)
                    # print(cm_soft)
                    # cm_hard=[3]
                    
                    name= "EnsembleC_MBI_all"+ option +"_"+ str_res[i] + ".txt"
                    file= open(name,"w")
                    file.write('Ensemble Classifier (Model_Based Integration) results: \n')
                    file.write(f'Dimension of Concatenated dataset: {dataset_concat.shape} \n')
                    file.write(f'Dimension of {names[0]} dataset:  {omics1.shape} \n')
                    file.write(f'Dimension of {names[1]} dataset:  {omics2.shape} \n')
                    file.write(f'Confusion Matrix Soft: \n {cm_soft} \n')
                    file.write(f'Accuracy Soft: \n {acc_soft} \n')
                    file.write(f'Roc Auc Soft : \n {roc_auc_soft} \n')
                    # file.write(f'Confusion Matrix Hard: \n {cm_hard} \n')
                    file.write(f'Accuracy Hard: \n {acc_hard} \n')
                    for name, data in importances_all[i].items():
                        # print(len(importances))
                        file.write(f'Variable Importance: \n {data} \n')
                    # file.write(f'Roc Auc Hard: \n {roc_auc_hard} \n')
                    file.write(f'Classification Report Soft: \n {cl_rs_soft} \n ')
                    file.write(f'Classification Report Hard: \n {cl_rs_hard} \n ')
                    
                    file.close()
                    