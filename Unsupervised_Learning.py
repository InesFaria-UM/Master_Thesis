# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 09:37:31 2021

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

# utils.install_packages("simfinapi")
# simfinapi=rpackages.importr("simfinapi")


# # utils.remove_packages("Rcpp")
# # utils.install_packages("C:/Users/Maria Ines/Desktop/Artigos/Tese/Rcpp_1.0.7.zip", type = "source")
# # Rcpp=rpackages.importr("Rcpp")
# utils.update_packages(ask= False)
# utils.install_packages("writexl")
# writexl=rpackages.importr("writexl")

# utils.install_packages("caret")
# caret=rpackages.importr("caret")
class Unsupervised_Learning:
    
    def __init__(self,train_test_datasets, omics, y_pred):
        self.omics= omics
        self.train_test_datasets=train_test_datasets
        self.y_pred=y_pred
    ''' Concatenation-based '''
    
    def MFA(self, type_): #Multiple Factor Analysis (MFA)
        
        robjects.r('''omic_without_y <- function(omic1, omic2, y_pred){
             omic1_new=omic1
             omic2_new=omic2
             omic1_new[,y_pred]=NULL
             omic2_new[y_pred]=NULL
             omics=list(omic1_new,omic2_new)
             return(omics)
            } ''')
    
    
        robjects.r(''' run_MFA <- function(dataset_concat,nvar1, nvar2,nvary, type, names, y_pred){
             install.packages("FactoMineR")
             library(FactoMineR)
             
             install.packages("factoextra")
             library("factoextra")
             
             res <- MFA(dataset_concat, group=c(nvar1,nvar2,nvary), type=type,
                        ncp=5, name.group=c(names[1], names[2],y_pred), num.group.sup=c(3))
             
             print(res)
             
              #The proportion of variances retained by the different dimensions
              #(axes) can be extracted using the function get_eigenvalue() 
              #[factoextra package] as follow:
                 
            eig.val <- get_eigenvalue(res)
            print(head(eig.val))
            
            #The function fviz_eig() or fviz_screeplot() [factoextra package]
            # can be used to draw the scree plot:
            
            png(file="MFA_scree.png")    
            print(fviz_screeplot(res))
            dev.off()
            
            #Graph of variables
            #Groups of variables
            #The function get_mfa_var() [in factoextra] is used to extract the results for 
            #groups of variables. This function returns a list containing the coordinates, 
            #the cos2 and the contribution of groups, as well as, the
            group <- get_mfa_var(res, "group")
            group
            
            #The different components can be accessed as follow:
            #Coordinates of groups
            print(head(group$coord))
            # Cos2: quality of representation on the factore map
            print(head(group$cos2))
            # Contributions to the  dimensions
            print(head(group$contrib))
            
            
            png(file="MFA_var.png") 
            print(fviz_mfa_var(res, "group"))
            dev.off()
           
            #red color = active groups of variables
            #green color = supplementary groups of variables
            #The plot above illustrates the correlation between groups and dimensions. 
            #Proteomics contributes more for Dim1 and microarrays for Dim2 
            #coordinates indicating a highest contribution to the second dimension.
            
            #To draw a bar plot of groups contribution to the dimensions, use the function 
            # fviz_contrib():
            # Contribution to the first dimension
            png(file="MFA_contrib_DIM1.png")
            
            print(fviz_contrib(res, "group", axes = 1))
            dev.off()
            
            png(file="MFA_contrib_DIM2.png")
            # Contribution to the second dimension
            print(fviz_contrib(res, "group", axes = 2))
            dev.off()
            
            
            #Quantitative variables
            #The function get_mfa_var() [in factoextra] is used to extract 
            # the results for quantitative variables. This function returns a 
            # list containing the coordinates,the cos2 and the contribution of 
            # variables:
            quanti.var <- get_mfa_var(res, "quanti.var")
            print(quanti.var) 
            
            #The different components can be accessed as follow:
            # Coordinates
            print(head(quanti.var$coord))
            # Cos2: quality on the factore map
            print(head(quanti.var$cos2))
            # Contributions to the dimensions
            print(head(quanti.var$contrib))
            
            
            
            #Correlation between quantitative variables and dimensions. The R code below 
            #plots quantitative variables colored by groups. The argument palette is used to
            #change group colors (see ?ggpubr::ggpar for more information about palette).
            # Supplementary quantitative variables are in dashed arrow and violet color. 
            #We use repel = TRUE, to avoid text overlapping.
            png(file="MFA_corr_coloredgroups.png")
            
            print(fviz_mfa_var(res, "quanti.var", palette = "jco", 
                          col.var.sup = "violet", repel = TRUE))
            
            dev.off()
            
            png(file="MFA_corr_points.png")
            #To make the plot more readable, we can use geom = c(“point”, “text”) instead of 
            #geom = c(“arrow”, “text”). We’ll change also the legend position from “right” to
            # “bottom”, using the argument legend = “bottom”:
            print(fviz_mfa_var(res, "quanti.var", palette = "jco", 
              col.var.sup = "violet", repel = TRUE,
              geom = c("point", "text"), legend = "bottom"))

            dev.off()
            
            
            # Contributions to dimension 1

            #The contribution of quantitative variables (in %) to the definition of the 
            #dimensions can be visualized using the function fviz_contrib() [factoextra package].
            # Variables are colored by groups. The R code below shows the top 20 variable 
            # categories contributing to the dimensions:
                
            png(file="MFA_contribution_DIM1.png")
            
            print(fviz_contrib(res, choice = "quanti.var", axes = 1, top = 20,
                          palette = "jco"))
            
            dev.off()
            
            png(file="MFA_contribution_DIM2.png")
            print(fviz_contrib(res, choice = "quanti.var", axes = 2, top = 20,
              palette = "jco"))
            dev.off()
            
            
            #The most contributing quantitative variables can be highlighted on the scatter 
            #plot using the argument col.var = “contrib”. This produces a gradient colors, 
            #which can be customized using the argument gradient.cols.
            png(file="MFA_most_contributed.png")
            
            
            print(fviz_mfa_var(res, "quanti.var", col.var = "contrib", 
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                          col.var.sup = "violet", repel = TRUE,
                          geom = c("point", "text"), title= "MFA most contributed colored by contributions"))
           
            dev.off()
            
              # Color by cos2 values: quality on the factor map

            #Similarly, you can highlight quantitative variables using their cos2 values 
            #representing the quality of representation on the factor map. If a variable 
            #is well represented by two dimensions, the sum of the cos2 is closed to one. 
            #For some of the row items, more than 2 dimensions might be required to perfectly
            #represent the data.
            
            png(file="MFA_most_contributed_coloredcos2.png")
            print(fviz_mfa_var(res, col.var = "cos2",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                          col.var.sup = "violet", repel = TRUE, title= "MFA most contributed qunatitative variables colored by cos2"))
            
            
            dev.off()             
            
            #Graph of individuals

            ind <- get_mfa_ind(res)
            print(ind)
            
            
            png(file="MFA_individuals.png")
            print(fviz_mfa_ind(res, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE))
            dev.off()
            
            png(file="MFA_individuals_coloredY.png")
            print(fviz_mfa_ind(res, 
              habillage = y_pred, # color by groups 
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              addEllipses = TRUE, ellipse.type = "confidence", 
              repel = TRUE # Avoid text overlapping
              )) 
            
            dev.off()
            
            png(file="MFA_factormap_Y.png")

            print(fviz_ellipses(res, c(y_pred), repel = TRUE))
            dev.off()
            
            png(file="MFA_factormap_all.png")

            print(fviz_ellipses(res, 1:2, geom = "point"))
            
            dev.off()
            
            #Graph of partial individuals

            #The results for individuals obtained from the analysis performed with a single
            # group are named partial individuals. In other words, an individual considered 
            #from the point of view of a single group is called partial individual.
            
            #In the default fviz_mfa_ind() plot, for a given individual, the point corresponds
            # to the mean individual or the center of gravity of the partial points of the 
            #individual. That is, the individual viewed by all groups of variables.
            
            #For a given individual, there are as many partial points as groups of variables.
            png(file="MFA_partial_individuals_ALL.png")

            print(fviz_mfa_ind(res, partial = "all")) 
            dev.off()
            
            png(file="MFA_partial_axes.png")
            
            #Graph of partial axes
            print(fviz_mfa_axes(res))
            
            dev.off()
            
            
            res_MFA=list(eig.val , group$coord, group$cos2, group$contrib,
                          quanti.var, quanti.var$coord, quanti.var$cos2, quanti.var$contrib)
            
            
             return (res_MFA) } ''')
    
        omics_=[]
        names=[]
        for name, lista in self.train_test_datasets.items():
            if name in self.omics:
                #confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        if len(omics_) ==2:
            omics1= omics_[0]
            omics2= omics_[1]
            omic_without_y= robjects.r['omic_without_y']
            omics_no_y=omic_without_y(omics1, omics2, self.y_pred)
            omic1_new=omics_no_y[0]
            omic2_new=omics_no_y[1]
            dataset_concat= pd.read_excel("dataset_concat.xlsx", index_col=0) 
            print(utils.head(dataset_concat))
            nvar_omic1=base.ncol(omic1_new) 
            nvar_omic2=base.ncol(omic2_new)
            nvar_y=1
            print("Number of variables "+ names[0]+ " has: ", nvar_omic1)
            print("Number of variables "+ names[1]+ " has: ", nvar_omic2)
            type__=robjects.StrVector(type_)
            print(type__.r_repr())
            names_=robjects.StrVector(names)
            print(names_.r_repr())
            run_MFA= robjects.r['run_MFA']
            MFA_res=run_MFA(dataset_concat,nvar_omic1,nvar_omic2,nvar_y, type__ ,names_,self.y_pred)
            eig_val=MFA_res[0]
            group_coord=MFA_res[1]
            group_cos2=MFA_res[2]
            group_contrib=MFA_res[3]
            quanti_var= MFA_res[4]
            quanti_coord= MFA_res[5]
            quanti_cos2=MFA_res[6]
            quanti_contrib=MFA_res[7]
            
            file= open("MFA_results.txt", "w")
            file.write(f'Eigen Value:\n {eig_val} \n')
            file.write(f'Group Coordenates: \n {group_coord} \n')
            file.write(f'Group Cos2:\n {group_cos2} \n')
            file.write(f'Group Contribution:\n {group_contrib} \n')
            file.write(f'Quantitative Variables: \n {quanti_var} \n')
            file.write(f'Quantitative Variables coordenates: \n: {quanti_coord} \n')
            file.write(f'Quantitative Variables cos2: \n {quanti_cos2} \n')
            file.write(f'Quantitative Variables contribution: \n {quanti_contrib} \n')
            file.close()
            
            
    def NEMO(self, num_clusters, num_neighbors, k ):
        
        
        robjects.r('''omic_without_y <- function(omic1, omic2, y_pred){
             omic1_new=omic1
             omic2_new=omic2
             omic1_new[,y_pred]=NULL
             omic2_new[y_pred]=NULL
             omics=list(omic1_new,omic2_new)
             return(omics)
            } ''')
        
        
        robjects.r('''run_NEMO <- function(omic1, omic2, dataset_concat, y_pred, num.clusters, num.neighbors, k){
            devtools::install_github('Shamir-Lab/NEMO/NEMO')
            library(NEMO)
            library(SNFtool)
            library(factoextra)
            
             
            # omics.list is a list of data frames, where in each data frame columns are samples and rows are features.
            # note that colnames for each data frame should be set.
            # Two toy datasets, omic1 and omic2, are also loaded whenever NEMO is loaded, but here we read them from a file.
            omic1_ = t(omic1)
            omic2_ = t(omic2)
            omics.list = list(omic1_, omic2_)
            
            # if(num.clusters=="None"){
            #         num.clusters = "NULL"
            #         }
            
            # if(num.neighbors=="None"){
            #         num.neighbors = NULL}
            
            if(k =="None"){
                    k=NA}
            
            if (num.clusters == "None" & num.neighbors== "None"){
                   clustering1 = nemo.clustering(omics.list)
                    # print(clustering)
                    
                    a= table(clustering1, dataset_concat[,y_pred])
                    print(table(clustering1, dataset_concat[,y_pred]))
                    
                    data= cbind(omic1,omic2)
                    object=list(data=data,cluster=clustering1)
                    #plot results of final k-means model
                    png(file="NEMO_clusters.png")
                    print(fviz_cluster(object))
                    dev.off()
              } else if (num.clusters=="None" & num.neighbors != "None" ){
                  
                  
                    # # the number of clusters and number of nearest neighbors used can also be given as input.
                    # # if they are not given, nemo decides the values itself.
                    clustering1=nemo.clustering(omics.list, num.neighbors=num.neighbors)
                    
                    a= table(clustering1, dataset_concat[,y_pred])
                    print(table(clustering1, dataset_concat[,y_pred]))
                    
                    data= cbind(omic1,omic2)
                    object=list(data=data,cluster=clustering1)
                    #plot results of final k-means model
                    png(file="NEMO_clusters.png")
                    print(fviz_cluster(object))
                    dev.off()
            
                  
                  } else if (num.clusters!="None" & num.neighbors == "None" ){
                  
                  
                    # the number of clusters and number of nearest neighbors used can also be given as input.
                    # if they are not given, nemo decides the values itself.
                    clustering1=nemo.clustering(omics.list, num.clusters=num.clusters)
                    
                    a= table(clustering1, dataset_concat[,y_pred])
                    print(table(clustering1, dataset_concat[,y_pred]))
                    
                    data= cbind(omic1,omic2)
                    object=list(data=data,cluster=clustering1)
                    #plot results of final k-means model
                    png(file="NEMO_clusters.png")
                    print(fviz_cluster(object))
                    dev.off()
        
                  }
            
            if (num.clusters!="None" & num.neighbors != "None" ){
                  
                  
                    # the number of clusters and number of nearest neighbors used can also be given as input.
                    # if they are not given, nemo decides the values itself.
                    clustering1=nemo.clustering(omics.list, num.clusters=num.clusters, num.neighbors= num.neighbors)
                    
                    a= table(clustering1, dataset_concat[,y_pred])
                    print(table(clustering1, dataset_concat[,y_pred]))
                    
                    data= cbind(omic1,omic2)
                    object=list(data=data,cluster=clustering1)
                    #plot results of final k-means model
                    png(file="NEMO_clusters.png")
                    print(fviz_cluster(object))
                    dev.off()
        
                  }
           
           
            
        
            
            # k can also be set to NA, in which case nemo chooses its value.
            # nemo.affinity.graph is the integrated affinity graph.
            
            
            affinity.graph = nemo.affinity.graph(omics.list, k = k)
            # print(affinity.graph)
            
            
            # ask nemo to estimate the number of clusters.
            num.clusters_ = nemo.num.clusters(affinity.graph)
            print(num.clusters_)
            
            
            # clustering is the cluster assignment vector, ordered by the columns of affinity.graph.
            
            
            
            clustering = spectralClustering(affinity.graph,2)
            names(clustering) = colnames(affinity.graph)
            # print(clustering)
                       
            
            print(table(clustering, dataset_concat[,y_pred]))
            
            
            data= cbind(omic1,omic2)
            object=list(data=data,cluster=clustering)
            #plot results of final k-means model
            png(file="NEMO_clusters_similarity_graph.png")
            print(fviz_cluster(object))
            dev.off()
            
            
            res_NEMO= list(table(clustering1, dataset_concat[,y_pred]), table(clustering, dataset_concat[,y_pred]))
            
            return(res_NEMO)}''')
            
            
        omics_=[]
        names=[]
        for name, lista in self.train_test_datasets.items():
            if name in self.omics:
                #confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        if len(omics_) ==2:
            omics1= omics_[0]
            omics2= omics_[1]
            omic_without_y= robjects.r['omic_without_y']
            omics_no_y=omic_without_y(omics1, omics2, self.y_pred)
            omic1_new=omics_no_y[0]
            omic2_new=omics_no_y[1]
            dataset_concat= pd.read_excel("dataset_concat.xlsx", index_col=0)
            run_NEMO=robjects.r["run_NEMO"]
            NEMO= run_NEMO(omic1_new, omic2_new, dataset_concat, self.y_pred, num_clusters= num_clusters,num_neighbors=num_neighbors, k=k)
            
            file= open("NEMO_results.txt", "w")
            file.write(f'NEMO clustering with only nemo.clustering():\n {NEMO[0]} \n')
            file.write(f'NEMO clustering with similarity graph: \n {NEMO[1]} \n')
            file.write('See image "NEMO_clusters.png" and "NEMO_clusters_similarity_graph.png" for the cluster plots.')
            file.close()
            
  
    def BCC(self, K, a=1,b=1,IndivAlpha = "FALSE", mu0 = list(),a0=list(),b0=list(),Concentration=1,NumDraws = 1000):
        '''
        Function: 
        
            BCC(X,K=2,a=1,b=1,IndivAlpha = FALSE, mu0 = list(),a0=list(),b0=list(),Concentration=1,NumDraws = 1000)

        Parameters
        ----------
            
        -X is a list of quantitative data sources, X = list(X1,...,XM),  Where each Xi is a D_i X N data matrix (D_i variables on N columns). 
        
        The columns should match for each data source.
        
        -K is the (maximum) number of clusters
        
        -a and b are hyperparameters for the Beta(a,b) prior distribution on alpha
        
        -IndivAlpha indicates whether the alpha should be separate for each data source ("TRUE" or "FALSE")
        
        -mu0 is a list of the prior mean vectors for each data source (i.e. mu0_i is the prior mean for the ith data source)
        
        -a0 and b0 are lists of the Gamma parameters (to determine variance) for each data source
        
        -Concentration is the Dirichlet concentration parameter for the overall cluster sized
        
        -NumDraws is the number of MCMC draws (NumDraws/2 is used as the "burn-in")

        Returns
        -------
        -Alpha: the average adherence (alpha).  If IndivAlpha = TRUE, Alpha[m] gives average adherence for data source m.
        
        -AlphaBounds:95% credible interval for Alpha
        
        -CBest: The "hard" overall clustering, as a binary matrix.  Cbest[i,j] =1 if sample i is in cluster j, 0 otherwise.
        
        -Lbest: A list of the separate clusterings.  Lbest[[m]][i,j] = 1 if sample i is in cluster j for source m, 0 otherwise.
        
        -AlphaVec: Vector of alpha values over MCMC draws, to assess mixing. 

        '''
        
        
        
        
        robjects.r('''omic_without_y <- function(omic1, omic2, y_pred){
             omic1_new=omic1
             omic2_new=omic2
             omic1_new[,y_pred]=NULL
             omic2_new[y_pred]=NULL
             omics=list(omic1_new,omic2_new)
             return(omics)
            } ''')
    
        robjects.r(''' run_BCC <- function(omic1, omic2,dataset_concat,y_pred, K ,a,b,IndivAlpha, mu0 ,a0,b0 ,Concentration,NumDraws){  
                install.packages("gtools")
                library(gtools)
                install.packages("factoextra")
                library("factoextra")
                source('BCC.r') 
                print(IndivAlpha)
                # print(head(t(omic1)))
                # print(head(omic2))
                print(dim(omic1))
                print(dim(omic2))
                # print((dataset_concat[,y_pred]))
                
                if (IndivAlpha =="FALSE"){
                        IndivAlpha = FALSE} 
                else if (IndivAlpha =="TRUE"){
                    IndivAlpha = TRUE
                            }
                print(IndivAlpha)
                if (length(mu0)==0 & length(a0)==0 & length(b0)==0){
                    # print("Here!")
                    # X= list(omic1=as.matrix(t(omic1)),omic2=as.matrix(t(omic2)))
                    # print(head(dataset_concat[,1:407]))
                    # print(head(dataset_concat[,408:814]))
                    X=list(t(omic1),t(omic2))
                    
                    Clusts = BCC(X,K=K,a=a,b=b,IndivAlpha = IndivAlpha,Concentration=Concentration ,NumDraws = NumDraws)
                    print("Clusts done!")
                    rownames(Clusts$Cbest)=rownames(omic1)
                    colnames(Clusts$Cbest)= seq(1:K)
                    # print(Clusts$Cbest) #gives overall clustering 
                    # print(Clusts$Lbest[[1]]) #gives clustering specific to X1
                    # print(table(t(Clusts$Cbest), dataset_concat$Berry))
                    clusters <-colnames(Clusts$Cbest)[which(Clusts$Cbest == 1, arr.ind=T)[, "col"]]
                    df= data.frame ("Clusters"=clusters, row.names =rownames(omic1) )
                    # print(head(df))
                   
                    
                    print(table(df$Clusters, dataset_concat[,y_pred]))
                        
                        }
                else if (length(a0)==0 & length(b0)==0 & length(mu0)!=0){
                     Clusts = BCC(list(t(omic1),t(omic2)),K=K,a=a,b=b,IndivAlpha = IndivAlpha,mu0=mu0,Concentration=Concentration ,NumDraws = NumDraws)
                    rownames(Clusts$Cbest)=rownames(omic1)
                    colnames(Clusts$Cbest)= seq(1:K)
                    # print(Clusts$Cbest) #gives overall clustering  
                    # print(Clusts$Lbest[[1]]) #gives clustering specific to X1
                    # print(table(t(Clusts$Cbest), dataset_concat$Berry))
                    clusters <-colnames(Clusts$Cbest)[which(Clusts$Cbest == 1, arr.ind=T)[, "col"]]
                    df= data.frame ("Clusters"=clusters, row.names =rownames(omic1) )
                    # print(df)
                   
                    
                    print(table(df$Clusters, dataset_concat[,y_pred]))
                    }
                
                else if(length(mu0)==0 & length(a0)!=0 & length(b0)!=0){
                     Clusts = BCC(list(t(omic1),t(omic2)),K=K,a=a,b=b,IndivAlpha = IndivAlpha,a0=a0,b0=b0,Concentration=Concentration ,NumDraws = NumDraws)
                    rownames(Clusts$Cbest)=rownames(omic1)
                    colnames(Clusts$Cbest)= seq(1:K)
                    # print(Clusts$Cbest) #gives overall clustering  
                    # print(Clusts$Lbest[[1]]) #gives clustering specific to X1
                    # print(table(Clusts, dataset_concat$Berry))
                    
                    clusters <-colnames(Clusts$Cbest)[which(Clusts$Cbest == 1, arr.ind=T)[, "col"]]
                   df= data.frame ("Clusters"=clusters, row.names =rownames(omic1) )
                    # print(df)
                   
                    
                    print(table(df$Clusters, dataset_concat[,y_pred]))  
                    }
                    
                
                
                
                
                
                  data= cbind(omic1,omic2)
                    object=list(data=data,cluster=df$Clusters)
                    #plot results of final k-means model
                    png(file="BCC_clusters.png")
                    print(fviz_cluster(object))
                    dev.off()
                
                   
                   return(table(df$Clusters, dataset_concat[,y_pred]))}''')
        omics_=[]
        names=[]
        for name, lista in self.train_test_datasets.items():
            if name in self.omics:
                #confirmar se têm o mesmo y
                # print(utils.head(lista[0]))
                omics_.append(lista[0])
                names.append(name)
        if len(omics_) ==2:
            omics1= omics_[0]
            omics2= omics_[1]
            # print(base.dim(omics1))
            omic_without_y= robjects.r['omic_without_y']
            omics_no_y=omic_without_y(omics1, omics2, self.y_pred)
            omic1_new=omics_no_y[0]
            omic2_new=omics_no_y[1]
            dataset_concat= pd.read_excel("dataset_concat.xlsx", index_col=0)
            
            K=K
            
            run_bcc= robjects.r["run_BCC"]
            BCC_res=run_bcc(omic1_new, omic2_new,dataset_concat,self.y_pred, K, a=a,b=b, IndivAlpha = IndivAlpha,mu0=mu0,a0=a0,b0=b0,Concentration=Concentration,NumDraws = NumDraws)
            
            file= open("BCC_results.txt", "w")
            file.write(f'BCC clustering:\n {BCC_res} \n')
            file.write('See image "BCC_clusters.png" for the cluster plot.')
            file.close()
            
            
            
            
            
            
            
    