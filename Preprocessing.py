# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 13:50:19 2021

@author: Maria Ines
"""
from sklearn.feature_selection import VarianceThreshold
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.0"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.1.0/bin"
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
utils.install_packages("dplyr")
dplyr = rpackages.importr("dplyr")

robjects.r('install.packages("BiocManager", '
   'repos="http://cran.r-project.org")')
robjects.r('BiocManager::install("genefilter")')

robjects.r('library(genefilter)')

utils.install_packages("nortest")
nortest = rpackages.importr("nortest")

stats = rpackages.importr("stats")
class Preprocessing:
    
    '''
    Preprocessing class allows preprocessing of the omics data that was previously read.
    The following parameters should be given:
        - data : dictionary with data that was previously read.
        -scaling: dictionary that provides information regarding the type of omics data
        (key) and if it is scaled (True of False as value).
    '''
    
    def __init__(self, data, variance_threshold = 0.3):
        self.data=data
        # self.data_wnas= None
        self.data_filter=None
        self.scaled=None
        self.pred_factor= None
        self.variance_threshold=variance_threshold
        
    def na_delete(self):
        for name, dataframe in self.data.items():
            if name.lower() != "metadata":
                 robjects.r(''' na_values <- function(file){
             exist_na= sum(sapply(sapply(file,is.na),sum))
             return (exist_na)}''')
                 na_values= robjects.r["na_values"]
                 has_na= na_values(dataframe)
                 if has_na != 0:
                    robjects.r('''delete_row_nas <- function(df){
                        df_wna <- df[complete.cases(df),]
                        return (df_wna)}''')
                        
                    delete_row_nas= robjects.r["delete_row_nas"]
                    df_w_nas= delete_row_nas(dataframe)
                    self.data.update({name:df_w_nas})
                    print(name, ": Rows with NAs were removed!")
                    print("Dimension of", name, ": ", base.dim(dataframe))
                    print("Dimension of", name, "without NAs: ", base.dim(df_w_nas))
                 elif has_na == 0:
                    print( name, ": Doesn't contain NAs.")
               #      print("Checking if contains rows of zero values...")
               #      robjects.r(''' row0_values <- function(file){
               # exist_rows0=sum( apply(file, 1, function(row) all(row ==0 )))
               # return (exist_rows0)}''')
               #      row0_values= robjects.r["row0_values"]
               #      has_zerorows= row0_values(dataframe)
               #      print(name, ": Has" , has_zerorows ," rows with only zero values.")
               #      if has_zerorows != 0:
               #           robjects.r('''remove_rows <- function(file){
               # # remove_0rows= file[rowSums(file[])!=0,]
               # remove_0rows=file[!!rowSums(abs(file[])),]

               # return (remove_0rows)}''')
               #           remove_rows= robjects.r["remove_rows"]
               #           remove_0rows= remove_rows(dataframe)
               #           self.data.update({name:remove_0rows})
               #           print(name, ": Rows with 0s were removed!")
               #           print("Dimension of", name, ": ", base.dim(dataframe))
               #           print("Dimension of", name, "without zeros: ", base.dim(remove_0rows))
        # print(self.da1ta["Metabolomics"].shape[1])
        return self.data
            
   
    
   
    # def remove_0rows(self):
    #     for name, dataframe in self.data.items():
    #         if name.lower() != "metadata":
    #             robjects.r('''remove_rows <- function(file){
    #          remove_0rows= file[rowSums(file[])>0,])
    #          return (remove_0rows)}''')
    #             remove_rows= robjects.r["remove_rows"]
    #             remove_0rows= remove_rows(dataframe)
    #             print(name, ": Rows with 0s were removed!")
    #             print("Dimension of", name, ": ", base.dim(dataframe))
    #             print("Dimension of", name, "without zeros: ", base.dim(remove_0rows))
        
    
    def filtering (self):
        robjects.r(''' filter_1 <- function(df, n_columns, start_col=1){
              df$Avg_score = rowMeans(df[,start_col:n_columns])
              df_filter_1 = filter(df,Avg_score>1)
              df_filter_1$Avg_score = NULL
              return (df_filter_1)}''' )
        
        robjects.r('''filter_median <- function(df){
                   
                   ## calculate median expression level
                   d= as.matrix(df)
                   cutoff <- median(d)

                   ## TRUE or FALSE for whether each gene is "expressed" in each sample
                   is_expressed <- d > cutoff

                   ## Identify genes expressed in more than 2 samples

                   keep <- rowSums(is_expressed) > 2


                   ## subset to just those expressed genes
                   data_filter <- d[keep,]
                   return (data_filter)}''')
            
        robjects.r('''filter_flat_patterns <- function(df, filter_median){
                   print(dim(df))
                   print(dim(filter_median))
                   sds= rowSds(filter_median)
                   m=median(sds)
                   df_r=df[sds>= 3*median(sds),]
                   return (df_r)}''')
               
        data_filter={}    
        for name, dataframe in self.data.items():
            if name.lower() == "metadata":
                data_filter[name]=dataframe
            
            if name.lower() == "transcriptomics":
                print("Dimension of", name," before filtering:", base.dim(self.data[name]))
                filter_1= robjects.r["filter_1"]
                # print(utils.str(dataframe))
                df_filter_1 = filter_1(dataframe,base.ncol(dataframe))
                print("Dimension df_filter_1:", base.dim(df_filter_1))
                filter_median=robjects.r["filter_median"]
                df_filter_median=filter_median(df_filter_1)
                
                filter_flat_patterns=robjects.r["filter_flat_patterns"]
                df_filter_fp=filter_flat_patterns(df_filter_1,df_filter_median)
                # print(utils.head(df_filter_fp))
                # print(df_filter_fp.shape)
                data_filter[name]=df_filter_fp
                print("Dimension of", name, "after filtering: ", base.dim(data_filter[name]))
            if name.lower() == "metabolomics":
                # print("Dimension of", name," before filtering:", base.dim(self.data[name]))
                # filter_1= robjects.r["filter_1"]
                # df_filter_1 = filter_1(dataframe,base.ncol(dataframe),start_col=2)
                # print("Dimension of", name, "after filtering: ", base.dim(df_filter_1))
                # data_filter[name]= df_filter_1
                if len(dataframe) < 300:
                    print("Number of features in metabolomics is already considered small.")
                    data_filter[name]= dataframe
                if len(dataframe) > 300: # should I use xcms? I need to try it first!
                    continue
                
            if name.lower() == "fluxomics":
                print("Dimension of", name," before filtering:", base.dim(self.data[name]))
                rownames= base.rownames(dataframe)
                from sklearn.feature_selection import VarianceThreshold
                dataframe_py=robjects.conversion.rpy2py(dataframe)
                selection = VarianceThreshold(threshold=self.variance_threshold)
                selection.fit(dataframe_py)
                
                dataframe_selected = pd.DataFrame(selection.transform(dataframe_py), columns=dataframe_py.columns[selection.get_support()])
                # print(dataframe_selected.head())
                # print(dataframe_selected.shape)
                dataframe_r=robjects.conversion.py2rpy(dataframe_selected)
                robjects.r('''change_rownames <- function(df,rownames){
                   
                  rownames(df)=rownames
                   return (df)}''')
                change_rownames=robjects.r["change_rownames"]  
                dataframe_r=change_rownames(dataframe_r,rownames)
                # print(utils.head(dataframe_r))
                data_filter[name]= dataframe_r
                # print("Dimension of", name," before filtering:", base.dim(self.data[name]))  
                # print("Dimension of", name, "after filtering: ", base.dim(data_filter[name]))
                print("Dimension of", name, "after filtering: ", base.dim(data_filter[name]))
                
        # print(len(data_filter))
        self.data_filter=data_filter
        return data_filter
    
    
    
    def normalization (self, scaled):
        self.scaled=scaled
        data_norm={}
        # robjects.r(''' scale <- function(df){
        #     df_scaled= scale(df)
        #     return (df_scaled)}''')
        
        
        for name, data in self.data_filter.items():
            if name.lower() =="metadata":
                data_norm[name]=data
            
            if name.lower() != "metadata":
                scal=self.scaled[name]
                
                if scal == False:
                    print( name, "doesn't follow a normal distribution. Scaling will be executed!")
                    # scale= robjects.r["scale"]
                    rownames= base.row_names(data)
                    colnames=base.colnames(data)
                    df_scaled=base.scale(data)
                    df_scaled_df=base.as_data_frame(df_scaled, row_names=rownames, col_names=base.as_vector(colnames))
                    # df_scaled_df= base.as_data_frame(df_scaled)
                    df_scaled_df=stats.setNames(df_scaled_df,colnames)
                    # print(utils.head(df_scaled_df))
                    data_norm[name]=df_scaled_df
                    
                if scal == True:
                    print(name," is already scaled. Scaling will not be executed!")
                    data_norm[name]=data
        # print(utils.head(self.data_filter["Metabolomics"]))
        # print(utils.head(data_norm["Metabolomics"]))
        return data_norm
                    
                
    def data_discretization(self, data_norm, col_index, cut, labels, name_variable):
        '''
        

        Parameters
        ----------
        cut : list of  two or more unique cut points or a 
        single number (greater than or equal to 2) giving the number of 
        intervals into which the data is to be cut.
        labels :list of labels for the levels of the resulting category.
        col_index : a numeric vector or string name which is to be converted to a factor by cutting.
        
        Returns
        -------
        Numeric Vector converted into categorical vector

        '''
        robjects.r('''change_col <- function(data, name_variable){
            colnames(data)[length(colnames(data))] = name_variable
            return (data)}''')
            
        robjects.r('''factor_pred <- function(data, col_index, cut, labels, name_variable){
            name_variable= cut(data[, col_index], cut, labels, right = FALSE)
            data = cbind(data, name_variable)
            # data= data.frame(metadata$`Time Point`, berry)
            return (data)}''')   
        for name, data in data_norm.items():
            if name.lower() == "metadata":
                breaks=robjects.IntVector(cut)
                label=robjects.StrVector(labels)
                
                # dados= data.rx(True,col_index)
                # pred_factor= base.cut(dados, breaks, label)
                # cbind=robjects.r["cbind"]
                # print(pred_factor)
                # pred_factor_f= robjects.FactorVector(pred_factor)
                # print(pred_factor_f)
                # data_metadata=base.cbind(data,pred_factor_f)
                # print(utils.head(data_metadata))
                # self.data=robjects.DataFrame(self.data)
                
                factor_pred=robjects.r["factor_pred"]
                data=factor_pred(data,col_index,breaks,label,name_variable)
                
                
                change_col= robjects.r["change_col"]
                data_metadata=change_col(data, name_variable)
                # base.colnames(data)[base.length(base.names(data))]<- name_variable
                # print(data_metadata)
                data_norm[name]=data_metadata
                
            # print(self.pred_factor)
            return data_norm
                
