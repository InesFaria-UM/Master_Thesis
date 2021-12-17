# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:09:54 2021

@author: Maria Ines
"""
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


class Annotation:
    
    '''A class that reads the annotations.The annotation file for transcriptomics,
    for example, should contain the ID, the UniprotKB e the Annotation of that gene
    
    Input:
        default - True or False whether the user wants the defualt file for annotation or not.
        
        annotation_file - If default is False the file that the user wants to read. 
        The file must have the header with this 3 columns: "ID", "UniprotKB" and "Annotation" 
        
    '''
    
    def __init__(self, default=True, annotation_file="C:/Users/Maria Ines/Desktop/Artigos/Tese/12864_2015_2115_MOESM5_ESM.xls"):
        self.default=default
        self. annotation_file=annotation_file
        
        
        
    def read_annotation(self, header=0, feature_col=0 , skip_rows=1):
        '''
        Default annotation is PN40024 Vitis vinifera reference genome and 
        annotation from http://plants.ensembl.org/Vitis_vinifera/Info/Index

        Parameters
        ----------
        header : INTEGER, Row (0-indexed) to use for the column labels of the parsed
                DataFrame. If a list of integers is passed those row positions will
                be combined into a ``MultiIndex``. Use None if there is no header. The default is 0.
        feature_col : INTEGER, ski. The default is 0.
        skip_rows : INTEGER, Line numbers to skip (0-indexed) or number of lines to skip (int) at the
    start of the file. The default is 1.

        Returns
        -------
        annotation : DATAFRAME
            Dataframe with the columns ID,UNIPROTKB and ANNOTATION to use in differential expression analysis.

        '''
        
        
        
        default="C:/Users/Maria Ines/Desktop/Artigos/Tese/12864_2015_2115_MOESM5_ESM.xls" 
        if self.default == True:
            df=pd.read_excel(default, index_col= feature_col, header=header, skiprows=skip_rows)
            df_r= robjects.conversion.py2rpy(df)
            columns= robjects.IntVector([1,2,4]) #The file must contain
            # base.colnames(df_r)[1]="ID"
            annotation= df_r.rx(True, columns)
            # base.colnames(annotation)[0]="ID"
            # print(utils.head(annotation))
            # print(base.names(annotation))
            
        elif self.default == False:
            df=pd.read_excel(self.annotation_file, index_col= feature_col, header=header, skiprows=skip_rows)
            robjects.r('''ID <- function(annotation){
                annotation$ID_MAISC= rownames(annotation)
                return (annotation)
                }''')
            df_r= robjects.conversion.py2rpy(df)
            print(utils.head(df_r))
            new_id=robjects.r["ID"]
            annotation_ID_col=new_id(df_r)
            print(utils.head(annotation_ID_col))
            columns= robjects.StrVector(["REF","Annotation"]) #The file must contain
            # base.colnames(df_r)[1]="ID"
            annotation= df_r.rx(True, columns)
            
        return annotation