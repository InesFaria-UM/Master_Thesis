# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:27:32 2021

@author: Maria Ines
"""
from Read_data import Read_data
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.0"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.1.0/bin"
os.environ["R_USER"]= r"C:\Users\Maria Ines\Anaconda3\Lib\site-packages\rpy2"
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import rpy2.robjects.numpy2ri as rpyn
from rpy2.robjects.lib import grid # install any dependency package if you get error like "module not found"
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

# utils.install_packages("nortest")
# nortest = rpackages.importr("nortest")

stats = rpackages.importr("stats")
graphics= rpackages.importr("graphics")

utils.install_packages("viridis")
viridis= rpackages.importr("viridis")

utils.install_packages("ggplot2")
ggplot2= rpackages.importr("ggplot2")

grdevices = importr('grDevices')

utils.install_packages("pheatmap")
pheatmap = importr('pheatmap')

stats = importr('stats')

utils.install_packages("ggrepel")
ggrepel= rpackages.importr("ggrepel")
# import rpy2.robjects.lib.ggplot2 ass ggplot2

robjects.r('BiocManager::install("limma")')

robjects.r('library(limma)')
utils.install_packages("limma")
limma= rpackages.importr("limma")

# utils.remove_packages("Rcpp")
# utils.install_packages("C:/Users/Maria Ines/Desktop/Artigos/Tese/Rcpp_1.0.7.zip", type = "source")
# Rcpp=rpackages.importr("Rcpp")



class Analysis:
    
    '''
    The Analysis class will execute the exploratory analysis of the data.
    The arguments that should be provided are:
        data_norm: Dictionary with data that was previously normalized.
        items: dictionary with the type o exploratory analysis the person wants to do and the column name or index.
    '''
    
    def __init__(self, data_norm_final, items):
        self.data_norm_final= data_norm_final
        self.items = items
        self.metadata= None
        
        
    def exploratory(self):
        
        robjects.r('''barplot <- function(dataf,nlevels,i){
            # install.packages("viridis")
            # library(viridis)
            p <- ggplot(dataf, aes(x=X, y=hline)) + geom_bar(stat="identity", fill= viridis::viridis(nlevels, begin=0.2, end = 0.8)) + xlab(" ") + ylab("Counts") + ggtitle(i)

            return (p) }''')
            
        robjects.r('''boxplot <- function(dataf,column){
            # install.packages("viridis")
            # library(viridis)
            p <- ggplot(data, aes(x=column)+ ggplot2.geom_boxplot(fill= viridis::viridis(1, begin=0.2, end = 0.8))) + xlab(" ") + ylab(column) 

            return (p) }''')
            
        robjects.r('''boxplot_ggplot <- function(data,X,Y,x,y){
            # install.packages("viridis")
            # library(viridis)
            print(summary(data[,"Berry Weight (g / berry)"]))
            p = ggplot(data, aes(x = X, y = Y, fill= X))+geom_boxplot(fill=viridis::viridis(2, begin=0.2, end = 0.8)) + ggtitle(y)
            #+ xlab(toString(x))+ ylab(toString(y))+ ggtitle(toString(y))) 
            return (p)}''')
            
        robjects.r('''plot_ggplot <- function(data,X,Y,x,y,color){
             
            p=ggplot(data, aes(x=X, y= Y, colour=data[,color])) + geom_point() +  ylab(y)+ xlab(x) + scale_colour_manual(values = c("#0073C2FF", "#EFC000FF"))
            
            
            return (p)}''')
            
            
        for name, data in self.data_norm_final.items():
            names=base.colnames(data)
            data= robjects.DataFrame(data)
            for a in base.range(base.length(names)):
                for type_ in self.items.keys():
                    if type_ == "boxplot":
                        column= self.items[type_]
                        for i in column:
                            if isinstance(i, tuple):
                                # print("Tuple")
                                # print(i)
                                #formula 2º is always the thing like species
                            
                                if i[0] in names:
                                    
                                        y=i[0]
                                        x= i[1]
                                        # y=data.rx(True,y)
                                        X=data.rx(True,x)
                                        Y=data.rx(True,y)
                                        # print(X)
                                        # print(y_vec)
                                        # robjects.r[''' boxplot <- function(data,x_vec,y_vec,column){
                                        #     p <-boxplot(x_vec~factor(y_vec),col=c("darkolivegreen2","darkolivegreen3"), main="Boxplot of "+ column, ylab=column)
                                        #     return (p) } ''']
                                        # graphics.boxplot(Formula(' x_vec~ factor(y_vec) '), col=c("darkolivegreen2","darkolivegreen3"), main="Boxplot of "+ column, ylab=column)
                                        grdevices = importr('grDevices')
            
                                        grdevices.png(file="boxplot_tuple.png", width=512, height=512)
                                        boxplot_ggplot= robjects.r["boxplot_ggplot"]
                                        print(boxplot_ggplot(data, X, Y, x,y))
                                        grdevices.dev_off()
                                    
                                
                                
                    
                            if not isinstance(i, tuple):
                                # print("Not tuple")
                                # print(i)
                                if i in names:
                                   vector=data.rx(True,i)
                                   # print(vector)
                                   grdevices = importr('grDevices')
           
                                   grdevices.png(file="boxplot.png", width=512, height=512)
                                   graphics.boxplot(vector, col="#404788FF",xlab="X", ylab= "Y", main= i)
                                   grdevices.dev_off()
                           
                            else:
                            
                               if i in names:
                                   # print(i)
                                   vector=data.rx(True,i)
                                   grdevices = importr('grDevices')
                       
                                   grdevices.png(file="boxplot"+ i + " .png", width=512, height=512)
                                   graphics.boxplot(vector)
                                   grdevices.dev_off()
                                   # x_dados= robjects.FloatVector(column)
                                   # graphics.boxplot(vector, col="#0073C2FF", main="Boxplot of "+ column, ylab=column)
                            
                        
                    
                    if type_ == "barplot":
                        column= self.items[type_]
                        
                        if isinstance(column, list):
                            for i in column:
                                # print(i)
                                vector=data.rx(True,i)
                                # print(vector)
                                as_factor=robjects.r["factor"]
                                levels=robjects.r["levels"]
                                factor_= as_factor(vector)
                                
                                # print(factor_)
                                # level=levels(factor_)
                                # nlevels=base.length(levels(factor_))
                                
                                factor=robjects.FactorVector(factor_)
                                # print(factor)
                                
                                
                                nlevels=factor.nlevels
                                # print(factor.nlevels)
                                levels_name=factor.levels
                                # print(factor.levels)
                                b=base.toString(i)
                                # print(b)
                                # print(robjects.StrVector(factor_vintage.levels))
                                
                                d = {"X": robjects.StrVector((levels_name)), "hline": base.as_vector(base.table(vector))}
                                dataf = robjects.DataFrame(d)
                                # hline=robjects.DataFrame(X=factor.levels, hline=base.as_vector(base.table(vector)))
                                print(dataf)
                                
                                grdevices = importr('grDevices')
        
                                
          
                                barplot=robjects.r['barplot']
                                grdevices.png(file="barplot_"+ i + ".png", width=512, height=512)
                                print(barplot(dataf,nlevels, i))
                                grdevices.dev_off()
                                
                                
                    if type_ =="plot":
                       column= self.items[type_]
                       for i in column:
                           if isinstance(i, tuple):
                                # print("Tuple")
                                # print(i)
                                #formula 2º is always the thing like species
                            
                                if i[0] in names:
                                        # print(base.colnames(data))
                                        y=i[1]
                                        x= i[0]
                                        color=i[2]
                                        
                                        # Y=data.rx(True,y)
                                        # print(x)
                                        X=data.rx(True,x)
                                        Y=data.rx(True,y)
                                        
                                        print(X)
                                        print(Y)
                                        # color=data.rx(True,color)
                                        # print(x_vec)
                                        # print(y_vec)
                                        # robjects.r[''' boxplot <- function(data,x_vec,y_vec,column){
                                        #     p <-boxplot(x_vec~factor(y_vec),col=c("darkolivegreen2","darkolivegreen3"), main="Boxplot of "+ column, ylab=column)
                                        #     return (p) } ''']
                                        # graphics.boxplot(Formula(' x_vec~ factor(y_vec) '), col=c("darkolivegreen2","darkolivegreen3"), main="Boxplot of "+ column, ylab=column)
                                        grdevices = importr('grDevices')
            
                                        grdevices.png(file="plot.png", width=512, height=512)
                                        plot_ggplot= robjects.r["plot_ggplot"]
                                        print(plot_ggplot(data,X,Y, x, y, color))
                                        grdevices.dev_off()
                                     
                                
                                
                    
            return print("Barplot and Boxplots were executed, see directory!")
                    
                
                
                
    def summary_data(self):
        for name, data in self.data_norm_final.items():
            if name.lower() != "metadata":
                print(name+ ":")
                print(utils.head(base.summary(data)))
        return
    
    
    
    def pheatmap(self, which_data, columns=None):
        data_samples = robjects.DataFrame({})
        # for name, data in self.data_norm.items():
        robjects.r('''pheatmap_simple <- function(data){
                corMatrix <- cor(data,use="c",method="spearman")
                # colnames(corMatrix )= rownames(data)
                pheatmap(corMatrix)
                return (pheatmap)}''')
            
        robjects.r('''pheatmap_data_samples <- function(data,data_samples){
                # data=data[complete.cases(data), ]
                corMatrix <- cor(data,use="c",method="spearman")
                # print(head(corMatrix))
                # print(rownames(data_samples))
                colnames(corMatrix)= rownames(data_samples)
                # data_samples= as.data.frame(data_samples)
                # rownames(data_samples) <- colnames(corMatrix)
                # print(dim(corMatrix))
                # print(data_samples)
                pheatmap(corMatrix,annotation_col=data_samples) 
                
                return (pheatmap)}''')
        
        
        if which_data == "ALL":
            if columns == None:
                for name, data in self.data_norm_final.items():
                    if name.lower() != "metadata":
                        # print(utils.head(self.data_norm[name]))
                        # print(base.rownames(data))
                        as_data_frame=robjects.r["as.data.frame"]
                        pheatmap_simple= robjects.r["pheatmap_simple"]
                        grdevices.png(file="pheatmap_ALL_None_"+ name.lower()+".png", width=512, height=512)
                        # pheatmap= rpackages.importr("pheatmap")
                        # pheatmap.pheatmap(stats.cor(as_data_frame(data), use="c", method= "spearman"))
                        print(pheatmap_simple(data))
                        # print(pheatmap)
                        grdevices.dev_off()
                   
            
            if columns != None:
                if len(columns)!=1:
                    column = robjects.StrVector(columns)
                    # print(column)
                    for name, data in self.data_norm_final.items():
                        data=robjects.DataFrame(data)
                        names=base.colnames(data)
                        # for a in base.range(base.length(names)):
                        if name.lower() == "metadata": 
                            for i in column:
                                for a in base.range(base.length(names)):
                                    if i in names:
                                        
                                        sample= data.rx(True,column)
                                        # print(sample)
                                        
                                        robjects.r(''' sampleInfo <- function(column,data,sample){
                                            for (i in column){
                                              sample[,i]=as.factor(sample[,i])     
                                              # a=as.factor(data[,i])
                                              # print(a)
                                              # sampleInfo=data.frame(column[i]=a)
                                            }
                                            
                                            # print(sample)
                                            return (sample)}''')
                                        sampleInfo_fun=robjects.r["sampleInfo"]
                                        sampleInfo=sampleInfo_fun(column,data,sample)
                                        # print(sampleInfo)
                        #                 as_data_frame=robjects.r["as.data.frame"]
                        #                 as_factor=robjects.r["factor"]
                        #                 factor_i=as_factor(sample[i])
                        #                 factor_i= robjects.FactorVector(factor_i)
                        #                 print(factor_i)
                        #                 data_samples_py=robjects.conversion.rpy2py(data_samples)
                        #                 if data_samples_py.empty is True: 
                        #                     data_samples = robjects.DataFrame({i: factor_i})
                        #                 else:
                        #                     data_samples=base.cbind(data_samples, factor_i)
                            
                        # data_samples=stats.setNames(data_samples, base.colnames(sample))
                        # rownames= robjects.r["rownames"]
                        # rownames_data=rownames(sample)
                        # data_samples=base.as_data_frame(data_samples, row_names=rownames_data)
                        # print(data_samples)
                        # # data_samples_df= robjects.DataFrame(data_samples)
                        # sampleInfo= data_samples.rx(True,column)
                        # print(sampleInfo)
                        if name.lower() != "metadata":
                            pheatmap_data_samples= robjects.r["pheatmap_data_samples"]
                           
                            grdevices.png(file="pheatmap_ALL_column_"+ name.lower()+".png", width=512, height=512)
                            pheatmap= pheatmap_data_samples(data, sampleInfo)
                            grdevices.dev_off()
                
                if len(columns)==1:  
                   column = robjects.StrVector(columns)
                   for name, data in self.data_norm_final.items():
                        data=robjects.DataFrame(data)
                        names=base.colnames(data)
                        # for a in base.range(base.length(names)):
                        if name.lower() == "metadata": 
                            for i in column:
                                for a in base.range(base.length(names)):
                                    if i in names:
                                       
                                      
                                        sample=data.rx(True,column)
                                        # print(sample)
                                        as_data_frame=robjects.r["as.data.frame"]
                                        as_factor=robjects.r["factor"]
                                        factor_i= as_factor(sample)
                                        # print(factor_i)
                                        robjects.r('''df <- function(data,column){
                                            y=factor(data[,column])
                                            y=as.data.frame(y)
                                            rownames(y)=rownames(data)
                                            colnames(y)=column
                                            return(y)
                                            
                                        
                                            }''')
                                        df=robjects.r["df"]
                                        sampleInfo=df(data,column)
                                        
                        if name.lower() != "metadata":
                            pheatmap_data_samples= robjects.r["pheatmap_data_samples"]
                           
                            grdevices.png(file="pheatmap_ALL_column_"+ name.lower()+".png", width=512, height=512)
                            print(pheatmap_data_samples(data, sampleInfo))
                            
                            grdevices.dev_off()
        if which_data != "ALL":
            if columns == None:
                for name, data in self.data_norm_final.items():
                    if name.lower == which_data.lower():
                        pheatmap_simple= robjects.r["pheatmap_simple"]
                        grdevices.png(file="pheatmap_"+ name.lower()+"_None.png", width=512, height=512)
                        pheatmap=pheatmap_simple(self.data_norm[name])
                        grdevices.dev_off()
                        # print(pheatmap)
        
            if columns != None:
                print(len(columns))
                # factors= []
                if len(columns)!= 1:
                    column = robjects.StrVector(columns)
                    # print(column)
                    for name, data in self.data_norm_final.items():
                        data=robjects.DataFrame(data)
                        names=base.colnames(data)
                        if name.lower() == "metadata": 
                            for i in column:
                                for a in base.range(base.length(names)):
                                    if i in names:
                                        
                                        sample= data.rx(True,column)
                                        as_data_frame=robjects.r["as.data.frame"]
                                        as_factor=robjects.r["factor"]
                                        factor_i=as_factor(sample[i])
                                        factor_i= robjects.FactorVector(factor_i)
    
                                        data_samples_py=robjects.conversion.rpy2py(data_samples)
                                        if data_samples_py.empty is True: 
                                            data_samples = robjects.DataFrame({i: factor_i})
                                        else:
                                            data_samples=base.cbind(data_samples, factor_i)
                            
                        data_samples=stats.setNames(data_samples, base.colnames(sample))
                        rownames= robjects.r["rownames"]
                        rownames_data=rownames(sample)
                        data_samples=base.as_data_frame(data_samples, row_names=rownames_data)
                        # print(data_samples)
                        data_samples_df= robjects.DataFrame(data_samples)
                        sampleInfo= data_samples_df.rx(True,column)
                        
                        pheatmap_data_samples= robjects.r["pheatmap_data_samples"]
                        grdevices.png(file="pheatmap_"+ which_data.lower()+"_columns.png", width=512, height=512)
                        for name, data in self.data_norm.items():
                            if name.lower() == which_data.lower():
                                print(pheatmap_data_samples(data, data_samples))
                                grdevices.dev_off()
                if len(columns)==1:
                    column = robjects.StrVector(columns)
                    # print(column)
                    for name, data in self.data_norm_final.items():
                        data=robjects.DataFrame(data)
                        names=base.colnames(data)
                        if name.lower() == "metadata": 
                            for i in column:
                                for a in base.range(base.length(names)):
                                    if i in names:
                                        
                                        sample=data.rx(True,column)
                                        # print(sample)
                                        as_data_frame=robjects.r["as.data.frame"]
                                        as_factor=robjects.r["factor"]
                                        factor_i= as_factor(sample)
                                        # print(factor_i)
                                        robjects.r('''df <- function(data,column){
                                            y=factor(data[,column])
                                            y=as.data.frame(y)
                                            rownames(y)=rownames(data)
                                            colnames(y)=column
                                            return(y)
                                            
                                        
                                            }''')
                                        df=robjects.r["df"]
                                        sampleInfo=df(data,column)
                            
                        # data_samples=stats.setNames(data_samples, base.colnames(sample))
                        # rownames= robjects.r["rownames"]
                        # rownames_data=rownames(sample)
                        # data_samples=base.as_data_frame(data_samples, row_names=rownames_data)
                        # print(data_samples)
                        # data_samples_df= robjects.DataFrame(data_samples)
                        # sampleInfo= data_samples_df.rx(True,column)
                        
                        pheatmap_data_samples= robjects.r["pheatmap_data_samples"]
                        grdevices.png(file="pheatmap_"+ which_data.lower()+"_columns.png", width=512, height=512)
                        for name, data in self.data_norm_final.items():
                            if name.lower() == which_data.lower():
                                pheatmap=pheatmap_data_samples(data, sampleInfo)
                                grdevices.dev_off()
        return print("Pheatmaps were created, see directory!")     
       
        
    def pca(self, which_data, variable):
        for name, data in self.data_norm_final.items():
            if name.lower() == "metadata":
                data=robjects.DataFrame(data)
                self.metadata = data
        
        # pcas= robjects.StrVector(pcas)
        # print(pcas[0])
        # print(pcas[1])
        robjects.r('''pca_simples <- function(metadata, data, variable){
            install.packages("ggrepel")
            library("ggrepel")
            install.packages("Rcpp")
            library(Rcpp)
            pca <- prcomp(t(data))
            summ=summary(pca)
            print(summ$importance[2,1:4])
            # print(factor(metadata[,variable]))
            nlevels= nlevels(factor(metadata[,variable]))
            # print(nlevels)
            # Join the PCs to the sample information
            print(cbind(metadata, pca$x) %>% 
            ggplot(aes(x =PC1, y=PC2, col=factor(metadata[,variable]),label=paste( variable))) + geom_point() + #geom_text_repel()+
            stat_ellipse()+scale_colour_manual(values =viridis::viridis(nlevels, begin=0.2, end = 0.8))+ ggtitle(paste0("Fig. PC1-PC2", variable))+
            xlab(paste0("PC1:", (summ$importance[2,1]*100), "% variance")) + 
            ylab(paste0("PC2:", (summ$importance[2,2]*100), "% variance"))) }''')
        
            
            
        robjects.r(''' pca_complex <- function(metadata, data, shape, color){
                   pca <- prcomp(t(data))
                   summ=summary(pca)
                   print(summ$importance[2,1:4])
                   
                   nlevels= nlevels(factor(metadata[,color]))
                   print(cbind(metadata, pca$x) %>% 
                   ggplot(aes(PC1, PC2, shape = factor(metadata[,shape]), colour = factor(metadata[,color]))) + geom_point(size=3) +
                   stat_ellipse()+
                   scale_shape_manual(values = c(0, 16)) + scale_colour_manual(values = viridis::viridis(nlevels, begin=0.2, end = 0.8)) +
                   xlab(paste0("PC1: ", (summ$importance[2,1]*100), "% variance")) + 
                   ylab(paste0("PC2: ", (summ$importance[2,2]*100), "% variance")))}''')   
            
            
            
            
        robjects.r(''' pca_color <- function(metadata, data, color){
                  
            pca <- prcomp(t(data))
                   summ=summary(pca)
                   print(summ$importance[2,1:4])
                   
                   nlevels= nlevels(factor(metadata[,color]))
                   print(cbind(metadata, pca$x) %>% 
                   ggplot(aes(PC1, PC2, colour = factor(metadata[,color]))) + geom_point(size=3) +
                   stat_ellipse()+
                   scale_colour_manual(values = viridis::viridis(nlevels, begin=0.2, end = 0.8)) +
                   xlab(paste0("PC1: ", (summ$importance[2,1]*100), "% variance")) + 
                   ylab(paste0("PC2: ", (summ$importance[2,2]*100), "% variance")))}''')
            
            
            
        robjects.r(''' pca_shape <- function(metadata, data, shape){
                   pca <- prcomp(t(data))
                   summ=summary(pca)
                   print(summ$importance[2,1:4])
                   
                   
                   print(cbind(metadata, pca$x) %>% 
                   ggplot(aes(PC1, PC2, shape = factor(metadata[,shape]),)) + geom_point(size=3) +
                   stat_ellipse()+
                   scale_shape_manual(values = c(0, 16))  +
                   xlab(paste0("PC1: ", (summ$importance[2,1]*100), "% variance")) + 
                   ylab(paste0("PC2: ", (summ$importance[2,2]*100), "% variance")))}''')
            
        if which_data == "ALL":
            for name, data in self.data_norm_final.items():
                data=robjects.DataFrame(data)
                
                if name.lower() != "metadata":
                     if not isinstance(variable, dict):
                        pca_simples= robjects.r["pca_simples"]
                        grdevices.png(file="pca_"+ name.lower()+ "_" + variable +".png", width=512, height=512)
                        pca= pca_simples(self.metadata, data, variable)
                        grdevices.dev_off()
                        
                     if isinstance(variable,dict):
                        for x, var in variable.items():
                            
                            if x =="shape":
                                shape= variable["shape"]
                            if x== "color":
                                color= variable["color"]
                        if shape == None:
                            pca_color=robjects.r["pca_color"]
                            grdevices.png(file="pca_color_"+ name.lower()+ ".png", width=512, height=512)
                            pca= pca_color(self.metadata, data, color)
                            grdevices.dev_off()
                        
                        if color == None:
                            pca_shape=robjects.r["pca_shape"]
                            grdevices.png(file="pca_shape_"+ name.lower()+ ".png", width=512, height=512)
                            pca= pca_shape(self.metadata, data, shape)
                            grdevices.dev_off()
                            
                        if color != None and shape != None:
                            pca_complex=robjects.r["pca_complex"]
                            grdevices.png(file="pca_complex_"+ name.lower()+ ".png", width=512, height=512)
                            pca= pca_complex(self.metadata, data, shape,color)
                            grdevices.dev_off()
                        
        
        if which_data != "ALL":
            for name, data in self.data_norm_final.items():
                data=robjects.DataFrame(data)
                if name.lower() == which_data.lower():
                    if not isinstance(variable, dict):
                        pca_simples= robjects.r["pca_simples"]
                        grdevices.png(file="pca_"+ which_data.lower()+ "_" + variable +".png", width=512, height=512)
                        pca= pca_simples(self.metadata, data, variable)
                        grdevices.dev_off()
                        
                    if isinstance(variable,dict):
                        for x, var in variable.items():
                            if x =="shape":
                                shape= variable["shape"]
                            if x== "color":
                                color= variable["color"]
                        pca_complex=robjects.r["pca_complex"]
                        grdevices.png(file="pca_complex_"+ name.lower() + "_.png", width=512, height=512)
                        pca= pca_complex(self.metadata, data, shape,color)
                        grdevices.dev_off()
        return print("PCAs were created, see directory!") 
    
    
    
    
    def differential_expression(self, annotation, which_data, y_pred, filter_results, default=True):
        robjects.r('''dif_exp <- function(annotation, metadata, data, y_pred){
                   
                    y_pred= metadata[,y_pred]
                    design_mv <- model.matrix(~0+y_pred)
                   

                    colnames(design_mv) <- c("PreV", "PostV")
                    # print(design_mv)
                    # print(head(data))

                     # fit_mv <- lmFit(data, design_mv) 
                     # # print(head(fit_mv$coefficients))

                     # contrasts_mv <- makeContrasts(PreV - PostV, levels=design_mv)


                     # fit2_mv <- contrasts.fit(fit_mv, contrasts_mv)
                     # fit2_mv <- eBayes(fit2_mv)
                     # # print(topTable(fit2_mv))
                    
                     # decideTests(fit2_mv)
                     # # print(table(decideTests(fit2_mv)))
                    
                    
                     aw <- arrayWeights(data,design_mv)
                    
                     fit_mv <- lmFit(data, design_mv,
                                   weights = aw)
                     contrasts_mv <- makeContrasts(PreV - PostV, levels=design_mv)
                     fit2_mv <- contrasts.fit(fit_mv, contrasts_mv)
                     fit2_mv <- eBayes(fit2_mv)
                     print(table(decideTests(fit2_mv)))

                     #Annotation
                    
                     t_mv= topTable(fit2_mv, sort.by = "P")
                     top_mv=rownames(t_mv)
                    
                     top10_mv=as.data.frame(annotation[top_mv,])
                     # print(top10_mc)
                    
                     top1000= topTable(fit2_mv,sort.by = "P", number=1000)
                     top_mv_1000 =rownames(top1000)
                    
                     top50= topTable(fit2_mv,sort.by = "P", number=50)
                     top_mv_50 =rownames(top50)
                    
                     top212= topTable(fit2_mv,sort.by = "P", number=212)
                     top_mv_212=rownames(top212)
                    
                     full_results_mv <- topTable(fit2_mv, number=Inf)
                     full_results_mv <- tibble::rownames_to_column(full_results_mv,"ID")

                     rownames(full_results_mv) <- full_results_mv$ID
                     topfull_mv=rownames(full_results_mv)
                    
                     full_results_mv["UniProtKB"]= annotation[topfull_mv,"UniProtKB"]
                     full_results_mv["Annotation"]=annotation[topfull_mv,"Annotation"]
                    
                      print(head(full_results_mv,30))
                      write.table(full_results_mv, "full_results_mv.csv", sep=",", col.names=TRUE, quote=FALSE, row.names=TRUE,dec = ".")
                  
                    # install.packages("xlsx")  
                    # library(xlsx)
                    # write.xlsx(x = full_results_mv, file = "full_results_mv.xlsx", row.names=TRUE, col.names=TRUE)
                  
                    
                    # fit_ <- eBayes(fit_mv) 
                  
                    # full_results_PreV <- topTable(fit_,coef="PreV",sort.by="P", number=Inf)
                    
                    #  full_results_PreV <- tibble::rownames_to_column(full_results_PreV,"ID")

                    #  rownames(full_results_PreV) <- full_results_PreV$ID
                    #  topfull_mv_PreV=rownames(full_results_PreV)
                    
                    #  full_results_PreV["UniProtKB"]= annotation[topfull_mv_PreV,"UniProtKB"]
                    #  full_results_PreV["Annotation"]=annotation[topfull_mv_PreV,"Annotation"]
                    
                    #   # print(head(full_results_PreV,30))
                    #   write.table(full_results_PreV, "full_results_PreV.csv", sep=",", col.names=TRUE, quote=FALSE, row.names=TRUE,dec = ".")
                    
                  
                    # full_results_PostV <- topTable(fit_,coef="PostV",sort.by="P", number=Inf)
                    
                    #  full_results_PostV <- tibble::rownames_to_column(full_results_PostV,"ID")

                    #  rownames(full_results_PostV) <- full_results_PostV$ID
                    #  topfull_mv_PostV=rownames(full_results_PostV)
                    
                    #  full_results_PostV["UniProtKB"]= annotation[topfull_mv_PostV,"UniProtKB"]
                    #  full_results_PostV["Annotation"]=annotation[topfull_mv_PostV,"Annotation"]
                    
                    #   # print(head(full_results_PostV,30))
                    #   write.table(full_results_PostV, "full_results_PostV.csv", sep=",", col.names=TRUE, quote=FALSE, row.names=TRUE,dec = ".")
                    
                  
                    
                  
                    
                  
                    
                  
                    
                     # Make sure you have ggplot2 loaded
                     ggplot(full_results_mv, aes(x = logFC, y=B)) + geom_point()
                    
                    
                     ## change according to your needs
                     p_cutoff <- 0.05
                     fc_cutoff <- 1
                    
                     full_results_mv %>% 
                       mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
                       ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
                    
                     p_cutoff <- 0.05
                     fc_cutoff <- 1
                     topN <- 20
                    
                    png(filename="full_results_mv.png")

 
                    print(full_results_mv %>% 
                     mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
                     mutate(Rank = 1:n(), Label = ifelse(Rank < topN, UniProtKB,"")) %>% 
                      ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point()) #+ geom_text_repel(col="black"))
                    
                    dev.off()


                     topN <- 20
                     ##
                     ids_of_interest <- mutate(full_results_mv, Rank = 1:n()) %>% 
                       filter(Rank < topN) %>% 
                       pull("ID")
                    
                     gene_names <- mutate(full_results_mv, Rank = 1:n()) %>% 
                       filter(Rank < topN) %>% 
                       pull(UniProtKB) 
                    
                    
                     ## Get the rows corresponding to ids_of_interest and all columns
                     gene_matrix <- data[ids_of_interest,]
                    
                     pheatmap_de= pheatmap(gene_matrix,
                               labels_row = annotation[ids_of_interest,1],
                               scale="row")
                   
                   
                   
                   return (pheatmap_de)}
                   ''')

    
        
        
        robjects.r('''top_mv <- function(annotation, metadata, data, y_pred, filter_results){
                   y_pred= metadata[,y_pred]
                 design_mv <- model.matrix(~0+y_pred)
                
 
                 colnames(design_mv) <- c("PreV", "PostV")
                 # print(design_mv)
                 # print(head(data))
 
                  # fit_mv <- lmFit(data, design_mv) 
                  # # print(head(fit_mv$coefficients))
 
                  # contrasts_mv <- makeContrasts(PreV - PostV, levels=design_mv)
 
 
                  # fit2_mv <- contrasts.fit(fit_mv, contrasts_mv)
                  # fit2_mv <- eBayes(fit2_mv)
                  # # print(topTable(fit2_mv))
                 
                  # decideTests(fit2_mv)
                  # # print(table(decideTests(fit2_mv)))
                 
                 
                  aw <- arrayWeights(data,design_mv)
                 
                  fit_mv <- lmFit(data, design_mv,
                                weights = aw)
                  contrasts_mv <- makeContrasts(PreV - PostV, levels=design_mv)
                  fit2_mv <- contrasts.fit(fit_mv, contrasts_mv)
                  fit2_mv <- eBayes(fit2_mv)
                  print(table(decideTests(fit2_mv)))
 
                  #Annotation
                 
                  t_mv= topTable(fit2_mv, number= filter_results)
                  top_mv=rownames(t_mv)
                 
                 
                  return (top_mv)}''')
        
            
            
            
        robjects.r('''dif_exp_CS2 <- function(annotation, metadata, data, y_pred){
               
                y_pred= metadata[,y_pred]
                design_mv <- model.matrix(~0+y_pred)
               

                colnames(design_mv) <- c("control", "drought")
                # print(design_mv)
                # print(head(data))

                 # fit_mv <- lmFit(data, design_mv) 
                 # # print(head(fit_mv$coefficients))

                 # contrasts_mv <- makeContrasts(control- drought, levels=design_mv)

                
                 aw <- arrayWeights(data,design_mv)
                
                 fit_mv <- lmFit(data, design_mv,
                               weights = aw)
                 contrasts_mv <- makeContrasts(control - drought, levels=design_mv)
                 fit2_mv <- contrasts.fit(fit_mv, contrasts_mv)
                 fit2_mv <- eBayes(fit2_mv)
                 print(table(decideTests(fit2_mv)))
                 
                 #Annotation
                
                 t_mv= topTable(fit2_mv, sort.by = "P")
                 top_mv=rownames(t_mv)
                
                 top10_mv=as.data.frame(annotation[top_mv,])
                 # print(top10_mc)
                
                 top1000= topTable(fit2_mv,sort.by = "P", number=1000)
                 top_mv_1000 =rownames(top1000)
                
                 top50= topTable(fit2_mv,sort.by = "P", number=50)
                 top_mv_50 =rownames(top50)
                
                 top212= topTable(fit2_mv,sort.by = "P", number=212)
                 top_mv_212=rownames(top212)
                
                 full_results_mv <- topTable(fit2_mv, number=Inf)
                 full_results_mv <- tibble::rownames_to_column(full_results_mv,"ID")
                 
                 print(head(full_results_mv))
                 
                 rownames(full_results_mv) <- full_results_mv$ID
                 topfull_mv=rownames(full_results_mv)
                
                print(head(topfull_mv))
                
                 full_results_mv["REF"]= annotation[topfull_mv,"REF"]
                 full_results_mv["Annotation"]=annotation[topfull_mv,"Annotation"]
                
                  print(head(full_results_mv,30))
                  write.table(full_results_mv, "full_results_mv_CS2.csv", sep=",", col.names=TRUE, quote=FALSE, row.names=TRUE,dec = ".")
              
                  # fit_ <- eBayes(fit_mv) 
                  
                  # full_results_Drought <- topTable(fit_,coef="Drought",sort.by="P", number=Inf)
                    
                  #     full_results_Drought <- tibble::rownames_to_column(full_results_Drought,"ID")

                  #     rownames(full_results_Drought) <- full_results_Drought$ID
                  #     topfull_mv_Drought=rownames(full_results_Drought)
                    
                  #     full_results_Drought["UniProtKB"]= annotation[topfull_mv_Drought,"UniProtKB"]
                  #     full_results_Drought["Annotation"]=annotation[topfull_mv_Drought,"Annotation"]
                    
                  #      print(head(full_results_Drought,30))
                  #     write.table(full_results_Drought, "full_results_Drought.csv", sep=",", col.names=TRUE, quote=FALSE, row.names=TRUE,dec = ".")
                    
              
                
              
                
                # install.packages("xlsx")  
                # library(xlsx)
                # write.xlsx(x = full_results_mv, file = "full_results_mv.xlsx", row.names=TRUE, col.names=TRUE)
              
                
                # fit_ <- eBayes(fit_mv) 
              
                # full_results_PreV <- topTable(fit_,coef="PreV",sort.by="P", number=Inf)
                
                #  full_results_PreV <- tibble::rownames_to_column(full_results_PreV,"ID")

                #  rownames(full_results_PreV) <- full_results_PreV$ID
                #  topfull_mv_PreV=rownames(full_results_PreV)
                
                #  full_results_PreV["UniProtKB"]= annotation[topfull_mv_PreV,"UniProtKB"]
                #  full_results_PreV["Annotation"]=annotation[topfull_mv_PreV,"Annotation"]
                
                #   # print(head(full_results_PreV,30))
                #   write.table(full_results_PreV, "full_results_PreV.csv", sep=",", col.names=TRUE, quote=FALSE, row.names=TRUE,dec = ".")
                
              
                # full_results_PostV <- topTable(fit_,coef="PostV",sort.by="P", number=Inf)
                
                #  full_results_PostV <- tibble::rownames_to_column(full_results_PostV,"ID")

                #  rownames(full_results_PostV) <- full_results_PostV$ID
                #  topfull_mv_PostV=rownames(full_results_PostV)
                
                #  full_results_PostV["UniProtKB"]= annotation[topfull_mv_PostV,"UniProtKB"]
                #  full_results_PostV["Annotation"]=annotation[topfull_mv_PostV,"Annotation"]
                
                #   # print(head(full_results_PostV,30))
                #   write.table(full_results_PostV, "full_results_PostV.csv", sep=",", col.names=TRUE, quote=FALSE, row.names=TRUE,dec = ".")
                
              
                
              
                
              
                
              
                
                 # Make sure you have ggplot2 loaded
                 ggplot(full_results_mv, aes(x = logFC, y=B)) + geom_point()
                
                
                 ## change according to your needs
                 p_cutoff <- 0.05
                 fc_cutoff <- 1
                
                 full_results_mv %>% 
                   mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
                   ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
                
                 p_cutoff <- 0.05
                 fc_cutoff <- 1
                 topN <- 20
                
                png(filename="full_results_mvCS2.png")

 
                print(full_results_mv %>% 
                 mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
                 mutate(Rank = 1:n(), Label = ifelse(Rank < topN, REF,"")) %>% 
                  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point()) #+ geom_text_repel(col="black"))
                
                dev.off()


                 topN <- 20
                 ##
                 ids_of_interest <- mutate(full_results_mv, Rank = 1:n()) %>% 
                   filter(Rank < topN) %>% 
                   pull("ID")
                
                 gene_names <- mutate(full_results_mv, Rank = 1:n()) %>% 
                   filter(Rank < topN) %>% 
                   pull(REF) 
                
                
                 ## Get the rows corresponding to ids_of_interest and all columns
                 gene_matrix <- data[ids_of_interest,]
                
                 pheatmap_de= pheatmap(gene_matrix,
                           labels_row = annotation[ids_of_interest,1],
                           scale="row")
               
               
               
               return (pheatmap_de)}
               ''')

    
            
            
            
            
            
        robjects.r('''top_CS2 <- function(annotation, metadata, data, y_pred, filter_results){
                   y_pred= metadata[,y_pred]
                   design_mv <- model.matrix(~0+y_pred)
                
 
                  colnames(design_mv) <- c("control", "drought")
                  # print(design_mv)
                  # print(head(data))
 
                  # fit_mv <- lmFit(data, design_mv) 
                  # # print(head(fit_mv$coefficients))
 
                  # contrasts_mv <- makeContrasts(control - drought, levels=design_mv)
 
 
                  # fit2_mv <- contrasts.fit(fit_mv, contrasts_mv)
                  # fit2_mv <- eBayes(fit2_mv)
                  # # print(topTable(fit2_mv))
                 
                  # decideTests(fit2_mv)
                  # # print(table(decideTests(fit2_mv)))
                 
                 
                  aw <- arrayWeights(data,design_mv)
                 
                  fit_mv <- lmFit(data, design_mv,
                                weights = aw)
                  contrasts_mv <- makeContrasts(control - drought, levels=design_mv)
                  fit2_mv <- contrasts.fit(fit_mv, contrasts_mv)
                  fit2_mv <- eBayes(fit2_mv)
                  print(table(decideTests(fit2_mv)))
 
                  #Annotation
                 
                  t_mv= topTable(fit2_mv, number= filter_results)
                  top_mv=rownames(t_mv)
                 
                 
                  return (top_mv)}''')    
                     
            
        top={}
        if default == True:           
            if which_data=="ALL":
            
               dif_exp= robjects.r["dif_exp"]
               for name, data in self.data_norm_final.items():
                   if name.lower() == "metadata":
                      self.metadata=data
                   else:
                       # if which_data== "ALL":
                          # for name, data in self.data_norm_final.items():
                          #     if name.lower() != "metadata":
                      print(name)
                      # print("Here!")
                      pheatmap_de = dif_exp(annotation, self.metadata, data, y_pred)
                      grdevices.png(file="pheatmap_de_"+ which_data.lower() + "_"+ name + ".png", width=512, height=512)
                      print(pheatmap_de)
                      grdevices.dev_off()
                      if filter_results is not None:
                         
                          top_mv= robjects.r["top_mv"]
                          top_mv=top_mv(annotation, self.metadata, data, y_pred, filter_results)
                          top[name]=top_mv
                      
                      if filter_results is None:
                          top[name]= None
                
                  
            if which_data != "ALL":
               dif_exp= robjects.r["dif_exp"]
               for name, data in self.data_norm_final.items():
                   if name.lower() == "metadata":
                      self.metadata=data
                   
                   top[name]=None
                   if name.lower() == which_data.lower():
                      pheatmap_de = dif_exp(annotation, self.metadata, data, y_pred)
                      grdevices.png(file="pheatmap_de_"+ which_data.lower() + "_"+ name + ".png", width=512, height=512)
                      print(pheatmap_de)
                      grdevices.dev_off()
                      if filter_results is not None:
                         
                          top_mv= robjects.r["top_mv"]
                          top_mv=top_mv(annotation, self.metadata, data, y_pred, filter_results)
                          top[name]=top_mv
                          
                      if filter_results is None:
                          top[name]= None
                        
        if default == False:
            ''' Aqui devo meter outras formas de anotação, como um package do R'''
            print("New_annotation")
        # print(utils.head(top_mv))
        
        
        
        if default =="CS2":
            # if which_data=="ALL":
            
            #    dif_exp= robjects.r["dif_exp_CS2"]
            #    for name, data in self.data_norm_final.items():
            #        if name.lower() == "metadata":
            #           self.metadata=data
            #        else:
            #            # if which_data== "ALL":
            #               # for name, data in self.data_norm_final.items():
            #               #     if name.lower() != "metadata":
            #           print(name)
            #           # print("Here!")
            #           # pheatmap_de = dif_exp(annotation, self.metadata, data, y_pred)
            #           # grdevices.png(file="pheatmap_de_"+ which_data.lower() + "_"+ name + ".png", width=512, height=512)
            #           # print(pheatmap_de)
            #           # grdevices.dev_off()
            #           if filter_results is not None:
                         
            #               top_cs2= robjects.r["top_CS2"]
            #               top_cs2=top_cs2(self.metadata, data, y_pred, filter_results)
            #               top[name]=top_cs2
                      
            #           if filter_results is None:
            #               top[name]= None
                
                  
            # if which_data != "ALL":
            #    dif_exp= robjects.r["dif_exp_CS2"]
            #    for name, data in self.data_norm_final.items():
            #         if name.lower() == "metadata":
            #            self.metadata=data
                    
            #         top[name]=None
            #         if name.lower() == which_data.lower():
            #            # pheatmap_de = dif_exp(annotation, self.metadata, data, y_pred)
            #            # grdevices.png(file="pheatmap_de_"+ which_data.lower() + "_"+ name + ".png", width=512, height=512)
            #            # print(pheatmap_de)
            #            # grdevices.dev_off()
            #            if filter_results is not None:
            #                # print("Here!")
            #                top_cs2= robjects.r["top_CS2"]
            #                top_cs2=top_cs2(self.metadata, data, y_pred, filter_results)
            #                top[name]=top_cs2
                           
            #            if filter_results is None:
            #                top[name]= None
                        
                     
                if which_data=="ALL":
                
                   dif_exp= robjects.r["dif_exp_CS2"]
                   for name, data in self.data_norm_final.items():
                       if name.lower() == "metadata":
                          self.metadata=data
                       else:
                           # if which_data== "ALL":
                              # for name, data in self.data_norm_final.items():
                              #     if name.lower() != "metadata":
                          print(name)
                          # print("Here!")
                          pheatmap_de = dif_exp(annotation, self.metadata, data, y_pred)
                          grdevices.png(file="pheatmap_de_"+ which_data.lower() + "_"+ name + "CS2.png", width=512, height=512)
                          print(pheatmap_de)
                          grdevices.dev_off()
                          if filter_results is not None:
                             
                              top_mv= robjects.r["top_CS2"]
                              top_mv=top_mv(annotation, self.metadata, data, y_pred, filter_results)
                              top[name]=top_mv
                          
                          if filter_results is None:
                              top[name]= None
                    
                      
                if which_data != "ALL":
                   dif_exp= robjects.r["dif_exp_CS2"]
                   for name, data in self.data_norm_final.items():
                       if name.lower() == "metadata":
                          self.metadata=data
                       
                       top[name]=None
                       if name.lower() == which_data.lower():
                          pheatmap_de = dif_exp(annotation, self.metadata, data, y_pred)
                          grdevices.png(file="pheatmap_de_"+ which_data.lower() + "_"+ name + "CS2.png", width=512, height=512)
                          print(pheatmap_de)
                          grdevices.dev_off()
                          if filter_results is not None:
                             
                              top_mv= robjects.r["top_CS2"]
                              top_mv=top_mv(annotation, self.metadata, data, y_pred, filter_results)
                              top[name]=top_mv
                              
                          if filter_results is None:
                              top[name]= None
        
        return top
            
                        
                      
          
                        
                       
                                   
                                
                                
                                
                            
                    

                        
            