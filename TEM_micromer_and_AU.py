# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 16:05:21 2021

@author: Usuario
"""

import numpy as np


import matplotlib.pyplot as plt

import pandas as pd
import os

from pylab import *
from scipy.optimize import curve_fit

import seaborn as sns

sns.set_theme()
sns.set_style("ticks")
sns.set_context("paper",font_scale=1.5)
sns.set_palette("tab10")


def gauss(x,mu,sigma,A):
    """
    Gaussian distribution function

    Parameters
    ----------
    x : Independent variable.
    mu : Mean of the distribution.
    sigma : Standard deviation of the distribution.
    A : Peak's Area.

    Returns
    -------
    Numpy array.

    """
    
    return A*exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    """
    Bimodal distribution function. Calculated as the sum of two gaussian distribution

    Parameters
    ----------
    x : Independent variable.
    mu1 : Mean of the first gaussian distribution .
    sigma1 : Standard deviation of the first gaussian distribution.
    A1 : Peak's Area of the first gaussian distribution.
    mu2 : Mean of the second gaussian distribution.
    sigma2 : Standard deviation of the second gaussian distribution.
    A2 : Peak's Area of the second gaussian distribution.
    
    Returns
    -------
    class 'numpy.ndarray'

    """
    
    
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)





def micromer_water():
    """
    Run plotting and fitting of the 25 nm Micromer-redF-Plain NPs dispersed in water.

    Returns
    -------
    None.

    """
    
    os.chdir(r"XXXX") #Choose directory containing the txt files with the measurements 
    
    
    size = [] #initialize size list
    
    #Read files and extract data
    
    for file in os.listdir("./"):
        if file[-4:]==".txt" and "Result" not in file:
            data = pd.read_table(file)
            print(data.head())
            for i in range(len(data.loc[:,"Length"])):
                size.append(data.loc[i,"Length"] )
    
    
    #Initialize plot
    
    fig, ax = plt.subplots()
    
    y,x,_=hist(size,20,alpha=1,label='data',edgecolor="k") #Plot data as histogram and extract the bins positions and count or fitting
    
    x=(x[1:]+x[:-1])/2 #Calculate bins' center point
    
    expected=(28,4,10) #Initial fitting parameters
    
    params,cov=curve_fit(gauss,x,y,expected) #Fit to gaussian distribution
    sigma=sqrt(diag(cov)) #Calculate errors of the fitting paramenters
    
    printing = pd.DataFrame(data={'params':params,'sigma':sigma},index=gauss.__code__.co_varnames[1:]) #Create DataFrame with the fitting results
        
    y1 = gauss(x,*params) #Obtain fitted gaussian distribution
    
    plt.plot(x,y1,"-",color = "red",linewidth=2) #Plot fitted gaussian distribution
    plt.xlabel("Diameter (nm)")
    plt.ylabel("Counts")
    
    plt.savefig(r"Saving_directory\Micromer_water_hist.tif",bbox_inches="tight",dpi=100) #Save image in the desired Saving directory
    
    print(printing) #Printfitting parameters

    return None


def micromer_AU():
    """
    Run plotting and fitting of the 25 nm Micromer-redF-Plain NPs dispersed in DMEM with measured protein corona.

    Returns
    -------
    None.

    """
    
    os.chdir(r"XXXX") #Choose directory containing the txt files with the measurements 
    
    
    size = [] #Initialize size list
    
    #Extract data and append to size list
    for file in os.listdir("./"):
        if file[-4:]==".txt" and "Result" not in file:
            data = pd.read_table(file)
            print(data.head())
            for i in range(len(data.loc[:,"Length"])):
                size.append(data.loc[i,"Length"] ) 
    
    
    #PInitialize plot
    fig, ax = plt.subplots()
    
    y,x,_=hist(size,20,alpha=1,label='data',edgecolor="k") #Plot data as histogram and extract the bins positions and count or fitting
    x=(x[1:]+x[:-1])/2   #Calculate bins' center point
    expected=(28,4,10) #Initial fitting parameters
    
    params,cov=curve_fit(gauss,x,y,expected)  #Fit to gaussian distribution
    sigma=sqrt(diag(cov)) #Calculate errors of the fitting paramenters
    
    printing = pd.DataFrame(data={'params':params,'sigma':sigma},index=gauss.__code__.co_varnames[1:]) #Create DataFrame with the fitting results
           
    
    y1 = gauss(x,*params) #Obtain fitted gaussian distribution
    
    plt.plot(x,y1,"-",color = "red",linewidth=2) #Plot fitted gaussian distribution
    
    plt.xlabel("Diameter (nm)")
    plt.ylabel("Counts")
    
    plt.savefig(r"Saving_directory\Micromer_AU_hist.tif",bbox_inches="tight",dpi=100) #Save image in the desired Saving directory
    print(printing)

def micromer_AU_core():
    
    
    """
    Run plotting and fitting of the 25 nm Micromer-redF-Plain NPs dispersed in DMEM with measured protein corona but measuring only the core.
   
    Returns
    -------
    None.
   
    """
    
    os.chdir(r"XXXX") #Choose directory containing the txt files with the measurements 
    
    
    size = [] #initialize size list
    
    #Read files and extract data
    for file in os.listdir("./"):
        if file[-4:]==".txt" and "Result" not in file:
            data = pd.read_table(file)
            print(data.head())
            for i in range(len(data.loc[:,"Length"])):
                size.append(data.loc[i,"Length"] )
    
    
    #Initialize plots
    fig, ax = plt.subplots()
    
    y,x,_=hist(size,20,alpha=1,label='data',edgecolor="k")
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
    expected=(28,4,10)
    params,cov=curve_fit(gauss,x,y,expected)
    sigma=sqrt(diag(cov))
    printing = pd.DataFrame(data={'params':params,'sigma':sigma},index=gauss.__code__.co_varnames[1:])
    mu1 = printing.loc["mu",:][0]
    sigma1 = printing.loc["sigma",:][0]
    A1 = printing.loc["A",:][0]
    y1 = gauss(x,mu1,sigma1,A1)
    
    plt.plot(x,y1,"-",color = "red",linewidth=2)
    plt.xlabel("Diameter (nm)")
    plt.ylabel("Counts")
    
    plt.savefig(r"C:\Users\Alberto\OneDrive - UAM\Tesis\Chapter 1\Figures\NPs\Micromer_AU_core_hist.tif",bbox_inches="tight",dpi=100)
    plt.show()
    plt.close()
    print(printing)


def prochimia_AU():
    
    """
    Run plotting and fitting of the 25 nm 2.5-Cy5-AuNPs dispersed in DMEM with measured protein corona.
   
    Returns
    -------
    None.
   
    """
    
    os.chdir(r"XXXX") #Choose directory containing the txt files with the measurements 
    
    size = []
    
    for file in os.listdir("./"):
        if file[-4:]==".txt" and "Result" not in file:
            data = pd.read_table(file)
            print(data.head())
            for i in range(len(data.loc[:,"Length"])):
                size.append(data.loc[i,"Length"] )
    
    fig, ax = plt.subplots()
    
    y,x,_=hist(size,30,alpha=1,label='data',edgecolor="k")
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
    expected=(30,7,10,50,10,10)
    params,cov=curve_fit(bimodal,x,y,expected)
    sigma=sqrt(diag(cov))
    
    printing = pd.DataFrame(data={'params':params,'sigma':sigma},index=bimodal.__code__.co_varnames[1:])
    mu1 = printing.loc["mu1",:][0]
    mu2 = printing.loc["mu2",:][0]
    sigma1 = printing.loc["sigma1",:][0]
    sigma2 = printing.loc["sigma2",:][0]
    A1 = printing.loc["A1",:][0]
    A2 = printing.loc["A2",:][0]
    y1 = gauss(x,mu1,sigma1,A1)
    plt.plot(x,y1,"-",color = "lime",linewidth=1)
    y2 = gauss(x,mu2,sigma2,A2)
    plt.plot(x,y2,"-",color = "lime",linewidth=1)
    print(printing)
    plot(x,bimodal(x,*params),color='red',lw=1.5,label='fit')
    
    
    plt.xlabel("Diameter (nm)")
    plt.ylabel("Counts")
    
    plt.savefig(r"C:\Users\Alberto\OneDrive - UAM\Tesis\Chapter 1\Figures\NPs\Prochimia_AU_hist.tif",bbox_inches="tight",dpi=100)
    print("Number of NPs measured: ", len(size))

def prochimia_AU_core():
    """
    Run plotting and fitting of the 25 nm 2.5-Cy5-AuNPs dispersed in DMEM with measured protein corona but measuring only the core.
   
    Returns
    -------
    None.
   
    """
    
    os.chdir(r"XXXX") #Choose directory containing the txt files with the measurements 
    
    size = []
    
    for file in os.listdir("./"):
        if file[-4:]==".txt" and "Result" not in file:
            data = pd.read_table(file)
            print(data.head())
            for i in range(len(data.loc[:,"Length"])):
                size.append(data.loc[i,"Length"] )
    
    fig, ax = plt.subplots()
    
    y,x,_=hist(size,30,alpha=1,label='data',edgecolor="k")
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
    expected=(12,7,10,25,10,10)
    params,cov=curve_fit(bimodal,x,y,expected)
    sigma=sqrt(diag(cov))
    
    printing = pd.DataFrame(data={'params':params,'sigma':sigma},index=bimodal.__code__.co_varnames[1:])
    mu1 = printing.loc["mu1",:][0]
    mu2 = printing.loc["mu2",:][0]
    sigma1 = printing.loc["sigma1",:][0]
    sigma2 = printing.loc["sigma2",:][0]
    A1 = printing.loc["A1",:][0]
    A2 = printing.loc["A2",:][0]
    y1 = gauss(x,mu1,sigma1,A1)
    plt.plot(x,y1,"-",color = "lime",linewidth=1)
    y2 = gauss(x,mu2,sigma2,A2)
    plt.plot(x,y2,"-",color = "lime",linewidth=1)
    print(printing)
    plot(x,bimodal(x,*params),color='red',lw=1.5,label='fit')
    
    
    plt.xlabel("Diameter (nm)")
    plt.ylabel("Counts")
    
    plt.savefig(r"C:\Users\Alberto\OneDrive - UAM\Tesis\Chapter 1\Figures\NPs\Prochimia_AU_core_hist.tif",bbox_inches="tight",dpi=100)
    print("Number of NPs measured: ", len(size))

prochimia_AU()
micromer_water()
micromer_AU()
prochimia_AU_core()