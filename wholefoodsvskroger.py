# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 00:45:40 2015

MTS Program that graphs difference in mahalanobis distances, as well as an MTS



@author: Naqib
"""
from __future__ import print_function
import numpy as np
import matplotlib as plt
import copy
import math




def standardize(healthydata,data):
    '''standardizes data, based off mean/std of 'healthydata'
    '''
    stand = (data - healthydata.mean(axis=0)) / healthydata.std(axis=0, ddof=1)
    return stand

def dividebyKvalue(mat):
    '''divides an array by the number of variables
    '''
    k = np.shape(mat)[1]
    x = mat/k
    return x

def getCor(healthydata):
    '''returns correlation coefficient of an array
    '''
    cor = np.corrcoef(np.matrix.transpose(healthydata))
    cor =np.asmatrix(cor)
    return cor
    
def getMD(healthydata,data2):
    '''
    calculates mahalanobis distance of paramater data2
    based on healthydata
    '''
    cor = getCor(healthydata) #gets correlation matrix
    data = standardize(healthydata,data2)#standardizes matrix
    inverse = np.linalg.inv(cor) #takes the inverse of the correlation matrix
    step4 = dividebyKvalue(data) #divides by k
    step5 = step4*inverse #explained based off excel sheet
    MD =  step5 *np.matrix.transpose(data) #step 6
    return np.diagonal(MD) #returns the diaganol of step 6 which is the MDs
    
def getOrth():
    '''
    Ask user for which orthogonal array they want to use
    '''
    name = raw_input("Which orthogonal array would you like to use?")
    ortho = np.genfromtxt(name, delimiter = ",")
    return ortho
    
def runMTS(healthydata, sickdata, ortho):
    '''
    runs a full mts
    plots on and off bar graph
    '''
    mahalanobisD= [] #empty list for all MD's
    for line in ortho: #iterates through ortho to determine on off
        rowstoremove = np.where(line==2)[0] #gets indices of vars to switch off
        s = copy.deepcopy(healthydata) #copies data to preserve original
        h = copy.deepcopy(sickdata) # same as before
        for value in reversed(rowstoremove): #reverse to not change indices
            s = np.delete(s, value, 1)#deletes off variables
            h = np.delete(h, value, 1)#same as before
        MD = getMD(s,h) #calculated MD based off run
        mahalanobisD.append(list(MD)) # MDs stored 
    TaguchiArray = np.reciprocal(np.square(mahalanobisD))# 1/md^2 
    StoN=[] #empty list for signal to noise
    for run in TaguchiArray: #iterates through each run
        sn = (-10)*math.log((0.1*sum(run)),10)
        StoN.append(sn)
    AverageOn = []#empty list for on 
    for row in np.transpose(ortho): # tranpos to determine when var was on
        x= 0#resets x
        rowstoaverage = np.where(row==1)[0] #determines where var was on
        for value in rowstoaverage: #sums values
            x+= StoN[value]
        AverageOn.append(x/len(rowstoaverage))#calculates average
    AverageOff = []#same process as above but for off values
    for row in np.transpose(ortho):
        x= 0
        rowstoaverage = np.where(row==2)[0]
        for value in rowstoaverage:
            x+= StoN[value]
        AverageOff.append(x/len(rowstoaverage))
#    fig, ax = plt.pyplot.subplots(figsize=(100,50))# creates a new plot
#    ind = range(len(AverageOn)) # x values are just the range of variables
#    ind2 = np.array(ind)+.35 #adds .35 to the range for to seperate bars 
#    On = ax.bar(ind,AverageOn,.35, color = 'r')#graphs on bars
#    Off = ax.bar(ind2,AverageOff, .35,color='c') # graphs off bars
#    ax.legend( (On[0], Off[0]), ('On', 'Off') )# creates leged
#    ax.set_xticks(ind2) # labels graph
#    ax.set_xticklabels(np.array(ind)+1) # aligns x labels
#    fig.set_dpi(400)
#    fig.savefig('taguchi')
#    plt.pyplot.show() #shows plot
#    
    ######
    fig, b = plt.pyplot.subplots(figsize=(30,15))# creates a new plot
    fig.suptitle('Change in Signal to Noise', fontsize=28)
    ind = range(len(AverageOn)) # x values are just the range of variables
    On = b.bar(ind,np.array(AverageOn)-np.array(AverageOff),.35, color = 'b')#graphs on bars
   # b.set_xticks(ind) # labels graph
    #b.set_xticklabels(np.array(ind)+1) # aligns x labels
       
    fig.set_dpi(400)
    fig.savefig('difference')
    plt.pyplot.show() #shows plot
    #####

        
def runMD(healthydata, sickdata):
    ''' 
    plots mahalanobis distances
    '''
    mdHealthy = getMD(healthydata,healthydata)
    mdSick= getMD(healthydata,sickdata)
    fig, lin = plt.pyplot.subplots(figsize=(15,10)) # creates new graph
    fig.suptitle('Mahalanobis Distances', fontsize=28)
    refgroup = lin.plot(mdHealthy,linewidth=2, color = 'g') 
    hgroup = lin.plot(mdSick,linewidth=2, color = 'b' )
    lin.legend( (refgroup[0],hgroup[0]),('reference group', 'outside group'))
    fig.set_dpi(400)
    fig.savefig('mahala')
    plt.pyplot.show()
    
    
def fillIN(aray, fill):
    ''' fills in the fill paramter with average values'''
    aray = np.transpose(aray)
    average = []
    for line in aray:
        y = list(line).count(fill)
    
        average.append((sum(line)-(fill*y))/ float((len (line)-y) ))
    x= 0
    avarray= []
    for line in aray:
        avarray.append([average[x] if n==fill else n for n in line])
        x+=1
        
        
    return np.transpose(np.array(avarray))
        

def main():
    
    #gets healthy data
    healthydata = np.genfromtxt ('whole_plums.csv', delimiter=",")
    #gets sick data
    sickdata = np.genfromtxt('kroger.csv', delimiter= ",")
    #gets ortho
    ortho = np.genfromtxt('ortho198.csv', delimiter = ",")
    #slicing out the first label rows
    runMD(fillIN(healthydata[1:,1:],2.5), fillIN(sickdata[1:,1:],2.5))
    runMTS(fillIN(healthydata[1:,1:],2.5), fillIN(sickdata[1:,1:],2.5), ortho)
    
#executes main 
if __name__ == '__main__':
    main()




