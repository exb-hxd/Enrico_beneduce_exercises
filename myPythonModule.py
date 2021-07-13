# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

import numpy as np
import re
from itertools import islice

class measurements:  #class for storing data from file with heading for every column on first line. nquantities dictionaries 
    def __init__(self, nquantities: int, nfields, headed_file_name, encoding = "utf8", consumeString = None, AddHeader = None):   ## with nfields entries each are created 
        self.nquantities = nquantities
        self.nfields = nfields
        
        self.data = []
        with open(headed_file_name, encoding=encoding ) as infile:
            
            print ("reading " + headed_file_name)
            
            
            if isinstance (consumeString, str) :
                for line in infile:
                    if re.search (consumeString, line ):
                        break
 
            
            #making out heading at the beginning of the readable part
            headers = ["\n"]
            titles = []
            
            
            
            while headers == ["\n"]:
                headers = infile.readline().split()
                print ("found headers: ", headers)
            numbers = []
            for h in headers:
                try:
                    numbers.append (float (h))
                except ValueError:  #this should happen with normal usage (headed file)
                    titles = headers
                    print (titles, " used as dictionary keys")
                    numbers = []
                    break
                    
            if len(numbers) == len(headers):
                if AddHeader == None:
                    titles = [("field" + str(ind)) for ind, n in enumerate (numbers)]
                    print ("no keys were found or specified: assigning automatic names")
                elif len (AddHeader) == len (numbers):
                    titles = AddHeader
                    print ("assigning key values ", AddHeader, " passed as argument")
                else:
                   raise ValueError( "AddHeader should have " + str (len (numbers)) + " elements, not " + len(AddHeader))

            
            
            
#           titles =  infile.readline()
            titles = np.array (titles)
            
            tempdata = np.loadtxt (infile, unpack=True)
            if len (numbers) == len (headers ): #also len(numbers) > 0 should be a sufficient contition here
                tempdata = np.insert (tempdata, 0, numbers, axis=1)

            if  isinstance (nfields, int):
                list_titles = titles.reshape ((self.nquantities, self.nfields))
                list_data = tempdata.reshape ((self.nquantities, self.nfields, -1))
            elif isinstance (nfields, tuple) :
                if len (self.nfields) != self.nquantities:
                    print ("warning: processing as ", len (self.nfields), "branches, not ", self.nquantities)
                donecols = 0
                list_titles = []
                list_data = []
                for n in nfields:
                    list_titles.append(titles [donecols : donecols + n])
                    list_data.append (tempdata [donecols : donecols + n])
                    donecols += n

            for names, branchData in zip(list_titles, list_data):
                tempdict = {}
                for field, leafData in zip (names, branchData):
                    tempdict[field] = leafData
                self.data.append (tempdict)
                
                
    def print_names (self):
        for ind, branch in enumerate (self.data):
            print ("branch " + str(ind) + ": ", ', '.join (branch), '\n')
                
            


def autocorrelationFromFile (filename, usecols:int, tlimit:int = 1000, readsize:int = 10000): 
    
    assert tlimit <= readsize
    
    Num1 = np.zeros(tlimit)
    Sum = 0.
    SumSq = 0.
    SaveFirst = np.zeros(tlimit)
    SaveLast = np.zeros(tlimit)
    
    ndata = 0
    heritage = np.zeros(tlimit)
    
    with open(filename) as file:
        while True:
            
            try:
                readvalues = np.loadtxt (islice (file, readsize), 
                                             usecols=usecols, unpack = True)
            except:
                break
            
            oldvalues = np.concatenate ((heritage, readvalues), axis=0)

            if readvalues.shape[0] == 0:
                break

            if ndata == 0:
                for t in range (1,tlimit):
                    Num1[t-1] += np.sum(readvalues[:-t] * readvalues[t:])
                    SaveFirst = readvalues[:tlimit]
            else:
                for t in range (1, tlimit):
                    Num1[t-1] += np.sum(oldvalues[-readvalues.shape[0]-t : -t] * readvalues, axis = 0)
                   
            Sum += np.sum(readvalues)
            SumSq += np.sum(readvalues**2)
            
            SaveLast = oldvalues[-tlimit:]
            
            ndata += readvalues.shape[0]
            heritage = readvalues [-tlimit:]
            
            print("processing " + str(ndata) + " values")
            
    correlations = np.empty(tlimit, dtype=float)
    for t in range(1, tlimit):
        correlations[t-1] = (Num1[t-1] - (Sum -  np.sum(SaveLast[-t:])) * (Sum - np.sum(SaveFirst[:t])) / float(ndata - t) )/ (SumSq - Sum**2 / ndata) * float(ndata)/float(ndata-t)

    return correlations
                         

def autocorrelation (array, tlimit):
    ndata = array.shape[0]
    Num1 = np.zeros(tlimit)
    Sum = np.sum (array)
    SumSq = np.sum (array**2)
    correlations = np.zeros(tlimit)

    
    for t in range (1,tlimit):
        Num1[t-1] += np.sum(array[:-t] * array[t:])
        correlations[t-1] = (Num1[t-1] - (Sum -  np.sum(array[-t:])) * (Sum - np.sum(array[:t])) / float(ndata - t) )/ (SumSq - Sum**2 / ndata) * float(ndata)/float(ndata-t)

    return correlations
            

def block_sigmas (array, sizes):
    
    sigmas = []

    for size in sizes :
        if (size == 0):
            continue
        if (size == 1):
            sigmas.append (0)
            continue
            
        values = 0
        if array.shape[0] % size == 0:
            values = np.average (np.reshape (array, (-1, size)), axis=1)
        else:
            values = np.average (np.reshape (array [0: - (array.shape[0] % size)], (-1, size)), axis=1)

        mean = np.average (values)
        sqMean = np.average (values ** 2)
        
        
        sigmas.append (np.sqrt ((sqMean - mean**2) / (array.shape [0] / (size-1) ) ) )
        
    return sigmas

                

class Converter:
    
    def __init__(self, ep, sig, m):
        self.kb = 1.38064852e-23
        self.epsilon = ep
        self.sigma = sig
        self.m = m
    
    def convert (self, x, quantity, conversion: float = 1):
        if quantity == "temp":
            return x * self.epsilon / self.kb * conversion
        if quantity == "E":
            return x * self.epsilon * conversion
        if quantity == "press":
            return x / self.sigma**3 * self.epsilon * conversion
        if quantity == "rho":
            return x * self.m / self.sigma**3 * conversion
        if quantity == "length":
            return x * self.sigma * conversion
            
      
            

            
