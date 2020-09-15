#import needed packages
from __future__ import division
import numpy as np
import random as rd
import pandas as pd
import numpy.random as rd

from classify import *
from find_resonance import *

folder = '/Users/talikhain/Desktop/ClassificationTest/'
N = 10 #number of clones

#a pandas dataframe with the barycentric aei elements of the best fit (use this to classify into some of the categories)
names_bary = pd.read_pickle('/Users/talikhain/Desktop/ClassificationTest/tno_bary.pkl')
num = len(names_bary)

for ind in range(num):
    tno = names_bary['tno'].values[ind]
    a = names_bary['ab'].values[ind]
    e = names_bary['eb'].values[ind]
    i = names_bary['incb'].values[ind]
    node = names_bary['lanb'].values[ind]
    peri = names_bary['aopb'].values[ind]
    M = names_bary['Mb'].values[ind]
    best_bary = [a, e, i, node, peri, M]
    print('tno #', ind, tno)
    
    classify(folder, tno, N, best_bary)
    find_resonance(folder, tno, N)
    