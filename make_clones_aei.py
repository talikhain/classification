# -*- coding: utf-8 -*-
"""
generate clones from .aei file
"""
#import needed packages
from __future__ import division
import numpy as np
import random as rd
import numpy.random as rd
import orbfit_tools_bary as bary
import orbfit_tools_helio as helio
import pandas as pd

#function to make clones
def make_clones(folder, tno, N):
    
    #start with barycentric orbfit
    abginfo = bary.orbfit_abg(folder + tno + '/' + tno + '.abg')
    jd0 = abginfo['jd0']
    elems_aei = abginfo['elems_aei']
    mu = bary.mean_anomaly(elems_aei, jd0)
    
    #read in from aei file
    aei_info = open(folder + tno + '/' + tno + '.aei').read().split()
    epoch = float(aei_info[8][:-1])
    a = float(aei_info[20])
    e = float(aei_info[21])
    inc = float(aei_info[22])
    lan = float(aei_info[23]) 
    aop = float(aei_info[24])
    top = float(aei_info[25])
    best_aei = [a, e, inc, lan, aop, top] #best fit
    
    #covariance matrix
    cov_aei = np.zeros([6,6])
    for i in xrange(6):
        for j in xrange(6):
            cov_aei[i,j] = float(aei_info[36 + 6*i + j])
    
    i = 0
    clones_aei = np.zeros([N+1, 6])
    clones = rd.multivariate_normal(best_aei, cov_aei, N)
    clones_aei[:N] = clones
    clones_aei[N] = best_aei
    
    #compute mean anomaly
    mu_list = []
    abginfo = helio.orbfit_abg(folder + tno + '/' + tno + '.abg')
    jd0 = abginfo['jd0']
    for i in xrange(N+1):
        clones_dict = {'a': clones_aei[i][0], 'e': clones_aei[i][1], 'i': clones_aei[i][2], 'lan': clones_aei[i][3], 'aop': clones_aei[i][4], 'top': clones_aei[i][5]}
        mu = helio.mean_anomaly(clones_dict, jd0)
        mu_list.append(mu)
    
    #function to make small.in file
    def prep_small(elems, mu_list, epoch, i):
        block = str(name)+' ep='+str(epoch)+'\n'+\
                ' '+str(elems[i][0])+' '+str(elems[i][1])+' '+str(elems[i][2])+' '+str(elems[i][4])+' '+str(elems[i][3])+' '+str(mu_list[i])+' 0 0 0 '+'\n'
        return block 
    
    #make small.in file
    header = ')O+_06 Small-body initial data  (WARNING: Do not delete this line!!)'+'\n'+\
    ') Lines beginning with ) are ignored.'+'\n'+\
    ')---------------------------------------------------------------------'+'\n'+\
    'style (Cartesian, Asteroidal, Cometary) = Ast'+'\n'+\
    ')---------------------------------------------------------------------'+'\n'

    #small.in file with best fit and N clones
    small_in = open(folder + tno + '/' + 'small.in','w+')
    small_in.write(header)
    for i in xrange(N):
        name = 't_'+str(i)
        clone = prep_small(clones_aei, mu_list, jd0, i)
        small_in.write(clone)
    
    name = 'best'
    best = prep_small(clones_aei, mu_list, jd0, N)
    small_in.write(best)
    
    small_in.close()


folder = '/Users/talikhain/Desktop/ClassificationTest/'
names = ['s12_good_3'] #pd.read_table('names.csv')
num = len(names)
N = 10 #number of clones

#make clones
for i in range(num):
    tno = names[i]#names.values[i][0]
    print(i, tno)
    make_clones(folder, tno, N)
