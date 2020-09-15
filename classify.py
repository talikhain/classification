#import needed packages
from __future__ import division
import numpy as np
import pandas as pd
import time as tm
import matplotlib as mpl
mpl.use('TkAgg') #use when running on a mac
import matplotlib.pyplot as plt
import pickle
from matplotlib import gridspec

##########################################
from matplotlib import rc
font = {'family': 'serif',
        'serif': ['Computer Modern'],
        'weight' : 'bold',
        'size'   : 12}

rc('font', **font)
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{amsmath}'
                          r'\boldmath')
##########################################


def classify(folder, tno, N, best_bary):
    print("...classifying...")
    headings = ['time', 'M', 'a', 'e', 'i', 'peri', 'node']
    columns = [0, 2, 3, 4, 5, 6, 7]
    Nep = pd.read_csv(folder + tno + '/' + 'NEPTUNE.aei', skiprows=4, delim_whitespace=True, usecols=columns, names=headings)
    data = pd.DataFrame()
    for i in range(N+1):
        if i == N:
            tmp = pd.read_csv(folder + tno + '/' + 'best.aei', skiprows=4, delim_whitespace=True, usecols=columns, names=headings)
            data = data.append({'time':tmp['time'].values, 'a':tmp['a'].values, 'e':tmp['e'].values, 'i':tmp['i'].values, 'node':tmp['node'].values, 'peri':tmp['peri'].values, 'M':tmp['M'].values}, ignore_index = True)
        else:
            tmp = pd.read_csv(folder + tno + '/' + 't_%s.aei' % (i), skiprows=4, delim_whitespace=True, usecols=columns, names=headings)
            data = data.append({'time':tmp['time'].values, 'a':tmp['a'].values, 'e':tmp['e'].values, 'i':tmp['i'].values, 'node':tmp['node'].values, 'peri':tmp['peri'].values, 'M':tmp['M'].values}, ignore_index = True)

    print("...plotting...")
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 8))

    plt.subplot(221)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['a'].values[i], alpha=0.75)
    plt.ylabel(r'\textbf{$a$ (AU)}', fontsize=16)

    plt.subplot(222)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['e'].values[i], alpha=0.75)
    plt.ylabel(r'\textbf{$e$}', fontsize=16)

    plt.subplot(223)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['i'].values[i], alpha=0.75)
    plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
    plt.ylabel(r'\textbf{$i$ (deg)}', fontsize=16)

    plt.subplot(224)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['a'].values[i]*(1 - data['e'].values[i]), alpha=0.75)
    plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
    plt.ylabel(r'\textbf{$q$ (AU)}', fontsize=16)
    fig.tight_layout()
    plt.savefig(folder + tno + '/' + 'aeiq.png', dpi=300)
    #plt.savefig(folder + tno + '/' + 'aeiq.pdf', dpi=300)
    
    ej = 0
    for i in range(N):
        if len(data['time'].values[i]) < len(Nep):
            ej = +1
    if ej > 0:
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 8))

        plt.subplot(221)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['a'].values[i], alpha=0.75)
        plt.ylabel(r'\textbf{$a$ (AU)}', fontsize=16)

        plt.subplot(222)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['e'].values[i], alpha=0.75)
        plt.ylabel(r'\textbf{$e$}', fontsize=16)

        plt.subplot(223)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['i'].values[i], alpha=0.75)
        plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
        plt.ylabel(r'\textbf{$i$ (deg)}', fontsize=16)

        plt.subplot(224)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['a'].values[i]*(1 - data['e'].values[i]), alpha=0.75)
        plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
        plt.ylabel(r'\textbf{$q$ (AU)}', fontsize=16)
        fig.tight_layout()
        plt.savefig(folder + tno + '/' + 'aeiq_wo_ej.png', dpi=300)
        
    print("...plotting with best fit...")
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 8))

    plt.subplot(221)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['a'].values[i], alpha=0.75, color='lightgrey')
    i = N
    plt.plot(data['time'].values[i]*1e-6, data['a'].values[i], alpha=0.75, color='blue')
    plt.ylabel(r'\textbf{$a$ (AU)}', fontsize=16)
    
    plt.subplot(222)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['e'].values[i], alpha=0.75, color='lightgrey')
    i = N
    plt.plot(data['time'].values[i]*1e-6, data['e'].values[i], alpha=0.75, color='blue')
    plt.ylabel(r'\textbf{$e$}', fontsize=16)

    plt.subplot(223)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['i'].values[i], alpha=0.75, color='lightgrey')
    i = N
    plt.plot(data['time'].values[i]*1e-6, data['i'].values[i], alpha=0.75, color='blue')
    plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
    plt.ylabel(r'\textbf{$i$ (deg)}', fontsize=16)

    plt.subplot(224)
    for i in range(N):
        plt.plot(data['time'].values[i]*1e-6, data['a'].values[i]*(1 - data['e'].values[i]), alpha=0.75, color='lightgrey')
    i = N
    plt.plot(data['time'].values[i]*1e-6, data['a'].values[i]*(1 - data['e'].values[i]), alpha=0.75, color='blue')
    plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
    plt.ylabel(r'\textbf{$q$ (AU)}', fontsize=16)
    fig.tight_layout()
    plt.savefig(folder + tno + '/' + 'aeiq_best.png', dpi=300)
    #plt.savefig(folder + tno + '/' + 'aeiq.pdf', dpi=300)
    
    
    ej = 0
    for i in range(N):
        if len(data['time'].values[i]) < len(Nep):
            ej = +1
    if ej > 0:
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 8))

        plt.subplot(221)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['a'].values[i], alpha=0.75, color='lightgrey')
        i = N
        if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['a'].values[i], alpha=0.75, color='blue')
        plt.ylabel(r'\textbf{$a$ (AU)}', fontsize=16)

        plt.subplot(222)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['e'].values[i], alpha=0.75, color='lightgrey')
        i = N
        if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['e'].values[i], alpha=0.75, color='blue')
        plt.ylabel(r'\textbf{$e$}', fontsize=16)

        plt.subplot(223)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['i'].values[i], alpha=0.75, color='lightgrey')
        i = N
        if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['i'].values[i], alpha=0.75, color='blue')
        plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
        plt.ylabel(r'\textbf{$i$ (deg)}', fontsize=16)

        plt.subplot(224)
        for i in range(N):
            if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['a'].values[i]*(1 - data['e'].values[i]), alpha=0.75, color='lightgrey')
        i = N
        if len(data['time'].values[i]) > len(Nep) - 3:
                plt.plot(data['time'].values[i]*1e-6, data['a'].values[i]*(1 - data['e'].values[i]), alpha=0.75, color='blue')
        plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)
        plt.ylabel(r'\textbf{$q$ (AU)}', fontsize=16)
        fig.tight_layout()
        plt.savefig(folder + tno + '/' + 'aeiq_best_wo_ej.png', dpi=300)
    
    a = best_bary[0]
    e = best_bary[1]
    i = best_bary[2]
    node = best_bary[3]
    peri = best_bary[4]
    M = best_bary[5]
    
    elements = open(folder + tno + '/' + 'elements.txt','w+')
    elements.write(str(a) + '\n')
    elements.write(str(e) + '\n')
    elements.write(str(i) + '\n')
    elements.write(str(peri) + '\n')
    elements.write(str(node) + '\n')
    elements.write(str(M))
    elements.close()
    
    flag = 0
    classification = open(folder + tno + '/' + 'classification.txt','w+')

    def is_comet(a, e, i):
        print("checking... comet")
        aJ = 5.204 #AU
        TJ = aJ/a + 2*np.sqrt(a/aJ * (1 - e**2))*np.cos(i) #Tisserand parameter
        q = a*(1 - e)
        return TJ < 3.05 and q < 7.35

    def is_trad_centaur(a, e):
        print("checking... traditional centaur")
        aN = 30 #AU
        return a < aN
    
    def is_scat_centaur(a, e):
        print("checking... scattering centaur")
        q = a*(1 - e)
        aN = 30 #AU
        return q < aN

    def is_oort(a):
        print("checking... oort")
        return a > 2000

    def is_scatter(data, N):
        print("checking... scatter")
        #different than gladman
        count = 0
        for i in range(N):
            delta_a = abs(data['a'].values[i] - data['a'].values[i][0]) 
            if max(delta_a)/a > (1.5/40):
                count +=1
                
        scattering = open(folder + tno + '/' + 'scattering.txt','w+')
        scattering.write(str(count))
        scattering.close()
    
        return count >= int(N/2)
    

    def is_detached(e):
        print("checking... detached")
        return e > 0.24


    if is_comet(a, e, i):
        classification.write("comet")
    elif is_trad_centaur(a, e):
        classification.write("traditional centaur")
    elif is_scat_centaur(a, e):
        classification.write("scattering centaur")
    elif is_oort(a):
        classification.write("inner oort cloud")
    elif is_scatter(data, N):
        classification.write("scattering")
    elif is_detached(e):
        classification.write("detached")
    else:
         classification.write("classical belt")

    classification.close()
    plt.close("all")