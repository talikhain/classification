# -*- coding: utf-8 -*-
"""
find resonance
"""
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
def find_resonance(folder, tno, N):
    def res_find(i, p, q, order, path, outer, c, d, inner, point_num):
        """
        arguments:
        i           index of TNO
        p           resonance coefficient
        q           resonance coefficient
        order       (int) r and s coefficients run from -order to order
        path        directory in which to save plots
        outer       dataframe with aei information of outer planet
        c           beginning of time interval, in number of points
        d           end of time interval, in number of points
        inner       dataframe with aei information of inner planet
        point_num   number of points in interval
        """

        #harmonics run-through
        R = np.arange(-1*order, order)  #harmonics preparation
        S = np.arange(-1*order, order)
        coeff_r = []
        coeff_s = []

        for r in R:
            for s in S:
                if p - q - r - s == 0:   #coefficients must add up to zero
                    coeff_r.append(r)
                    coeff_s.append(s)

        res_sum_pts = np.zeros((len(coeff_r), 4)) #res_sum (# grid squares in res), col_sum (# of intervals in res), r, s
        for u in range(len(coeff_r)):
            r = coeff_r[u]
            s = coeff_s[u]

            #resonance angle construction
            phi = p*(inner['peri'].values[c:d] + inner['node'].values[c:d] + inner['M'].values[c:d]) - q*(outer['peri'].values[i][c:d] + outer['node'].values[i][c:d] + outer['M'].values[i][c:d]) - r*(inner['peri'].values[c:d] + inner['node'].values[c:d]) - s*(outer['peri'].values[i][c:d] + outer['node'].values[i][c:d])

            #overlay a grid on the plot (phi vs. time)
            #j runs from time[c] to time[d], 20 columns
            #k runs from 0 to 360

            k_number = 18
            j_number = 20*int(np.ceil(abs(d-c)/point_num))

            scaling = int(abs(d-c)/j_number) #total number of points/j_number
            threshold = 1

            #approx rectangles:
            k_length = int(360/k_number)

            flags = np.zeros([k_number, j_number])

            #loop through grid:
            for j in range(0, j_number):
                for k in range(0, k_number):
                    #find number of points in one rectangle
                    array = phi[int(j*scaling):int((j+1)*scaling)]%360
                    pt_number = ((k*k_length <= array) & (array < (k+1)*k_length)).sum()

                    if pt_number > threshold:
                        flags[k][j] = 0
                    else:
                        flags[k][j] = 1

            #resonance check: number of white rectangles
            flag_col_sum = np.sum(flags, axis = 0)
            
            if any(m > 1 for m in flag_col_sum): 

                col_sum = sum(flag_col_sum > 1)*scaling #at least two squares must pass to be a resonance; stores # of points
                res_sum = sum(flag_col_sum) #number of passing squares

                #throw away false positive
                j_index = np.argmax(flag_col_sum)
                edge_ind = np.nonzero(flags[:, j_index])[0]
                res_flag = 0
                for ind_check in range(len(edge_ind)):
                    if flags[(edge_ind[ind_check] + 1)%k_number][j_index] < 1: #make sure that empty squares are mostly adjacent (one gap allowed)
                        res_flag +=1

                if res_flag > 1:    
                    res_sum = 0
                    col_sum = 0

                res_sum_pts[u][0] = res_sum #used to find best r & s
                res_sum_pts[u][1] = col_sum #used to find fraction of time in resonance
                res_sum_pts[u][2], res_sum_pts[u][3] = r, s
        
        return res_sum_pts, j_number, k_number

    def res_scan(N_max, res_width, i, T_av, order, path, outer, c, d, inner, resonance_table, point_num):
        """
        arguments:
        N_max            how many high-order resonances we consider
        res_width        considered period ratio range (ex: res_width = 0.25 --> T_av - 0.25, T_av + 0.25)
        i                index of simulation/TNO/etc
        T_av             average period ratio over given time interval
        order            (int) r and s coefficients run from -order to order
        path             directory in which to save plots
        outer            dataframe with aei information of outer planet
        c                beginning of time interval, in number of points
        d                end of time interval, in number of points
        inner            dataframe with aei information of inner planet
        resonance_table  list with important information about successful resonances
        point_num        number of points in interval
        """

        ratios = np.empty([N_max**2, 3])  #array: p, q, p/q
        for n in range(N_max):
            ratios[n*N_max:(n+1)*N_max,0] = np.arange(1, N_max+1)
            ratios[n*N_max:(n+1)*N_max,1] = np.ones(N_max)*(n+1)
            ratios[n*N_max:(n+1)*N_max,2] = ratios[n*N_max:(n+1)*N_max, 0]/ratios[n*N_max:(n+1)*N_max, 1]

        u, indices = np.unique(ratios[:,2], return_index=True) #removing repeating rows
        newratios = ratios[indices,:]

        N_res = len(newratios[:,2]) #number of resonances we consider

        T_max = T_av + res_width
        T_min = T_av - res_width

        for k in range(N_res):
            if newratios[k][2] <= T_max and newratios[k][2] >= T_min:
                p = min(newratios[k][0], newratios[k][1]) #need to be changed based on outer/inner
                q = max(newratios[k][0], newratios[k][1])

                res_sum_pts, j_number, k_number = res_find(i, p, q, order, path, outer, c, d, inner, point_num)

                #save the best r & s combination for a given p & q
                if sum(res_sum_pts[:,0] > 0) > 0:

                    u = np.argmax(res_sum_pts[:,0])
                    r, s = res_sum_pts[u][2], res_sum_pts[u][3]
                    
                    if res_sum_pts[u][0] > j_number*k_number/9:
                        
                        resonance_table.append([i, p, q, max(res_sum_pts[:,0]), max(res_sum_pts[:,1])])
                        if min(len(inner['time'].values), len(outer['time'].values[i])) < d:
                            d = min(len(inner['time'].values), len(outer['time'].values[i]))

                        phi = p*(inner['peri'].values[c:d] + inner['node'].values[c:d] + inner['M'].values[c:d]) - q*(outer['peri'].values[i][c:d] + outer['node'].values[i][c:d] + outer['M'].values[i][c:d]) - r*(inner['peri'].values[c:d] + inner['node'].values[c:d]) - s*(outer['peri'].values[i][c:d] + outer['node'].values[i][c:d])

                        plt.plot(outer['time'].values[i][c:d]*1e-6, phi%360, linestyle='None', marker='.', markersize=0.9)
                        plt.xlim(outer['time'].values[i][c]*1e-6, outer['time'].values[i][d-1]*1e-6)
                        plt.ylim(0, 360)
                        plt.ylabel(r'\textbf{$\phi_{res}$ (deg)}', fontsize=16)
                        plt.xlabel(r'\textbf{time (Myr)}', fontsize=16)

                        if int(r) == 0:
                            plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} - %s\varpi_{TNO}$}' % (int(p), int(q), int(s)))
                            if int(s) == 0:
                                plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO}$}' % (int(p), int(q)))
                            if int(s) < 0:
                                plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} + %s\varpi_{TNO}$}' % (int(p), int(q), abs(int(s))))

                        if int(r) != 0:
                            plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} - %s\varpi_{N} - %s\varpi_{TNO}$}' % (int(p), int(q), int(r), int(s)))
                            if int(s) == 0:
                                plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} - %s\varpi_{N}$}' % (int(p), int(q), int(r)))
                                if int(r) < 0:
                                    plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} + %s\varpi_{N}$}' % (int(p), int(q), abs(int(r))))

                            if int(r) < 0:
                                plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} + %s\varpi_{N} - %s\varpi_{TNO}$}' % (int(p), int(q), abs(int(r)), int(s)))
                            if int(s) < 0:
                                plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} - %s\varpi_{N} + %s\varpi_{TNO}$}' % (int(p), int(q), int(r), abs(int(s))))
                            if int(r) < 0 and int(s) < 0:
                                plt.title(r'\textbf{$\phi = %s\lambda_{N} - %s\lambda_{TNO} + %s\varpi_{N} + %s\varpi_{TNO}$}' % (int(p), int(q), abs(int(r)), abs(int(s))))

                        plt.savefig('%s/%s_%s_%s_%s_%s.png' % (path, i, int(p), int(q), u, c), dpi=100)
                        plt.gcf().clear()

    def res_av(chunk_length, chunk_num, point_num, inner, outer):
        """
        arguments:
        chunk_length   length of time interval over which period ratio is averaged
        chunk_num      number of chunks in total integration time
        point_num      number of points per chunk
        i              index of simulation/TNO/etc
        inner          dataframe with aei information of inner planet
        outer          dataframe with aei information of outer planet
        """

        T_av = []

        for m in range(chunk_num):
            T_av.append(np.mean((outer['a'].values[i][m*point_num:(m+1)*point_num]/inner['a'].values[m*point_num:(m+1)*point_num])**1.5))

        error = 0.01 #0.01 #may need adjusting
        k = 0
        l = 0
        T_res = []

        c_pts = []
        d_pts = []
        while k < len(T_av) - 1:
            T = []
            l_range = []
            while abs(T_av[l] - T_av[l+1]) < error:  #collect steps that are "close"

                l_range.append(l)
                T.append(T_av[l])
                l +=1
                if l+1 == len(T_av):
                    break

            T.append(T_av[l])
            l_range.append(l)

            if len(T) > 0:  #average over the "close" steps, 
                T_res.append(np.mean(T))  #save the average 
                c_pts.append(point_num*l_range[0]) #save start and stop intervals
                d_pts.append(point_num*(l_range[-1]+1))

            l +=1
            k = l

        return T_res, c_pts, d_pts

######################################################

    headings = ['time', 'M', 'a', 'e', 'i', 'peri', 'node']
    columns = [0, 2, 3, 4, 5, 6, 7]
    Nep = pd.read_csv(folder + tno + '/' + 'NEPTUNE.aei', skiprows=4, delim_whitespace=True, usecols=columns, names=headings)
    data = pd.DataFrame()
    for i in range(N):
        tmp = pd.read_csv(folder + tno + '/' + 't_%s.aei' % (i), skiprows=4, delim_whitespace=True, usecols=columns, names=headings)
        data = data.append({'time':tmp['time'].values, 'a':tmp['a'].values, 'e':tmp['e'].values, 'i':tmp['i'].values, 'node':tmp['node'].values, 'peri':tmp['peri'].values, 'M':tmp['M'].values}, ignore_index = True)

    t0 = tm.time()
    N_max = 25
    res_width = 0.2 
    order = 25
    path = folder + tno + '/' + 'res' 

    resonance_table = []  

    ejected_flag = 0
    for i in range(N):
        print('index:', i)
        if len(data['time'].values[i]) < len(Nep):
            print('index:', i, 'is skipped')
            ejected_flag +=1
            continue
        chunk_length = 0.5e7
        chunk_num = abs(int(np.ceil(data['time'].values[i][-1]/(chunk_length)))) #number of chunks, rounded up
        pt_time_conv = Nep['time'].values[-1]/len(Nep)
        T_chunk = chunk_length/pt_time_conv
        point_num = abs(int(T_chunk)) #number of points in chunk, generally 5000 is a good number

        T_res, c_pts, d_pts = res_av(chunk_length, chunk_num, point_num, Nep, data)

        for f in range(len(c_pts)):
            if d_pts[f] == c_pts[f]:
                d_pts[f] = c_pts[f+1] 

        for j in range(len(T_res)):
            T_av = T_res[j]
            c = c_pts[j]
            d = d_pts[j]
            res_scan(N_max, res_width, i, T_av, order, path, data, c, d, Nep, resonance_table, point_num)
    
    if ejected_flag > 0:
        N = N - ejected_flag
        
    ejected = open(folder + tno + '/' + 'ejected.txt','w+')
    ejected.write(str(ejected_flag))
    ejected.close()

    t1 = tm.time()
    print("resonance identification took t =", (t1-t0)/60, "minutes")

    ################ end of identification ############
    print("...plotting summary plots...")

    with open(path + '/' + 'resonance_table', 'wb') as f:
        pickle.dump(resonance_table, f)

    total_points = 0
    for i in range(N):
        total_points += len(data['time'].values[i])

    total_res_sum = 0
    for z in range(len(resonance_table)):
        total_res_sum += resonance_table[z][4]

    print(round(total_res_sum/total_points*100, 3), '% of time in resonance')
    resonance = open(folder + tno + '/' + 'resonance.txt','w+')
    if round(total_res_sum/total_points*100, 3) > 0.00001:
        k = 0
        summed_res_table = [[resonance_table[k][1], resonance_table[k][2], resonance_table[k][4]]]
        for k in range(1, len(resonance_table)):  #take a p and a q
            flag = 0
            for h in range(len(summed_res_table)): #go through the rows of the master list
                if resonance_table[k][1] == summed_res_table[h][0] and resonance_table[k][2] == summed_res_table[h][1]: #if we find a match
                    summed_res_table[h][2] += resonance_table[k][4] #add
                else:
                    flag +=1
            if flag == len(summed_res_table):
                summed_res_table.append([resonance_table[k][1], resonance_table[k][2], resonance_table[k][4]])

        sorted_res = np.zeros([len(summed_res_table), 4])
        for i in range(len(summed_res_table)):
            sorted_res[i][0] = summed_res_table[i][0]/summed_res_table[i][1]
            sorted_res[i][1], sorted_res[i][2], sorted_res[i][3] =  summed_res_table[i][0], summed_res_table[i][1], summed_res_table[i][2]

        ind = np.argsort(sorted_res[:, 0])
        sorted_res = sorted_res[ind.ravel(),:]

        fig = plt.figure(figsize=(6, 4)) #12, 4
        gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) #5, 1

        ax0 = plt.subplot(gs[0])
        objects = []
        res = []
        for l in range(len(sorted_res)):
            objects.append(r'\textbf{%s:%s}' % (int(sorted_res[l][1]), int(sorted_res[l][2])))
            res.append(sorted_res[l][3]/total_points)

        y_pos = np.arange(len(objects))

        barlist = plt.bar(y_pos, res, width=0.8, align='center', alpha=0.5, color='royalblue')
        plt.xticks(y_pos, objects)
        plt.ylabel(r'\textbf{fraction of time}', fontsize=16)
        plt.axis('tight')
        fig.tight_layout()

        ax1 = plt.subplot(gs[1])
        objects2 = []
        res2 = []
        objects2.append(r'\textbf{res}')
        res2.append(sum(res))

        objects2.append(r'\textbf{non-res}')
        res2.append(1 - sum(res))

        y_pos2 = np.arange(len(objects2))

        barlist = plt.bar(y_pos2, res2, width=0.8, align='center', alpha=0.5, color='royalblue')
        barlist[-1].set_color('salmon')
        plt.xticks(y_pos2, objects2)
        plt.axis('tight')
        fig.tight_layout()
        plt.savefig(path + '/' + 'res_summary.png', dpi=300)
        plt.close("all")

        
        resonance.write("resonant"+'\n')
        resonance.write(str(round(sum(res)*100, 3))+'\n')
        resonance.write(str(round(max(res)*100, 3))+'\n') 
        resonance.write('%s:%s' % (int(sorted_res[np.argmax(res)][1]), int(sorted_res[np.argmax(res)][2])))
    resonance.close()
