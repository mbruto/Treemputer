#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-


import argparse
import sys
import re
import os
import subprocess
import numpy
from multiprocessing import Pool
from Amalphy_Script_Tree import readTree, getLeavesNames, getLeavesNumNames, distanceFrom

import ctypes
iaf_lib = ctypes.cdll.LoadLibrary("Library_iaf.so")

#------------------
#

def process_input (sub_options) :

    sub_list = []

    for opt in sub_options:
        
        for file in sorted(os.listdir(opt)) :
            
            sub_list.append(opt+"/"+file)
    
    return sub_list


#------------------

def phylip_to_ndarray(sub_phy):

    with open(sub_phy) as phy_file:
        strain_num = int(phy_file.readline())
        n = 0
        sub_id_list = []
        
        sub_ndarray = numpy.full(shape=(strain_num,strain_num), fill_value=-99.0, dtype=float)
        for a in range(0, int(strain_num)):
            sub_ndarray[a][a] = 0.0
        
        for line in phy_file:
        
            if not line.strip():
                continue
        
            temp_tab = re.split("\s+", line)
            sub_id_list.append(temp_tab[0])
            
            for j in range(1, int(strain_num)+1):
                sub_ndarray[n, j-1] = temp_tab[j]
                
            n += 1
    
    return sub_ndarray, sub_id_list

#------------------
#

def nd_array_to_phylip(sub_ndarray, sub_id_list):

    output = ""

    output += "\t"+str(len(sub_id_list))+"\n"

    for m in range(0, len(sub_id_list)) :
        name = sub_id_list[m] + ' ' * (10 - len(sub_id_list[m]))
        output += str(name)+" "
        for n in range(0, len(sub_id_list)):
            output += str(round(sub_ndarray[m][n], 6))+" "
        output += "\n"
    
    return output

#------------------
#

def tree_2_taxa_list(sub_file) :

    taxa_list = []
    sub_tree_num = 0

    with open(sub_file) as sub_tree_file :
        for line in sub_tree_file :
        
            if not line.rstrip() :
                continue
        
            Tree_obj = readTree(line.rstrip())
            for tax in getLeavesNames(Tree_obj) :
                if tax not in taxa_list :
                    taxa_list.append(tax)
            
            sub_tree_num += 1
    
    return taxa_list, sub_tree_num

#------------------
#
'''
def mean_dist_calc ( sub_a, sub_b, sub_list, sub_tree, sub_array ) :

    print(str(sub_a)+" "+str(sub_b))

    total_dist = 0
    count = 0
    mean_dist = -99.0

    cmd_grep = "grep \'"+sub_list[sub_a]+":\' "+str(sub_tree)+" |grep \'"+sub_list[sub_b]+":\'"
    result = subprocess.getoutput(cmd_grep)
    
    if len(result.split("\n")) == 1 :
        return sub_array
    
    else :
        for tree in result.split("\n") :
            Tree_obj = readTree(tree)
            dict_uto = getLeavesNumNames(Tree_obj)
            total_dist += distanceFrom(Tree_obj, dict_uto[sub_list[sub_a]], dict_uto[sub_list[sub_b]])
            count += 1
        
        sub_array[sub_a][sub_b] = total_dist / count
        sub_array[sub_b][sub_a] = total_dist / count
        return sub_array
'''
#------------------
#

'''
def tree_2_ndarray(sub_tree_file, sub_threads) :

    sub_tax_list = tree_2_taxa_list(sub_tree_file)
    print(sub_tax_list)
    
    #Creation of distance matrix ndarray
    sub_ndarray = numpy.full(shape=(len(sub_tax_list),len(sub_tax_list)), fill_value=-99.0, dtype=float)
    
    for i in range(0, len(sub_tax_list)) :
        print(i)
        

    return sub_ndarray, sub_tax_list
'''

#------------------
#Ameliorable !!!

def tree_2_dist (sub_tree_file) :

    count = 1
    sub_dict = {}
    it = 0.05
    thresh = 0

    sub_tax_list, tree_num = tree_2_taxa_list(sub_tree_file)
    
    sub_ndarray = numpy.full(shape=(len(sub_tax_list),len(sub_tax_list)), fill_value=-99.0, dtype=float)
    
    sys.stderr.write("Importing trees and creating a distance matrix with mean values ...\n")
    sys.stderr.write("\r %a %%" % thresh)
    sys.stderr.flush()
    thresh += 5
    
    with open(sub_tree_file) as tree_file :
        for line in tree_file :
        
            if not line.rstrip() :
                continue
            
            p = (count / tree_num) * 100
            if p > thresh :
                sys.stderr.write("\r %a %%" % thresh)
                sys.stderr.flush()
                thresh += 5

            sub_id_dict = {}
            Tree_obj = ""
            
            Tree_obj = readTree(line.rstrip())
            sub_id_dict = getLeavesNumNames(Tree_obj)
            list_keys = list(sub_id_dict.keys())

            for k in range(0, len(list_keys)) :
                i = sub_tax_list.index(list_keys[k])
                if i not in sub_dict :
                    sub_dict.update( { i : {} } )
                for l in range(k+1, len(list_keys)) :
                    j = sub_tax_list.index(list_keys[l])
                    if j not in sub_dict :
                        sub_dict.update( { j : {} } )
                    if i not in sub_dict[j] :
                        sub_dict[j].update( { i : { 'total_dist' : 0, 'count' : 0 } } )
                    if j not in sub_dict[i] :
                        sub_dict[i].update( { j : { 'total_dist' : 0, 'count' : 0 } } )
            
                    dist = distanceFrom(Tree_obj, sub_id_dict[list_keys[k]], sub_id_dict[list_keys[l]])
                    sub_dict[i][j]['total_dist'] += dist
                    sub_dict[i][j]['count'] += 1
                    sub_dict[j][i]['total_dist'] += dist
                    sub_dict[j][i]['count'] += 1
            
            count += 1
            
    thresh += 5
    sys.stderr.write("\r %a %%\n" % thresh)
    sys.stderr.flush()
    
    for i in range(0, len(sub_tax_list)):

        for j in range(0, len(sub_tax_list)):
        
            if i == j:
                sub_ndarray[i, j] = 0.0
            
            elif j in sub_dict[i]:
                sub_ndarray[i, j] = round(sub_dict[i][j]["total_dist"] / sub_dict[i][j]["count"], 8)
                sub_ndarray[j, i] = round(sub_dict[i][j]["total_dist"] / sub_dict[i][j]["count"], 8)
            
            elif i in sub_dict[j]:
                sub_ndarray[i, j] = round(sub_dict[j][i]["total_dist"] / sub_dict[j][i]["count"], 8)
                sub_ndarray[j, i] = round(sub_dict[j][i]["total_dist"] / sub_dict[j][i]["count"], 8)
    
    return sub_ndarray, sub_tax_list
    
#------------------
#Ameliorable !!!
    
def tree_2_dist_var (sub_tree_file) :

    count = 1
    sub_dict = {}
    it = 0.05
    thresh = 0
    temp_list_var = []

    sub_tax_list, tree_num = tree_2_taxa_list(sub_tree_file)
    
    sub_ndarray = numpy.full(shape=(len(sub_tax_list),len(sub_tax_list)), fill_value=-99.0, dtype=float)
    
    sys.stderr.write("Importing trees and creating a distance matrix with variance values ...\n")
    sys.stderr.write("\r %a %%" % thresh)
    sys.stderr.flush()
    thresh += 5
    
    with open(sub_tree_file) as tree_file :
        for line in tree_file :
        
            if not line.rstrip() :
                continue
            
            p = (count / tree_num) * 100
            if p > thresh :
                sys.stderr.write("\r %a %%" % thresh)
                sys.stderr.flush()
                thresh += 5

            sub_id_dict = {}
            Tree_obj = ""
            
            Tree_obj = readTree(line.rstrip())
            sub_id_dict = getLeavesNumNames(Tree_obj)
            list_keys = list(sub_id_dict.keys())

            for k in range(0, len(list_keys)) :
                i = sub_tax_list.index(list_keys[k])
                if i not in sub_dict :
                    sub_dict.update( { i : {} } )
                for l in range(k+1, len(list_keys)) :
                    j = sub_tax_list.index(list_keys[l])
                    if j not in sub_dict :
                        sub_dict.update( { j : {} } )
                    if i not in sub_dict[j] :
                        sub_dict[j].update( { i : { 'total_dist' : "" } } )
                    if j not in sub_dict[i] :
                        sub_dict[i].update( { j : { 'total_dist' : "" } } )
            
                    dist = distanceFrom(Tree_obj, sub_id_dict[list_keys[k]], sub_id_dict[list_keys[l]])
                    sub_dict[i][j]['total_dist'] += "XXX"+str(dist)
                    sub_dict[j][i]['total_dist'] += "XXX"+str(dist)
            
            count += 1
            
    thresh += 5
    sys.stderr.write("\r %a %%\n" % thresh)
    sys.stderr.flush()
    
    for i in range(0, len(sub_tax_list)):

        for j in range(0, len(sub_tax_list)):
        
            if i == j:
                sub_ndarray[i, j] = 0.0
            
            elif j in sub_dict[i]:
            
                temp_list_var = sub_dict[i][j]["total_dist"].split("XXX")
                temp_list_var.pop(0)
                
                if len(temp_list_var) > 1 :
                    sub_ndarray[i, j] = round(numpy.var(list(map(float, temp_list_var))), 8)
                    sub_ndarray[j, i] = round(numpy.var(list(map(float, temp_list_var))), 8)
                else :
                    sub_ndarray[i, j] = temp_list_var[0]
                    sub_ndarray[j, i] = temp_list_var[0]
                
            elif i in sub_dict[j]:
                temp_list_var = sub_dict[i][j]["total_dist"].split("XXX")
                temp_list_var.pop(0)
                
                if len(temp_list_var) > 1 :
                    sub_ndarray[i, j] = round(numpy.var(list(map(float, temp_list_var))), 8)
                    sub_ndarray[j, i] = round(numpy.var(list(map(float, temp_list_var))), 8)
                else :
                    sub_ndarray[i, j] = temp_list_var[0]
                    sub_ndarray[j, i] = temp_list_var[0]
    
    return sub_ndarray
    
#------------------
#

#Python implementation of quartet distance estimation
'''
def quartet_distance_imputation(sub_nd_array, sub_matrix_size):

    for i in range(0, sub_matrix_size):
    
        print(i)

        for j in range(i+1, sub_matrix_size):
    
            if sub_nd_array[i][j] != -99.0:
                continue
            
            k_list = []
            l_list = []
            opti_dist = 0.0
            
            for x in range(0, sub_matrix_size) :
                if sub_nd_array[i][x] != -99.0 and sub_nd_array[j][x] != -99.0 :
                    
                    k_list.append((x, sub_nd_array[i][x]))
                    l_list.append((x, sub_nd_array[j][x]))
            
            k_list = sorted(k_list, key=ValueSort)
            l_list = sorted(l_list, key=ValueSort)
            
            
            for k in range(0, len(k_list)) :
            
                k_index = k_list[k][0]
            
                for l in range(0, len(l_list)) :
                
                    l_index = l_list[l][0]
                    
                    if(k_index == l_index) :
                        continue
                    
                    if sub_nd_array[k_index][l_index] == -99.0 :
                        continue
                    
                    if sub_nd_array[k_index][l_index] < sub_nd_array[i][k_index] and sub_nd_array[k_index][l_index] < sub_nd_array[i][l_index] and sub_nd_array[k_index][l_index] < sub_nd_array[j][k_index] and sub_nd_array[k_index][l_index] < sub_nd_array[j][l_index] :
                        continue
                    
                    d_1 = 0
                    d_2 = 0
                    d_3 = 0
                    d_4 = 0
                    d_5 = 0
                    
                    if sub_nd_array[i][k_index] > sub_nd_array[i][l_index] and sub_nd_array[j][l_index] > sub_nd_array[j][k_index] and sub_nd_array[i][k_index] < sub_nd_array[k_index][l_index] and sub_nd_array[j][l_index] < sub_nd_array[k_index][l_index] :
                        
                        d_1 = 0.5 * (sub_nd_array[i][l_index] + sub_nd_array[i][k_index] - sub_nd_array[k_index][l_index]);
                        d_2 = 0.5 * (sub_nd_array[i][l_index] + sub_nd_array[k_index][l_index] - sub_nd_array[i][k_index]);
                        d_3 = 0.5 * (sub_nd_array[j][l_index] + sub_nd_array[j][k_index] - sub_nd_array[k_index][l_index]);
                        d_4 = 0.5 * (sub_nd_array[j][k_index] + sub_nd_array[k_index][l_index] - sub_nd_array[j][l_index]);
                        d_5 = 0.5 * (sub_nd_array[i][k_index] + sub_nd_array[j][l_index] - sub_nd_array[i][l_index] - sub_nd_array[j][k_index]);
                        
                        if d_1 < 0 or d_2 < 0 or d_3 < 0 or d_4 < 0 or d_5 < 0 :
                            continue
                        
                        if d_5 < d_1 or d_5 < d_2 or d_5 < d_3 or d_5 < d_4  :
                            continue
                        
                        opti_dist = sub_nd_array[j][l_index] + sub_nd_array[i][k_index] - sub_nd_array[k_index][l_index];
                        
                        if opti_dist > 0.0 :
                            sub_nd_array[i, j] = opti_dist
                            sub_nd_array[j, i] = opti_dist
                        else :
                            opti_dist = 0.0
                    
                    
                    
                    elif sub_nd_array[i][k_index] < sub_nd_array[i][l_index] and sub_nd_array[j][l_index] < sub_nd_array[j][k_index] and sub_nd_array[i][l_index] < sub_nd_array[k_index][l_index] and sub_nd_array[j][k_index] < sub_nd_array[k_index][l_index] :
                        
                        d_1 = 0.5 * (sub_nd_array[i][k_index] + sub_nd_array[i][l_index] - sub_nd_array[k_index][l_index]);
                        d_2 = 0.5 * (sub_nd_array[i][k_index] + sub_nd_array[k_index][l_index] - sub_nd_array[i][l_index]);
                        d_3 = 0.5 * (sub_nd_array[j][k_index] + sub_nd_array[j][l_index] - sub_nd_array[k_index][l_index]);
                        d_4 = 0.5 * (sub_nd_array[j][l_index] + sub_nd_array[k_index][l_index] - sub_nd_array[j][k_index]);
                        d_5 = 0.5 * (sub_nd_array[i][l_index] + sub_nd_array[j][k_index] - sub_nd_array[i][k_index] - sub_nd_array[j][l_index]);
                        
                        if d_1 < 0 or d_2 < 0 or d_3 < 0 or d_4 < 0 or d_5 < 0 :
                            continue
                        
                        if d_5 < d_1 or d_5 < d_2 or d_5 < d_3 or d_5 < d_4  :
                            continue
                        
                        opti_dist = sub_nd_array[j][k_index] + sub_nd_array[i][l_index] - sub_nd_array[k_index][l_index];
                        
                        if opti_dist > 0.0 :
                            sub_nd_array[i, j] = opti_dist
                            sub_nd_array[j, i] = opti_dist
                        else :
                            opti_dist = 0.0
                    
                    if opti_dist != 0 :
                        break
                
                if opti_dist != 0 :
                    break
                    
    return sub_nd_array
'''
            
            
#------------------
#

def quartet_distance_imputation(sub_nd_array, sub_matrix_size):

    func_iaf = iaf_lib.distance_estimation_quartet
    func_iaf.restype = ctypes.c_double
    dist_est = ctypes.c_double()

    for i in range(0, sub_matrix_size):

        for j in range(i+1, sub_matrix_size):
    
            if sub_nd_array[i][j] != -99.0:
                continue
            
            dist_est = func_iaf(i, j, sub_matrix_size, sub_nd_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            sub_nd_array[i, j] = dist_est
            sub_nd_array[j, i] = dist_est

            #print(str(i)+" "+str(j)+" "+str(dist_est))
    
    return sub_nd_array


#------------------
#

def quartet_distance_imputation_mean(sub_nd_array, sub_var_array, sub_matrix_size):

    func_iaf = iaf_lib.distance_estimation_quartet_mean
    func_iaf.restype = ctuple
    dist_est = ctuple

    for i in range(0, sub_matrix_size):

        for j in range(i+1, sub_matrix_size):
    
            if sub_nd_array[i][j] != -99.0:
                continue
            
            dist_est = func_iaf(i, j, sub_matrix_size, sub_nd_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 10)
            
            sub_nd_array[i, j] = dist_est.value_1
            sub_nd_array[j, i] = dist_est.value_1
            
            sub_var_array[i, j] = dist_est.value_2
            sub_var_array[j, i] = dist_est.value_2
    
    return sub_nd_array, sub_var_array

#------------------
#
'''
def min_distance_imputation(sub_nd_array, sub_matrix_size) :

    for i in range(0, sub_matrix_size):
    
        print(i)

        for j in range(i+1, sub_matrix_size):
        
            test_dist = 0
            opti_dist = 0.0
            Min_dist = 10000.0
            
            if sub_nd_array[i][j] == -99.0:
                
                for k in range(0, sub_matrix_size):
                
                    if k == i or k == j:
                        continue
                
                    #All entries except ij must be known
                    if sub_nd_array[i][k] == -99.0 or sub_nd_array[j][k] == -99.0:
                        continue
                    
                    for l in range(k+1, sub_matrix_size):
                        
                        if l == i or l == j :
                            continue
                        
                        #All entries except ij must be known
                        if sub_nd_array[i][l] == -99.0 or sub_nd_array[j][l] == -99.0 or sub_nd_array[k][l] == -99.0 :
                            continue
                        
                        test_dist = sub_nd_array[i][k] + sub_nd_array[i][l] + sub_nd_array[j][k] + sub_nd_array[j][l] + sub_nd_array[k][l]
                        
                        if test_dist < Min_dist :
                            
                            if abs(sub_nd_array[i][k] + sub_nd_array[j][l] - sub_nd_array[i][l] + sub_nd_array[j][k]) < 0.00001:
                                continue

                            elif sub_nd_array[i][k] + sub_nd_array[j][l] > sub_nd_array[i][l] + sub_nd_array[j][k] :
                                if sub_nd_array[i][k] + sub_nd_array[j][l] - sub_nd_array[k][l] > 0 :
                                    Min_dist = test_dist;
                                    opti_dist = sub_nd_array[i][k] + sub_nd_array[j][l] - sub_nd_array[k][l]
                                    
                            
                            elif sub_nd_array[i][k] + sub_nd_array[j][l] < sub_nd_array[i][l] + sub_nd_array[j][k] :
                                if sub_nd_array[j][k] + sub_nd_array[i][l] - sub_nd_array[k][l] > 0 :
                                    Min_dist = test_dist;
                                    opti_dist = sub_nd_array[j][k] + sub_nd_array[i][l] - sub_nd_array[k][l]
                                    
                
                if opti_dist > 0.0 :
                    sub_nd_array[i, j] = opti_dist
                    sub_nd_array[j, i] = opti_dist
    
    return sub_nd_array
'''

#------------------
#

def min_distance_imputation(sub_nd_array, sub_matrix_size) :

    func_iaf = iaf_lib.min_dist_imputation
    func_iaf.restype = ctypes.c_double
    dist_est = ctypes.c_double()

    for i in range(0, sub_matrix_size):

        for j in range(i+1, sub_matrix_size):
    
            if sub_nd_array[i][j] != -99.0:
                continue
            
            dist_est = func_iaf(i, j, sub_matrix_size, sub_nd_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            sub_nd_array[i, j] = dist_est
            sub_nd_array[j, i] = dist_est
            
            print(dist_est)
            exit()
    
    return sub_nd_array

#------------------
#

def discard_unknown_leaf (sub_nd_array, sub_tax_name_list) :


    while len(numpy.where(sub_nd_array == -99.0)[0]) > 0 :

        x = numpy.where(sub_nd_array == -99.0)[0]
        dict_count = dict((i, list(x).count(i)) for i in x)
        key_to_del = max(dict_count, key=dict_count.get)
        print(key_to_del)
        
        if dict_count[key_to_del] > 1 :
            sub_nd_array = numpy.delete(sub_nd_array, key_to_del, axis=0)
            sub_nd_array = numpy.delete(sub_nd_array, key_to_del, axis=1)
            sub_tax_name_list.pop(key_to_del)
        
        
        
    return sub_nd_array, sub_tax_name_list
    
#------------------
#

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Impute missing distances from sparse matrxix')
    parser.add_argument('-t', '--tree_list', required=False, help="File containing one tree per line (newick)")
    parser.add_argument('-m', '--matrix', required=False, help="Matrix in phylip format with unknown values coded as \'-99.0\'")
    parser.add_argument('-d', '--method', required=True, help="Method of imputation to use : \'min\' for slower matrix exploration ; \'quartet\' for faster imputation of missing distances; 'both' for using fast then slow imputation of missing distances")
    parser.add_argument('-o', '--output', required=True, help="Output file name; /!\ the extension \".phy\" will be added at the end")
    #parser.add_argument('-T', '--threads', required=True, help="Number of threads to use")
    args = parser.parse_args()
    
    dict_dist = {}
    count = 1
    prev_unknown = -1

    if args.tree_list :
        nd_array_dist, tax_name_list = tree_2_dist(args.tree_list)
        phy_array_mean = nd_array_to_phylip(nd_array_dist, tax_name_list)
        #MAT = open(args.output+".mean_dist.phy", "w")
        #MAT.write(phy_array_mean)
        #MAT.close()
        #phy_array_mean = ""
        #sys.stderr.write("Mean distance matrix was written in "+str(args.output)+".mean_dist.phy for future use ...\n")
        #sys.stderr.flush()
        
        '''
        nd_array_var = tree_2_dist_var(args.tree_list)
        phy_array_var = nd_array_to_phylip(nd_array_var, tax_name_list)
        MAT = open(args.output+".var_dist.phy", "w")
        MAT.write(phy_array_var)
        MAT.close()
        phy_array_var = ""
        sys.stderr.write("Variance distance matrix was written in "+str(args.output)+".var_dist.phy for future use ...\n")
        sys.stderr.flush()
        '''

    elif args.matrix :
    
        nd_array_dist, tax_name_list = phylip_to_ndarray(args.matrix)
    
    else :
        sys.stderr.write("Please provide a valid input file; either trees (-t | --tree_list) or matrix (-m | --matrix). See help for further informations\n")
        exit()

    temp_log = open("log_imput.txt", "w")
    temp_log.close()

    if args.method == 'min' :
        print("Number of unknow, values in the matrix before imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
        nd_array_dist = min_distance_imputation(nd_array_dist, len(tax_name_list) )
        print("Number of unknow, values in the matrix after imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
    elif args.method == 'quartet' :
        print("Number of unknow, values in the matrix before imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
        
        while len(numpy.where(nd_array_dist == -99.0)[0]) > 0 or len(numpy.where(nd_array_dist == -99.0)[0]) != prev_unknown:
            nd_array_dist = quartet_distance_imputation(nd_array_dist, len(tax_name_list) )
            print("Number of unknow, values in the matrix after imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
            prev_unknown = len(numpy.where(nd_array_dist == -99.0)[0])
        
        print("Number of unknow, values in the matrix after imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
    elif args.method == 'both' :
        print("Number of unknow, values in the matrix before quartet imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
        nd_array_dist = quartet_distance_imputation(nd_array_dist, len(tax_name_list) )
        print("Number of unknow, values in the matrix after quartet imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
        nd_array_dist = min_distance_imputation(nd_array_dist, len(tax_name_list) )
        print("Number of unknow, values in the matrix after min imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
        
    else :
        print("ERROR : please provide a valid method to impute missing distances => either 'min' or 'quartet'")
        exit()
    
    #nd_array_dist, tax_name_list = discard_unknown_leaf(nd_array_dist, tax_name_list)
    
    
    phy_array_dist = nd_array_to_phylip(nd_array_dist, tax_name_list)
    outfile = open(str(args.output)+".mean.phy", "w")
    outfile.write(phy_array_dist)
    outfile.close()
    
    '''
    phy_array_var = nd_array_to_phylip(nd_array_var, tax_name_list)
    outfile = open(str(args.output)+".var.phy", "w")
    outfile.write(phy_array_var)
    outfile.close()
    '''
    
