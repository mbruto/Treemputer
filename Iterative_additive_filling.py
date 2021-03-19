#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-


import argparse
import sys
import re
import os
import subprocess
import numpy
import ctypes
from multiprocessing import Pool

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from Amalphy.Amalphy_Script_Tree import readTree, getLeavesNames, getLeavesNumNames, distanceFrom
iaf_lib = ctypes.cdll.LoadLibrary(os.path.join(os.path.dirname(os.path.realpath(__file__)), "Library_iaf.so"))


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

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Impute missing distances from sparse matrxix')
    parser.add_argument('-t', '--tree_list', required=False, help="File containing one tree per line (newick)")
    parser.add_argument('-m', '--matrix', required=False, help="Matrix in phylip format with unknown values coded as \'-99.0\'")
    parser.add_argument('-o', '--output', required=True, help="Output file name; /!\ the extension \".phy\" will be added at the end")
    args = parser.parse_args()
    
    dict_dist = {}
    count = 1
    prev_unknown = -1

    if args.tree_list :
        nd_array_dist, tax_name_list = tree_2_dist(args.tree_list)
        phy_array_mean = nd_array_to_phylip(nd_array_dist, tax_name_list)
        MAT = open(args.output+".mean_dist.phy", "w")
        MAT.write(phy_array_mean)
        MAT.close()
        phy_array_mean = ""
        sys.stderr.write("Mean distance matrix was written in "+str(args.output)+".mean_dist.phy for future use ...\n")
        sys.stderr.flush()

    elif args.matrix :
    
        nd_array_dist, tax_name_list = phylip_to_ndarray(args.matrix)
    
    else :
        sys.stderr.write("Please provide a valid input file; either trees (-t | --tree_list) or matrix (-m | --matrix). See help for further informations\n")
        exit()

    temp_log = open("log_imput.txt", "w")
    temp_log.close()

    print("Number of unknow, values in the matrix before imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
        
    while len(numpy.where(nd_array_dist == -99.0)[0]) > 0 or len(numpy.where(nd_array_dist == -99.0)[0]) != prev_unknown:
        nd_array_dist = quartet_distance_imputation(nd_array_dist, len(tax_name_list) )
        print("Number of unknow, values in the matrix after imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
        prev_unknown = len(numpy.where(nd_array_dist == -99.0)[0])
        
    print("Number of unknow, values in the matrix after imputation : "+str(len(numpy.where(nd_array_dist == -99.0)[0])))
    
    phy_array_dist = nd_array_to_phylip(nd_array_dist, tax_name_list)
    outfile = open(str(args.output)+".mean.phy", "w")
    outfile.write(phy_array_dist)
    outfile.close()
