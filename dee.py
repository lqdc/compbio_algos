#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Author: Roman Sinayev
Date: Thu May  3 04:35:00 EDT 2012
Description: 
DEE algorithm for singles and pairs using Goldstein elimination criteria
"""
from numpy import *
from pylab import mlab
from time import time
#import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

def symmetrize(a):
    return a + a.T - diag(a.diagonal())

def get_digonal_block(matrix, i, current_sub_height):
    start = 0
    if i > 0:
        for block_height in current_sub_height[0:i]:
            start += block_height
            end = start + current_sub_height[i]
    else:
        end = current_sub_height[0]
    return matrix[start:end, start:end]

def get_pairs_block(matrix, i, j, current_sub_height):
    start1 = 0
    start2 = 0
    if i > 0:
        for block_height in current_sub_height[0:i]:
            start1 += block_height
            end1 = start1 + current_sub_height[i]
    else:
        end1 = current_sub_height[0]
    if j > 0:
        for block_height in current_sub_height[0:j]:
            start2 += block_height
            end2 = start2 + current_sub_height[j]
    else:
            end2 = current_sub_height[0]
    return matrix[start1:end1, start2:end2]
    
def delete_stuff(sym_matrix, eliminated_rows):
    sym_matrix = delete(sym_matrix, eliminated_rows, 0)
    sym_matrix = delete(sym_matrix, eliminated_rows, 1)
    return sym_matrix

def eliminate_singles(sym_matrix, positions, current_sub_height, height_now , eliminated, pairs_binary):
    #print 'begin singles'
    temp_sub_height = current_sub_height.copy()
    print current_sub_height

    while eliminated:
        height_now.append(current_sub_height)
        eliminated = False
        pairs_matrices = []
        mask_matrices = []
        eliminated_rows = []
        total_min = 0
        absolute_row_index = -1
        # absolute_row2_index = 0
        
        for block_row_i in xrange(positions):
            diag_block = get_digonal_block(sym_matrix, block_row_i , current_sub_height)
            for block_column_j in xrange(positions):
                if block_column_j != block_row_i:
                    pairs_matrices.append(get_pairs_block(sym_matrix, block_row_i, block_column_j, current_sub_height))
                    mask_matrices.append(get_pairs_block(pairs_binary, block_row_i, block_column_j, current_sub_height))
            for s in xrange(current_sub_height[block_row_i]):
                absolute_row_index += 1
                s_value = diag_block[s, s] 
                for t in xrange(current_sub_height[block_row_i]):
                    if s != t:
                        total_min=0
                        t_value = diag_block[t,t]
                        for m_index, pairs_matrix in enumerate(pairs_matrices):
                            min_row = pairs_matrix[s, :] - pairs_matrix[t, :]
                            min_row_mask_s = mask_matrices[m_index][s, :]
                            min_row_zeros = min_row[where(min_row_mask_s == 0)[0]]
                            if len(min_row_zeros) > 0:
                                total_min+=min_row_zeros.min()
                            else:
                                total_min=999999
                        total_sum = 0
                        total_sum = s_value - t_value + total_min
                        if total_sum > 0:
                            eliminated_rows.append(absolute_row_index)
                            temp_sub_height[block_row_i] -= 1
                            eliminated = True
                            break
            pairs_matrices = []
            mask_matrices = []
        if eliminated:
            sym_matrix = delete_stuff(sym_matrix, eliminated_rows)
            pairs_binary = delete_stuff(pairs_binary, eliminated_rows)
            eliminated_rows = []
            current_sub_height = temp_sub_height.copy()
            print current_sub_height

    return (sym_matrix, current_sub_height, eliminated, pairs_binary, height_now)


def get_height_until_j(current_sub_height, j):
    if j == 0:
        return 0
    else:
        total_height = 0
        for i in xrange(j):
            total_height+= current_sub_height[i]
        return total_height

def sum_pairs(sym_matrix, block_row_i, block_row_j, current_sub_height, s, t, u, v, positions):
    current_interaction = 0
    total_min = 0 #begin pairs summation
    for block_column_k in xrange(positions):
        if block_column_k != block_row_i and block_column_k != block_row_j:
            pair_interactions = zeros(current_sub_height[block_column_k])
            if block_row_i > block_column_k:
                pair_block_i_k = get_pairs_block(sym_matrix, block_column_k, block_row_i, current_sub_height)
            else:
                pair_block_i_k = get_pairs_block(sym_matrix, block_row_i, block_column_k, current_sub_height)
            if block_row_j > block_column_k:
                pair_block_j_k = get_pairs_block(sym_matrix, block_column_k, block_row_j, current_sub_height)
            else:
                pair_block_j_k = get_pairs_block(sym_matrix, block_row_j, block_column_k, current_sub_height)
            for w in xrange(current_sub_height[block_column_k]):
                if block_row_i > block_column_k:
                    current_interaction = pair_block_i_k[w , s] - pair_block_i_k[w, u]
                else:
                    current_interaction = pair_block_i_k[s , w] - pair_block_i_k[u, w]
                if block_row_j > block_column_k:
                    current_interaction += pair_block_j_k[w, t] - pair_block_j_k[w, v]
                else:
                    current_interaction += pair_block_j_k[t, w] - pair_block_j_k[v, w]
                pair_interactions[w] = current_interaction
            pairs_min = pair_interactions.min()
            total_min += pairs_min
    return total_min


def create_uv_values(positions, current_sub_height):
    u_list = []
    v_list = []
    for position in current_sub_height:
        u_list.append(arange(position))
        v_list.append(arange(position))
    for i in xrange(positions):
        random.shuffle(u_list[i])
        random.shuffle(v_list[i])
        u_list[i] = list(u_list[i])
        v_list[i] = list(v_list[i])
    return (u_list, v_list)

def eliminate_doubles(sym_matrix, positions, current_sub_height, eliminated, pairs_binary):
    print 'perform pairs'
    absolute_row_index = 0
    absolute_column_index = 0
    eliminated = False
    total_columns_eliminated = 0
    break_out_of_u = False
    
    for block_row_i in range(positions):
        diag_block_i = get_digonal_block(sym_matrix, block_row_i , current_sub_height)
        height_until_i = get_height_until_j(current_sub_height, block_row_i)
        for block_row_j in range(positions):
            if block_row_j > block_row_i:#reciprocal pair as i-j = j-i
                diag_block_j = get_digonal_block(sym_matrix, block_row_j , current_sub_height)
                if block_row_i > block_row_j:
                    pair_block_i_j = get_pairs_block(sym_matrix, block_row_j, block_row_i, current_sub_height)
                else:# i < j
                    pair_block_i_j = get_pairs_block(sym_matrix, block_row_i, block_row_j, current_sub_height)
                height_until_j = get_height_until_j(current_sub_height, block_row_j)
                for s in range(current_sub_height[block_row_i]):
                    if total_columns_eliminated > 20000:
                        break
                    absolute_row_index  = height_until_i + s
                    s_value = diag_block_i[s, s]
                    for t in range(current_sub_height[block_row_j]):
                        absolute_column_index = height_until_j + t
                        if pairs_binary[absolute_row_index][absolute_column_index] == 1:
                            continue
                        t_value = diag_block_j[t, t]
                        if block_row_i > block_row_j:
                            s_t_pair = pair_block_i_j[t, s]
                        else:
                            s_t_pair = pair_block_i_j[s, t]
                        for u in range(current_sub_height[block_row_i]):
                            if u != s:
                                u_value = diag_block_i[u, u]
                                for v in range(current_sub_height[block_row_j]):
                                    if v != t:
                                        v_value = diag_block_j[v,v]
                                        if block_row_i > block_row_j:
                                            u_v_pair = pair_block_i_j[v, u]
                                        else:
                                            u_v_pair = pair_block_i_j[u, v]
                                        total_min = sum_pairs(sym_matrix, block_row_i, block_row_j, current_sub_height, s, t, u, v, positions)
                                        total_energy = s_value + t_value + s_t_pair - u_value - v_value - u_v_pair + total_min
                                        if total_energy > 0:
                                            pairs_binary[absolute_row_index][absolute_column_index] = 1
                                            eliminated= True
                                            break_out_of_u = True
                                            total_columns_eliminated += 1
                                            break
                            if break_out_of_u: 
                                break_out_of_u = False
                                break
    pairs_binary = triu(pairs_binary)
    pairs_binary = symmetrize(pairs_binary)
    print "eliminated %d pairs" % total_columns_eliminated
    total_columns_eliminated = 0
    return (eliminated, pairs_binary)

def evaluate_energy(sym_matrix, indices, current_sub_height, positions, pairs_combinations):
    total_energy = 0
    for i in range(positions):
        total_energy += get_digonal_block(sym_matrix, i, current_sub_height)[indices[i]][indices[i]]
    for combination in pairs_combinations:
        total_energy += get_pairs_block(sym_matrix, combination[0], combination[1], current_sub_height)\
        [indices[combination[0]]][indices[combination[1]]]
    return total_energy

def enumerate_remaining(sym_matrix, current_sub_height, positions):
    indices = zeros((positions),dtype=np.int)
    some_other_index = 0
    pairs_combinations = combinations(arange((positions),dtype=np.int),2)
    best_energy = 9999
    best_index = []
    best_absolute_loc = []
    abs_index = 0
    while True:
        current_energy = evaluate_energy(sym_matrix, indices, current_sub_height, positions, pairs_combinations)
        if current_energy < best_energy:
            best_energy = current_energy
            best_index = indices.copy()
        for position in range(positions):
            indices[position] += 1
            if indices[position] < current_sub_height[position]:
                break
            indices[position] = 0
            some_other_index = position
        if some_other_index == positions-1:
            break
    for i in range(positions):
        height_until_j = get_height_until_j(current_sub_height, i)
        for j in range(current_sub_height[i]):
            abs_index = height_until_j+j
            if best_index[i] == j:
                best_absolute_loc.append(abs_index)
    
    print "energy of the final conformation is %f" % best_energy
    return sym_matrix[best_absolute_loc][:,best_absolute_loc]
    
if __name__ == '__main__':

    # pair_small = zeros((405,405))
    # f = open('./pair_small.dat')
    # positions = 3

    pair_small = zeros((1530,1530))
    positions = 5
    try:
        f = open('./pair.dat')
    except:
        print "could not find the pair.dat file"

    pair_small = mlab.load(f, pair_small)
    f.close()
    print 'loaded'

    pairs_binary = zeros(pair_small.shape, dtype=int16)
    sym_matrix = symmetrize(pair_small)
    sub_matrix_length =  sym_matrix[:,0].size / positions
    current_sub_height = zeros((positions), dtype=int16)
    current_sub_height.fill(sub_matrix_length)
    temp_sub_height = current_sub_height.copy()
    eliminated = True
    small_u_list = []
    small_v_list = []
    doubles_loops = 0
    height_now = []
    
    while eliminated:
        time1 = time()
        sym_matrix, current_sub_height, eliminated, pairs_binary, height_now = \
            eliminate_singles(sym_matrix, positions, current_sub_height, height_now, eliminated, pairs_binary)
        time2 = time()
        total_time = time2-time1
#        print 'time for singles is %f' % total_time
        u_list, v_list = create_uv_values(positions, current_sub_height)
        
        time1 = time()
        eliminated, pairs_binary = eliminate_doubles(sym_matrix, positions, current_sub_height, eliminated, pairs_binary)
        time2 = time()
        total_time = time2 - time1
#        print 'time for doubles is %f' % total_time
#    print sym_matrix.shape
    print "Enumerating remaining fragments..."
    sym_matrix = enumerate_remaining(sym_matrix, current_sub_height, positions)
    print "Done."
    print "Final Symmetric Matrix"
    print sym_matrix
    
    #save the final matrix
    f = open('final_matrix.csv', 'w')
    mlab.save(f, sym_matrix)
    f.close()
    
    #plot stuff
    x = arange(len(height_now))
    transpose = [list(c) for c in zip(*height_now)]
    for ind, i in enumerate(transpose):
        plt.plot(x,i,linewidth=2.0,label="position %d" % ind)
        plt.legend()
    plt.title('Large Matrix')
    plt.xlabel('Singles Iterations')
    plt.ylabel('Rotamers Left')
    plt.show()