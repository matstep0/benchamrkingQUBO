#!/usr/bin/env python
import numpy as np
import random

#Generate triangle matrix (n) with density parameter d (0,1).
n = 20 #size of array
d = 0.7

def random_tril_matrix(n, d):
    matrix = np.tril(numpy.random.rand(n,n))
    x,y = np.nonzero(matrix)
    ind = list(zip(x,y))
    ind = random.choices(ind, k = int((1-d)*len(ind))) 
    for (i,j) in ind:
        matrix[i,j] = 0
    return matrix


def random_qubo():

    # Generating model with parameters so QUBO problem with matrix Mij

    spiny=[Spin("s"+str(i)) for i in range(1, n+1)]
    H = 0
    for i in range(0, n):
        for j in range(0, n):
            H += matrix[i,j]*spiny[i]*spiny[j]
    model = H.compile()
    qubo, offset = model.to_qubo()

    return qubo, offset

def anneal():
    sampler = neal.SimulatedAnnealingSampler()
    bqm = model.to_bqm()
    sampleset = sampler.sample(bqm, num_reads=10)
    decoded_samples = model.decode_sampleset(sampleset)
    best_sample = min(decoded_samples, key=lambda x: x.energy)
    print(best_sample.sample) 

