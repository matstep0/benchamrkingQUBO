#TODOLIST
#seed do losowania
#funckja extract problem która wczytuje problem do pliku
#funckja read problem która wczytuje problem z pliku
# przepisać na sparce matrix


#!/usr/bin/env python
import numpy as np
import random
import pyqubo
import neal
from scipy import sparse

class  Problem_Generator:
    """This class generate problem desinged to be solved by QUBO.
    Square matrix A and vector b and binary vector x can be generated, 
    so Ax=b where x is binary vector. 
    Random A traingle matrix with given density d (0,1) and random binary array x
    are generated and than b is calculated as b=Ax
    Than  vector x is exact solution and sum_i (A_ij x_j -b_j)^2 is then minimized.
    
    Args: 
        size: size of square matrix 
        density: density of matrix A that will be sampled
    """
    def __init__(self,size,density=1):
        self.size=size
        self.density=density
        self.A=None      
        self.b=None
        self.x=None    
        self.H=0
        self.gen_sol=None
        self.qubo=None
        self.offset=None
    def __gen_A(self):
        "Generate random instance of traingle matrix with given density"
        self.A=sparse.random(self.size,self.size,density=self.density)
        self.A=sparse.tril(self.A)
        return self.A

    
    def __gen_x(self): 
        """Generate solution for problem"""
        return sparse.csr_matrix(np.random.randint(2,size=(self.size,1)))
    
    def __gen_b(self):
        """Calculating b based on generated A and x, b=Ax"""
        return self.A.multiply(self.x)    

    def get_A(self):
        return self.A
    
    def get_x(self):
        return self.x
    
    def get_b(self):
        return self.b
    
    def generate_problem(self):
        """Generate problem. A set of equation written in matrix form Ax=b
        
        Return: 
           ::class:: 'numpy.matrix' """
           
        self.A=self.__gen_A()
        self.x=self.__gen_x()
        self.b=self.__gen_b()
        return self.A,self.x,self.b

    def Binary_Hamiltonian(self):
        """Generate Hamiltonian in with binary varables.
        Return:
            <class 'cpp_pyqubo.Add'>"""
        from pyqubo import Binary
        self.binar=np.array([Binary("x"+str(i)) for i in range(1,self.size+1)] )
        self.H=0
        for line, bi in zip(self.A, self.b):
            self.H+=(np.sum(np.multiply(line,self.binar))-bi)**2
            #print(line)
            #print(self.binar)
            #print(np.multiply(line,self.binar))
        return self.H
    
    def compile(self):
        """Compile Hamiltonian"""
        self.H=self.H.compile()
        return self.H      
    
    def to_qubo(self):
        """Calculate and return QUBO model coefficiencts and offset"""
        self.qubo , self.offset= self.H.to_qubo()
        return self.qubo, self.offset
    def anneal(self):
        """Simulate algorith of quantum annealing from deave-neal package 
        Return:
            Python dictionary as solution, 
            use method solution to extract it to binary array"""
        sampler = neal.SimulatedAnnealingSampler()
        bqm = self.H.to_bqm()
        sampleset = sampler.sample(bqm, num_reads=1)
        decoded_samples = self.H.decode_sampleset(sampleset)
        best_sample = min(decoded_samples, key=lambda x: x.energy)
        self.gen_sol=best_sample.sample 
        return self.gen_sol
    
    def solution(self):
        """Gives calculated solution as numpy array
        Return:
            class:numpy.array"""
        res=[]
        for i in range(1,self.size+1):
            res.append(self.gen_sol['x'+str(i)])
        return np.array(res) 
    
    def exact_solution(self):
        """Give solution as 1-dim numpy array"""
        return np.ndarray.flatten(self.x)
    def cost(self,t):
        """Claculate square error of given solution.
        solution must by in numpy array form"""
        #print(self.get_A())
        #print(t.transpose())
        #print(np.matmul(self.get_A(),t.transpose()))
        #print(self.get_b())
        #print(np.matmul(self.get_A(),t.transpose()) - self.get_b().flatten() ) 
        return np.sum((np.matmul(self.get_A(),t.transpose()) - self.get_b().flatten())**2)

"""Code for checking performance of code"""
import time

def calculate_execution_time(res,sizes,densities):
    for s in sizes:
        for den in densities:
            start = time.process_time()
            Problem_Generator(s,den).generate_problem()    
            end=time.process_time()
            res[s,den]=round(end-start,4)
    return res
def pretty_print(dic,sizes,densities):
    print("x    ",end=' ')
    for den in densities:
        print(den,end='    ')
    print()
    for s in sizes:
        print(s,end=' ')
        for den in densities:
            print(dic[s,den],end=' ')
        print()
def save_result():
    print("TU bedize funkcja")

sizes=[10,20,30,40,50]
densities=[ 0.2,0.4,0.6,0.8,1]
dic={}
dic=calculate_execution_time(dic,sizes,densities)
pretty_print(dic,sizes,densities)


"""Execution code validity for testing"""
"""
x=Problem_Generator(4,0.8) #create class parameters size-20 density=0.8
x.generate_problem()       #sample problem Ax=b
print(x.get_A(),x.get_x(),x.get_b())
x.Binary_Hamiltonian()     #create hamiltonianian for generated problem
x.compile()             #compile model
qubo, offset = x.to_qubo() #get QUBO problem 
print(qubo)
planted_solution=x.exact_solution()  #get solution as numpy array
x.anneal()             # run simulated annealing to obtain solution 
generated_solution=x.solution() #get this calculated solution as numpy array
print(planted_solution) 
print(generated_solution)
error_planted=x.cost(planted_solution)         #calculate square error
error_generated=x.cost(generated_solution)
print(error_planted)
print(error_generated)
"""

"""
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

"""