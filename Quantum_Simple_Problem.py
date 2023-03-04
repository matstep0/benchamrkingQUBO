#!/usr/bin/env python
import numpy as np
import random
import pyqubo
import neal


class  Problem_Generator:
    """This class generate problem desinged to be solved by QUBO.
    Matrix A and vector b and binary vector x can be generated, so Ax=b where x is binary vector. 
    Than  vector x is exact solution and sum_i (A_ij x_j -b_j)^2 is then minimized.
    """
    def __init__(self,size,density):
        self.size=size
        self.density=density
        self.A=None      
        self.b=None
        self.x=None    
        self.H=0
        self.gen_sol=None
        
    def __gen_A(self):
        "Generate random instance of traingle matrix with given density"
        temp_mat = np.tril(np.random.rand(self.size,self.size))
        x,y = np.nonzero(temp_mat)
        ind = list(zip(x,y))
        ind = random.choices(ind, k = int((1-self.density)*len(ind))) 
        for (i,j) in ind:
            temp_mat[i,j] = 0   #poprawić
        return temp_mat
    
    def __gen_x(self): 
        #Getting exact solution to minimizing problem
        return np.random.randint(2,size=(self.size,1))
    
    def __gen_b(self):
        #Gettin b based on generated A and x b=Ax
        return np.matmul(self.A,self.x)    

    def get_A(self):
        return self.A
    
    def get_x(self):
        return self.x
    
    def get_b(self):
        return self.b
    
    def generate_problem(self):
        self.A=self.__gen_A()
        self.x=self.__gen_x()
        self.b=self.__gen_b()
        return

    def Ising_Hamiltonian(self):
        """Generate and return Hamiltonian"""
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
        self.H=self.H.compile()
        return        
    
    def to_qubo(self):
       return self.H.to_qubo()
    
    def anneal(self):
        sampler = neal.SimulatedAnnealingSampler()
        bqm = self.H.to_bqm()
        sampleset = sampler.sample(bqm, num_reads=10)
        decoded_samples = self.H.decode_sampleset(sampleset)
        best_sample = min(decoded_samples, key=lambda x: x.energy)
        self.gen_sol=best_sample.sample 
        return self.gen_sol
    
    def solution(self):
        """Gives generated solution as list"""
        print(type(self.gen_sol))
        res=[]
        for i in range(1,self.size+1):
            res.append(self.gen_sol['x'+str(i)])
        return np.array(res) 
    
    def real_solution(self):
        """Give solution as 1-dim numpy array"""
        return np.ndarray.flatten(self.x)
    
 

"""Execution code for testing"""

x=Problem_Generator(3,0.4)
x.generate_problem()
x.Ising_Hamiltonian() ##H wychodzi jako tablica z jednym elementem...
x.compile()
qubo, offset = x.to_qubo()
print(qubo)
planted_solution=x.real_solution()
x.anneal()
generated_solution=x.solution()
print(planted_solution)
print(generated_solution)


error_planted=np.sum((np.matmul(x.get_A(),x.get_x()) - x.get_b())**2)    #square error
error_generated=np.sum((np.matmul(x.get_A(),generated_solution.transpose()) - x.get_b())**2)

print(error_planted)
print(error_generated)

#qubo, offset = model.to_qubo()
#print(qubo)


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

