
from typing import Type
import numpy as np
 
import pyqubo


class Problem3R3X():
    """Generates planted solution instance of 3regular3xorsat
    Each variable belong to 3 equations and each equation have exactly 3 variables
    so XOR(x_i1, x_i2, x_i3) = b_i 

    Able to convert it to Ising model coefficients
    """
    def __init__(self, n=10):
        self.eqn_num=n
        self.solution=np.zeros(self.eqn_num, dtype=bool)
        self.variables=np.ndarray(shape=(3,self.eqn_num), dtype=int)
        self.variables_t=self.variables.transpose()
        self.b=np.zeros(self.eqn_num, dtype=bool)
        self.spins_=pyqubo.Array.create('s',shape=(1,self.eqn_num),vartype='SPIN')
        self.model_=None
        
    def generate(self) -> None: #generate random instance of 3R3X
        for i in range (0,3):
            self.variables[i]=np.random.permutation(self.eqn_num)
        self.variables_t=self.variables.transpose()
        self.solution=np.random.randint(0,2, size=self.eqn_num, dtype=bool)
        x=self.solution[[self.variables_t]]
        self.b=np.remainder(np.sum(x,axis=1),2)
        self.__makeHam()
        return 
    
    def __G3X(self, num: int, b: bool) : # declare  class 'cpp_pyqubo.Base'
        ind=self.variables_t[num]
        xksy=self.solution[ind]
        x1,x2,x3=xksy[0],xksy[1],xksy[2]
        s1=pyqubo.Spin('s'+str(ind[0]))
        s2=pyqubo.Spin('s'+str(ind[1]))
        s3=pyqubo.Spin('s'+str(ind[2]))
        sa=pyqubo.Spin('s'+str(self.eqn_num+num))
        if(b==0):
            if(x1==0 and x2==0 and x3==0):
                (h,hbar,J,Jbar)=(-1,-2,1,2)
            else: 
                (h,hbar,J,Jbar)=(-1,2,1,-2)
        else: 
            if(x1==1 and x2==1 and x3==1):
                (h,hbar,J,Jbar)=(-1,2,-2,-1)
            else:
                (h,hbar,J,Jbar)=(2,1,-1,-2)
        H=h*(s1+s2+s3) +hbar*sa +J*(s1*s2+s2*s3+s3*s1)+Jbar*sa*(s1+s2+s3)
        return H
            
    def __makeHam(self) -> None:
        H=0
        for i in range(0,self.eqn_num):
            H+=self.__G3X(i,self.b[i])
        self.model_=H.compile()
        return  
            
    def to_ising(self):
        return self.model_.to_ising()    
    def to_qubo(self):
        return self.model_.to_qubo()
    def to_bqm(self):
        return self.model_.to_bqm()
    
    def get_solution(self):
        return self.solution

    def is_solution(self,arr: np.array):
        tmp=arr[[self.variables_t]]
        calculated_b=np.remainder(np.sum(tmp,axis=1),2)
        print(self.b-calculated_b)


def convert_solution_to_array(n, sol):
        """Convert solution from dict to array"""
        res = []
        for i in range(1, n+1):
            res.append(sol['s'+str(i)])
        return np.array(res) 
   
def check(n):
    import neal
    x=Problem3R3X(n)
    x.generate()
    model=x.model_
    bqm=x.to_bqm()
    sa = neal.SimulatedAnnealingSampler()
    sampleset = sa.sample(bqm, num_reads=10)
    samples = model.decode_sampleset(sampleset)
    best_sample = min(samples, key=lambda s: s.energy)
    sol_calc=convert_solution_to_array(n, best_sample.sample)
    print(sol_calc) 
    print(x.get_solution())
    x.is_solution(sol_calc)

check(3)