import numpy as np

class InitialConditions:
    
    def __init__(self, model, size):
        
        self.model=model
        self.size=size
    
    def around_steady_state(self):
        
#        seed_int = np.random.randint(0,500)
        seed_int = 10000
        np.random.seed(seed_int)
        print(seed_int)

        
        from scipy.integrate import odeint
        total_t = 500;dt=0.01
        rn=np.random.rand()
        y0 = [self.model.sand_t*rn,self.model.sand_t*(1-rn)]
        t=np.arange(0,total_t,dt)
        sol = odeint(self.model.reaction, y0, t, args=(self.model.pheromone_ampl,))        
        uss,vss=sol[-1]
        
        R1 = np.random.rand(self.size,self.size);
        R2 = np.random.rand(self.size,self.size);
        Sp=uss*R1
        Sf=vss*R2#(1-R)
        
        return Sp,Sf
    
    def random(self):
        
        seed_int = 10000#np.random.randint(0,500)
        print(seed_int)
        np.random.seed(seed_int)        
        
        Sp = self.model.sand_t*np.random.rand(self.size,self.size)
        Sf = self.model.sand_t - Sp
        
        return Sp,Sf
        
