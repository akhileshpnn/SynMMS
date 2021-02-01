import numpy as np

class InitialConditions:
    
    def __init__(self, model, size):
        
        self.model=model
        self.size=size
    
    def around_steady_state(self):
        
#        seed_int = np.random.randint(0,500)
        seed_int = 10000
        np.random.seed(seed_int)
        print('random number seed:'+str(seed_int))

        
        from scipy.integrate import odeint
        total_t = 500;dt=0.01
        rn=np.random.rand()
        y0 = [self.model.c2*rn,self.model.c2*(1-rn)]
        t=np.arange(0,total_t,dt)
        sol = odeint(self.model.reaction, y0, t, args=(self.model.rho,))        
        uss,vss=sol[-1]
        
        R1 = np.random.rand(self.size,self.size);
        R2 = np.random.rand(self.size,self.size);
        Up=uss*R1
        Vf=vss*R2#(1-R)
        
        return Up,Vf
    
    def random(self):
        
        seed_int = 10000#np.random.randint(0,500)
        print(seed_int)
        np.random.seed(seed_int)   
        print('random number seed:'+str(seed_int))
        
        Up = self.model.c2*np.random.rand(self.size,self.size)
        Vf = self.model.c2 - Up
        
        return Up,Vf
        
