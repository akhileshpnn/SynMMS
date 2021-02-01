import numpy as np

class InitialConditions:
    
    def __init__(self, model, size):
        
        self.model=model
        self.size=size
    
    def around_steady_state(self):
        
#        seed_int = 33#np.random.randint(0,500)
#        print(seed_int)
#        np.random.seed(seed_int)
        
        from scipy.integrate import odeint
        total_t = 500;dt=0.01
        rn1=np.random.rand()/2; rn2=np.random.rand()/2
        y0 = [self.model.c2*rn1,self.model.c2*rn2,self.model.c2*(1-(rn1+rn2))]
        t=np.arange(0,total_t,dt)
        sol = odeint(self.model.reaction, y0, t, args=(self.model.rho,))        
        umss,vmss,vcss=sol[-1]
        
        R1 = np.random.rand(self.size,self.size);
        R2 = np.random.rand(self.size,self.size);
        R3 = np.random.rand(self.size,self.size);
        Up=umss*R1
        Vfm=vmss*R2
        Vfc=vcss*R3
        
#        Sfc=self.model.sand_t-(Sp+Sfm)
        
        return Up,Vfm,Vfc
    
    def random(self):
        
#        seed_int = np.random.randint(0,500)
#        print(seed_int)
#        np.random.seed(seed_int)        
        
#        Sp = self.model.sand_t/(self.size*self.size)*np.random.rand(self.size,self.size)/2
#        Sfm = self.model.sand_t/(self.size*self.size)*np.random.rand(self.size,self.size)/2
#        Sfc = self.model.sand_t/(self.size*self.size) - (Sp+Sfm)
        
        Up = self.model.c2*np.random.rand(self.size,self.size)/2
        Vfm = self.model.c2*np.random.rand(self.size,self.size)/2
        Vfc = self.model.c2 - (Up+Vfm)
        
        return Up,Vfm,Vfc
        
