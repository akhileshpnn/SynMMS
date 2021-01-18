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
        y0 = [self.model.sand_t*rn1,self.model.sand_t*rn2,self.model.sand_t*(1-(rn1+rn2))]
        t=np.arange(0,total_t,dt)
        sol = odeint(self.model.reaction, y0, t, args=(self.model.pheromone_ampl,))        
        umss,vmss,vcss=sol[-1]
        
        R1 = np.random.rand(self.size,self.size);
        R2 = np.random.rand(self.size,self.size);
        R3 = np.random.rand(self.size,self.size);
        Sp=umss*R1
        Sfm=vmss*R2
        Sfc=vcss*R3
        
#        Sfc=self.model.sand_t-(Sp+Sfm)
        
        return Sp,Sfm,Sfc
    
    def random(self):
        
#        seed_int = np.random.randint(0,500)
#        print(seed_int)
#        np.random.seed(seed_int)        
        
#        Sp = self.model.sand_t/(self.size*self.size)*np.random.rand(self.size,self.size)/2
#        Sfm = self.model.sand_t/(self.size*self.size)*np.random.rand(self.size,self.size)/2
#        Sfc = self.model.sand_t/(self.size*self.size) - (Sp+Sfm)
        
        Sp = self.model.sand_t*np.random.rand(self.size,self.size)/2
        Sfm = self.model.sand_t*np.random.rand(self.size,self.size)/2
        Sfc = self.model.sand_t - (Sp+Sfm)
        
        return Sp,Sfm,Sfc
        
