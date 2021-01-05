import numpy as np

class Model:
    def print_model_name(self):
        print(type(self))
        
    def reaction(self, u, v,u2,v2, pher):
        return


class SynMMS(Model):
    c1 = 1000 #Mt
    c2 = 2750 # Aur

    gamma1 = 2; omega1 = .5; k01 = 0.067 #MT
    gamma2 = 0.5; omega2 = 5; k02 = 0.067 # Aurora
    
    k3=1.0 # Aur->MT
    k4=0.5 # MT -> Aur

#    k4=10 # MT -> Aur Supp Fig. 9g   
#    k3=0;k4=0 # MT subsystem

    du1 = .1; dv1 = 40;
    du2 = .1; dv2 = 40;

    
    def reaction(self, u1, v1, u2, v2, pher_mt, pher_aurora):
        
        fu1v1 =(self.k01 + self.k3*u2+ pher_mt*self.gamma1*u1*u1)*v1 - self.omega1*u1
        gu1v1=-fu1v1
        fu2v2 =(self.k02 + self.k4*u1+ pher_aurora*self.gamma2*u2*u2)*v2 - self.omega2*u2
        gu2v2=-fu2v2

        return fu1v1, gu1v1, fu2v2, gu2v2
