import numpy as np


class SynMMS():
    
    #c1 = total amount of MTs; gamma1 = SIC; c2 =total amount of Aurora; gamma2 = CC; k3= causal link from Aurora -> MT; k4= causal link from MT -> Aurora 

    c1=1300; gamma1 = 1; omega1 = .5; k0 = 0.067; c2 = 0; gamma2 = 0; omega2 = 0; k3=0; k4=0 # Fig. 1g
    c1=1100; gamma1 = 1; omega1 = .5; k0 = 0.067; c2 = 0; gamma2 = 0; omega2 = 0; k3=0; k4=0 # Supplementary Fig. 1h
    c1=2700; gamma1 = 1; omega1 = .5; k0 = 0.067; c2 = 0; gamma2 = 0; omega2 = 0; k3=0; k4=0 # Supplementary Fig. 1i

    c1=1000; gamma1 = 1; omega1 = .5; k0 = 0.067; c2 = 2750; gamma2 = 0.5; omega2 = 5; k3=1.0; k4=0.5 # Fig. 8b, Supplementary Fig. 9d,f (SIC>CC)
    c1=1000; gamma1 = 1; omega1 = .5; k0 = 0.067; c2 = 2750; gamma2 = 5.0; omega2 = 5; k3=1.0; k4=0.5 # Fig. 8c, Supplementary Fig. 9e  (SIC<CC)
    c1=1000; gamma1 = 1; omega1 = .5; k0 = 0.067; c2 = 2750; gamma2 = 0.5; omega2 = 5; k3=1.0; k4=10 # Supplementary Fig. 9g (SIC>CC)

    Du1 = .1; Dv1 = 40; # diffusion contants for MTs in protrusions and free MTs
    Du2 = .1; Dv2 = 40; # diffusion contants for Aurora in clusters and free Aurora monomer
    
    def reaction(self, u1, v1, u2, v2, rho):
        
        fu1v1 =(self.k0 + self.k3*u2+ self.gamma1*u1*u1)*v1 - self.omega1*u1
        gu1v1=-fu1v1
        
        fu2v2 =(self.k0 + self.k4*u1+ rho*self.gamma2*u2*u2)*v2 - self.omega2*u2
        gu2v2=-fu2v2

        return fu1v1, gu1v1, fu2v2, gu2v2
