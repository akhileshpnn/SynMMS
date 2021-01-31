class AuroraSubSystem:
    
    c2=5;rho = 1;gamma2 = 10;omega2 = 5;k02 = 0.067 
    Du2=1;Dv2=50;
    
    def reaction(self, y, t, rho):
        u2,v2 = y    
        fuv =(self.k02 + rho*self.gamma2*u2*u2)*v2 - self.omega2*u2
        guv=-fuv
        return [fuv, guv]
