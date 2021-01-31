class AuroraSubSystem3var:
      
    c2=8;rho=1;gamma2 = 5;omega2 = 5;k02 = 0.067
    kon=0.1;koff=0.2
    Du2=1;Dv2=1;Dw2=50;
    
    def reaction(self, y, t, rho):
        u2,v2,w2 = y    
        f =(self.k02 + rho*self.gamma2*u2*u2)*v2 - self.omega2*u2
        g=-f+self.kon*w2-self.koff*v2
        h=-self.kon*w2+self.koff*v2
        return [f, g, h]
