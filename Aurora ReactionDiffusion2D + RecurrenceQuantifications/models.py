class AuroraSubSystem:
    
    sand_t=5;pheromone_ampl = 1;gamma = 10;omega = 5;k0 = 0.067 
    Dp=1;Df=50;
    
    def reaction(self, y, t, pher_aurora):
        u,v = y    
        fuv =(self.k0 + pher_aurora*self.gamma*u*u)*v - self.omega*u
        guv=-fuv
        return [fuv, guv]
