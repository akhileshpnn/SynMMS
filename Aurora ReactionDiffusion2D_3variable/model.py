class AuroraSubSystem3var:
      
    sand_t=8; pheromone_ampl = 1;gamma = 5;omega = 5;k0 = 0.1
    kon=0.1;koff=0.2;eta=0.01
    Dm1=1;Dm2=1;Dc=50;
    
    def reaction(self, y, t, pher_aurora):
        um,vm,vc = y    
        fuv =(self.k0 + pher_aurora*self.gamma*um*um)*vm - self.omega*um
        guv=-fuv+self.kon*vc-self.koff*vm
        huv=-self.kon*vc+self.koff*vm
        return [fuv, guv, huv]
