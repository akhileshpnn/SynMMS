import numpy as np
from ReactionDiffusion1D_main import ReactionDiffusion1D as rd1d


class TimeVaryingPheromone(): 
        
     
    def constant_pheromone(self,X,const):
        Pher = np.ones((len(X),))*const
        return Pher

    def gaussian_pheromone(self,X,mean, sigma, const):
        gussian_ampl=0.5
        Pher = gussian_ampl*np.exp(-0.5*((X-mean)/sigma)**2) + const
        return Pher
    
    
    def add_pheromone_mt(self,X, t):

        if 0<=t<=rd1d.t_total:
            Pher_mt = self.constant_pheromone(X,1)
        return Pher_mt
    
    def add_pheromone_aur(self,X, t):
        
        if len(rd1d.f)==0: # initial state
            Pher_aur = self.constant_pheromone(X,1)
            
        elif len(rd1d.f)==2:  # 1 stimulus          
            [f1,f2]=rd1d.f       
            t_beg1=rd1d.t_total*f1;t_end1=rd1d.t_total*f2;            
            if 0<=t<=t_beg1: 
                Pher_aur = self.constant_pheromone(X,1)
            elif  t_beg1<=t<= t_end1:
                Pher_aur = self.gaussian_pheromone(X,65,5,1)
                
        elif len(rd1d.f)==4: # 2 stimuli
            
            [f1,f2,f3,f4]=rd1d.f
        
            t_beg1=rd1d.t_total*f1;t_end1=rd1d.t_total*f2;
            t_beg2=rd1d.t_total*f3;t_end2=rd1d.t_total*f4;
            
            if 0<=t<=t_beg1:
                Pher_aur = self.constant_pheromone(1)
            elif  t_beg1<=t<= t_end1:
                Pher_aur = self.gaussian_pheromone(48,3.5,1)
            elif t_beg2 <t<= t_end2:
                Pher_aur = self.gaussian_pheromone(70,3.5,1)

        return Pher_aur
            
