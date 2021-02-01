import numpy as np
from reaction_diffusion1D_main import ReactionDiffusion1D as rd1d


class TimeVaryingCue(): 
             
    def constant_cue(self,X,const):
        Light = np.ones((len(X),))*const
        return Light

    def gaussian_cue(self,X,mean, sigma, const):
        gussian_ampl=0.5
        Light = gussian_ampl*np.exp(-0.5*((X-mean)/sigma)**2) + const
        return Light

    
    def add_light_aurora(self,X, t):
        
        if len(rd1d.f)==0: # no light cue genrating initial state
            Light_aur = self.constant_cue(X,1)
            
        elif len(rd1d.f)==2:  # single stimulus case         
            [f1,f2]=rd1d.f       
            t_beg1=rd1d.N*f1;t_end1=rd1d.N*f2;            
            if 0<=t<=t_beg1: 
                Light_aur = self.constant_cue(X,1)
            elif  t_beg1<=t<= t_end1:
                Light_aur = self.gaussian_cue(X,65,5,1)
                
        elif len(rd1d.f)==4: # multi-stimulus case, here with 2 signals
            
            [f1,f2,f3,f4]=rd1d.f
        
            t_beg1=rd1d.N*f1;t_end1=rd1d.N*f2;
            t_beg2=rd1d.N*f3;t_end2=rd1d.N*f4;
            
            if 0<=t<=t_beg1:
                Light_aur = self.constant_cue(1)
            elif  t_beg1<=t<= t_end1:
                Light_aur = self.gaussian_cue(48,3.5,1)
            elif t_beg2 <t<= t_end2:
                Light_aur = self.gaussian_cue(70,3.5,1)

        return Light_aur
            
