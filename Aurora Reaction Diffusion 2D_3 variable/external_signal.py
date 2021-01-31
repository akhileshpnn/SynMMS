import numpy as np

class Uniform:
    
    def add_pheremone(self, ampl, size):
        Pher = ampl*np.ones((size,size))
        return Pher
