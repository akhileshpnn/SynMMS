import numpy as np

class Uniform:
    
    def add_light(self, ampl, size):
        light = ampl*np.ones((size,size))
        return light
