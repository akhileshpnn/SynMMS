class PeriodicBoundaryConditions:
    
    def update_boundaries(self, Z):
        for U in Z:
            U[0, :] = U[-1, :]
            U[:, 0] = U[:, -1]


class NeumannBoundaryConditions:
    
    def update_boundaries(self, Z):
        for U in Z:
            U[0, :] = U[1, :]
            U[-1, :] = U[-2, :]
            U[:, 0] = U[:, 1]
            U[:, -1] = U[:, -2]
