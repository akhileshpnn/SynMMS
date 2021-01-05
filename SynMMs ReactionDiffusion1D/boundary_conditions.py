class BoundaryConditions:
    def update_boundaries(self):
        return


class PeriodicBoundaryConditions(BoundaryConditions):
    def update_boundaries(self, Z):
        for field in Z.dtype.names:
            Z[field][0] = Z[field][-2]
            Z[field][-1] = Z[field][1]


class ReflectingBoundaryConditions(BoundaryConditions):
    def update_boundaries(self, Z):
        for field in Z.dtype.names:
            Z[field][1] += Z[field][0]
            Z[field][-2] += Z[field][-1]
            Z[field][0] = 0
            Z[field][-1] = 0


class NeumannBoundaryConditions(BoundaryConditions):
    def update_boundaries(self, Z):
        for field in Z.dtype.names:
            Z[field][0] = Z[field][1]
            Z[field][-1] = Z[field][-2]

