import numpy as np

class Field:
    def __init__(self, grid):
        # Number of cubic cells in x-, y- and z-directions
        nx = grid.nx
        ny = grid.ny
        nz = grid.nz

        # electric fields; located on cube edges
        # Ex oriented in x+ direction -> one more value than cell count in y- and z-direction
        self.Ex = np.zeros(nx, ny+1, nz+1)
        # Ey oriented in y+ direction -> one more value than cell count in x- and z-direction
        self.Ey = np.zeros(nx+1, ny, nz+1)
        # Ez oriented in z+ direction -> one more value than cell count in x- and y-direction
        self.Ez = np.zeros(nx+1, ny+1, nz)

        # magnetic fields; located on cube faces
        # Bx oriented in x+ direction -> one more value than cell count in x-direction
        self.Bx = np.zeros(nx+1, ny, nz)
        # By oriented in y+ direction -> one more value than cell count in y-direction
        self.By = np.zeros(nx, ny+1, nz)
        # Bz oriented in z+ direction -> one more value than cell count in z-direction
        self.Bz = np.zeros(nx, ny, nz+1)

        # TODO: material parameters