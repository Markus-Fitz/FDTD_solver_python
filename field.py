import numpy as np

class Field:
    def __init__(self, grid):
        # Number of cubic cells in x-, y- and z-directions
        nx = grid.nx
        ny = grid.ny
        nz = grid.nz

        # global cell-size in x-, y- and z-direction
        self.dx = grid.d
        self.dy = grid.d
        self.dz = grid.d

        # number of cells in x-, y- and z-direction; TODO: this does not belong here, modify logic so that it is not needed
        self.nx = grid.nx
        self.ny = grid.ny
        self.nz = grid.nz

        # electric fields; located on cube edges
        # Ex oriented in x+ direction -> one more value than cell count in y- and z-direction
        self.Ex = np.zeros((nx, ny+1, nz+1))
        # Ex export field for colocated visualization
        self.Ex_exp = np.zeros((nx, ny, nz))
        # Ey oriented in y+ direction -> one more value than cell count in x- and z-direction
        self.Ey = np.zeros((nx+1, ny, nz+1))
        # Ey export field for colocated visualization
        self.Ey_exp = np.zeros((nx, ny, nz))
        # Ez oriented in z+ direction -> one more value than cell count in x- and y-direction
        self.Ez = np.zeros((nx+1, ny+1, nz))
        # Ez export field for colocated visualization
        self.Ez_exp = np.zeros((nx, ny, nz))

        # magnetic fields; located on cube faces
        # Bx oriented in x+ direction -> one more value than cell count in x-direction
        self.Bx = np.zeros((nx+1, ny, nz))
        # Bx export field for colocated visualization
        self.Bx_exp = np.zeros((nx, ny, nz))
        # By oriented in y+ direction -> one more value than cell count in y-direction
        self.By = np.zeros((nx, ny+1, nz))
        # By export field for colocated visualization
        self.By_exp = np.zeros((nx, ny, nz))
        # Bz oriented in z+ direction -> one more value than cell count in z-direction
        self.Bz = np.zeros((nx, ny, nz+1))
        # Bz export field for colocated visualization
        self.Bz_exp = np.zeros((nx, ny, nz))

        # TODO: material parameters