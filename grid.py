class Grid:
    def __init__(
        self, nx, ny, nz, d
    ):  # number of cells centers in x, y, and z-direction, side-length of cubic cells d in m
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.d = d
