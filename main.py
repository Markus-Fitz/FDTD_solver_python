try:
    import cupy as np
except ImportError:
    import numpy as np

# import grid and field classes
from grid import Grid
from field import Field

# define domain