try:
    import cupy as np
except ImportError:
    import numpy as np

# import grid and field classes
from grid import Grid
from field import Field
from vtk_export import export_data

# define domain
domain = Grid(51, 51, 51, 0.005) # domain size is 25cm x 25cm x 25cm
fields = Field(domain)

# define global parameters
#vacuum permittivity
eps0 = 8.8542*10e-12

#vacuum permeability
mu0 = 1.2566*10e-6

#end-time of simulation in seconds
tend = 0.00000001

#desired courant number, has to be lower than one
C = 0.1

#calculation of timestep from courant number
dt = C*domain.d*np.sqrt(eps0*mu0)

#timesteps
tsteps = int(np.floor(tend/dt))

print("dt is set to ", dt, ". The number of timesteps is ", tsteps)

print("periods of source simulated:", tend * 1000000000)

# definition of update functions
def update_Bx(fields):
    for i in range(fields.nx + 1):
        for j in range(fields.ny):
            for k in range(fields.nz):
                # values of Bx at position (i-1/2)*dx, j*dy, k*dz are evaluated;
                # calculation is shifted half an index in negative x-direction!
                B_x_new = fields.Bx[i, j, k] - dt*((fields.Ey[i, j, k+1] - fields.Ey[i, j, k])/fields.dz - (fields.Ez[i, j+1, k] - fields.Ez[i, j, k])/fields.dy)
                
                fields.Bx[i, j, k] = B_x_new

def update_By(fields):
    for i in range(fields.nx):
        for j in range(fields.ny+1):
            for k in range(fields.nz):
                # values at position i*dx, (j+1/2)*dy, k*dz are evaluated;
                # calculation is shifted half an index in negative y-direction!
                B_y_new = fields.By[i, j, k] - dt*((fields.Ez[i+1, j, k] - fields.Ez[i, j, k])/fields.dx - (fields.Ex[i, j, k+1] - fields.Ex[i, j, k])/fields.dz)

                fields.By[i, j, k] = B_y_new
		
def update_Bz(fields):
    for i in range(fields.nx):
        for j in range(fields.ny):
            for k in range(fields.nz+1):
                # NOTE: values at position i*dx, j*dy, (k+1/2)*dz are evaluated;
                # calculation is shifted half an index in negative z-direction!
                B_z_new = fields.Bz[i, j, k] - dt*((fields.Ex[i, j+1, k] - fields.Ex[i, j, k])/fields.dy - (fields.Ey[i+1, j, k] - fields.Ey[i, j, k])/fields.dx) 
                
                fields.Bz[i, j, k] = B_z_new
		
def update_Ex(fields):
    for i in range(fields.nx):
        for j in range(1, fields.ny):
            for k in range(1, fields.nz):
                # NOTE: Only the inner edges can be evaluated!
                # Calculation is shifted half an index in negative y- and z-dir!
                E_x_new = fields.Ex[i, j, k] + (dt/(mu0*eps0))*((fields.By[i, j, k] - fields.By[i, j, k-1])/fields.dz - (fields.Bz[i, j, k] - fields.Bz[i, j-1, k])/fields.dy)# - (dt/eps0)*J_x[i, j, k]

                fields.Ex[i, j, k] = E_x_new
		
def update_Ey(fields):
    for i in range(1, fields.nx):
        for j in range(fields.ny):
            for k in range(1, fields.nz):
                # NOTE: Only the inner edges can be evaluated!
                # Calculation is shifted half an index in negative x- and z-dir!
                E_y_new = fields.Ey[i, j, k] + (dt/(mu0*eps0))*((fields.Bz[i, j, k] - fields.Bz[i-1, j, k])/fields.dx - (fields.Bx[i, j, k] - fields.Bx[i, j, k-1])/fields.dz)# - (dt/eps0)*J_y[i, j, k]
                
                fields.Ey[i, j, k] = E_y_new
		
def update_Ez(fields):
    for i in range(1, fields.nx):
        for j in range(1, fields.ny):
            for k in range(fields.nz):
                # NOTE: Only the inner edges can be evaluated!
                # Calculation is shifted half an index in negative x- and y-dir!
                E_z_new = fields.Ez[i, j, k] + (dt/(mu0*eps0))*((fields.Bx[i, j, k] - fields.Bx[i, j-1, k])/fields.dy - (fields.By[i, j, k] - fields.By[i-1, j, k])/fields.dx)# - (dt/eps0)*J_z[i, j, k]
                
                fields.Ez[i, j, k] = E_z_new

# interpolation calculations for E field exports
def calc_Ex_exp(fields):
    for i in range(fields.nx):
        for j in range(fields.ny):
            for k in range(fields.nz):
                fields.Ex_exp[i, j, k] = (fields.Ex[i, j, k]+fields.Ex[i, j+1, k]+fields.Ex[i, j, k+1]+fields.Ex[i, j+1, k+1])/4.0

def calc_Ey_exp(fields):
    for i in range(fields.nx):
        for j in range(fields.ny):
            for k in range(fields.nz):
                fields.Ey_exp[i, j, k] = (fields.Ey[i, j, k]+fields.Ey[i+1, j, k]+fields.Ey[i, j, k+1]+fields.Ey[i+1, j, k+1])/4.0

def calc_Ez_exp(fields):
    for i in range(fields.nx):
        for j in range(fields.ny):
            for k in range(fields.nz):
                fields.Ez_exp[i, j, k] = (fields.Ez[i, j, k]+fields.Ez[i+1, j, k]+fields.Ez[i, j+1, k]+fields.Ez[i+1, j+1, k])/4.0

# interpolation calculations for B field exports

def calc_Bx_exp(fields):
    for i in range(fields.nx):
        for j in range(fields.ny):
            for k in range(fields.nz):
                fields.Bx_exp[i, j, k] = (fields.Bx[i, j, k] + fields.Bx[i+1, j, k])/2.0

def calc_By_exp(fields):
    for i in range(fields.nx):
        for j in range(fields.ny):
            for k in range(fields.nz):
                fields.By_exp[i, j, k] = (fields.By[i, j, k] + fields.By[i, j+1, k])/2.0

def calc_Bz_exp(fields):
    for i in range(fields.nx):
        for j in range(fields.ny):
            for k in range(fields.nz):
                fields.Bz_exp[i, j, k] = (fields.Bz[i, j, k] + fields.Bz[i, j, k+1])/2.0

def apply_BCs_E(fields):
    # BCs for Ex fields
    # Ex BCs for z = 0 and z = nz+1
    for i in range(fields.nx):
        for j in range(fields.ny+1):
            fields.Ex[i, j, 0] = fields.Ex[i, j, 1]
            fields.Ex[i, j, fields.nz] = fields.Ex[i, j, fields.nz-1]

    # Ex BCs for y = 0 and y = ny+1
    for i in range(fields.nx):
        for k in range(fields.nz+1):
            fields.Ex[i, 0, k] = fields.Ex[i, 1, k]
            fields.Ex[i, fields.ny, k] = fields.Ex[i, fields.ny-1, k]

    # BCs for Ey fields
    # Ey BCs for z = 0 and z = nz+1
    for i in range(fields.nx+1):
        for j in range(fields.ny):
            fields.Ey[i, j, 0] = fields.Ey[i, j, 1]
            fields.Ey[i, j, fields.nz] = fields.Ey[i, j, fields.nz-1]

    # Ey BCs for x = 0 and x = nx+1
    for j in range(fields.ny):
        for k in range(fields.nz+1):
            fields.Ey[0, j, k] = fields.Ey[1, j, k]
            fields.Ey[fields.nx, j, k] = fields.Ey[fields.nx-1, j, k]

    # BCs for Ez fields
    # Ez BCs for y = 0 and y = ny+1
    for i in range(fields.nx+1):
        for k in range(fields.nz):
            fields.Ez[i, 0, k] = fields.Ez[i, 1, k]
            fields.Ez[i, fields.ny, k] = fields.Ez[i, fields.ny-1, k]

    # Ez BCs for x = 0 and x = nx+1
    for j in range(fields.ny+1):
        for k in range(fields.nz):
            fields.Ez[0, j, k] = fields.Ez[1, j, k]
            fields.Ez[fields.nx, j, k] = fields.Ez[fields.nx-1, j, k]

def apply_BCs_B(fields):
    # BCs for Bx
    for j in range(fields.ny):
        for k in range(fields.nz):
            fields.Bx[0, j, k] = fields.Bx[1, j, k]
            fields.Bx[fields.nx, j, k] = fields.Bx[fields.nx-1, j, k]

    # BCs for By
    for i in range(fields.nx):
        for k in range(fields.nz):
            fields.By[i, 0, k] = fields.By[i, 1, k]
            fields.By[i, fields.ny, k] = fields.By[i, fields.ny-1, k]

    # BCs for Bz
    for i in range(fields.nx):
        for j in range(fields.ny):
            fields.Bz[i, j, 0] = fields.Bz[i, j, 1]
            fields.Bz[i, j, fields.nz] = fields.Bz[i, j, fields.nz-1]

for t in range(tsteps):

    # Ey-source in the middle of the domain
    fields.Ey[10, 10, 10] = 1.0*np.sin(2*np.pi*1000000000*dt*t)
    fields.Ey[11, 10, 10] = 1.0*np.sin(2*np.pi*1000000000*dt*t)
    fields.Ey[10, 10, 11] = 1.0*np.sin(2*np.pi*1000000000*dt*t)
    fields.Ey[11, 10, 11] = 1.0*np.sin(2*np.pi*1000000000*dt*t)

    update_Ex(fields)
    update_Ey(fields)
    update_Ez(fields)

    apply_BCs_E(fields)

    # interpolation of Ex, Ey and Ez values for Paraview export
    calc_Ex_exp(fields)
    calc_Ey_exp(fields)
    calc_Ez_exp(fields)

    update_Bx(fields)
    update_By(fields)
    update_Bz(fields)

    # apply_BCs_B(fields)

    # interpolation of Bx, By and Bz values for Paraview export
    calc_Bx_exp(fields)
    calc_By_exp(fields)
    calc_Bz_exp(fields)

    filename = "timestep " + str(t)

    export_data(fields, domain, filename)