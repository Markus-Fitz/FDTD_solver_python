# Finite difference time domain solver for Maxwell's equations.

This project's aim is to create a simple and light-weight numerical solver for Maxwell's equations.

This solver is an educational project first, which preferres explicit expression over elegance and performance. It will be optimized as much as it can be in this context, though other implementations in different languages may follow for that.

The default case is an electric field-source in the bottom-left-front corner of a domain which is 25cm x 25cm x 25cm in size. The source radiates at 1V/m field strength with a frequency of 1 GHz.

The following images show the evolution of the fields after 100, 200, 300, 400 and 500 timesteps.

<img width="621" height="623" alt="Screenshot From 2026-02-14 16-10-37" src="https://github.com/user-attachments/assets/3608f08c-7cb0-486d-89f7-63c81ed38e97" />
<img width="621" height="623" alt="Screenshot From 2026-02-14 16-10-26" src="https://github.com/user-attachments/assets/16694bbb-8970-46ea-99e0-cd2ff42e1c3d" />
<img width="621" height="623" alt="Screenshot From 2026-02-14 16-10-11" src="https://github.com/user-attachments/assets/0d2bc6fb-fd3e-448e-b018-4a69ffc93406" />
<img width="621" height="623" alt="Screenshot From 2026-02-14 16-09-57" src="https://github.com/user-attachments/assets/08a83bf6-8b1b-424f-aa77-a8e5e319bd60" />
<img width="621" height="623" alt="Screenshot From 2026-02-14 16-09-44" src="https://github.com/user-attachments/assets/e32c8225-8418-45ee-9bb9-dd6482374a4c" />

## Maxwell's equations

$$
\nabla \cdot E = \frac{\rho}{\epsilon_0}
$$
$$
\nabla \cdot B = 0
$$
$$
\nabla \times E = - \frac{\partial B}{\partial t}
$$
$$
\nabla \times B = \mu_0 (J + \epsilon_0 \frac{\partial E}{\partial t})
$$

## Yee grid

The Yee grid FDTD method is based on a staggered grid with the electric and magnetic fields stored at different places within each cell.

The electric fields are separated into the x-, y- and z-components and stored on the edges of each cubic cell. The directions of the respective field is aligned with the cell edge on which it is placed.

The magnetic fields are stored in the middle of each cell face. The magnetic fields are placed on the cell faces for which the field direction and the normal vector of the cell face align. The By field components are stored in the middle of the cell-faces in y+ and y- directions of the cubic cell.

## Discretized operators

The functions specifying the time-evolution for the electric and magnetic field components require the evaluation of the cross-product of both the electric and the magnetic fields.

To come up with the discretized partial differential equations, the curl of a vector field

$$
\nabla \times A
$$

can first be expressed in three separate equations

$$
(\nabla \times A)_x = \frac{\partial A_y}{\partial z} - \frac{\partial A_z}{\partial y}
$$
$$
(\nabla \times A)_y = \frac{\partial A_z}{\partial x} - \frac{\partial A_x}{\partial z}
$$
$$
(\nabla \times A)_z = \frac{\partial A_x}{\partial y} - \frac{\partial A_y}{\partial x}
$$

. Each partial differential acting on a field $F$ at a point specified by the index i, j, and k for the three spacial directions (from now on noted as a superscript on the corresponding field) can be expressed by the approximations

$$
\frac{\partial F}{\partial x} |_{i, j, k} \approx \frac{F^{i+\frac{1}{2}, j, k} - F^{i-\frac{1}{2}, j, k}}{\Delta x}
$$
$$
\frac{\partial F}{\partial y} |_{i, j, k} \approx \frac{F^{i, j+\frac{1}{2}, k} - F^{i, j-\frac{1}{2}, k}}{\Delta y}
$$
$$
\frac{\partial F}{\partial x} |_{i, j, k} \approx \frac{F^{i, j, k+\frac{1}{2}} - F^{i, j, k-\frac{1}{2}}}{\Delta z}
$$

, where the $\Delta x$, $\Delta y$ and $\Delta z$ describe the length of the cells in each direction (which is equal to the distance in between each field quantity).

With this approximation, the cross product of the curl operator acting on a vector field can be expressed as

$$
(\nabla \times A)_x^{i, j, k} \approx \frac{A_y^{i, j, k+\frac{1}{2}} - A_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{A_z^{i, j+\frac{1}{2}, k} - A_z^{i, j-\frac{1}{2}, k}}{\Delta y}
$$
$$
(\nabla \times A)_y^{i, j, k} \approx \frac{A_z^{i+\frac{1}{2}, j, k} - A_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{A_x^{i, j, k+\frac{1}{2}} - A_x^{i, j, k-\frac{1}{2}}}{\Delta z}
$$
$$
(\nabla \times A)_z^{i, j, k} \approx \frac{A_x^{i, j+\frac{1}{2}, k} - A_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{A_y^{i+\frac{1}{2}, j, k} - A_y^{i-\frac{1}{2}, j, k}}{\Delta x}
$$

. The nature of the staggered Yee grid is very convenient for evaluating these approximations to the partial differential equations. As is seen in the approximations above, the field values shifted by half a cell distance is needed to evaluate the central differences. This fits perfectly with the staggered grid proposed by Yee in his paper "Numerical Solution of Initial Boundary Value Problems Involving Maxwell's Equations in Isotropic Media". The values of $E_x$, $E_y$, $E_z$ as well as $B_x$, $B_y$ and $B_z$ are placed in such a way that the central difference scheme described above comes naturally.

## Time evolution of fields

Assuming that the initial fields adhere to all four equations, the induction law of Faraday (equation nr. 3) and Ampère's law (equation nr. 4) describe the time-evolution of both the electric and magnetic fields. The time-step becomes relevant now, and is noted by the time-index $n$, used to describe the $n$-th timestep $t_n$.

The time-derivative $\frac{1}{\partial t}$ can be approximated in a similar way to the spacial derivatives. For a field $F(t)$, the approximation

$$
\frac{\partial F}{\partial t} \approx \frac{F(t_{n+1}) - F(t_n)}{\Delta t}
$$

is used, with $\Delta t$ being one time-step in the discretized time domain.

With this approximation, the Induction law of Faraday,

$$
\nabla \times E = - \frac{\partial B}{\partial t}
$$

, can be rewritten as the three equations

$$
\frac{E_y^{i, j, k+\frac{1}{2}} - E_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{E_z^{i, j+\frac{1}{2}, k} - E_z^{i, j-\frac{1}{2}, k}}{\Delta y} = - \frac{B_x^{i, j, k}(t_{n+1}) - B_x^{i, j, k}(t_n)}{\Delta t}
$$
$$
\frac{E_z^{i+\frac{1}{2}, j, k} - E_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{E_x^{i, j, k+\frac{1}{2}} - E_x^{i, j, k-\frac{1}{2}}}{\Delta z} = - \frac{B_y^{i, j, k}(t_{n+1}) - B_y^{i, j, k}(t_n)}{\Delta t}
$$
$$
\frac{E_x^{i, j+\frac{1}{2}, k} - E_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{E_y^{i+\frac{1}{2}, j, k} - E_y^{i-\frac{1}{2}, j, k}}{\Delta x} = - \frac{B_z^{i, j, k}(t_{n+1}) - B_z^{i, j, k}(t_n)}{\Delta t}
$$

. This directly leads to an expression for the three vector fields $B_x$, $B_y$ and $B_z$ at the time-step $t_{n+1}$, based only on information we have at time-step $t_n$. The following three equations therefore are our update-equations for the B-field at the next timestep

$$
B_x^{i, j, k}(t_{n+1}) = B_z^{i, j, k}(t_n) - \Delta t \left( \frac{E_y^{i, j, k+\frac{1}{2}} - E_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{E_z^{i, j+\frac{1}{2}, k} - E_z^{i, j-\frac{1}{2}, k}}{\Delta y} \right)
$$
$$
B_y^{i, j, k}(t_{n+1}) = B_y^{i, j, k}(t_n) - \Delta t \left( \frac{E_z^{i+\frac{1}{2}, j, k} - E_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{E_x^{i, j, k+\frac{1}{2}} - E_x^{i, j, k-\frac{1}{2}}}{\Delta z} \right)
$$
$$
B_z^{i, j, k}(t_{n+1}) = B_z^{i, j, k}(t_n) - \Delta t \left( \frac{E_x^{i, j+\frac{1}{2}, k} - E_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{E_y^{i+\frac{1}{2}, j, k} - E_y^{i-\frac{1}{2}, j, k}}{\Delta x} \right)
$$

. Similarly, Ampère's law can be approximated in the following three equations

$$
\frac{B_y^{i, j, k+\frac{1}{2}} - B_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{B_z^{i, j+\frac{1}{2}, k} - B_z^{i, j-\frac{1}{2}, k}}{\Delta y} = \mu_0 \left( J_x^{i, j, k} + \epsilon_0 \frac{E_x^{i, j, k}(t_{n+1}) - E_x^{i, j, k}(t_n)}{\Delta t} \right)
$$
$$
\frac{B_z^{i+\frac{1}{2}, j, k} - B_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{B_x^{i, j, k+\frac{1}{2}} - B_x^{i, j, k-\frac{1}{2}}}{\Delta z} = \mu_0 \left( J_y^{i, j, k} + \epsilon_0 \frac{E_y^{i, j, k}(t_{n+1}) - E_y^{i, j, k}(t_n)}{\Delta t} \right)
$$
$$
\frac{B_x^{i, j+\frac{1}{2}, k} - B_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{B_y^{i+\frac{1}{2}, j, k} - B_y^{i-\frac{1}{2}, j, k}}{\Delta x} = \mu_0 \left( J_z^{i, j, k} + \epsilon_0 \frac{E_z^{i, j, k}(t_{n+1}) - E_z^{i, j, k}(t_n)}{\Delta t} \right)
$$

, leading to the update equations for the E-fields

$$
E_x^{i, j, k}(t_{n+1}) = E_x^{i, j, k}(t_n) + \frac{\Delta t}{\epsilon_0} \left( \frac{1}{\mu_0} \left( \frac{B_y^{i, j, k+\frac{1}{2}} - B_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{B_z^{i, j+\frac{1}{2}, k} - B_z^{i, j-\frac{1}{2}, k}}{\Delta y} \right) - J_x^{i, j, k} \right)
$$

$$
E_y^{i, j, k}(t_{n+1}) = E_y^{i, j, k}(t_n) + \frac{\Delta t}{\epsilon_0} \left( \frac{1}{\mu_0} \left( \frac{B_z^{i+\frac{1}{2}, j, k} - B_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{B_x^{i, j, k+\frac{1}{2}} - B_x^{i, j, k-\frac{1}{2}}}{\Delta z} \right) - J_y^{i, j, k} \right)
$$

$$
E_z^{i, j, k}(t_{n+1}) = E_z^{i, j, k}(t_n) +  \frac{\Delta t}{\epsilon_0} \left( \frac{1}{\mu_0} \left( \frac{B_x^{i, j+\frac{1}{2}, k} - B_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{B_y^{i+\frac{1}{2}, j, k} - B_y^{i-\frac{1}{2}, j, k}}{\Delta x} \right) - J_z^{i, j, k} \right)
$$

. The discretization and structure of the grid becomes important for the implementation, so that step is next.

## Yee grid and implementation

The Yee grid proposed in 1966 is based on a rectilinear (in this case even more simplified - cubic) grid. The E-fields are contained in the edges of the cubes, with each E-field component fixed to the edge in its corresponding directions. The B-fields are fixed on the centers of the cell faces and are perpendicular to the cell face. 

The index notation for the coordinates of the values - $i$, $j$ and $k$ for the x-, y- and z-directions is discussed more closely in the following to avoid confusion with the half-index nature of the staggered grid.

The indices in each direction starts at 0. The indices should be thought of as indexing the cell centers, not the field values. The field values are anchored to the cells, but always shifted by half of the cell-size in one (B-fields) or two (E-fields) of the three cardinal directions. The following illustration demonstrates the staggered grid and the places at which the field values are known / evaluated.

![[Pasted image 20260213085030.png]]

$E_x^{0, 0, 0}$ is shifted half a cell size in negative y-direction and half a cell-size in negative z-direction compared to the center of the cell denoted by the indices $(0, 0, 0)$. This modifies the update-equation for $E_x$ to:

$$
	E_x^{i, j-\frac{1}{2}, k-\frac{1}{2}}(t_{n+1}) = E_x^{i, j-\frac{1}{2}, k-\frac{1}{2}}(t_n) + \frac{\Delta t}{\epsilon_0 \mu_0} \left( \frac{B_y^{i, j-\frac{1}{2}, k} - B_y^{i, j-\frac{1}{2}, k-1}}{\Delta z} - \frac{B_z^{i, j, k-\frac{1}{2}} - B_z^{i, j-1, k-\frac{1}{2}}}{\Delta y} \right) - \frac{\Delta t}{\epsilon_0}J_x^{i, j-\frac{1}{2}, k-\frac{1}{2}}
$$

.

$E_y^{0, 0, 0}$ is shifted half a cell size in negative x-direction and half a cell-size in negative z-direction compared to the center of the cell denoted by the indices $(0, 0, 0)$. This modifies the update-equation for $E_y$ to:

$$
E_y^{i-\frac{1}{2}, j, k-\frac{1}{2}}(t_{n+1}) = E_y^{i-\frac{1}{2}, j, k-\frac{1}{2}}(t_n) + \frac{\Delta t}{\epsilon_0 \mu_0} \left( \frac{B_z^{i, j, k-\frac{1}{2}} - B_z^{i-1, j, k-\frac{1}{2}}}{\Delta x} - \frac{B_x^{i-\frac{1}{2}, j, k} - B_x^{i-\frac{1}{2}, j, k-1}}{\Delta z} \right) - \frac{\Delta t}{\epsilon_0} J_y^{i-\frac{1}{2}, j, k-\frac{1}{2}}
$$

.

$E_z^{0, 0, 0}$ is shifted half a cell size in negative x-direction and half a cell-size in negative y-direction compared to the center of the cell denoted by the indices $(0, 0, 0)$. This modifies the update-equation for $E_z$ to:

$$
E_z^{i-\frac{1}{2}, j-\frac{1}{2}, k}(t_{n+1}) = E_z^{i-\frac{1}{2}, j-\frac{1}{2}, k}(t_n) +  \frac{\Delta t}{\epsilon_0 \mu_0} \left( \frac{B_x^{i-\frac{1}{2}, j, k} - B_x^{i-\frac{1}{2}, j-1, k}}{\Delta y} - \frac{B_y^{i, j-\frac{1}{2}, k} - B_y^{i-1, j-\frac{1}{2}, k}}{\Delta x} \right) - \frac{\Delta t}{\epsilon_0}J_z^{i-\frac{1}{2}, j-\frac{1}{2}, k}
$$

. Note in the adjusted equations, that the current J is defined at the edges parallel to its direction - at the same points as the electric fields in the same direction.

$B_x^{0, 0, 0}$ is shifted half a cell size in negative x-direction compared to the center of the cell denoted by the indices $(0, 0, 0)$. This modifies the update-equation for $B_x$ to:

$$
B_x^{i-\frac{1}{2}, j, k}(t_{n+1}) = B_x^{i-\frac{1}{2}, j, k}(t_n) - \Delta t \left( \frac{E_y^{i-\frac{1}{2}, j, k+\frac{1}{2}} - E_y^{i-\frac{1}{2}, j, k-\frac{1}{2}}}{\Delta z} - \frac{E_z^{i-\frac{1}{2}, j+\frac{1}{2}, k} - E_z^{i-\frac{1}{2}, j-\frac{1}{2}, k}}{\Delta y} \right)
$$

.

$B_y^{0, 0, 0}$ is shifted half a cell size in negative y-direction compared to the center of the cell denoted by the indices $(0, 0, 0)$. This modifies the update equation for $B_y$ to:

$$
B_y^{i, j-\frac{1}{2}, k}(t_{n+1}) = B_y^{i, j-\frac{1}{2}, k}(t_n) - \Delta t \left( \frac{E_z^{i+\frac{1}{2}, j-\frac{1}{2}, k} - E_z^{i-\frac{1}{2}, j-\frac{1}{2}, k}}{\Delta x} - \frac{E_x^{i, j-\frac{1}{2}, k+\frac{1}{2}} - E_x^{i, j-\frac{1}{2}, k-\frac{1}{2}}}{\Delta z} \right)
$$

.

$B_z^{0, 0, 0}$ is shifted half a cell size in negative z-direction compared to the center of the cell denoted by the indices $(0, 0, 0)$. This modifies the update equation for $B_z$ to:

$$
B_z^{i, j, k-\frac{1}{2}}(t_{n+1}) = B_z^{i, j, k-\frac{1}{2}}(t_n) - \Delta t \left( \frac{E_x^{i, j+\frac{1}{2}, k-\frac{1}{2}} - E_x^{i, j-\frac{1}{2}, k-\frac{1}{2}}}{\Delta y} - \frac{E_y^{i+\frac{1}{2}, j, k-\frac{1}{2}} - E_y^{i-\frac{1}{2}, j, k-\frac{1}{2}}}{\Delta x} \right)
$$

.
At this point the Yee-grid seems very convenient, as the values in the equations above are directly accessible without averaging values over multiple cells etc. They can be directly calculated from values in the stored arrays.

Note though, that the E-fields can only be evaluated for the inner edges. They would need B-field values outside of the simulated domain. The B-field values can however be evaluated in all cells.

For the dimension of the field arrays, this means that the size has to be the following for a volume with $n_x$ cells in x-direction, $n_y$ cells in y-direction and $n_z$ cells in z-direction:

$E_x[n_x, n_y+1, n_z+1]$
$E_y[n_x+1, n_y, n_z+1]$
$E_z[n_x+1, n_y+1, n_z]$

$B_x[n_x+1, n_y, n_z]$
$B_y[n_x, n_y+1, n_z]$
$B_z[n_x, n_y, n_z+1]$

. In Python, these arrays are implemented in six different numpy-arrays. The update functions derived above can then be implemented as for-loops iterating over all cell indices. The update-functions for the electric and magnetic fields are split into the three cardinal directions as the index is one more than the cell-domain for one or two directions depending on which field and which direction is in question.

The update equation for the $B_x$ field becomes:

``` python

for i in range(n_x + 1):
	for j in range(n_y):
		for k in range(n_z):
			# NOTE: values at position (i+1/2)*dx, j*dy, k*dz are evaluated;
			# calculation is shifted half an index in negative x-direction!
			B_x_new = B_x[i, j, k] + dt*((E_y[i, j, k+1] - E_y[i, j, k])/dz - (E_z[i, j+1, k] - E_z[i, j, k])/dy)
			B_x[i, j, k] = B_x_new
```

. For the $B_y$ field it becomes:

``` python

for i in range(n_x):
	for j in range(n_y+1):
		for k in range(n_z):
			# NOTE: values at position i*dx, (j+1/2)*dy, k*dz are evaluated;
			# calculation is shifted half an index in negative y-direction!
			B_y_new = B_y[i, j, k] + dt*((E_z[i+1, j, k] - E_z[i, j, k])/dx - (E_x[i, j, k+1] - E_x[i, j, k])/dz)
			B_y[i, j, k] = B_y_new
```

. And finally, for the $B_z$ field it becomes:

``` python

for i in range(n_x):
	for j in range(n_y):
		for k in range(n_z+1):
			# NOTE: values at position i*dx, j*dy, (k+1/2)*dz are evaluated;
			# calculation is shifted half an index in negative z-direction!
			B_z_new = B_z[i, j, k] + dt*((E_x[i, j+1, k] - E_x[i, j, k])/dy - (E_y[i+1, j, k] - E_y[i, j, k])/dx) 
			B_z[i, j, k] = B_z_new
```

. The same is done for the electric fields.
The update equation for the $E_x$ field is:

``` python

for i in range(n_x):
	for j in range(1, n_y-1):
		for k in range(1, n_z-1):
			# NOTE: Only the inner edges can be evaluated!
			# Calculation is shifted half an index in negative y- and z-dir!
			E_x_new = E_x[i, j, k] + (dt/(mu0*eps0))*((B_y[i, j, k] - B_y[i, j, k-1])/dz - (B_z[i, j, k] - B_z[i, j-1, k])/dy) - (dt/eps0)*J_x[i, j, k]
			E_x[i, j, k] = E_x_new
```

.The update equation for the $E_y$ field is:

``` python

for i in range(1, n_x-1):
	for j in range(n_y):
		for k in range(1, n_z-1):
			# NOTE: Only the inner edges can be evaluated!
			# Calculation is shifted half an index in negative x- and z-dir!
			E_y_new = E_y[i, j, k] + (dt/(mu0*eps0))*((B_z[i, j, k] - B_z[i-1, j, k])/dx - (B_x[i, j, k] - B_x[i, j, k-1])/dz) - (dt/eps0)*J_y[i, j, k]
			E_y[i, j, k] = E_y_new
```

. And finally, the update equation for the $E_z$ field is:

``` python

for i in range(1, n_x-1):
	for j in range(1, n_y-1):
		for k in range(n_z):
			# NOTE: Only the inner edges can be evaluated!
			# Calculation is shifted half an index in negative x- and y-dir!
			E_z_new = E_z[i, j, k] + (dt/(mu0*eps0))*((B_x[i, j, k] - B_x[i, j-1, k])/dy - (B_y[i, j, k] - B_y[i-1, j, k])/dx) - (dt/eps0)*J_z[i, j, k]
			E_z[i, j, k] = E_z_new
```

.
