# **cs_sia**

# a simple matlab and python shallow-ice model to run on the MITgcm cube-sphere grid

### Equations solved

The functions solve the isothermal shallow-ice equations:

$$ \frac{\partial H}{\partial t} + \nabla\cdot\vec{q} = \dot{a}. $$

Here $H$ is ice thickness, i.e. $H = z_s - z_b$ where $z_s$ and $z_b$ are surface and bed elevation, resp. $\dot{a}$ is the *surface mass balance*, the amount of accumulation or ablation of ice at a location. $\vec{q}$ is depth-integrated mass flux. According to the *shallow-ice* approximation (SIA), this can be written

$$ q = -\frac{2A}{n+2} {\tau}_d^n H^2 $$

where ${\tau}_d$ is the *driving stress* $\rho_i g H \alpha$, with $\alpha$ the surface slope; $A$ is the *Glen's Law* coefficient and $n$ is taken to be 3. (See Cuffey and Paterson 2011, Ch 8 for detail on how this equation arises.) The above does not consider basal sliding -- but there are common parameterisations used.

Altogether this gives

$$ \frac{\partial H}{\partial t} = \nabla\cdot(D\nabla(H)) + \nabla\cdot(D\nabla(z_b)) + \dot{a} \quad (1), $$

$$ D = \frac{2A}{n+2} (\rho_i g)^n |\nabla z_s|^{n-1} H^{n+2}. $$

The form of $q$ arises from several assumptions:
- the ice sheet is a *power-law shear thinning* viscous rheology with $\dot{\varepsilon} = \tau^n$,
- the only type of deformation that is important is *vertical shearing*,
- there is no movement at the ice bed.

### Numerical Scheme

We start by integrating Equation (1) over a single computational cell ($C_{ij}$) to get:

$$ \frac{\partial H_{i,j}}{\partial t} = \frac{1}{|C_{ij}|}\int_{\partial C_{ij}} D \nabla H \cdot \vec{n} ds + \frac{1}{|C_{ij}|}\int_{\partial C_{ij}} D \nabla z_b \cdot \vec{n} ds + \frac{1}{|C_{ij}|} \dot{a}. $$

Now the 1st term on the right hand side can be discretised as

$$ 
-\frac{1}{|C_{ij}|\Delta x_{c,i-1/2,j}} \Delta y_{i-1/2,j} D_{i-1/2,j}(H_{ij}-H_{i-1,j}) +  \frac{1}{|C_{ij}|\Delta x_{c,i+1/2,j}} \Delta y_{i+1/2,j} D_{i+1/2,j}(H_{i+1,j}-H_{ij}) $$ 
$$ -  \frac{1}{|C_{ij}|\Delta y_{c,i,j-1/2}} \Delta x_{i,j-1/2} D_{i,j-1/2}(H_{ij}-H_{i,j-1}) + \frac{1}{|C_{ij}|\Delta y_{c,i,j+1/2}} \Delta x_{i,j+1/2} D_{i,j+1/2}(H_{i,j+1}-H_{ij}),
$$

where for instance $\Delta x_{i,j-1/2}$ is the width of the "south" (bottom) edge of the cell and $\Delta y_{i-1/2,j}$ of the "west" edge. Meanwhile $\Delta x_{c,i+1/2,j}$ is the distance from the cell center to the centre of the cell's "eastern" neighbor, and similarly for other distances. $D_{i-1/2,j}$ is the function $D$ evaluated at the western edge of the cell. (Actually $D$ is found in the centre of each cell in the code, and then averaged to cell edges.) $|C_{ij}|$ is cell area. the second term on the right hand side is evaluated similarly.

The left hand side is approximated by 

$$ \frac{H_{ij}^{(k+1)}-H_{ij}^{(k)}}{\Delta t} $$

where the $k+1$ and $k$ superscript indicate time step $k$ and $k+1$. This raises the question of whether the $H$ terms on the right hand side are at time level $k$ or $k+1$. If the former, this is an explicit method, which is subject to instability when time steps are too small -- and is not recommended. 

A "semi" implicity method is where the $H$ terms in the right hand side, such as $(H_{ij}-H_{i-1,j})$ and similar, are at level $k+1$, i.e. $(H_{ij}^{(k+1)}-H_{i-1,j}^{(k+1)})$ and the same for similar terms -- but the $D_{i-1/2,j}$ and similar, which depend on thickness, are found with $H^{(k)}$, as is the mass balance $\dot{a}$ (if applicable). 

Multiplying through by $\Delta t$, we can then rearrange so all terms involving thickness at time $(k+1)$ on the left hand side.. and all else on the right. This forms a *linear system* 

$$ \boldsymbol{B} H^{(k+1)} = r $$

to be solved for the solution vector $H$ at time step $k+1$, where $B$ and $r$ are a matrix and vector, respectively, the elements of which can be found at time $(k)$. 

For instance, the $i$ th diagonal entry of $B$, $B_{ii}$, is given by 

$$ B_{ii} = 1 - \frac{1}{|C_{ij}|\Delta x_{c,i-1/2,j}} \Delta y_{i-1/2,j} D_{i-1/2,j} - \frac{1}{|C_{ij}|\Delta x_{c,i+1/2,j}} \Delta y_{i+1/2,j} D_{i+1/2,j} $$

$$ - \frac{1}{|C_{ij}|\Delta y_{c,i,j-1/2}} \Delta x_{i,j-1/2} D_{i,j-1/2} - \frac{1}{|C_{ij}|\Delta y_{c,i,j+1/2}} \Delta x_{i,j+1/2} D_{i,j+1/2}$$

### Lateral boundary conditions

There is nothing too complex about the boundary. A growing ice sheet will advance across an ice-free continent. An exception is a continental edge. We assume that ice thickness/surface in an adjacent ocean cell is always zero. Thus ice is allowed to flow into the ocean (and be lost to the sea) but will not expand.

It is to be discussed whether this flux needs to be tracked for mass conservation. (Though there is no point in doing so if surface mass balance is not also conservative in terms of the ice-ocean-atmosphere system.)

### Numerical Scheme -- MATLAB details

The difficulty in implementing such a scheme in MATLAB or Python is the notion of the 'east' and 'west' and 'north' and 'south' neighbors. Each grid cell has them.. but they cannot be found by incrementing rows and colums in an array. Rather, we transform all grid variables into 1D vectors, and adopt a consistent numbering of the cells. We access east, west, south and north positions through arrays of indices stored in the array `dof_matrix`.

To solve the matrix, we do not actually form it, but store its coefficients in a 5-column matrix `Acols`. A user-defined function `cg_func` calculates the *action* of the matrix on a vector, using `dof_matrix`. The conjugate gradient method with a black-box operator then solves the linear system.

The $\dot{a}$ field is generated with the `smb` function -- a very simple latitude- and height-dependent function.

### Numerical Scheme -- Python details

TBD

## References

Cuffey, Kurt M., and William Stanley Bryce Paterson. The physics of glaciers. Academic Press, 2010.


