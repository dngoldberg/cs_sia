# **cs_sia**

# a simple matlab and python shallow-ice model to run on the MITgcm cube-sphere grid

### Equations solved

The functions solve the isothermal shallow-ice equations:

$$ \frac{\partial H}{\partial t} + \nabla\cdot\vec{q} = \dot{a}. $$

Here $H$ is ice thickness, i.e. $H = z_s - z_b$ where $z_s$ is surface and bed elevation, resp. $\dot{a}$ is the \textit{surface mass balance}, the amount of accumulation or ablation of ice at a location. $\vec{q}$ is depth-integrated mass flux. According to the \textit{shallow-ice} approximation (SIA), this can be written

$$ q = -\frac{2A}{n+2} {\tau}_d^n H^2 $$

where ${\tau}_d$ is the *driving stress* $\rho_i g H \alpha$, with $\alpha$ the surface slope; $A$ is the *Glen's Law* coefficient and $n$ is taken to be 3. (See Cuffey and Paterson 2011, Ch 8.)

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
-\frac{1}{|C_{ij}|\Delta x_{i-1/2,j}} \Delta y_{i-1/2,j} D_{i-1/2,j}(H_{ij}-H_{i-1,j}) +  \frac{1}{|C_{ij}|\Delta x_{i+1/2,j}} \Delta y_{i+1/2,j} D_{i+1/2,j}(H_{i+1,j}-H_{ij}) $$ 
$$ -  \frac{1}{|C_{ij}|\Delta y_{i,j-1/2}} \Delta x_{i,j-1/2} D_{i,j-1/2}(H_{ij}-H_{i,j-1}) + \frac{1}{|C_{ij}|\Delta y_{i,j+1/2}} \Delta x_{i,j+1/2} D_{i,j+1/2}(H_{i,j+1}-H_{ij})
$$
