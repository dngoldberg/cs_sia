# **cs_sia**

# a simple matlab and python shallow-ice model to run on the MITgcm cube-sphere grid

### Equations solved

The functions solve the isothermal shallow-ice equations:

$$ \frac{\partial H}{\partial t} + \nabla\cdot\vec{q} = \dot{a}. $$

Here $H$ is ice thickness, i.e. $H = z_s - z_b$ where $z_s$ is surface and bed elevation, resp. $\dot{a}$ is the \textit{surface mass balance}, the amount of accumulation or ablation of ice at a location. $\vec{q}$ is depth-integrated mass flux. According to the \textit{shallow-ice} approximation (SIA), this can be written

$$ q = -\frac{2A}{n+2} {\tau}_d^n H^2 $$

where ${\tau}_d$ is the *driving stress* $\rho_i g H \alpha$, with $\alpha$ the surface slope; $A$ is the *Glen's Law* coefficient and $n$ is taken to be 3. (See Cuffey and Paterson 2011, Ch 8.)

Altogether this gives

$$ \frac{\partial H}{\partial t} = \nabla\cdot(D\nabla(H+z_b)) + \dot{a}, $$

$$ D = \frac{2A}{n+2} (\rho_i g)^n |\nabla z_s|^{n-1} H^{n+2}. $$

The form of $q$ arises from several assumptions:
- the ice sheet is a *power-law shear thinning* viscous rheology with $\dot{\varepsilon} = \tau^n$,
- the only type of deformation that is important is *vertical shearing*,
- there is no movement at the ice bed.
