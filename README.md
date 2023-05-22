# Einstein_IFE
This repository contains example LAMMPS scripts for computing Interfacial Free Energies (IFEs) via Einstein crystals. The general idea is to use an Einstein crystal as a reference state to which the solid bulk and slab may be easily transformed into. As the free energy of an Einstein crystal does not depend upon the cartesian positions of the constituent atoms an explicit reconstruction is avoided.



$$\gamma_{Interface} = \gamma_{Liquid} + \frac{\Delta F_{Bulk}^{Ein.} - \Delta F_{Slab}^{Ein.}}{2A}$$

Where 

Transformation of the solid bulk and slab is achieved via [Thermodynamic Integration](#Thermodynamic-Integration).

The IFE of a liquid is more efficiently calculated via the [Kirkwood-Buff method](#The-Kirkwood-Buff-Method).

Full details available at [DOI: 10.1063/5.0095130](https://doi.org/10.1063/5.0095130).

## Thermodynamic Integration

The Free Energy difference between system **A** and **B** is given by:

$$\Delta F_A^B = \int_{\lambda=0}^{\lambda=1} \left< \frac{\partial H(\lambda)}{\partial\lambda} \right>_\lambda d\lambda$$

Where $\lambda$ is a control parameter which smoothly changes between system **A** and system **B**:

$$H(\lambda) = \left( 1 - \lambda \right) H_A + \lambda H_B$$

Rather than use $\lambda$ directly we use a sigmoid function of $\lambda$ to improve convergence as $\lambda$ approaches 0 or 1:

$$f(\lambda) = \lambda^5 \left( 70 \lambda^4 - 315 \lambda^3 + 540 \lambda^2 - 420 \lambda + 126 \right)$$

The derivative $\frac{\partial H(\lambda)}{\partial \lambda}$ is then given by the chain rule:

$$\frac{\partial H(\lambda)}{\partial \lambda} = \frac{\partial H(f(\lambda))}{\partial f(\lambda)} \frac{\partial f(\lambda)}{\partial \lambda}$$

Where $\frac{\partial f(\lambda)}{\partial \lambda}$ is given analytically by:

$$\frac{\partial f(\lambda)}{\partial \lambda} = 630 \left( \lambda^2 - \lambda \right)^4$$

The derivative $\frac{\partial H(f(\lambda))}{\partial f(\lambda)}$ is calculated numerically via a central differencing scheme:

$$\frac{\partial H(f(\lambda))}{\partial f(\lambda)} \approx \frac{H(f(\lambda) + \delta(\lambda)) - H(f(\lambda) - \delta(\lambda))}{2 \delta(\lambda)}$$

We set $\delta$ to be a function of $\lambda$ to avoid numerical issues which may occur for small values of $\lambda$ when using a fixed $\delta$. Setting $\delta$ to be 1% the value of $f(\lambda)$ generally works well:

$$\delta(\lambda) = 0.01 \times f(\lambda)$$

## The Kirkwood-Buff Method

$$\gamma_{Liquid} = \frac{1}{2} \int_{0}^{L_{x}} [P_{xx} - 0.5(P_{yy}+P_{zz})] \,dx$$

Where $\gamma_{Liquid}$ is the surface tension of the liquid, $L_x$ is the length of the simulation cell in the $x$ direction and $P_{xx}$, $P_{yy}$ and $P_{zz}$ are the diagonal components of the pressure tensor. The integral is performed numerically by dividing the simulation cell into thin slices over the perpendicular direction $x$ and computing the average pressure tensor in each slice. The prefactor of $1/2$ accounts for the presence of two liquid/vacuum interfaces.

Note the KB method is applicable only to fluid interfaces and cannot be used for solids.

## A Worked Example

