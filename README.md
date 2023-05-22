# Einstein_IFE
This repository contains example LAMMPS scripts for computing Interfacial Free Energies (IFEs) via Einstein crystals. The general idea is to use an Einstein crystal as a reference state into which the solid bulk and slab may be easily transformed. As the free energy of an Einstein crystal does not depend upon the cartesian positions of the constituent atoms, an explicit reconstruction is avoided.

The thermodynamic pathway to create the interface may be visualised as generating a vacuum gap in a liquid (creating 2 liquid/vacuum interfaces), conversion of bulk material into an Einstein crystal, and conversion of the Einstein crystal into a slab within the vacuum gap in contact with the 2 liquid interfaces. The total free energy of this process is divided by twice the area of the simulation cell to obtain an IFE. In practice, the slab is converted into an Einstein crystal and the free energy of this part of the process is negated. The free energy to create the liquid/vacuum interfaces also simplifies:

$$\gamma_{Interface} = \gamma_{Liquid} + \frac{\Delta F_{Bulk}^{Ein.} - \Delta F_{Slab}^{Ein.}}{2A}$$

Where $\gamma_{Interface}$ is the IFE of the interface, $\gamma_{Liquid}$ is the IFE (surface tension) of the liquid with vacuum, $\Delta F_{Bulk}^{Ein.}$ is the free energy required to transform the solid bulk into an Einstein crystal (of the same stoichiometry as the slab), $\Delta F_{Slab}^{Ein.}$ is the free energy required to transform the slab into an Einstein crystal and $A$ is the area of one interface.

Transformation of the solid bulk and slab to an Einstein crystal is achieved via [Thermodynamic Integration](#Thermodynamic-Integration).

The IFE of a liquid is more efficiently calculated via the [Kirkwood-Buff method](#The-Kirkwood-Buff-Method).

A [worked example](#A-Worked-Example) is available below.

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






<!--- The total kinetic energy is our calculations is constant so we can use the internal energy, $U$ instead of $H$. --->







We set $\delta$ to be a function of $\lambda$ to avoid numerical issues which may occur for small values of $\lambda$ when using a fixed $\delta$. Setting $\delta$ to be 1% the value of $f(\lambda)$ generally works well:

$$\delta(\lambda) = 0.01 \times f(\lambda)$$

The numerical integration may be performed with any quadrature rule; however there are two main considerations. Using a fixed number of calculation points and a simple quadrature scheme like the the Trapezoidal rule or Simpson's rule will allow the calculation points to be easily taskfarmed on large HPC machines. Alternatively, an adaptive quadrature rule may result in improved convergence with fewer calculation points, at the expense of not knowing which points to include *a priori*. An adaptive method which suggests batches of calculation points is likely to be most efficient.

## The Kirkwood-Buff Method

In a liquid the IFE is equivalent to the surface tension, as described by the Shuttleworth equation. Therefore the IFE can be calculated much more efficiently via methods other than Thermodynamic Integration. The method of Kirkwood and Buff relates the surface tension to fluctuations in the local pressure tensor:

$$\gamma_{Liquid} = \frac{1}{2} \int_{0}^{L_{x}} [P_{xx} - 0.5(P_{yy}+P_{zz})] \,dx$$

Where $\gamma_{Liquid}$ is the surface tension of the liquid, $L_x$ is the length of the simulation cell in the $x$ direction and $P_{xx}$, $P_{yy}$ and $P_{zz}$ are the diagonal components of the pressure tensor. The integral is performed numerically by dividing the simulation cell into thin slices over the perpendicular direction $x$ and computing the average pressure tensor in each slice. The prefactor of $1/2$ accounts for the presence of two liquid/vacuum interfaces in a slab configuration.

Note that the Kirkwood-Buff method is applicable only to fluid interfaces and **CANNOT** be used for interfaces with solids.

## A Worked Example

In this example the IFE of the NaCl {100} surface with pure water is calculated. The files may be found in [examples/NaCl_water_example/](examples/NaCl_water_example/). These example directories make extensive use of symlinks which enables efficient use of storage space and also helps to ensure consistency between calculations.

The worked example provided here uses a very minimal number of Thermodynamic Integration points (values of $\lambda$). Full publication quality calculations **WILL** require many additional points. Alongside the free energy integrals calculated in this example, fully converged values using 128 points and Romberg's method will be presented in paranthesis.

The three components required for the calculation of IFEs are:

1. NaCl Bulk to Einstein Crystal
2. NaCl Slab to Einstein Crystal
3. Water/Vacuum Surface Tension

### NaCl Bulk to Einstein Crystal

The first task is to compute the free energy of transforming NaCl bulk into an Einstein crystal. This is performed in 5 stages:

1. Calculate average lattice vectors
2. Calculate enthalpy
3. Activate harmonic wells
4. Deactivate interactions
5. Compute free energy per formula unit

#### Calculate average lattice vectors
















#### Calculate enthalpy

#### Activate harmonic wells

The minima of the harmonic well for each atom is set to be at it's initial position. Therefore the same `data.lmp` file **MUST** be used for **ALL** Thermodynamic Integration calculations to maintain consistency.

| $\lambda$ | $f(\lambda)$| $\frac{\partial f(\lambda)}{\partial \lambda}$ | $\delta(\lambda)$ | $H(f(\lambda) - \delta(\lambda))$ | $H(f(\lambda))$ | $H(f(\lambda) + \delta(\lambda))$ | $\frac{\partial H(\lambda)}{\partial \lambda}$ |
| --- | ---| --- | --- | --- | --- | --- | --- |
0.0000 | 0.0 | 0.0 | 0.0 |  |  |  | 0.0 |
0.1250 | 0.00248228 | 0.090159774 | 2.48228E-05 | -56241.20199 | -56241.08048 | -56240.95891 | 441.4538144 |
0.2500 | 0.048927307 | 0.778656006 | 0.000489273 | -56161.86652 | -56160.34391 | -56158.82123 | 2423.224422 |
0.3750 | 0.216618016 | 1.901015639 | 0.00216618 | -56021.43951 | -56018.02223 | -56014.60486 | 2999.006927 |
0.5000 | 0.5 | 2.4609375 | 0.005 | -55923.06645 | -55918.6806 | -55914.29455 | 2158.710143 |
0.6250 | 0.783381984 | 1.901015639 | 0.00783382 | -55874.67204 | -55869.91453 | -55865.15704 | 1154.492179 |
0.7500 | 0.951072693 | 0.778656006 | 0.009510727 | -55855.71546 | -55850.83414 | -55845.95246 | 399.6548102 |
0.8750 | 0.99751772 | 0.090159774 | 0.009975177 | -55851.28348 | -55846.37467 | -55841.46615 | 44.36653954 |
1.0000 | 1 | 0.0 | 0.01 |  |  |  | 0.0 |

#### Deactivate interactions

#### Compute free energy per formula unit

The free energy of transforming bulk NaCl into an Einstein crystal is given by the free energy of activating harmonic wells plus the free energy of deactivating the interactions. Dividing by the number of NaCl formula units gives a value per formular unit that may be scaled to match any slab calculation.











This value only needs to be calculated once and can be re-used for all NaCl IFE calculations.

### NaCl Slab to Einstein Crystal

In the slab calculations the correction of [Balleneggar *et al.*](https://doi.org/10.1063/1.3216473) is applied. This correction is only strictly required for slabs with a dipole or charge, but has minimal computational overhead and is also included here as an example.

#### Calculate average lattice vectors

#### Calculate enthalpy

#### Activate harmonic wells

#### Deactivate interactions

### Water/Vacuum Surface Tension

The files may be found in [examples/KB_water_tension/](examples/KB_water_tension/).

## Additional Examples

A significantly more complex potential file for CaSO<sub>4</sub>.xH<sub>2</sub>O/solution interfaces is included in [scripts/](scripts/).


















