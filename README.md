# Einstein_IFE
This repository contains example LAMMPS scripts for computing Interfacial Free Energies (IFEs) via Einstein crystals. All scripts are tested to work with the 23 June 2022 stable release of [LAMMPS](https://github.com/lammps/lammps/tree/stable_23Jun2022).

The general idea of this approach is to use an Einstein crystal as a reference state into which the solid bulk and slab may be easily transformed. As the free energy of an Einstein crystal does not depend upon the cartesian positions of the constituent atoms, an explicit reconstruction of the material is avoided.

The thermodynamic pathway to create the interface may be described as:

1. Generate a vacuum gap in the liquid (creating 2 liquid/vacuum interfaces)
2. Convert the solid bulk material into an Einstein crystal
3. Convert the Einstein crystal into a solid slab (within the vacuum gap, and in contact with the 2 liquid interfaces)

The total free energy of this process is divided by twice the area of the simulation cell to obtain an IFE. In practice, the slab is converted into an Einstein crystal and the free energy of this part of the process is negated. The free energy required to create the liquid/vacuum interfaces also simplifies:

$$\gamma_{Interface} = \gamma_{Liquid} + \frac{\Delta F_{Bulk}^{Ein.} - \Delta F_{Slab}^{Ein.}}{2A}$$

Where $\gamma_{Interface}$ is the IFE of the interface, $\gamma_{Liquid}$ is the IFE (surface tension) of the liquid with vacuum, $\Delta F_{Bulk}^{Ein.}$ is the free energy required to transform the solid bulk into an Einstein crystal (of the same stoichiometry as the slab), $\Delta F_{Slab}^{Ein.}$ is the free energy required to transform the slab into an Einstein crystal and $A$ is the area of one interface.

Transformation of the solid bulk and slab to an Einstein crystal is achieved via [Thermodynamic Integration](#Thermodynamic-Integration).

The IFE of a liquid is more efficiently calculated via the [Kirkwood-Buff method](#The-Kirkwood-Buff-Method).

A [worked example](#A-Worked-Example) is available below.

Full details available at [DOI: 10.1063/5.0095130](https://doi.org/10.1063/5.0095130).

## Theoretical Background

### Thermodynamic Integration

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

### The Kirkwood-Buff Method

In a liquid the IFE is equivalent to the surface tension, as described by the [Shuttleworth equation](https://doi.org/10.1088/0370-1298/63/5/302). Therefore the IFE can be calculated much more efficiently via methods other than Thermodynamic Integration. The method of Kirkwood and Buff relates the surface tension to fluctuations in the local pressure tensor:

$$\gamma_{Liquid} = \frac{1}{2} \int_{0}^{L_{x}} [P_{xx} - 0.5(P_{yy}+P_{zz})] \,dx$$

Where $\gamma_{Liquid}$ is the surface tension of the liquid, $L_x$ is the length of the simulation cell in the $x$ direction and $P_{xx}$, $P_{yy}$ and $P_{zz}$ are the diagonal components of the pressure tensor. The integral is performed numerically by dividing the simulation cell into thin slices over the perpendicular direction $x$ and computing the average pressure tensor in each slice. The prefactor of $1/2$ accounts for the presence of two liquid/vacuum interfaces in a slab configuration.

Note that the Kirkwood-Buff method is applicable only to fluid interfaces and **CANNOT** be used for interfaces with solids.

## The Scripts

These scripts provide the core functionality of the method and can be found in [scripts/](scripts/).

### input_lattequi.lmp

The `input_lattequi.lmp` script reads a `data.lmp` file. The simulation cell is expanded and an NPT simulation using the specified settings and produces a `lattequi_data.lmp` file with the expanded cell and average lattice vectors. These settings are available at the top of the script:

```
variable seed equal          94275     # random number seed
variable T equal             300.0     # temperature (kelvin)
variable P equal               0.0     # pressure (bar)
variable isotropy equal          1     # isotropy flag (1 = isotropic, 2 = anisotropic, 3 = full triclinic,
                                       # 4 = triclinic with xy and xz fixed, 5 = triclinic with x, xy and xz fixed)

variable dt equal            0.001     # timestep (femtoseconds)
variable screen equal         1000     # screen and writing output frequency (steps)

variable rep equal               1     # replication flag (1 = TRUE)
variable repx equal             12     # replication factor for x
variable repy equal             12     # replication factor for y
variable repz equal             12     # replication factor for z

variable Nlatt equal        600000     # lattice averaging time (steps)
variable ldelay equal       100000     # delay before starting lattice averaging (steps)
variable lsample equal         100     # lattice vector sampling frequency (steps)

variable traj equal              0     # 1 = print traj files
```

### input_run.lmp

The `input_run.lmp` script reads a `data.lmp` file.  An NVT equilibration simulation is performed followed by an NVT production simulation that calculates the average potential energy alongside a compressed trajectory (`prod_traj.lmp.gz`) sampled at the same points. These settings are available at the top of the script:

```
variable seed equal          94275     # random number seed
variable T equal             300.0     # temperature (kelvin)
variable dt equal            0.001     # timestep (femtoseconds)
variable screen equal        10000     # screen output frequency (steps)
variable Esample equal        1000     # energy sampling and traj output frequency (steps)

variable rep equal               0     # replication flag (1 = TRUE)
variable repx equal              1     # replication factor for x
variable repy equal              1     # replication factor for y
variable repz equal              1     # replication factor for z

variable Nequi equal        100000     # equilibration time (steps)
variable Nprod equal        500000     # production time (steps)

variable traj equal              1     # 1 = print traj files
```

### input_rerun.lmp

The `input_run.lmp` script reads a `data.lmp` file and a `prod_traj.lmp.gz` trajectory file. The energy is re-computed at each frame of the trajectory based on the current `potential.lmp`. The re-computed potential energy is averaged and printed at the end. These settings are available at the top of the script:

```
variable Esample equal        1000     # energy sampling frequency (steps)
variable screen equal        10000     # screen and writing output frequency (steps)
```

### potential.lmp

The `potential.lmp` script is the most complex.

At the top are the settings to control the Thermodynamic Integration:

```
# Einstein settings
variable ein_lambda equal 0.0000
variable ein_delta equal 0.0
#	variable ein_delta equal 0.01
#	variable ein_delta equal -0.01

# Potential settings
variable pot_lambda equal 1.0000
variable pot_delta equal 0.0
#	variable pot_delta equal 0.01
#	variable pot_delta equal -0.01
```

The `ein_lambda` variable controls the strength of the harmonic wells while the `pot_lambda` variable controls the strength of the scaled interactions. Together they are used to transform a solid into an Einstein crystal in a two stage process.

The `ein_delta` and `pot_delta` variables control the value of $\delta(\lambda)$ for use in the central differencing scheme. These variables are given as a fraction of $f(\lambda)$ (see [Thermodynamic Integration](#Thermodynamic-Integration)).
 
Below this section are additional settings:

```
# einstein group
group einstein_group type 1 2

# einstein spring constant
variable ein_spring equal 10.0

# Additional potential settings are located further down.
# That section requires extensive modification for each system and so is impractical to include here.
```

The `einstein_group` variable is used to indicate which atom types are to be transformed into an Einstein crystal and the `ein_spring` variable indicates the spring constant to be used.

Further down the script the interaction potential of the system is given. The internal solid-solid interactions and the solid-liquid interactions must be scaled by $f(\lambda) + \delta(\lambda)$ while the liquid-liquid interactions remain unchanged. This is achieved by scaling the appropriate interactions by the `pot_eff` variable.

A basic `potential.lmp` file for the NaCl/water interfaces is located at [scripts/potential_NaCl.lmp](scripts/potential_NaCl.lmp).

A significantly more complex `potential.lmp` file for CaSO<sub>4</sub>.xH<sub>2</sub>O/solution interfaces is included in [scripts/potential_CaSO4.lmp](scripts/potential_CaSO4.lmp).

### KB.lmp







### slab_correction.lmp

This script implements the slab dipole correction given by [Balleneggar *et al.*](https://doi.org/10.1063/1.3216473). This correction is the same as implemented in LAMMPS by the `kspace_modify slab` keyword, but extended to non-orthorhombic boxes and any slab orentation. The only option is:

```
# set slab orientation (1 = yz plane, 2 = xz plane, 3 = xy plane)
variable orient equal 1
```

It is also prudent to implement walls in the `input.lmp` script to prevent translation of species across the boundary. For a slab aligned in the yz plane this would take the form:

```
fix x_walls all wall/lj93 xlo EDGE 0.001 3.5 3.0 xhi EDGE 0.001 3.5 3.0 units box pbc yes	# wall repulsion for slabs
include slab_correction.lmp									# slab dipole correction
```

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

In this calculation a cubic unit cell of NaCl is read from `data.lmp` and expanded 12x in each dimension. The simulation cell is then run under the chosen conditions in an isotropic NPT ensemble for 0.6 ns. The first 0.1 ns is discarded and then the lattice vectors are recorded every 100 fs. At the end of the simulation the average lattice vectors are calculated and a new `lattequi_data.lmp` is written with the average lattice vectors imposed.

The files to perform this calculation are located in [examples/NaCl_water_example/1_bulk/1_lattice_equilibration](examples/NaCl_water_example/1_bulk/1_lattice_equilibration).

#### Calculate enthalpy

This stage is optional but useful if the interfacial enthalpy is desired at a later date.

First, the `lattequi_data.lmp` file is copied from the previous stage and renamed `data.lmp`. 












The files to perform this calculation are located in [examples/NaCl_water_example/1_bulk/2_enthalpy](examples/NaCl_water_example/1_bulk/2_enthalpy).

#### Activate harmonic wells

In this stage multiple calculations are set up, one for each Thermodynamic Integration point of interest.

The `prod_data.lmp` is copied in from the previous calculation and named `data.lmp`. The same `data.lmp` **MUST** be used for **ALL** Thermodynamic Integration calculations to maintain consistency. This is because the scripts set the minima of the harmonic well for each atom is set to be at it's initial position, and drift of the harmonic well between points will result in an inconsistent pathway. Here the consistency has been assured by the use of symlinks in each sub-directory which point to the same `data.lmp`.

Each sub-directory is labelled with the value of $\lambda$ that has been set, e.g. `lambda_0.125/`. All sub-directories are identical save for the `potential.lmp` file which has a different value set for the variable `ein_lambda` which controls the point along the thermodynamic pathway. The value of the other important variable `pot_lambda` is kept at 1.0 for this stage indicating the interactions are still at full strength.

Running LAMMPS in each sub-directory performs a 100 ps NVT equilibration phase, followed by a 500 ps NVT production phase.

During the production phase the potential energy is sampled every 1 ps and a trajectory is recorded at the same time. The average potential energy is printed at the end of the simulation and the trajectory at the same sample points is written to `prod_traj.lmp.gz`.













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
















