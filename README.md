# Thermo4PFM

Library to compute equilibrium compositions and equal chemical potential
compositions in binary and ternary alloys.
This library supports thermodynamic data formulated according to CALPHAD
(Computer Coupling of Phase Diagrams and Thermochemistry) models, as
well as dilute binary alloys.

One major aim is to solve the Kim-Kim-Suzuki (KKS) equations to obtain
local compositions given a local mixture of two phases.
These quantities are then used to evaluate the phase-field model driving
force (right-hand side of time-evolution equation), as well as the
energy density.

In the CALPHAD approach, the Gibbs energy of each individual phase is
defined, and the model parameters are collected in a thermodynamic
database.
For multicomponent solution phases, the Gibbs energy is modeled as the
sum of 3 components:

* the contribution from the mechanical mixing of the pure components
* is the ideal mixing contribution
* the excess Gibbs energy of mixing due to non-ideal interactions

The values of coefficients used for parameterizing these energies are
stored in JSON (JavaScript Object Notation) files.
These files are parsed with the Boost Property Tree library.

Current functionalities include:

* KKS equations solver for binary (CALPHAD and dilute) and ternary
  alloys (CALPHAD)
* equilibrium compositions solvers for binary and ternary alloys
  (CALPHAD)

Solvers are based on internal Newton solver.
Specific solvers are derived classes of base class "Newton" and
implement specific Jacobian and right-hand side calculations.
Cramer's rule is used to solve the linear system at each Newton
iteration.
Solvers are templated on the dimension of the system of equations to be
solved.

Issues related to compositions possibly taking values outside of
$[0, 1]$ during iterative solve are mostly handled by using continuous
extensions of $x \log(x)$ and $(1-x) \log(1-x)$ functions.

## Dependencies

* [BOOST](http://www.boost.org)

## References

J.-L. Fattebert, S. DeWitt, A. Perron, J. Turner,
"Thermo4PFM: Facilitating Phase-field simulations of alloys with thermodynamic driving forces",
_Comput. Phys. Comm._, __288__ (2023), 108739

S.G. Kim, W.T. Kim, T. Suzuki,
"Phase-field model for binary alloys",
_Phys. Rev. E_ __60__:6 (1999) 7186.

J.-L. Fattebert, M. E. Wickett, P. E. A. Turchi,
"Phase-field modeling of coring during solidification of Au-Ni alloy
using quaternions and CALPHAD input",
_Acta Mater._, __62__ (2014), 89-104
