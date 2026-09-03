Introduction  {#mainpage}
---------

This documentation describes the HFBTHO solver.
HFBTHO is an implementation of the nuclear energy density functional theory
(DFT), where the energy of the nucleus is obtained by integration over space of
some phenomenological energy density, which is itself a functional of the
intrinsic neutron and proton densities. In HFBTHO, the energy density derives
either from the zero-range Skyrme or finite range Gogny effective two-body
interaction, or is parametrized directly as a functional form of the particle
densities as for, e.g., the SeaLL1 functional. Nuclear superfluidity is treated
at the Hartree-Fock-Bogoliubov (HFB) approximation; in the case of nuclei with
an integer number of particles (either protons, neutrons, or both), the
blocking prescription at the equal filling approximation to determine the
energy of the nucleus. The code also allows solving the HFB equation at finite
temperature. In the case of Skyrme functionals, symmetry restoration techniques
are available: projection on particle number can be performed either in the
canonical basis or in the Bogoliubov basis; angular momentum and parity
projection are performed in the Bogoliubov basis and can be combined with
particle number projection. Projection is only implemented after variation.
Axial symmetry of the nuclear shape and time-reversal invariance are built-in
symmetries.

The program expands the neutron and proton wave functions on the cylindrical
harmonic oscillator (HO) basis. Expansion coefficients are determined by block
diagonalization of the HFB matrix, which is defined as the functional
derivative of the energy with respect to the neutron and proton densities.
Since the HFB matrix depends on the densities, and therefore on its own
eigenfunctions, the problem to solve is highly non-linear and requires an
iterative method (self-consistency). Suitable constraints are used to obtain
states with specific configurations. The constraints current available can be
applied on the mass multipole moments or on the expectation value of the
Gaussian neck operator.

