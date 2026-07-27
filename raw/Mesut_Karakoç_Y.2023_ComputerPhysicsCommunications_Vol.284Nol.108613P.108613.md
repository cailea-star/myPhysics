---
[source_json]: Mesut_Karakoç_Y.2023_ComputerPhysicsCommunications_Vol.284Nol.108613P.108613.json
[doi]: 10.1016/j.cpc.2022.108613
[first_author]: Mesut, Karakoç
[affiliations]: Department of Physics, Faculty of Science, Akdeniz University, TR 07070, Antalya, Türkiye
---

# BiFold: A Python code for the calculation of double-folded (bifold) potentials with density-in/dependent nucleon-nucleon interactions

### Abstract

BiFold calculates the density-dependent (DDM3Yn, BDM3Yn, CDM3Yn) or independent double-folded potentials between two colliding spherical nuclei. It is written in a Python package form to give the ability to use the potentials directly in a nuclear reaction/structure code. In addition to using Woods-Saxon/Fermi or Gaussian functions, the code also allows the definition of nuclear matter densities using pre-calculated densities in a data file. The manuscript provides an overview of the double folding model and the use of the code. Program Title: BiFold CPC Library link to program files:https://doi.org/10.17632/rmdx86p9dg.1 Developer's repository link:https://github.com/mkarakoc/BiFold Code Ocean capsule:https://codeocean.com/capsule/5286231 Licensing provisions: GPLv3 Programming language: Python 3.x Nature of problem: BiFold calculates the real part of the nuclear potential between two colliding spherical nuclei by integrating a density-independent/dependent nucleon-nucleon (NN) interaction [1,2,3] over the nuclear matter densities of the two nuclei. The code is based on M3Y Reid/Paris NN interactions [1,2,3] by default, but it is possible to define custom NN interactions when necessary. Solution method: The code uses the Fourier transform method in spherical coordinates to calculate the potential. The method simplifies the sixfold integration [1] and makes the calculation significantly faster. The integration is done by default using Simpson's integration method, but Filon's integration method is also available. [1]G.R. Satchler, W.G. Love, Phys. Rep. 55 (1979) 183.[2]A. Kobos, B. Brown, P. Hodgson, G. Satchler, A. Budzanowski, Nucl. Phys. A 384 (1) (1982) 65–87.[3]D.T. Khoa, W. von Oertzen, H.G. Bohlen, S. Ohkubo, J. Phys. G, Nucl. Part. Phys. 34 (3) (2007) R111–R164.




### Motivation

##### quotation-01

```tags
[claim-type]: motivation
[tags]: density_dependent
[source]: section 1.
```

```quote
Many studies (see the references in the present work) have used these DD treatments of NN interactions, but there are very few published codes [7, 8] to be able to reproduce the results of these works.
```

##### quotation-02

```tags
[claim-type]: background
[tags]: nuclear_matter_density, spherical
[source]: section 1.
```

```quote
The usual DF potential between two spherical nuclei is constructed by integrating over an effective nucleon-nucleon (NN) interaction with nuclear matter densities representing nucleons of both nuclei.
```

##### quotation-03

```tags
[claim-type]: background
[tags]: density_dependent
[source]: section 1.
```

```quote
The main reason for density dependence (DD) at NN interaction is the Pauli principle effects in the nuclear medium of both colliding nuclei.
```




### Methods

##### quotation-04

```tags
[claim-type]: background
[tags]: double_folding_potential, m3y_interaction
[source]: section 2.2.
```

```quote
Although BiFold can use a wide range of density-independent NN interactions in the DF potential calculations, M3Y type interactions [10, 9, 11] are defined by default in the code since they are perhaps the most widely used ones.
```

##### quotation-05

```tags
[claim-type]: definition
[tags]: fourier_transform
[source]: section 2.2, Eqs. (15)–(17).
```

```quote
Now, the exchange potential with Fourier transforms in spherical coordinates will take the form [1, 17]:
```

```math
U^{EX}(R)=4\pi\int_0^\infty G(R,s)j_0(K(R)s/M)v^{EX}(s)s^2\,ds.
```

##### quotation-06

```tags
[claim-type]: method
[tags]: double_folding_potential, m3y_interaction, fourier_transform
[source]: section 4.1.
```

```quote
The reason for choosing this distribution is to obtain an exact analytical double-folded potential since the density-independent NN effective interaction (M3Y-Reid [...] also has a similar mathematical form. Then this is easily achieved by using the Fourier transform techniques [...] for the double folding integral [...], but this time with an analytical integration.
```

##### quotation-07

```tags
[claim-type]: assumption
[tags]: double_folding_potential, frozen_density_approximation
[source]: section 2.3, Eq. (19).
```

```quote
The overlapping density for the direct part and the ZR exchange part of the folded potential has been approximated in most of the folding potential calculations [2, 33, 34, 27, 9, 35, 29, 32] as
```

```math
\rho=\rho_p(\vec r_p)+\rho_t(\vec r_t).
```




### Results

##### quotation-08

```tags
[claim-type]: result
[tags]: double_folding_potential
[source]: section 3, Fig. 3.
```

```quote
The output of the BiFold calculation shown in Fig. 3 gives individual information about the potentials and the functions used in the calculations. [...] The output also contains printout of the calculated potential UR versus radial distance R between two nuclei.
```

##### quotation-09

```tags
[claim-type]: result
[tags]: double_folding_potential
[source]: Fig. 4.
```

```quote
It is the double folding potential [dash-dot] with direct [solid] and exchange [dash] parts for the α + 40Ca elastic scattering.
```

##### quotation-10

```tags
[claim-type]: definition
[tags]: mean_relative_error, double_folding_potential
[source]: section 4, Eq. (24), Table 4.
```

```quote
This formula defines a mean relative error (mre) [8] for comparing BiFold calculations with reference calculations where U A and U B are the reference and the present double folding potentials [...], respectively. [...] The results of the mre calculations are in Table 4, and the details of the test cases are in the following sections.
```

##### quotation-11

```tags
[claim-type]: comparison
[tags]: mean_relative_error, double_folding_potential
[source]: section 4.1, Fig. 5, Table 4.
```

```quote
Both computations agree very well, as supported by the ξ² values in Table 4 and shown in Fig. 5a. The ξ² values are almost the same for both integration methods.
```

##### quotation-12

```tags
[claim-type]: result
[tags]: double_folding_potential
[source]: Fig. 6.
```

```quote
The full circles, up triangles, and down triangles are the calculations of the total, direct and exchange parts of the double folding potentials using DFPOT [7], respectively. The solid lines are the computations using BiFold.
```

##### quotation-13

```tags
[claim-type]: result
[tags]: double_folding_potential
[source]: Fig. 7.
```

```quote
The full circles, up triangles, and down triangles are the calculations of the total, direct and exchange parts of the double folding potentials from Khoa et al. [35], respectively. The solid lines are the computations using BiFold.
```




### Meanings

##### quotation-14

```tags
[claim-type]: innovation
[tags]: m3y_interaction
[source]: section 1.
```

```quote
For example, DFPOT [7] cannot calculate the potentials with the DDM3Y, BDM3Yn [...] and CDM3Yn [...] type density-dependent interactions [3]. While DFMSPH [8] can calculate many of those, it does not support BDM3Y2, BDM3Y3 [9], and the first version of DDM3Y [2]. [...] Therefore, the code BiFold will help the community in these manners.
```

##### quotation-15

```tags
[claim-type]: comparison
[tags]: mean_relative_error, double_folding_potential
[source]: section 4.3, Fig. 7, Table 4.
```

```quote
As can be seen from the ξ² values in Table 4 and the potentials in Fig. 7, BiFold is also consistent with the final reference work.
```





### Secondary Citations

##### quotation-16

```tags
[claim-type]: background
[tags]: double_folding_potential
[source]: section 1.
[ref]: [1] G.R. Satchler and W.G. Love, Phys. Rep. 55 (1979) 183.
[doi]: 10.1016/0370-1573(79)90081-4
```

```quote
The folded potential model [1] is a well-known model for describing the mean-field nuclear interaction between two colliding nuclei.
```

##### quotation-17

```tags
[claim-type]: background
[tags]: frozen_density_approximation
[source]: section 1.
[ref]: [1] G.R. Satchler and W.G. Love, Phys. Rep. 55 (1979) 183; [2] A. Kobos et al., Nucl. Phys. A 384 (1982) 65–87; [3] D.T. Khoa et al., J. Phys. G 34 (2007) R111–R164.
[doi]: [1] 10.1016/0370-1573(79)90081-4; [2] 10.1016/0375-9474(82)90305-0; [3] 10.1088/0954-3899/34/3/R01
```

```quote
The DD of NN interaction has several treatments in the literature, but the code BiFold is built on the treatments of Satchler and Love [1], Kobos et al. [2], and Khoa et al. [3]. All these three treatments have the frozen density approximation [...]
```

##### quotation-18

```tags
[claim-type]: background
[tags]: m3y_interaction
[source]: section 2.2.
[ref]: [12] R.V. Reid, Ann. Phys. 50 (1968) 411–448; [13] M. Lacombe et al., Phys. Rev. C 21 (1980) 861–873; [14] A. Anantaraman et al., Nucl. Phys. A 398 (1983) 269.
[doi]: [12] 10.1016/0003-4916(68)90126-7; [13] 10.1103/PhysRevC.21.861; [14] 10.1016/0375-9474(83)90487-6
```

```quote
These interactions are called M3Y-Reid [12] and M3Y-Paris [13] effective interactions. The former is derived from the solution of the Bethe-Goldstone equation with Reid [12] soft-core interaction [...] to obtain G-matrix. The latter is derived from a more fundamental Paris NN potential [13] to generate all components of the effective interaction [14].
```

##### quotation-19

```tags
[claim-type]: definition
[tags]: mean_relative_error, double_folding_potential
[source]: section 4, Eq. (24).
[ref]: [8] Gontchar et al., Comput. Phys. Commun. 181 (2010) 168.
[doi]: 10.1016/j.cpc.2009.09.007
```

```quote
This formula defines a mean relative error (mre) [8] for comparing BiFold calculations with reference calculations where U A and U B are the reference and the present double folding potentials with the radial distance Ri, respectively.
```


### Gaps

- `[section]: Results | [item]: coverage | [target]: Figure 1 | [reason]: method schematic; no current-work output quantity`
- `[section]: Results | [item]: coverage | [target]: Table 1 | [reason]: interaction support table; no current-work output quantity`
- `[section]: Results | [item]: coverage | [target]: Figure 2 | [reason]: code tree; no current-work output quantity`
- `[section]: Results | [item]: coverage | [target]: Table 2 | [reason]: function inventory; no current-work output quantity`
- `[section]: Results | [item]: coverage | [target]: Table 3 | [reason]: density-function inventory; no current-work output quantity`
