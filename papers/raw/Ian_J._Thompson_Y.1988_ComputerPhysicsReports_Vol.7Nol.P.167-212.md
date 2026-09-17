---
[source_json]: Ian_J._Thompson_Y.1988_ComputerPhysicsReports_Vol.7Nol.P.167-212.json
[doi]: 10.1016/0167-7977(88)90005-6
[first_author]: Ian J., Thompson
[affiliations]: Department of Engineering Mathematics, Bristol University, Bristol, BS8 1TR, UK, Department of Engineering Mathematics, Bristol University, Bristol BS8 1TR, UK
---

# Coupled reaction channels calculations in nuclear physics

### Abstract

This paper describes the components and methods of a comprehensive code for coupled reaction channels calculations in nuclear physics. Procedures are described which are common to the modelling of reactions induced by light and medium-mass ions, and which are sufficient to calculate the effects of successive processes to any order.




### Motivation

##### quotation-01

```tags
[claim-type]: background
[tags]: nuclear_reaction
[source]: section 1, p. 170.
```

> Inelastic excitations may occur, for example when one or both of the nuclei are deformed or deformable, with the result that higher-energy states of the nuclei may become populated.

##### quotation-02

```tags
[claim-type]: motivation
[tags]: coupled_channels_method
[source]: section 1, p. 170.
```

> One-step processes have been traditionally described with the Distorted Wave Born Approximation (DWBA), and although second-order DWBA expressions can be written down and computed, I shall be mainly concerned with coupled-channels formalisms, in order to predict the effects of multi-step processes to any or all orders.




### Methods

##### quotation-03

```tags
[claim-type]: definition
[tags]: coupled_channels_method
[source]: section 2, p. 171, equation (2.1).
```

> The coupled reaction channels (CRC) model of direct reactions in nuclear physics proceeds by constructing a model of the system wave function, and solving Schrödinger’s equation as accurately as possible within that model space. The model used here projects the complete wave function $\overline{\Psi}$ onto a product $\phi_i\equiv\phi_{ip}\times\phi_{it}$ of projectile and target states with a wave function $\psi_i(\mathbf R_i)$ describing their relative motion:

```math
P\overline{\Psi}\equiv\Psi
=\sum_{i=1}^{N}\phi_i\psi_i(\mathbf R_i),
\qquad
\phi_i\equiv\phi_{ip}\times\phi_{it}.
\tag{2.1}
```

##### quotation-04

```tags
[claim-type]: method
[tags]: coupled_channels_method
[source]: section 2, p. 172, equation (2.2).
```

> For a complete Hamiltonian $\overline{\mathcal H}$ and total energy $E$, Schrödinger’s equation $[\overline{\mathcal H}-E]\overline{\Psi}=0$ becomes $[\mathcal H-E]\Psi=0$ in the model space with [14]

```math
\mathcal H
=
P\overline{\mathcal H}P
-
P\overline{\mathcal H}Q
\frac{1}{Q\overline{\mathcal H}Q-E-i\epsilon}
Q\overline{\mathcal H}P,
\qquad
Q\equiv1-P.
\tag{2.2}
```

##### quotation-05

```tags
[claim-type]: definition
[tags]: coupled_channels_method
[source]: section 2, p. 172, equation (2.3).
```

> If $E_i$ is the asymptotic kinetic energy in the $i$th channel, then the channel-projected Hamiltonian $H_i$ satisfies [...] and will be composed of a kinetic energy term and a diagonal optical potential.

```math
H_i-E_i
=
\langle\phi_i|\mathcal H-E|\phi_i\rangle.
\tag{2.3}
```

##### quotation-06

```tags
[claim-type]: definition
[tags]: coupled_channels_method, coupling_interaction
[source]: section 2, p. 172, equation (2.4) and following sentence.
```

> The ‘interaction potential’ $V_i$ is then defined to be everything in $\mathcal H$ not included in $H_i$, so [...] This construction gives $V_i$ which have vanishing diagonal matrix elements $\langle\phi_i|V_i|\phi_i\rangle=0$.

```math
H_i-E_i+V_i=\mathcal H-E.
\tag{2.4}
```

##### quotation-07

```tags
[claim-type]: method
[tags]: coupled_channels_method
[source]: section 2.1, p. 172, equation (2.5).
```

> If we take the model Schrödinger’s equation $[\mathcal H-E]\Psi=0$, and project separately onto the different basis states $\phi_i$, we derive the set of equations [...] which couple together the unknown wave functions $\psi_i(\mathbf R_i)$.

```math
[E_i-H_i]\psi_i(\mathbf R_i)
=
\sum_{j\ne i}
\langle\phi_i|\mathcal H-E|\phi_j\rangle
\psi_j(\mathbf R_j).
\tag{2.5}
```

##### quotation-08

```tags
[claim-type]: definition
[tags]: coupled_channels_method, coupling_interaction
[source]: section 2.1, pp. 172–173, equation (2.7) and following definitions.
```

> The matrix element $\langle\phi_i|\mathcal H-E|\phi_j\rangle$ has two different forms, depending on whether we expand [...] where $V_{ij}^{post}\equiv\langle\phi_i|V_i|\phi_j\rangle$, $V_{ij}^{prior}\equiv\langle\phi_i|V_j|\phi_j\rangle$ and $K_{ij}\equiv\langle\phi_i|\phi_j\rangle$.

```math
\begin{aligned}
\langle\phi_i|\mathcal H-E|\phi_j\rangle
&=V_{ij}^{\mathrm{post}}+[H_i-E_i]K_{ij}
&&\text{(post)}\\
&=V_{ij}^{\mathrm{prior}}+K_{ij}[H_j-E_j]
&&\text{(prior)}.
\end{aligned}
\tag{2.7}
```

##### quotation-09

```tags
[claim-type]: method
[tags]: coupled_channels_method, coupling_interaction
[source]: section 2.2, pp. 173–174, equations (2.8)–(2.12), initial-channel equation, and definition of chi_N.
```

> This approximation always feeds flux ‘forwards’ in the sequence $1\to2\to\cdots\to N$ neglecting the back couplings. [...] If the post interaction had been used for all the couplings to this last channel, then

```math
[E_i-H_i]\psi_i(\mathbf R_i)
=
\sum_{j=1}^{i-1}
\langle\phi_i|\mathcal H-E|\phi_j\rangle
\psi_j(\mathbf R_j).
\tag{2.8}
```

Editorial note: In equation (2.8), the printed \(H_j\) on the left-hand side is corrected to \(H_i\), consistent with the channel index and the following equations.

```math
[E_1-H_1]\psi_1(\mathbf R_1)=0.
```

```math
[E_2-H_2]\psi_2(\mathbf R_2)
=
\langle\phi_2|\mathcal H-E|\phi_1\rangle
\psi_1(\mathbf R_1).
\tag{2.9}
```

```math
[E_N-H_N]\psi_N(\mathbf R_N)
=
\sum_{j=1}^{N-1}
\langle\phi_N|\mathcal H-E|\phi_j\rangle
\psi_j(\mathbf R_j).
\tag{2.10}
```

```math
[E_N-H_N]\psi_N
=
\sum_{j=1}^{N-1}
\langle\phi_N|V_N|\phi_j\rangle\psi_j
+
[H_N-E_N]
\sum_{j=1}^{N-1}
\langle\phi_N|\phi_j\rangle\psi_j.
\tag{2.11}
```

```math
[E_N-H_N]\chi_N(\mathbf R_N)
=
\sum_{j=1}^{N-1}V_{Nj}^{\mathrm{post}}\psi_j.
\tag{2.12}
```

```math
\chi_N(\mathbf R_N)
=
\psi_N+
\sum_{j=1}^{N-1}
\langle\phi_N|\phi_j\rangle\psi_j
\equiv
\langle\phi_N|\Psi\rangle.
```

##### quotation-10

```tags
[claim-type]: method
[tags]: coupled_channels_method, coupling_interaction
[source]: section 2.3, pp. 174–175, equations (2.13)–(2.17).
```

> There are schemes available which can iterate all channels with an arbitrary choice of post or prior interactions for all the couplings. [...] The following iterative scheme [18] ($n=1,2,\ldots$) on convergence then solves the CRC equations (2.5):

```math
\theta_{ij}=
\begin{cases}
1,&\text{post on the }j\to i\text{ coupling},\\
0,&\text{prior on the }j\to i\text{ coupling}.
\end{cases}
\tag{2.13}
```

```math
\psi_i^{(0)}
=
\delta(i,i_0)\psi_{\mathrm{elastic}},
\qquad
\Delta S_i^{(0)}=\Delta\psi_i^{(0)}=0.
\tag{2.14}
```

```math
[H_i-E_i]\chi_i^{(n)}+S_i^{(n-1)}=0.
\tag{2.15}
```

```math
S_i^{(n-1)}
=
\sum_j
\left[
\theta_{ij}V_{ij}^{\mathrm{post}}
+
(1-\theta_{ij})V_{ij}^{\mathrm{prior}}
\right]\psi_j^{(n-1)}
-
\Delta S_i^{(n-1)}.
\tag{2.16}
```

```math
\begin{aligned}
\Delta\psi_i^{(n)}
&=
\sum_j\theta_{ij}
\langle\phi_i|\phi_j\rangle
\psi_j^{(n-1)},\\
\Delta S_i^{(n)}
&=
\sum_j(1-\theta_{ij})
\langle\phi_i|\phi_j\rangle
\left[
S_j^{(n-1)}
+
[H_j-E_j]\Delta\psi_j^{(n)}
\right],\\
\psi_i^{(n)}
&=
\chi_i^{(n)}-\Delta\psi_i^{(n)}.
\end{aligned}
\tag{2.17}
```




### Results

...




### Meanings

...





### Secondary Citations

##### quotation-xx

```tags
[claim-type]: ...
[tags]: ...
[source]: ...
[ref]: [ab] reference index infromation from the original paper
[doi]: doi number of the reference or unknown.
```

> This is a quotation with secondary citations from the original paper[ab], which is tagged with the above tags.


### Gaps

- `[section]: [section name] | [item]: [claim-type|coverage|quotation] | [target]: [exact requirement or expanded target] | [reason]: [why no valid evidence exists, up to 10 words]`
