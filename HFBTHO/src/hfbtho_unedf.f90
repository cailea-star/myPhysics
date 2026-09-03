!***********************************************************************
!
!    Copyright (c) 2021, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Lead developer: Nicolas Schunck, schunck1@llnl.gov
!    HFBTHO
!    -----
!      LLNL-CODE-826901 All rights reserved.
!      LLNL-CODE-728299 All rights reserved.
!      LLNL-CODE-573953 All rights reserved.
!
!      Copyright 2021, P. Marevic, N. Schunck, E. Ney, R. Navarro Perez,
!                      M. Verriere, J. O'Neal
!      Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang,
!                      J. Sarich
!      Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                      N. Michel, J. Sarich, S. Wild
!      Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!
!    This file is part of HFBTHO.
!
!    HFBTHO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HFBTHO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************

! ==================================================================== !
!                                                                      !
!                             UNEDF PACKAGE                            !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides the main \theCode DFT solver. It includes
!> routines for the calculation and diagonalization of the HFB matrix;
!> the definition of densities on the quadrature mesh and in configuration
!> space; the self-consistent loop the calculation of expectation values of
!> observables. It also includes a package of a few routines to perform
!> particle number projection.
!>
!>  @author
!>    Markus Kortelainen, Mario Stoitsov, Rodrigo Navarro Perez, Nicolas Schunck
!----------------------------------------------------------------------
Module UNEDF

  Use HFBTHO_utilities

  Implicit None

  Character(16), Private :: Version='17'

  Public calculate_U_parameters,default_UNEDF_NAMELIST,read_UNEDF_NAMELIST,set_all_gaussians, &
         set_functional_parameters,print_functional_parameters,t_from_C,Vexternal

  Private gogny_force,skforce,calculate_natural_units,calculate_C_from_NM,calculate_finite_range_NM,P_SNM_FR, &
          calculate_NM_properties,find_NM_RHOC,C_from_t

  Logical, Public :: use_charge_density, use_cm_cor,use_DME3N_terms,   &
                     use_j2terms,use_full_cm_cor,use_INM,use_Namelist, &
                     Print_Namelist,finite_range,hb0_charge_dependent, &
                     force_is_DME,use_3N_couplings,override_3N_couplings, &
                     coulomb_gaussian, TDDFT_compatibility
  Logical, Public :: harmonic_trap = .True.
  Integer(ipr), Public :: DMEorder,DMElda,use_TMR_pairing
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorho,Urhotau,UrhoDrho,Unablarho  ! ph amplitudes
  Real(pr), Public, Dimension(0:3,0:7) :: UJnablarho,UrhonablaJ,UJJ,UJabJba
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorhopr                           ! pp amplitudes
  Real(pr), Public, Dimension(0:1) :: UEnonstdr,UFnonstdr,URnonstdr           ! Other amplitudes
  Real(pr), Public :: hbzero,sigma,e2charg,CExPar                             ! hbr^2/2m, DD sigma, e^2 charge, coul.exch.
  Real(pr), Public :: hbzeron,hbzerop                                         ! hbr^2/2m_n,hbr^2/2m_p
  Real(pr), Public, Dimension(0:1) :: Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CJbar,CpV0,CpV1 ! Time-even coupling constants (and pairing)
  Real(pr), Public, Dimension(0:1) :: Cs,Cds,Csj,CsDs,CsdJ,cT,CsNabs,CF            ! Time-odd coupling constants
  Real(pr), Public, Dimension(0:1) :: Cnrho,CJdr
  Real(pr), Public :: E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM
  Real(pr), Public :: CHrho                                                   ! Crho(0) from the Hartree term in NM
  Real(pr), Public :: mass_pion,gA,fpi,c1,c3,c4,cd,ce,LambdaX
  Real(pr), Public :: t0s,t0a,drs,dra,ts,ta,t3alp,t3al0,t3alm,t324,alp,alm,wla0, &
                      wla1,ta7,ta8,tb7,tb8,tv1,tv2,tv3,tv4,tv5,tv6,ts1,ts2,t4o3
  Real(pr), Public :: t0_pub,t1_pub,t2_pub,t3_pub,x0_pub,x1_pub,x2_pub,x3_pub,b4_pub,b4p_pub,te_pub,to_pub
  ! Gogny parameters
  Logical, Private, Save :: is_gogny_force
  Real(pr), Public, Allocatable, Save :: mu_g(:), W_g(:), B_g(:), H_g(:), M_g(:), mu_g_all(:)
  Integer(ipr), Public, Save :: n_g=0, n_g_all=0
  ! Coulomb by Gaussians parameters
  Integer(ipr), Public, Parameter :: n_g_coul = 9
  Real(pr), Dimension(1:n_g_coul), Public, Save :: mu_g_coul,V_g_coul
  ! Bulgac functionals
  Logical, Public :: is_NEDF
  Real(pr), Public, Dimension(0:2) :: a_NEDF,b_NEDF,c_NEDF
  Real(pr), Public :: eta_NEDF , W0_NEDF !
  Real(pr), Public, Save :: alfa_dme

  Real(pr), Private, Dimension(0:1), Save :: nuCrho,nuCdrho,nuCtau,nuCrDr  ! basic coupling constants in natural units
  Real(pr), Private, Dimension(0:1), Save :: nuCrdJ,nuCJ,nuCpV0,nuCpV1     !
  Real(pr), Private, Save :: t0,t1,t2,t3,x0,x1,x2,x3,b4,b4p,te,to
  Real(pr), Private, Save :: nuLambda,nufpi                                ! parameters associated to natural units
  Real(pr), Private, Save :: eps                                           ! dmelda: 0=Kf-LDA, 1=CB-LDA
  Real(pr), Private, Save :: kfconst,CK                                    ! (3Pi^2/2)^(1/3)
  Real(pr), Parameter, Private :: mevfm=197.30_pr
  Real(pr), Private, Save :: rho(0:1),tau(0:1),nrho2(0:1),lrho(0:1)
  Real(pr), Private, Save :: arhorho,brhorho,arhodrho,brhodrho,arhotau,brhotau,ajj,bjj,adrdr,bdrdr
  Real(pr), Private, Save :: darhorho,dbrhorho,darhodrho,dbrhodrho,darhotau,dbrhotau,dajj,dbjj,dadrdr,dbdrdr
  Real(pr), Private, Save :: ddarhodrho,ddbrhodrho,ddarhotau,ddbrhotau,ddarhorho,ddbrhorho
  Real(pr), Private, Save :: hrho0rho0,hrho1rho1,hdr0dr0,hdr1dr1,hrho0Drho0,hrho1Drho0, &
                             hrho1Drho1,hrho0tau0,hrho1tau0,hrho1tau1,hJ0dr0,hrho0DJ0,hJ1dr1,hrho1DJ1, &
                             hJ0dr1,hrho1DJ0,hJ1dr0,hJ0J0,hJ0J1,hJ1J1
  Real(pr), Private, Save :: dhrho0rho0,dhrho1rho1,dhdr0dr0,dhdr1dr1,dhrho0Drho0, &
                             dhrho1Drho0,dhrho1Drho1,dhrho0tau0,dhrho1tau0,dhrho1tau1,dhJ0dr0,dhrho0DJ0, &
                             dhJ1dr1,dhrho1DJ1,dhJ0dr1,dhrho1DJ0,dhJ1dr0,dhJ0J0,dhJ0J1,dhJ1J1
  Real(pr), Private, Save :: ddhrho0rho0,ddhrho1rho1,ddhrho0Drho0,ddhrho1Drho0, &
                             ddhrho1Drho1,ddhrho0tau0,ddhrho1tau0,ddhrho1tau1
  Real(pr), Private :: ac2,ac3,acoord
  Parameter(acoord=0.50_pr,ac2=4.0_pr*(acoord**2-acoord+0.50_pr),ac3=2.0_pr*(acoord**2-acoord+0.50_pr))
  Character(Len=30) :: FunctionalName

  Namelist /UNEDF_NAMELIST/ FunctionalName, DMEorder, DMElda, use_INM, hbzero, use_TMR_pairing,   &
                            Crho, Cdrho, Ctau, CrDr, CrdJ, CJ, sigma, CpV0, CpV1, e2charg,        &
                            E_NM, K_NM, SMASS_NM, RHO_NM, ASS_NM, LASS_NM, VMASS_NM,              &
                            mass_pion, gA, fpi, c1, c3, c4, cd, ce, LambdaX,                            &
                            use_cm_cor, use_charge_density, use_DME3N_terms, use_j2terms, CExPar, &
                            Print_Namelist, use_INM, use_cm_cor,       &
                            finite_range, hbzeron, hbzerop, alfa_dme, force_is_dme
  Namelist /GOGNY_NAMELIST/ mu_g, W_g, B_g, H_g, M_g

Contains
  !=======================================================================
  !> This routine computes the amplitudes \f$ U^{m}_{tt'} \f$ amplitude, where m runs over the bilinear
  !> combinations \f$ \{ \rho_t\rho_{t'}, \rho_t \tau_{t'} \gras{J}_t \gras{J}_{t'}, \rho_t \Delta\rho_{t'},
  !> \rho_t \gras{\nabla}\gras{J}_{t'} \} \f$ and \f$ t,t' = 0,1 \f$; see \cite stoitsov2010microscopically.
  !> The arrays \f$ \mathtt{Urhorho(0:3,0:7)} \f$ containing the U amplitudes have have the following
  !> storing convention:
  !> * First index
  !>   - i = 0 -> (t,t') = (0,0)
  !>   - i = 1 -> (t,t') = (1,1)
  !>   - i = 2 -> (t,t') = (0,1)
  !>   - i = 3 -> (t,t') = (1,0)
  !> * Second index
  !>   - j = 0 -> U
  !>   - j = 1 -> \f$ \displaystyle\frac{\partial U}{\partial\rho_0} \f$
  !>   - j = 2 -> \f$ \displaystyle\frac{\partial U}{\partial\rho_1} \f$
  !>   - j = 3 -> \f$ \displaystyle\frac{\partial^2 U}{\partial\rho_0 \partial\rho_0} \f$
  !>   - j = 4 -> \f$ \displaystyle\frac{\partial^2 U}{\partial\rho_1 \partial\rho_1} \f$
  !>   - j = 5 -> \f$ \displaystyle\frac{\partial^2 U}{\partial\rho_0 \partial\rho_1} \f$
  !>   - j = 6 -> \f$ \displaystyle\frac{\partial U}{\partial \tau_0} \f$
  !>   - j = 7 -> \f$ \displaystyle\frac{\partial U}{\partial \partial\Delta \rho_0} \f$
  !>
  !>  @param[in] rho0_in - Isoscalar density \f$ \rho_0 \f$ at the current point
  !>  @param[in] rho1_in - Isovector density \f$ \rho_1 \f$ at the current point
  !>  @param[in] tau0_in - Isoscalar kinetic density \f$ \tau_0 \f$ at the current point
  !>  @param[in] tau1_in - Isovector kinetic density \f$ \tau_1 \f$ at the current point
  !>  @param[in] laprho0 - Laplacian of the isoscalar density \f$ \Delta\rho_0 \f$ at the current point
  !>  @param[in] laprho1 - Laplacian of the isovector density \f$ \Delta\rho_1 \f$ at the current point
  !>  @param[in] nablarho0s - Gradient of the isoscalar density \f$ \gras{\nabla}\rho_0 \f$ at the current point
  !>  @param[in] nablarho1s - Gradient of the isovector density \f$ \gras{\nabla}\rho_1 \f$ at the current point
  !=======================================================================
  Subroutine calculate_U_parameters(rho0_in,rho1_in,tau0_in,tau1_in,laprho0,laprho1,nablarho0s,nablarho1s)
    Real(pr), Intent(In) :: rho0_in,rho1_in,tau0_in,tau1_in
    Real(pr), Intent(In), Optional :: nablarho0s,nablarho1s,laprho0,laprho1

    Integer(ipr) :: t
    Real(pr) :: u,du,ddu,dtu,dlu
    Real(pr) :: ph,ucut,ucut3n

    ucut=0.1_pr; ucut3n=0.6_pr

    rho(0)=rho0_in; rho(1)=rho1_in;
    tau(0)=tau0_in; tau(1)=tau1_in;

    lrho=0.0_pr; nrho2=0.0_pr;
    If(Present(laprho0)) lrho(0)=laprho0
    If(Present(laprho1)) lrho(1)=laprho1
    If(Present(nablarho0s)) nrho2(0)=nablarho0s
    If(Present(nablarho1s)) nrho2(1)=nablarho1s

    arhorho  = 0.0_pr; darhorho  = 0.0_pr; ddarhorho  = 0.0_pr
    brhorho  = 0.0_pr; dbrhorho  = 0.0_pr; ddbrhorho  = 0.0_pr
    arhodrho = 0.0_pr; darhodrho = 0.0_pr; ddarhodrho = 0.0_pr
    brhodrho = 0.0_pr; dbrhodrho = 0.0_pr; ddbrhodrho = 0.0_pr
    arhotau  = 0.0_pr; darhotau  = 0.0_pr; ddarhotau  = 0.0_pr
    brhotau  = 0.0_pr; dbrhotau  = 0.0_pr; ddbrhotau  = 0.0_pr
    adrdr = 0.0_pr; dadrdr = 0.0_pr
    bdrdr = 0.0_pr; dbdrdr = 0.0_pr
    ajj = 0.0_pr; dajj = 0.0_pr
    bjj = 0.0_pr; dbjj = 0.0_pr

    hrho0rho0=0.0_pr; hrho1rho1=0.0_pr; hdr0dr0=0.0_pr; hdr1dr1=0.0_pr
    hrho0Drho0=0.0_pr; hrho1Drho0=0.0_pr; hrho1Drho1=0.0_pr
    hrho0tau0=0.0_pr; hrho1tau0=0.0_pr; hrho1tau1=0.0_pr
    hJ0dr0=0.0_pr; hrho0DJ0=0.0_pr; hJ1dr1=0.0_pr; hrho1DJ1=0.0_pr
    hJ0dr1=0.0_pr; hrho1DJ0=0.0_pr; hJ1dr0=0.0_pr
    hJ0J0=0.0_pr; hJ0J1=0.0_pr; hJ1J1=0.0_pr
    dhrho0rho0=0.0_pr; dhrho1rho1=0.0_pr; dhdr0dr0=0.0_pr; dhdr1dr1=0.0_pr
    dhrho0Drho0=0.0_pr; dhrho1Drho0=0.0_pr; dhrho1Drho1=0.0_pr
    dhrho0tau0=0.0_pr; dhrho1tau0=0.0_pr; dhrho1tau1=0.0_pr
    dhJ0dr0=0.0_pr; dhrho0DJ0=0.0_pr; dhJ1dr1=0.0_pr; dhrho1DJ1=0.0_pr
    dhJ0dr1=0.0_pr; dhrho1DJ0=0.0_pr; dhJ1dr0=0.0_pr
    dhJ0J0=0.0_pr; dhJ0J1=0.0_pr; dhJ1J1=0.0_pr
    ddhrho0rho0=0.0_pr; ddhrho1rho1=0.0_pr
    ddhrho0Drho0=0.0_pr; ddhrho1Drho0=0.0_pr; ddhrho1Drho1=0.0_pr
    ddhrho0tau0=0.0_pr; ddhrho1tau0=0.0_pr; ddhrho1tau1=0.0_pr

    u=0.0_pr; du=0.0_pr; ddu=0.0_pr; dtu=0.0_pr; dlu=0.0_pr

    Urhorho   = 0.0_pr; Urhotau    = 0.0_pr; UrhoDrho = 0.0_pr; Unablarho = 0.0_pr
    UJnablarho= 0.0_pr; UrhonablaJ = 0.0_pr
    Urhorhopr = 0.0_pr; UJJ        = 0.0_pr; UJabJba  = 0.0_pr
    UEnonstdr = 0.0_pr; UFnonstdr = 0.0_pr;  URnonstdr= 0.0_pr

    ! Notations for Uamplitudes(0:3,0:7)
    ! t for Uamplitudes(t,*)
    ! 0 -> 0,0
    ! 1 -> 1,1
    ! 2 -> 0,1
    ! 3 -> 1,0
    ! n for Uamplitudes(*,n)
    ! 0 -> U
    ! 1 -> dU/dRHO_0
    ! 2 -> dU/dRHO_1
    ! 3 -> d2U/(dRHO_0*dRHO_0)
    ! 4 -> d2U/(dRHO_1*dRHO_1)
    ! 5 -> d2U/(dRHO_0*dRHO_1)
    ! 6 -> dU/d(TAU_0)
    ! 7 -> dU/d(Delta RHO_0)
    !
    !! 2N terms
    Do t=0,1
       ph=1.0_pr
       If(t == 1) ph=-1.0_pr
       Urhorho(t,0)    = Crho(t) + Cdrho(t)*rho(0)**sigma + 0.50_pr*(arhorho+ph*brhorho)*mevfm
       Urhotau(t,0)    = Ctau(t) + 0.50_pr*(arhotau+ph*brhotau)*mevfm        !! These two determine the
       UrhoDrho(t,0)   = Crdr(t) + ac2*0.50_pr*(arhoDrho+ph*brhoDrho)*mevfm !! effective mass (when recoupling to p-n)??
       UJJ(t,0)        = CJ(t)   + 0.50_pr*(ajj+ph*bjj)*mevfm
       Unablarho(t,0)  = Cnrho(t)+ 0.50_pr*(adrdr+ph*bdrdr)*mevfm
       UrhonablaJ(t,0) = Crdj(t)
       UJnablarho(t,0) = Cjdr(t)

       Urhorho(t,1)  = sigma*Cdrho(t)*(rho(0)**sigma)/(rho(0)+eps) + 0.50_pr*(darhorho+ph*dbrhorho)*du*mevfm
       Urhotau(t,1)  = 0.50_pr*(darhotau+ph*dbrhotau)*du*mevfm
       UrhoDrho(t,1) = ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*du*mevfm
       UJJ(t,1)      = 0.50_pr*(dajj+ph*dbjj)*du*mevfm
       Unablarho(t,1)= 0.50_pr*(dadrdr+ph*dbdrdr)*du*mevfm

       Urhorho(t,6)   = 0.50_pr*(darhorho+ph*dbrhorho)*dtu*mevfm
       Urhotau(t,6)   = 0.50_pr*(darhotau+ph*dbrhotau)*dtu*mevfm
       UrhoDrho(t,6)  = ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dtu*mevfm
       UJJ(t,6)       = 0.50_pr*(dajj+ph*dbjj)*dtu*mevfm
       Unablarho(t,6) = 0.50_pr*(dadrdr+ph*dbdrdr)*dtu*mevfm

       Urhorho(t,7)   = 0.50_pr*(darhorho+ph*dbrhorho)*dlu*mevfm
       Urhotau(t,7)   = 0.50_pr*(darhotau+ph*dbrhotau)*dlu*mevfm
       UrhoDrho(t,7)  = ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dlu*mevfm
       UJJ(t,7)       = 0.50_pr*(dajj+ph*dbjj)*dlu*mevfm
       Unablarho(t,7) = 0.50_pr*(dadrdr+ph*dbdrdr)*dlu*mevfm

       Urhorho(t,3) = sigma*(sigma-1.0_pr)*Cdrho(t)*(rho(0)**sigma)/(rho(0)**2+eps) &
                    + 0.50_pr*(darhorho+ph*dbrhorho)*ddu*mevfm + 0.50_pr*(ddarhorho+ph*ddbrhorho)*du*du*mevfm
       Urhotau(t,3) = 0.50_pr*(darhotau+ph*dbrhotau)*ddu*mevfm + 0.50_pr*(ddarhotau+ph*ddbrhotau)*du*du*mevfm
       UrhoDrho(t,3)= ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*ddu*mevfm + ac2*0.50_pr*(ddarhoDrho+ph*ddbrhoDrho)*du*du*mevfm

    End Do
    Urhorhopr(0,0) =  CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr) + CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(1,0) =  CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr) + CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(2,0) = (CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr) - CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr))*2.0_pr
    Urhorhopr(0,1) =(-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(1,1) =(-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(2,1) = 2.0_pr*(-CpV0(0)*CpV1(0)+CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr = Urhorhopr/16.0_pr

    !! 3N terms
    If (use_DME3N_terms) Then

       Urhorho(0,0)   = Urhorho(0,0)   + hrho0rho0*rho(0)
       Urhorho(1,0)   = Urhorho(1,0)   + hrho1rho1*rho(0)
       Unablarho(0,0) = Unablarho(0,0) + hdr0dr0*rho(0)
       Unablarho(1,0) = Unablarho(1,0) + hdr1dr1*rho(0)
       UrhoDrho(0,0)  = UrhoDrho(0,0)  + hrho0Drho0*rho(0)*ac3
       UrhoDrho(1,0)  = UrhoDrho(1,0)  + hrho1Drho1*rho(0)*ac3
       UrhoDrho(3,0)  = UrhoDrho(3,0)  + hrho1Drho0*rho(1)*ac3
       Urhotau(0,0)   = Urhotau(0,0)   + hrho0tau0*rho(0)
       Urhotau(1,0)   = Urhotau(1,0)   + hrho1tau1*rho(0)
       Urhotau(3,0)   = Urhotau(3,0)   + hrho1tau0*rho(1)
       UJnablarho(0,0)= UJnablarho(0,0)+ hJ0dr0*rho(0)
       UJnablarho(1,0)= UJnablarho(1,0)+ hJ1dr1*rho(0)
       UJnablarho(2,0)= UJnablarho(2,0)+ hJ0dr1*rho(1)
       UJnablarho(3,0)= UJnablarho(3,0)+ hJ1dr0*rho(1)
       UrhonablaJ(0,0)= UrhonablaJ(0,0)+ hrho0dJ0*rho(0)
       UrhonablaJ(1,0)= UrhonablaJ(1,0)+ hrho1dJ1*rho(0)
       UrhonablaJ(3,0)= UrhonablaJ(3,0)+ hrho1dJ0*rho(1)
       UJJ(0,0)       = UJJ(0,0)       + hJ0J0*rho(0)
       UJJ(1,0)       = UJJ(1,0)       + hJ1J1*rho(0)
       UJJ(2,0)       = UJJ(2,0)       + hJ0J1*rho(1)

       Urhorho(0,1)    = Urhorho(0,1)    + dhrho0rho0*du*rho(0) + hrho0rho0
       Urhorho(1,1)    = Urhorho(1,1)    + dhrho1rho1*du*rho(0) + hrho1rho1
       Unablarho(0,1)  = Unablarho(0,1)  + dhdr0dr0*du*rho(0) + hdr0dr0
       Unablarho(1,1)  = Unablarho(1,1)  + dhdr1dr1*du*rho(0) + hdr1dr1
       UrhoDrho(0,1)   = UrhoDrho(0,1)   + (dhrho0Drho0*du*rho(0) + hrho0Drho0)*ac3
       UrhoDrho(1,1)   = UrhoDrho(1,1)   + (dhrho1Drho1*du*rho(0) + hrho1Drho1)*ac3
       UrhoDrho(3,1)   = UrhoDrho(3,1)   + dhrho1Drho0*du*rho(1)*ac3
       Urhotau(0,1)    = Urhotau(0,1)    + dhrho0tau0*du*rho(0)   + hrho0tau0
       Urhotau(1,1)    = Urhotau(1,1)    + dhrho1tau1*du*rho(0)   + hrho1tau1
       Urhotau(3,1)    = Urhotau(3,1)    + dhrho1tau0*du*rho(1)
       UJnablarho(0,1) = UJnablarho(0,1) + dhJ0dr0*du*rho(0) + hJ0dr0
       UJnablarho(1,1) = UJnablarho(1,1) + dhJ1dr1*du*rho(0) + hJ1dr1
       UJnablarho(2,1) = UJnablarho(2,1) + dhJ0dr1*du*rho(1)
       UJnablarho(3,1) = UJnablarho(3,1) + dhJ1dr0*du*rho(1)
       UrhonablaJ(0,1) = UrhonablaJ(0,1) + dhrho0dJ0*du*rho(0) + hrho0dJ0
       UrhonablaJ(1,1) = UrhonablaJ(1,1) + dhrho1dJ1*du*rho(0) + hrho1dJ1
       UrhonablaJ(3,1) = UrhonablaJ(3,1) + dhrho1dJ0*du*rho(1)
       UJJ(0,1)        = UJJ(0,1)        + dhJ0J0*du*rho(0) + hJ0J0
       UJJ(1,1)        = UJJ(1,1)        + dhJ1J1*du*rho(0) + hJ1J1
       UJJ(2,1)        = UJJ(2,1)        + dhJ0J1*du*rho(1)

       Urhorho(0,6)    = Urhorho(0,6)    + dhrho0rho0*dtu*rho(0)
       Urhorho(1,6)    = Urhorho(1,6)    + dhrho1rho1*dtu*rho(0)
       Unablarho(0,6)  = Unablarho(0,6)  + dhdr0dr0*dtu*rho(0)
       Unablarho(1,6)  = Unablarho(1,6)  + dhdr1dr1*dtu*rho(0)
       UrhoDrho(0,6)   = UrhoDrho(0,6)   + dhrho0Drho0*dtu*rho(0)*ac3
       UrhoDrho(1,6)   = UrhoDrho(1,6)   + dhrho1Drho1*dtu*rho(0)*ac3
       UrhoDrho(3,6)   = UrhoDrho(3,6)   + dhrho1Drho0*dtu*rho(1)*ac3
       Urhotau(0,6)    = Urhotau(0,6)    + dhrho0tau0*dtu*rho(0)
       Urhotau(1,6)    = Urhotau(1,6)    + dhrho1tau1*dtu*rho(0)
       Urhotau(3,6)    = Urhotau(3,6)    + dhrho1tau0*dtu*rho(1)
       UJnablarho(0,6) = UJnablarho(0,6) + dhJ0dr0*dtu*rho(0)
       UJnablarho(1,6) = UJnablarho(1,6) + dhJ1dr1*dtu*rho(0)
       UJnablarho(2,6) = UJnablarho(2,6) + dhJ0dr1*dtu*rho(1)
       UJnablarho(3,6) = UJnablarho(3,6) + dhJ1dr0*dtu*rho(1)
       UrhonablaJ(0,6) = UrhonablaJ(0,6) + dhrho0dJ0*dtu*rho(0)
       UrhonablaJ(1,6) = UrhonablaJ(1,6) + dhrho1dJ1*dtu*rho(0)
       UrhonablaJ(3,6) = UrhonablaJ(3,6) + dhrho1dJ0*dtu*rho(1)
       UJJ(0,6)        = UJJ(0,6)        + dhJ0J0*dtu*rho(0)
       UJJ(1,6)        = UJJ(1,6)        + dhJ1J1*dtu*rho(0)
       UJJ(2,6)        = UJJ(2,6)        + dhJ0J1*dtu*rho(1)

       Urhorho(0,7)    = Urhorho(0,7)    + dhrho0rho0*dlu*rho(0)
       Urhorho(1,7)    = Urhorho(1,7)    + dhrho1rho1*dlu*rho(0)
       Unablarho(0,7)  = Unablarho(0,7)  + dhdr0dr0*dlu*rho(0)
       Unablarho(1,7)  = Unablarho(1,7)  + dhdr1dr1*dlu*rho(0)
       UrhoDrho(0,7)   = UrhoDrho(0,7)   + dhrho0Drho0*dlu*rho(0)*ac3
       UrhoDrho(1,7)   = UrhoDrho(1,7)   + dhrho1Drho1*dlu*rho(0)*ac3
       UrhoDrho(3,7)   = UrhoDrho(3,7)   + dhrho1Drho0*dlu*rho(1)*ac3
       Urhotau(0,7)    = Urhotau(0,7)    + dhrho0tau0*dlu*rho(0)
       Urhotau(1,7)    = Urhotau(1,7)    + dhrho1tau1*dlu*rho(0)
       Urhotau(3,7)    = Urhotau(3,7)    + dhrho1tau0*dlu*rho(1)
       UJnablarho(0,7) = UJnablarho(0,7) + dhJ0dr0*dlu*rho(0)
       UJnablarho(1,7) = UJnablarho(1,7) + dhJ1dr1*dlu*rho(0)
       UJnablarho(2,7) = UJnablarho(2,7) + dhJ0dr1*dlu*rho(1)
       UJnablarho(3,7) = UJnablarho(3,7) + dhJ1dr0*dlu*rho(1)
       UrhonablaJ(0,7) = UrhonablaJ(0,7) + dhrho0dJ0*dlu*rho(0)
       UrhonablaJ(1,7) = UrhonablaJ(1,7) + dhrho1dJ1*dlu*rho(0)
       UrhonablaJ(3,7) = UrhonablaJ(3,7) + dhrho1dJ0*dlu*rho(1)
       UJJ(0,7)        = UJJ(0,7)        + dhJ0J0*dlu*rho(0)
       UJJ(1,7)        = UJJ(1,7)        + dhJ1J1*dlu*rho(0)
       UJJ(2,7)        = UJJ(2,7)        + dhJ0J1*dlu*rho(1)

       Urhorho(0,3)  = Urhorho(0,3)  + 2.0_pr*dhrho0rho0*du + ddhrho0rho0*du*du*rho(0) + dhrho0rho0*ddu*rho(0)
       Urhorho(1,3)  = Urhorho(1,3)  + 2.0_pr*dhrho1rho1*du + ddhrho1rho1*du*du*rho(0) + dhrho1rho1*ddu*rho(0)
       UrhoDrho(0,3) = UrhoDrho(0,3) +(2.0_pr*dhrho0Drho0*du + ddhrho0Drho0*du*du*rho(0) + dhrho0Drho0*ddu*rho(0))*ac3
       UrhoDrho(1,3) = UrhoDrho(1,3) +(2.0_pr*dhrho1Drho1*du + ddhrho1Drho1*du*du*rho(0) + dhrho1Drho1*ddu*rho(0))*ac3
       Urhotau(0,3)  = Urhotau(0,3)  + 2.0_pr*dhrho0tau0*du + ddhrho0tau0*du*du*rho(0) + dhrho0tau0*ddu*rho(0)
       Urhotau(1,3)  = Urhotau(1,3)  + 2.0_pr*dhrho1tau1*du + ddhrho1tau1*du*du*rho(0) + dhrho1tau1*ddu*rho(0)
       UrhoDrho(3,3) = UrhoDrho(3,3) +(ddhrho1Drho0*du*du*rho(1) + dhrho1Drho0*ddu*rho(1))*ac3
       Urhotau(3,3)  = Urhotau(3,3)  + ddhrho1tau0*du*du*rho(1) + dhrho1tau0*ddu*rho(1)

       UrhoDrho(3,2)   = UrhoDrho(3,2)   + hrho1Drho0*ac3
       Urhotau(3,2)    = Urhotau(3,2)    + hrho1tau0
       UJnablarho(3,2) = UJnablarho(3,2) + hJ1dr0
       UJnablarho(2,2) = UJnablarho(2,2) + hJ0dr1
       UrhonablaJ(3,2) = UrhonablaJ(3,2) + hrho1dJ0
       UJJ(2,2)        = UJJ(2,2)        + hJ0J1

       UrhoDrho(3,5)=UrhoDrho(3,5) +dhrho1Drho0*du*ac3
       Urhotau(3,5)=Urhotau(3,5) +dhrho1tau0*du

    End If   !! if(use3Nterm..

    UEnonstdr=0.0_pr; UFnonstdr=0.0_pr; URnonstdr=0.0_pr

    If(is_nedf) Then
       Do t = 0,1
          Urhorho(t,0) =   a_NEDF(t)/(rho(0)**(1/three)+eps) + b_NEDF(t) + c_NEDF(t)*rho(0)**( 1/three)
          Urhorho(t,1) = (-a_NEDF(t)/(rho(0)**(4/three)+eps)             + c_NEDF(t)/(rho(0)**(2/three)+eps))/three
       End Do
       Urhorho(2,0) =  ( a_NEDF(2)/(rho(0)**(10/three)+eps) + &
                         b_NEDF(2)/(rho(0)**(3)+eps) + &
                         c_NEDF(2)/(rho(0)**( 8/three)+eps))*rho(1)**3
       Urhorho(2,1) =  (-10*a_NEDF(2)/(rho(0)**(13/three)+eps) - &
                          9*b_NEDF(2)/(rho(0)**(4)+eps) - &
                          8*c_NEDF(2)/(rho(0)**(11/three)+eps))*rho(1)**3/three
       Urhorho(2,2) =3*( a_NEDF(2)/(rho(0)**(10/three)+eps) + &
                         b_NEDF(2)/(rho(0)**(3)+eps) + &
                         c_NEDF(2)/(rho(0)**( 8/three)+eps))*rho(1)**2
    End If

    If(.Not.use_j2terms) Then
       UJJ= zero
       UJabJba = zero
    End If

  End Subroutine calculate_U_parameters
  !=======================================================================
  !> This routine presets default values for Skyrme-like functionals.
  !=======================================================================
  Subroutine default_UNEDF_NAMELIST()

    ! Parameters
    eps     = Spacing(one)
    kfconst = (1.5_pr*PI**2)**(one/three)
    CK      = three/five*kfconst**2
    ! Some default values for all cases
    Print_Namelist=.False.
    FunctionalName="Bla-Bla"
    ! Kind of the functional
    use_INM                = .False.
    use_DME3N_terms        = .False.
    use_charge_density     = .False.
    use_cm_cor             = .False.
    use_full_cm_cor        = .False.
    use_j2terms            = .False.
    hb0_charge_dependent   = .False.
    finite_range           = .False.
    force_is_DME           = .False.
    override_3N_couplings  = .False.
    alfa_dme               =  zero
    is_nedf                = .False.
    TDDFT_compatibility    = .False.
    coulomb_gaussian       = .False.
    use_TMR_pairing        =  0
    DMEorder               = -1
    DMElda                 =  0
    ! Coupling constants: ph channel (UNEDF1 values)
    Crho(0)=   -779.373008720865300_pr; Crho(1)  =   287.7221315832867960_pr
    CDrho(0)=   891.477890442349690_pr; CDrho(1) =  -200.5877743178848790_pr
    Ctau(0)=  -0.989915057807676746_pr; Ctau(1)  =   -33.6320970701835549_pr
    CrDr(0)=   -45.1351310222373030_pr; CrDr(1)  =  -145.3821679080569990_pr
    CrdJ(0)=   -74.0263331764599002_pr; CrdJ(1)  =   -35.6582611147916992_pr
    CJ(0)=      0.00000000000000000_pr; CJ(1)    =     0.0000000000000000_pr
    ! Coupling constants: pp channel (UNEDF1 values)
    CpV0(0)=   -186.065399575123990_pr; CpV0(1)  =  -206.5795938901059970_pr
    CpV1(0)=   0.500000000000000000_pr; CpV1(1)  =     0.5000000000000000_pr
    ! Various
    Cnrho    =    0.000000000000000_pr; CJdr     =     0.0000000000000000_pr
    sigma    =    0.2700180115027076_pr;
    hbzero   =   20.7355300000000007_pr;
    e2charg  =    1.4399784085965135_pr ; CExPar = 1.0_pr
    ! DME
    mass_pion=   138.03_pr/197.3_pr; fpi = 92.4_pr/197.3_pr; gA = 1.29_pr
    c1 =    -0.81_pr/1000.0_pr*197.3_pr
    c3 =    -3.40_pr/1000.0_pr*197.3_pr
    c4 =     3.40_pr/1000.0_pr*197.3_pr
    cd = -2062.00_pr/1000.0_pr
    ce =  -625.00_pr/1000.0_pr
    ! Natural units
    LambdaX  = 700.0_pr/197.3_pr; nuLambda = 700.0_pr; nufpi = 93.0_pr
    ! Nuclear matter (UNEDF1 values)
    E_NM     = -15.80000000000000000_pr; RHO_NM   =  0.158706769332587039_pr
    K_NM     = 220.00000000000002800_pr; SMASS_NM =  0.992423332283363990_pr
    ASS_NM   =  28.98678905777210700_pr; LASS_NM  = 40.004790480413632300_pr
    VMASS_NM =   1.24983857423226996_pr
    ! NEDF
    a_NEDF   = zero
    b_NEDF   = zero
    c_NEDF   = zero
    eta_NEDF = zero
    W0_NEDF  = zero

  End Subroutine default_UNEDF_NAMELIST
  !===============================================================================================
  ! This routine defines the coupling constants for the energy functional.
  !===============================================================================================
  Subroutine read_UNEDF_NAMELIST(fname,noForces,filename_nml)
    Character(Len=30), Intent(Inout) :: fname
    Integer(ipr), Intent(Out) :: noForces
    Character(Len=256),  Intent(In), Optional :: filename_nml

    Character(Len=30) :: inforcedname
    Integer(ipr) :: ios,lnamelist=16
    ! TODO: Make this length a program-wide constant
    Character(Len=256) :: filename

    If(present(filename_nml)) Then
       filename = filename_nml
    Else
       filename = 'UNEDF_NAMELIST.dat'
    End If

    ! parameters
    eps     = Spacing(1.0_pr)
    kfconst =(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK      = 3.0_pr/5.0_pr*kfconst**2
    n_g     = 0 ! initial number of Gaussians

    use_Namelist=.True.
    Do
       !---------------------------------------------------------------------
       ! Some default values for all cases
       !---------------------------------------------------------------------
       Call default_UNEDF_NAMELIST()
       !---------------------------------------------------------------------
       ! Select the functional: start with interaction
       !---------------------------------------------------------------------
       noForces=0 ! No forces to start with
       Call skforce(fname,noForces)

       If (noForces == 1) Then
           inforcedname='FORCE'
           use_Namelist=.False.
       Else
           FUNCTIONAL: Select Case (Trim(fname))
           Case ('READ')
              inforcedname='READ'
              use_Namelist=.False.
              ! Set to the same
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.599763_pr
              CrDr(1)  = -143.935074_pr
              CpV0(0)  = -234.380010_pr
              CpV0(1)  = -260.437001_pr
              CrdJ(0)  =  -73.946388_pr
              CrdJ(1)  =  -51.912548_pr
              CJ(0)    =    0.000000_pr
              CJ(1)    =    0.000000_pr
              CExPar   =    1.000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156247_pr
              E_NM     =  -15.800000_pr
              P_NM     =    0.000000_pr
              K_NM     =  244.668204_pr
              ASS_NM   =   28.668204_pr
              LASS_NM  =   40.109081_pr
              SMASS_NM =    1.067970_pr
              VMASS_NM =    1.249838574232270_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('RAND')
              inforcedname='RAND'
              use_Namelist=.False.
              ! Set to the same
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.599763_pr
              CrDr(1)  = -143.935074_pr
              CpV0(0)  = -234.380010_pr
              CpV0(1)  = -260.437001_pr
              CrdJ(0)  =  -73.946388_pr
              CrdJ(1)  =  -51.912548_pr
              CJ(0)    =    0.000000_pr
              CJ(1)    =    0.000000_pr
              CExPar   =    1.000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156247_pr
              E_NM     =  -15.800000_pr
              P_NM     =    0.000000_pr
              K_NM     =  244.668204_pr
              ASS_NM   =   28.668204_pr
              LASS_NM  =   40.109081_pr
              SMASS_NM =    1.067970_pr
              VMASS_NM =    1.249838574232270_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('FITS')
              inforcedname='FITS'
              use_Namelist=.True.
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              P_NM     =    0.000000000000000_pr
              VMASS_NM =    1.249838000000000_pr
              !
              force_is_dme   = .False.
              finite_range   = .False.
              use_cm_cor     = .False.
              use_j2terms    = .True.
              CExPar         = 1.0000000000000000_pr
              use_INM        = .True.
              VMASS_NM       =    0.9916634575083850_pr
              RHO_NM         =    0.1545314206232500_pr
              E_NM           =  -14.9428151261574005_pr
              K_NM           =  254.6493849869999906_pr
              ASS_NM         =   28.2555644286320984_pr
              LASS_NM        =   64.9345887215435056_pr
              SMASS_NM       =    0.8026870082347940_pr
              CrDr(0)        =  -40.1069828003495985_pr
              CrDr(1)        =  -35.8476436441168005_pr
              CpV0(0)        = -233.9624836153689955_pr
              CpV0(1)        = -284.1118846887910081_pr
              CrdJ(0)        =  -74.0943831060370002_pr
              CrdJ(1)        =   10.3610837684245993_pr
              CJ(0)          =  -81.3789370214482943_pr
              CJ(1)          =  -66.9270694667199990_pr
           Case ('UNE0')
              inforcedname='UNE0'
              use_Namelist=.False.
              ! kind of the functional
              use_INM    = .True.
              use_cm_cor = .True.
              ! Surface coefficients
              CrDr(0)  =  -55.260600000000000_pr
              CrDr(1)  =  -55.622600000000000_pr
              CpV0(0)  = -170.374000000000000_pr
              CpV0(1)  = -199.202000000000000_pr
              CrdJ(0)  =  -79.530800000000000_pr
              CrdJ(1)  =   45.630200000000000_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.160526000000000_pr
              E_NM     =  -16.055900000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  230.000000000000000_pr
              ASS_NM   =   30.542900000000000_pr
              LASS_NM  =   45.080400000000000_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838000000000_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('UNE1')
              inforcedname='UNE1'
              use_Namelist=.False.
              ! kind of the functional
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.135131022237300_pr
              CrDr(1)  = -145.382167908057000_pr
              CpV0(0)  = -186.065399575124000_pr
              CpV0(1)  = -206.579593890106000_pr
              CrdJ(0)  =  -74.026333176459900_pr
              CrdJ(1)  =  -35.658261114791700_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.158706769332587_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  220.000000000000000_pr
              ASS_NM   =   28.986789057772100_pr
              LASS_NM  =   40.004790480413600_pr
              SMASS_NM =    0.992423332283364_pr
              VMASS_NM =    1.249838574232270_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('UNE2')
              inforcedname='UNE2'
              use_Namelist=.False.
              ! kind of the functional
              use_INM     = .True.
              use_j2terms = .True.
              ! Surface coefficients
              CrDr(0)  =  -46.831409147060600_pr
              CrDr(1)  = -113.163790795259000_pr
              CpV0(0)  = -208.889001962571000_pr
              CpV0(1)  = -230.329984038628000_pr
              CrdJ(0)  =  -64.308862415783800_pr
              CrdJ(1)  =  -38.650194685135500_pr
              CJ(0)    =  -54.433363597372100_pr
              CJ(1)    =  -65.903031044593800_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156310622197074_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  239.929568022437000_pr
              ASS_NM   =   29.131006470773700_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    1.073763804147980_pr
              VMASS_NM =    1.249838574232270_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('HFB1')
              ! J. Phys. G: Nucl. Part. Phys. 42 (2015) 034024
              inforcedname='HFB1'
              use_Namelist=.False.
              ! kind of the functional
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.599763_pr
              CrDr(1)  = -143.935074_pr
              CpV0(0)  = -234.380010_pr
              CpV0(1)  = -260.437001_pr
              CrdJ(0)  =  -73.946388_pr
              CrdJ(1)  =  -51.912548_pr
              CJ(0)    =    0.000000_pr
              CJ(1)    =    0.000000_pr
              CExPar   =    1.000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156247_pr
              E_NM     =  -15.800000_pr
              P_NM     =    0.000000_pr
              K_NM     =  244.668204_pr
              ASS_NM   =   28.668204_pr
              LASS_NM  =   40.109081_pr
              SMASS_NM =    1.067970_pr
              VMASS_NM =    1.249838574232270_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('N0LO')
              inforcedname='N0LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 0
              ! Surface coefficients
              CrDr(0)  =  -67.437_pr
              CrDr(1)  =   21.551_pr
              CpV0(0)  = -241.203_pr
              CpV0(1)  = -252.818_pr
              CrdJ(0)  =  -95.451_pr
              CrdJ(1)  =  -65.906_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('N1LO')
              inforcedname='N1LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 1
              ! Surface coefficients
              CrDr(0)  =  -63.996_pr
              CrDr(1)  =    -9.276_pr
              CpV0(0)  = -241.484_pr
              CpV0(1)  = -252.222
              CrdJ(0)  =  -95.463_pr
              CrdJ(1)  =  -60.800_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           Case ('N2LO')
              inforcedname='N2LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .True.
              DMEorder        = 2
              ! Surface coefficients
              CrDr(0)  = -197.132_pr
              CrDr(1)  =  -12.503_pr
              CpV0(0)  = -272.164_pr
              CpV0(1)  = -193.188_pr
              CrdJ(0)  = -193.188_pr
              CrdJ(1)  =   37.790_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
              If(use_3N_couplings) Then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              End If
           !---------------------------------------------------------------------
           !NEDF Nuclear energy density functional
           !---------------------------------------------------------------------
           Case ('NEDF')
              inforcedname='NEDF'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,31.4703220539_pr, -64.0583220539_pr]
              b_NEDF=[-674.797522134_pr,59.2294137863_pr, 500.1281083477_pr]
              c_NEDF=[ 804.984434780_pr, 0._pr          ,-695.8644347800_pr]
              eta_NEDF = 90.300256558980_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('SeaLL1')
              inforcedname='SeaLL1'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              TDDFT_compatibility = .True.
              hb0_charge_dependent = .False.
              a_NEDF=[   0._pr         ,  64.2474102072_pr, -96.8354102072_pr]
              b_NEDF=[-684.524043779_pr, 119.8621469590_pr, 449.221896820_pr]
              c_NEDF=[ 827.262878410_pr,-256.4927039210_pr,-461.650174489_pr]
              eta_NEDF = 81.3917529003_pr
              W0_NEDF  = 73.5210618422_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF1')
              inforcedname='NEDF1'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,30.1396291921_pr, -62.7276291921_pr]
              b_NEDF=[-672.784776452_pr,56.7249539338_pr, 500.6198225180_pr]
              c_NEDF=[ 802.203802846_pr, 0._pr          ,-693.0838028460_pr]
              eta_NEDF = 105.629528244_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF2')
              inforcedname='NEDF2'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,30.3723779104_pr, -62.9603779104_pr]
              b_NEDF=[-672.212514042_pr,57.1630037930_pr, 499.6095102490_pr]
              c_NEDF=[ 801.413215568_pr, 0._pr          ,-692.2932155680_pr]
              eta_NEDF = 96.1410827861_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF3')
              inforcedname='NEDF3'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,30.3259301650_pr, -62.913930165_pr]
              b_NEDF=[-671.889100224_pr,57.0755857892_pr, 499.373514435_pr]
              c_NEDF=[ 800.966415558_pr, 0._pr          ,-691.846415558_pr]
              eta_NEDF = 91.896020403_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF4')
              inforcedname='NEDF4'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,31.5030085475_pr, -64.0910085475_pr]
              b_NEDF=[-672.625018435_pr,59.2909321226_pr, 497.8940863130_pr]
              c_NEDF=[ 801.983095259_pr, 0._pr          ,-692.8630952590_pr]
              eta_NEDF = 87.2919112779_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF5')
              inforcedname='NEDF5'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0.0_pr        , 43.8215204792_pr,  -76.4095204792_pr]
              b_NEDF=[-669.239299724_pr, 80.7199308819_pr,  473.0793688420_pr]
              c_NEDF=[ 797.305684868_pr,-92.7644849555_pr, -595.4211999120_pr]
              eta_NEDF = 82.95333556289553_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case default
              inforcedname=fname
              use_Namelist=.True.
           End Select FUNCTIONAL
       End If
       !---------------------------------------------------------------------
       ! Exit loop condition
       !---------------------------------------------------------------------
       If(.Not.use_Namelist) Exit
       !---------------------------------------------------------------------
       ! Read namelists
       !---------------------------------------------------------------------
       Open(UNIT=lnamelist,file=filename,status='OLD',delim='APOSTROPHE')
       Read(UNIT=lnamelist,NML=UNEDF_NAMELIST,iostat=ios)
       If(ios /= 0) Then
          ! WRong entry within UNEDF_NAMELIST.DAT file
          Write(*,'(1X,/,A)') 'ATTENTION: WRONG INPUT!'
          Write(*,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(*,*) 'INSIDE THE UNEDF_NAMELIST.DAT FILE IS WRONG.'
          Write(*,*) 'PLEASE CORRECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN read_UNEDF_NAMELIST'
       End If
       Close(lnamelist)
       If(Trim(FunctionalName) == Trim(inforcedname)) Exit
    End Do
    !---------------------------------------------------------------------
    ! See what the namelists modified
    !---------------------------------------------------------------------
    INFORCED_FUNCTIONAL: Select Case (Trim(inforcedname))
    Case ("FORCE")
       FunctionalName='FORCE'
    Case ("UNE0")
       FunctionalName='UNE0'
    Case ("UNE1")
       FunctionalName='UNE1'
    Case ("UNE2")
       FunctionalName='UNE2'
    Case ("HFB1")
       FunctionalName='HFB1'
    Case ("N0LO")
       FunctionalName='N0LO'
    Case ("N1LO")
       FunctionalName='N1LO'
    Case ("N2LO")
       FunctionalName='N2LO'
    Case ("NEDF")
       FunctionalName='NEDF'
    Case ("SeaLL1")
       FunctionalName='SeaLL1'
    Case ("NEDF1")
       FunctionalName='NEDF1'
    Case ("NEDF2")
       FunctionalName='NEDF2'
    Case ("NEDF3")
       FunctionalName='NEDF3'
    Case ("NEDF4")
       FunctionalName='NEDF4'
    Case ("NEDF5")
       FunctionalName='NEDF5'
    Case ("READ")
       FunctionalName='READ'
    Case ("RAND")
       FunctionalName='RAND'
    Case ("FITS")
       FunctionalName='FITS'
    Case ("DME_LO")
       FunctionalName='DME_LO'
    Case ("DME_NLO")
       FunctionalName='DME_NLO'
    Case ("DME_N2LO")
       FunctionalName='DME_N2LO'
    Case ("DME_NLOD")
       FunctionalName='DME_NLOD'
    Case ("DME_N2LOD")
       FunctionalName='DME_N2LOD'
    Case ("REG_LO")
       FunctionalName='REG_LO'
    Case ("REG_NLO")
       FunctionalName='REG_NLO'
    Case ("REG_N2LO")
       FunctionalName='REG_N2LO'
    Case ("REG_NLOD")
       FunctionalName='REG_NLOD'
    Case ("REG_N2LOD")
       FunctionalName='REG_N2LOD'
    Case default
       ! Missing entry within hfbtho_NAMELIST.dat file
       If(Trim(FunctionalName) /= Trim(inforcedname)) Then
          Write(*,'(1X,/,A)') 'ATTENTION: MISSING INPUT!'
          Write(*,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(*,*) 'IS MISSING INSIDE THE UNEDF_NAMELIST.DAT FILE.'
          Write(*,*) 'PLEASE CORECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN SET_FUNCTIONAL_PARAMETERS'
       End If
    End Select INFORCED_FUNCTIONAL

  End Subroutine read_UNEDF_NAMELIST
  !=======================================================================
  !> This routine sets up the parameters of the Gogny force.
  !=======================================================================
  Subroutine gogny_force(fname)
    Character(Len=30), Intent(Inout) :: fname

    is_gogny_force = .False.

    INTERACTION: Select Case (Trim(fname))
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1')
       is_gogny_force = .True.
       n_g = 2
       If(Allocated(mu_g)) Deallocate(mu_g,W_g,B_g,H_g,M_g)
       Allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -402.40_pr, -21.30_pr ]
       B_g = [ -100.00_pr, -11.77_pr ]
       H_g = [ -496.20_pr,  37.27_pr ]
       M_g = [  -23.56_pr, -68.81_pr ]
    !---------------------------------------------------------------------
    ! CPC 63, 365 (1991)
    !---------------------------------------------------------------------
    Case ('D1S')
       is_gogny_force = .True.
       n_g = 2
       If(Allocated(mu_g)) Deallocate(mu_g,W_g,B_g,H_g,M_g)
       Allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -1720.30_pr,  103.64_pr ]
       B_g = [  1300.00_pr, -163.48_pr ]
       H_g = [ -1813.53_pr,  162.81_pr ]
       M_g = [  1397.60_pr, -223.93_pr ]
       !W_g = [ -1720.30_pr,  103.639_pr ]
       !B_g = [  1300.00_pr, -163.483_pr ]
       !H_g = [ -1813.53_pr,  162.812_pr ]
       !M_g = [  1397.60_pr, -223.934_pr ]
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1p')
       is_gogny_force = .True.
       n_g = 2
       If(Allocated(mu_g)) Deallocate(mu_g,W_g,B_g,H_g,M_g)
       Allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -402.40_pr, -21.30_pr ]
       B_g = [ -100.00_pr, -11.77_pr ]
       H_g = [ -496.20_pr,  37.27_pr ]
       M_g = [  -23.56_pr, -68.81_pr ]
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('D1N')
       is_gogny_force = .True.
       n_g = 2
       If(Allocated(mu_g)) Deallocate(mu_g,W_g,B_g,H_g,M_g)
       Allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.8_pr, 1.2_pr ]
       W_g = [ -2047.61_pr,  293.02_pr ]
       B_g = [  1700.00_pr, -300.78_pr ]
       H_g = [ -2414.93_pr,  414.59_pr ]
       M_g = [  1519.35_pr, -316.84_pr ]
    !---------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------
    Case default
        !Write(6,'("No Gogny interaction defined in routine skforce()")')
    End Select INTERACTION

  End Subroutine gogny_force
  !=======================================================================
  !> This routine defines all gaussian ranges (Gogny and/or Coulomb).
  !=======================================================================
  Subroutine set_all_gaussians(icoul)
    Integer(ipr), Intent(In) :: icoul

    If(Allocated(mu_g_all)) Deallocate(mu_g_all)
    If(icoul < 0) Then
       coulomb_gaussian = .True.
       n_g_all = n_g + n_g_coul
       Allocate(mu_g_all(1:n_g_all))
       If(finite_range) mu_g_all(1:n_g) = mu_g(1:n_g)
       mu_g_all(n_g+1:n_g_all) = mu_g_coul(1:n_g_coul)
    Else If(finite_range) Then
       n_g_all = n_g
       Allocate(mu_g_all(1:n_g_all))
       mu_g_all = mu_g
    End If

  End Subroutine set_all_gaussians
  !=======================================================================
  !> This routine sets up the parameters of the pairing and Skyrme force and some useful combinations
  !> of those.
  !=======================================================================
  Subroutine skforce(fname,noForces)
    Integer(ipr), Intent(Inout) :: noForces
    Character(Len=30), Intent(Inout) :: fname

    Real(pr) :: wls,TA7,TA8
    Real(pr) :: zero,one,two,three,four,five,six,seven,eight,nine
    Real(pr) :: half,pp16,pp24

    zero = 0.0_pr; one = 1.0_pr; two = 2.0_pr; three = 3.0_pr; four = 4.0_pr
    five = 5.0_pr; six = 6.0_pr; seven = 7.0_pr; eight = 8.0_pr; nine = 9.0_pr
    half = 0.5_pr; pp16 = 16.0_pr; pp24 = 24.0_pr
    !
    ! Default for all forces if not modified
    hbzero = 1.0d0/0.04823_pr ! DMSHB0=1/hbzero
    sigma = one
    t0 = zero; x0 = zero
    t1 = zero; x1 = zero
    t2 = zero; x2 = one
    t3 = zero; x3 = one
    wls= zero; b4 = wls/two; b4p=wls/two
    te = zero; to = zero
    CExPar=1.0_pr

    noForces=0 ! No forces at all

    INTERACTION: Select Case (Trim(fname))
    !---------------------------------------------------------------------
    ! SIII, Beiner et al., Nucl. Phys. A 238, 29 (1975)
    !---------------------------------------------------------------------
    Case ('SIII')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73533_pr
        t0 = -.1128750d+04; x0 = +0.4500000_pr
        t1 = +.3950000d+03; x1 = +0.0000000_pr
        t2 = -.9500000d+02; x2 = +0.0000000_pr
        t3 = +.1400000d+05; x3 = +1.0000000_pr
        wls= +.1200000d+03; sigma = one
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKM*, Bartel et al., Nucl. Phys. A 386, 79 (1982)
    !---------------------------------------------------------------------
    Case ('SKM*')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73_pr
        t0 = -.2645000d+04; x0 = +.0900000_pr
        t1 = +.4100000d+03; x1 = +.0000000_pr
        t2 = -.1350000d+03; x2 = +.0000000_pr
        t3 = +.1559500d+05; x3 = +.0000000_pr
        wls= +.1300000d+03; sigma = one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! Modified SKM* for TDDFT calculations
    !   - center of mass correction is switched off
    !   - hbar^2/2m is set to the default value in TDDFT codes
    !---------------------------------------------------------------------
    Case ('SKM*mod')
        ! ph-Force
        noForces=1
        use_cm_cor = .False.
        TDDFT_compatibility = .True.
        hbzero = 20.7355180041_pr
        t0 = -.2645000d+04; x0 = +.0900000_pr
        t1 = +.4100000d+03; x1 = +.0000000_pr
        t2 = -.1350000d+03; x2 = +.0000000_pr
        t3 = +.1559500d+05; x3 = +.0000000_pr
        wls= +.1300000d+03; sigma = one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKP, Dobaczewski et al., Nucl. Phys. A 422, 103 (1984)
    !---------------------------------------------------------------------
    Case ('SKP')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.730_pr
        t0 =-0.2931696d+04; x0 = 0.2921515_pr
        t1 = 0.3206182d+03; x1 = 0.6531765_pr
        t2 =-0.3374091d+03; x2 =-0.5373230_pr
        t3 = 0.1870896d+05; x3 = 0.1810269_pr
        wls= 0.1000000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SLY4, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY4')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.2488913d+04; x0 = 0.8340000_pr
        t1 = 0.4868180d+03; x1 =-0.3440000_pr
        t2 =-0.5463950d+03; x2 =-1.0000000_pr
        t3 = 0.1377700d+05; x3 = 1.3540000_pr
        wls= 0.1230000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -325.2500_pr, -340.0625_pr /) ! HFB
    !---------------------------------------------------------------------
    ! Modified SLy4 for TDDFT calculations
    !   - center of mass correction is switched off
    !   - hbar^2/2m is set to the default value in TDDFT codes
    !---------------------------------------------------------------------
    Case ('SLY4mod')
        ! ph-Force
        noForces=1
        use_cm_cor = .False.
        TDDFT_compatibility = .True.
        hbzero = 20.7355180041_pr
        t0 =-0.2488913d+04; x0 = 0.8340000_pr
        t1 = 0.4868180d+03; x1 =-0.3440000_pr
        t2 =-0.5463950d+03; x2 =-1.0000000_pr
        t3 = 0.1377700d+05; x3 = 1.3540000_pr
        wls= 0.1230000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -370.000_pr, -370.000_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY5, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY5')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73553_pr
        t0 =-0.2483450d+04; x0 = 0.7760000_pr
        t1 = 0.4842300d+03; x1 =-0.3170000_pr
        t2 =-0.5566900d+03; x2 =-1.0000000_pr
        t3 = 0.1375700d+05; x3 = 1.2630000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY6, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY6')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2479500d+04; x0 = 0.8250000_pr
        t1 = 0.4621800d+03; x1 =-0.4650000_pr
        t2 =-0.4486100d+03; x2 =-1.0000000_pr
        t3 = 0.1367300d+05; x3 = 1.3550000_pr
        wls= 0.1220000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY7, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY7')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_j2terms     = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2480800d+04; x0 = 0.8480000_pr
        t1 = 0.4612900d+03; x1 =-0.4920000_pr
        t2 =-0.4339300d+03; x2 =-1.0000000_pr
        t3 = 0.1366900d+05; x3 = 1.3930000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SKI3, P.G.-Reinhard et al. Nucl. Phys. A584, 467  (1995)
    !---------------------------------------------------------------------
    Case ('SKI3')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.7525d0
        t0 =-0.176288d+04; x0 = 0.30830_pr
        t1 = 0.561608d+03; x1 =-1.17220_pr
        t2 =-0.227090d+03; x2 =-1.09070_pr
        t3 = 0.810620d+04; x3 = 1.29260_pr
        sigma=one/four
        b4 = 94.254_pr; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -357.2324_pr, -388.5625_pr /)
    !---------------------------------------------------------------------
    ! SKO, P.-G. Reinhard et al. Phys. Rev. C 60, 014316 (1999)
    !---------------------------------------------------------------------
    Case ('SKO')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.21036530d+04; x0 = -0.2107010_pr
        t1 = 0.30335200d+03; x1 = -2.8107520_pr
        t2 = 0.79167400d+03; x2 = -1.4615950_pr
        t3 = 0.13553252d+05; x3 = -0.4298810_pr
        wls= 0.12300000d+03; sigma=one/four
        b4 = 0.17657800d+03; b4p=-0.1987490d+03
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! SKO', P.-G. Reinhard et al. Phys. Rev. C 60, 014316 (1999)
    !---------------------------------------------------------------------
    Case ('SKOP')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        use_j2terms     = .True.
        hbzero = 20.735530_pr
        t0 =-0.20994190d+04; x0 = -0.029503_pr
        t1 = 0.30153100d+03; x1 = -1.325732_pr
        t2 = 0.15478100d+03; x2 = -2.323439_pr
        t3 = 0.13526464d+05; x3 = -0.147404_pr
        wls= 0.12300000d+03; sigma=one/four
        b4 = 0.14389500d+03; b4p=-0.828888d+02
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -253.753488522_pr, -274.679033802_pr /) ! See NB above
    !---------------------------------------------------------------------
    ! SKX, A.Brown Phys.Rev. C 58, 220 (1998)
    !---------------------------------------------------------------------
    Case ('SKX')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73_pr
        t0 = -1445.300_pr; x0 = 0.340_pr
        t1 =   246.900_pr; x1 = 0.580_pr
        t2 =  -131.800_pr; x2 = 0.127_pr
        t3 = 12103.900_pr; x3 = 0.030_pr
        sigma=one/two
        b4 = 0.0743d+03; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! Gogny D1, Decharge et al. PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1350.00_pr; x3 = one;
        wls= 115.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    ! Gogny D1S, Berger et al. CPC 63, 365 (1991)
    !---------------------------------------------------------------------
    Case ('D1S')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        !hbzero = 20.73667552957479_pr
        hbzero = 20.73667622931579050281_pr
        t3 = 6.0_pr*1390.600_pr; x3 = one;
        wls= 130.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    ! Gogny D1', Decharge et al. PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1p')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1350.00_pr; x3 = one;
        wls= 130.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('D1N')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1609.50_pr; x3 = one;
        wls= 115.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    !T0X0 A very simple skyrme functional
    !---------------------------------------------------------------------
    Case ('T0X0')
        noForces=1
        use_cm_cor  = .True.
        t0 =  0.0_pr; x0 = 0.0_pr;
        ! t0 =  -1128.75_pr; x0 = 0.45_pr;
        hbzero = 20.73667622931579_pr
        ! pp-Forces
        CpV1= zero
        CpV0= zero !No pairing here.
    !---------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------
    Case default
        ! Write(6,'("No Skyrme interaction defined in routine skforce()")')
    End Select INTERACTION

    If(noForces == 1) Then
       ! obtain coupling constants
       Call C_from_t()
       ! Frequent combinations entering the energy
       tv1   =  t0*(one+half*x0)*half;    tv2 = t0*(x0+half)*half
       tv3   =  t3*(one+half*x3)/12.0_pr; tv4 = t3*(x3+half)/12.0_pr
       tv5   = (t1*(one+half*x1)+t2*(one+half*x2))/four
       tv6   = (t2*(half+x2)-t1*(half+x1))/four
       ts1   = (t2*(one+half*x2)-three*t1*(one+half*x1))/pp16
       ts2   = (t1*(half+x1)*three+t2*(half+x2))/pp16
       t4o3  =  four/three; t324 = t3/pp24
       ! Frequent combinations entering the potential
       t0s   =  t0*(one-x0)*half; t0a = t0*(one+x0*half)
       drs   = (t2*(one+x2)-t1*(one-x1))*three/pp16
       dra   = (t2*(one+half*x2)-three*t1*(one+half*x1))/eight
       ts    = (t1*(one-x1) + three*t2*(one+x2))/eight
       ta    = (t1*(one+half*x1) + t2*(one+half*x2))/four
       t3alp = t3*(two+sigma)*(two+x3)/pp24
       t3al0 = t3*(x3+half)/six; t3alm = t3*sigma*(one+two*x3)/pp24
       alp   = one + sigma; alm = sigma - one
       wla0  = CrdJ(0)+CrdJ(1); wla1  = CrdJ(0)-CrdJ(1);
       TA7   = zero; TA8 = zero
       If(use_j2terms) Then
          TA7=(T1*(ONE-X1)-T2*(ONE+X2))/eight + five*to/four
          TA8=-(T1*X1+T2*X2)/four             + five*(te+to)/four
       End If
       TB7 = TA7; TB8 = TA8*half
       If(use_3N_couplings) Then
          use_3N_couplings = .False.
          override_3N_couplings = .True.
       End If
    End If

  End Subroutine skforce
  !=======================================================================
  !> This routine defines the parameters of the energy functional based on nuclear matter properties.
  !=======================================================================
  Subroutine set_functional_parameters(fname,npr)
    Character(Len=30), Intent(Inout) :: fname !< - Code name of the EDF (SkM*, SLy4, D1S, etc.)
    Integer(ipr), Dimension(3), Intent(In) :: npr !< - Number of particles: 1 neutrons, 2 protons, 3 total

    Integer(ipr), Parameter :: lin=15

    ! parameters
    FunctionalName = fname
    eps = Spacing(1.0_pr)
    kfconst = (1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK = 3.0_pr/5.0_pr*kfconst**2
    nuLambda = 700.0_pr ; nufpi = 93.0_pr

    ! exact Hartree CHrho from INM
    CHrho=0.0_pr

    If(use_INM) Then
       Call calculate_C_from_NM(E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,npr)
    Else
       Crho(0)=Crho(0)+CHrho
    End If
    Call calculate_NM_properties()

    Crho(0)=Crho(0)-CHrho

    Call calculate_natural_units()

  End Subroutine set_functional_parameters
  !=======================================================================
  !> This routine prints the parameters of the energy functional.
  !=======================================================================
  Subroutine print_functional_parameters(fname)
    Use HFBTHO_utilities, Only: lout,lfile

    Character(Len=30), Optional, Intent(In) :: fname !< - Code name of the EDF (SkM*, SLy4, D1S, etc.)

    Integer(ipr) :: iw

    Do iw=lout,lfile
       Write(iw,'("  ---------------------------------------")')
       Write(iw,'("           UNEDF Module Version: ",a)') Trim(Version)
       Write(iw,'("         M.Kortelainen & M.Stoitsov ")')
       Write(iw,'("  ---------------------------------------")')

       Write(iw,*)
       If(Present(fname)) Then
          Write(iw,'(2x,a," functional")') Trim(fname)
       Else
          Write(iw,'(2x,a," functional")') Trim(FunctionalName)
       End If
       Write(iw,'("  ----------------------------------------")')
       If(.not.is_nedf) Then
          Write(iw,'("  Crho(0)= ",g26.18,"; Crho(1)= ",g26.18)') Crho
          Write(iw,'("  CDrho(0)=",g26.18,"; CDrho(1)=",g26.18)') CDrho
          Write(iw,'("  Ctau(0)= ",g26.18,"; Ctau(1)= ",g26.18)') Ctau
          Write(iw,'("  CrDr(0)= ",g26.18,"; CrDr(1)= ",g26.18)') Crdr
          Write(iw,'("  CrdJ(0)= ",g26.18,"; CrdJ(1)= ",g26.18)') CrdJ
          Write(iw,'("  CJ(0)=   ",g26.18,"; CJ(1)=   ",g26.18)') CJ
          Write(iw,'("  CpV0(0)= ",g26.18,"; CpV0(1)= ",g26.18)') CpV0
          Write(iw,'("  CpV1(0)= ",g26.18,"; CpV1(1)= ",g26.18)') CpV1
          Write(iw,'("  sigma=   ",g26.18,"; hbzero=  ",g26.18)') sigma,hbzero
          Write(iw,'("  functional has DME couplings: ",L1)') force_is_dme
          If(override_3N_couplings) Then
             Write(iw,'("  functional ",a," has no 3N force, ignoring namelist input")') trim(FunctionalName)
          Else
             Write(iw,'("  use 3-Nucleon DME couplings: ",L1)') use_3n_couplings
          End If
       Else
          Write(iw,'("  a0=      ",g26.18,"; a1=      ",g26.18,"; a2=      ",g26.18)') a_NEDF
          Write(iw,'("  b0=      ",g26.18,"; b1=      ",g26.18,"; b2=      ",g26.18)') b_NEDF
          Write(iw,'("  c0=      ",g26.18,"; c1=      ",g26.18,"; c2=      ",g26.18)') c_NEDF
          Write(iw,'("  eta=     ",g26.18,"; W0=      ",g26.18)') eta_NEDF,W0_NEDF
          Write(iw,'("  CpV0(0)= ",g26.18,"; CpV0(1)= ",g26.18)') CpV0
          Write(iw,'("  CpV1(0)= ",g26.18,"; CpV1(1)= ",g26.18)') CpV1
          Write(iw,'("  hbzeron= ",g26.18,"; hbzerop= ",g26.18)') hbzeron,hbzerop
       End If
       ! Finite-range force (Gogny force)
       If(finite_range) Then
          Write(iw,*)
          Write(iw,'("  Finite-range potential")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("  mu=",10f26.18)') mu_g
          Write(iw,'("  W=",10f26.18)') W_g
          Write(iw,'("  B=",10f26.18)') B_g
          Write(iw,'("  H=",10f26.18)') H_g
          Write(iw,'("  M=",10f26.18)') M_g
       End If
       ! Common features
       Write(iw,'("  e^2 chrg=",g26.18,"; CExPar=  ",g26.18)') e2charg,CExPar
       Write(iw,'("  c.m. correction: ",L1,", chr. density in direct Coul: ",L1)') use_cm_cor,use_charge_density
       Write(iw,'("  use tensor terms: ",L1)') use_j2terms
       ! Natural units
       If(.not.is_nedf .And. .Not.finite_range) Then
          Write(iw,*)
          Write(iw,'("  Coupling constants in natural units")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("  Crho_nu(0)= ",g26.18,"; Crho_nu(1)= ",g26.18)') nuCrho
          Write(iw,'("  CDrho_nu(0)=",g26.18,"; CDrho_nu(1)=",g26.18)') nuCDrho
          Write(iw,'("  Ctau_nu(0)= ",g26.18,"; Ctau_nu(1)= ",g26.18)') nuCtau
          Write(iw,'("  CrDr_nu(0)= ",g26.18,"; CrDr_nu(1)= ",g26.18)') nuCrdr
          Write(iw,'("  CrdJ_nu(0)= ",g26.18,"; CrdJ_nu(1)= ",g26.18)') nuCrdJ
          Write(iw,'("  CJ_nu(0)=   ",g26.18,"; CJ_nu(1)=   ",g26.18)') nuCJ
          Write(iw,'("  CpV0_nu(0)= ",g26.18,"; CpV0_nu(1)= ",g26.18)') nuCpV0
          Write(iw,'("  CpV1_nu(0)= ",g26.18,"; CpV1_nu(1)= ",g26.18)') nuCpV1
          Write(iw,'("  fpi_nu=     ",g26.18,"; Lambda_nu=  ",g26.18)') nufpi,nuLambda
       End If
       ! DME
       If(dmeorder >= 0) Then
          Write(iw,*)
          Write(iw,'("  DME parameters")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("       gA=",f12.6," mpi [1/fm]=",f12.6," fpi [1/fm]=",f12.6)') gA,mass_pion,fpi
          Write(iw,'("  c1 [fm]=",f12.6,"    c3 [fm]=",f12.6,"    c4 [fm]=",f12.6)') c1,c3,c4
          Write(iw,'("       cd=",f12.6,"         ce=",f12.6," LamX[1/fm]=",f12.6)') cd,ce,LambdaX
          Write(iw,'("  ->CHrho=",f12.6)') CHrho
          If(dmeorder >= 2) Write(iw,'("  use 3N terms: ",L1)') use_DME3N_terms
       End If
       Write(iw,*)
       Write(iw,'("  Nuclear matter properties")')
       Write(iw,'("  ----------------------------------------")')
       Write(iw,'("  E_NM=    ",g26.18,"; K_NM=     ",g26.18)') E_NM,K_NM
       Write(iw,'("  P_NM=    ",g26.18,"; RHO_NM=   ",g26.18)') P_NM,RHO_NM
       Write(iw,'("  ASS_NM=  ",g26.18,"; LASS_NM=  ",g26.18)') ASS_NM,LASS_NM
       Write(iw,'("  SMASS_NM=",g26.18,"; VMASS_NM= ",g26.18)') SMASS_NM,VMASS_NM

       If(.Not. is_gogny_force) Call t_from_C()

       ! (t,x) parametrization of the Skyrme functional
       If(.not.is_nedf) Then
          Write(iw,*)
          Write(iw,'("  Associated (t,x)-coupling constants")')
          Write(iw,'("  ----------------------------------------")')
          If(.Not.finite_range .Or. force_is_dme) Then
             Write(iw,'("  t0=    ",g26.18,"; x0=     ",g26.18)') t0,x0
             Write(iw,'("  t1=    ",g26.18,"; x1=     ",g26.18)') t1,x1
             Write(iw,'("  t2=    ",g26.18,"; x2=     ",g26.18)') t2,x2
         End If
          Write(iw,'("  t3=    ",g26.18,"; x3=     ",g26.18)') t3,x3
          Write(iw,'("  b4=    ",g26.18,"; b4p=    ",g26.18)') b4,b4p
          Write(iw,'("  te=    ",g26.18,"; to=     ",g26.18)') te,to
          Write(iw,'("  sigma= ",g26.18,"; hbzero= ",g26.18)') sigma,hbzero
       End If

       If(Print_Namelist) Then
          Write(iw,*)
          SELECTED_FUNCTIONAL: Select Case (Trim(FunctionalName))
          Case ("UNEDF","SKYRME")
                Write(iw,'("NAMELIST CONTENT (cannot be modified for functional names UNEDF,SKYRME)")')
                Write(iw,'("-----------------------------------------------------------------------")')
          Case ("FITS")
                Write(iw,'("NAMELIST CONTENT (advanced usage: modify all but not C-, NM-, and more...)")')
                Write(iw,'("--------------------------------------------------------------------------")')
          Case default
                Write(iw,'("NAMELIST CONTENT (copy/past to UNEDF_NAMELIST.DAT and modify)")')
                Write(iw,'("-------------------------------------------------------------")')
          End Select SELECTED_FUNCTIONAL
          Write(*,'(" !NB: FUNCTIONALNAME should be always in quotations")')
          Write(*,UNEDF_NAMELIST)
       End If
    End Do

  End Subroutine print_functional_parameters
  !=======================================================================
  !> This routine calculates the coupling constants in natural units
  !=======================================================================
  Subroutine calculate_natural_units()

    nuCrho  = Crho *(nufpi**2)/(mevfm**3)
    nuCdrho = Cdrho*(nufpi**2)*((nuLambda*nufpi*nufpi)**sigma)/(mevfm**(3.0_pr*(1.0_pr+sigma)))
    nuCtau  = Ctau*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrDr  = CrDr*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrdJ  = CrdJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCJ    = CJ  *((nufpi*nuLambda)**2)/(mevfm**5)
    nuCpV0  = CpV0*( nufpi**2)/(mevfm**3)
    nuCpV1  = CpV1*( nufpi**4)*nuLambda/(mevfm**6)

  End Subroutine calculate_natural_units
  !=======================================================================
  !> This routine calculates the volume coupling constants (and the density dependency \f$ \sigma \f$) from
  !> nuclear matter properties
  !>  @param[in] E - Energy per nucleon of nuclear matter
  !>  @param[in] K - Incompressibility parameter
  !>  @param[in] SMASS -Inverse of the scalar effective mass \f$ 1/m_{s}^{*} \f$
  !>  @param[in] RHO - Saturation density of nuclear matter
  !>  @param[in] ASS - Asymmetry coefficient \f$ a_{\rm sym} \f$ of nuclear matter
  !>  @param[in] LASS - Slope of the energy per nucleon curve \f$ L_{\rm sym} \f$
  !>  @param[in] VMASS - Vector effective mass
  !>  @param[in] npr - Number of particles: 1 neutrons, 2 protons, 3 total
  !>  @param[in] sigma_NM - Optional: density dependency. If not present, it is computed from nuclear matter
  !=======================================================================
  Subroutine calculate_C_from_NM(E,K,SMASS,RHO,ASS,LASS,VMASS,npr,sigma_NM)
    Real(pr), Intent(In) :: E,K,SMASS,RHO,ASS,LASS,VMASS
    Real(pr), Intent(In), Optional :: sigma_NM
    Integer(ipr), Dimension(3), Intent(In) :: npr

    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho2,rho0,rho1
    Real(pr) :: E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr

    tauc=CK*RHO**c23; u=(kfconst/mass_pion)*RHO**c13; rho2=rho**2

    ! Initial guess for isovector density based on current number of neutrons and protons. For DME functionals,
    ! the routine calculate_U_parameters() is called at each iteration with the proper values of the densities
    rho0=RHO; rho1=Real(npr(1)-npr(2),Kind=pr)/Real(npr(1)+npr(2),Kind=pr)*RHO
    Call calculate_U_parameters(rho0,rho1,tauc*rho0,tauc*rho1,0.0_pr,0.0_pr)
    If(finite_range) Then
       Call calculate_finite_range_NM(rho,E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr)
    Else
       E_fr=0_pr; P_fr=0_pr; K_fr=0_pr; SMASS_fr=0_pr
       ASS_fr=0_pr; LASS_fr=0_pr; KA_fr=0._pr; VMASS_fr=0_pr
    End If

    aRho0Rho0   = 0.50_pr*(aRhoRho   + bRhoRho)*mevfm
    aRho1Rho1   = 0.50_pr*(aRhoRho   - bRhoRho)*mevfm
    aRho0Tau0   = 0.50_pr*(aRhoTau   + bRhoTau)*mevfm
    aRho1Tau1   = 0.50_pr*(aRhoTau   - bRhoTau)*mevfm
    daRho0Rho0  = 0.50_pr*(daRhoRho  + dbRhoRho)*mevfm
    daRho1Rho1  = 0.50_pr*(daRhoRho  - dbRhoRho)*mevfm
    daRho0Tau0  = 0.50_pr*(daRhoTau  + dbRhoTau)*mevfm
    daRho1Tau1  = 0.50_pr*(daRhoTau  - dbRhoTau)*mevfm
    ddaRho0Rho0 = 0.50_pr*(ddaRhoRho + ddbRhoRho)*mevfm
    ddaRho1Rho1 = 0.50_pr*(ddaRhoRho - ddbRhoRho)*mevfm
    ddaRho0Tau0 = 0.50_pr*(ddaRhoTau + ddbRhoTau)*mevfm
    ddaRho1Tau1 = 0.50_pr*(ddaRhoTau - ddbRhoTau)*mevfm

    ! set/calculate sigma
    If(Present(sigma_NM)) Then
       sigma=sigma_NM
    Else
       sigma = ((1.0_pr/3.0_pr)*(-(K-K_fr)+tauc*hbzero*(-3.0_pr+4.0_pr*(SMASS-SMASS_fr))-9.0_pr*(E-E_fr+P_fr/rho) &
             + 9.0_pr*RHO2*hRho0Rho0 + 21.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+5.0_pr*tauc*daRho0Tau0 &
             + 7.0_pr*RHO*dhRho0Rho0+11.0_pr*tauc*RHO*dhRho0Tau0+u*ddaRho0Rho0 &
             + u*tauc*ddaRho0Tau0+u*RHO*ddhRho0Rho0+u*tauc*RHO*ddhRho0Tau0))) &
             / (tauc*hbzero*(-3.0_pr+2.0_pr*(SMASS-SMASS_fr))+3.0_pr*(E-E_fr+P_fr/rho)+3.0_pr*RHO2*hRho0Rho0 &
             + 3.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+tauc*daRho0Tau0 &
             +  RHO*dhRho0Rho0+tauc*RHO*dhRho0Tau0))
    End If

    ! Isoscalar
    Crho(0) = (c13*(tauc*hbzero*(-3.0_pr+(2.0_pr-3.0_pr*sigma)*(SMASS-SMASS_fr)) &
            + 3.0_pr*(1.0_pr+sigma)*(E-E_fr)+3.0_pr*P_fr/rho-3.0_pr*sigma*RHO*aRho0Rho0 &
            + 3.0_pr*(1.0_pr-sigma)*RHO2*hRho0Rho0+3.0_pr*tauc*RHO2*hRho0Tau0 &
            + u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
            + tauc*RHO*dhRho0Tau0)))/(sigma*RHO)
    Cdrho(0) = (c13*RHO**(-1.0_pr-sigma)*(tauc*hbzero*(3.0_pr-2.0_pr*(SMASS-SMASS_fr))&
             - 3.0_pr*(E-E_fr+P_fr/rho)-3.0_pr*RHO**2*hRho0Rho0-3.0_pr*tauc*RHO2*hRho0Tau0&
             - u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
             + tauc*RHO*dhRho0Tau0)))/sigma
    Ctau(0) = (hbzero*(SMASS-SMASS_fr-1.0_pr)-RHO*(aRho0Tau0+RHO*hRho0Tau0))/RHO

    ! Isovectors
    Crho(1) = (27.0_pr*(ASS-ASS_fr)*(1.0_pr+sigma)-9.0_pr*(LASS-LASS_fr) &
            + 5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*(VMASS-VMASS_fr)+3.0_pr*sigma*(-4.0_pr+3.0_pr*(VMASS-VMASS_fr))) &
            + 20.0_pr*tauc*(2.0_pr-3.0_pr*sigma)*RHO*aRho0Tau0 &
            + RHO*(-27.0_pr*sigma*aRho1Rho1+5.0_pr*tauc*(11.0_pr-12.0_pr*sigma)*RHO*hRho0Tau0 &
            - 27.0_pr*(-1.0_pr+sigma)*RHO*hRho1Rho1+9.0_pr*tauc*(5.0_pr-3.0_pr*sigma)*RHO*hRho1Tau0 &
            + 45.0_pr*tauc*RHO*hRho1Tau1+40.0_pr*tauc*Ctau(0)-60.0_pr*tauc*sigma*Ctau(0) &
            + 5.0_pr*u*tauc*daRho0Tau0+9.0_pr*u*daRho1Rho1+15.0_pr*u*tauc*daRho1Tau1 &
            + 5.0_pr*u*tauc*RHO*dhRho0Tau0+9.0_pr*u*RHO*dhRho1Rho1+9.0_pr*u*tauc*RHO*dhRho1Tau0 &
            + 15.0_pr*u*tauc*RHO*dhRho1Tau1))/(27.0_pr*sigma*RHO)
    Cdrho(1) = -(RHO**(-1.0_pr-sigma)*(27.0_pr*(ASS-ASS_fr)-9.0_pr*(LASS-LASS_fr) &
             + 5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*(VMASS-VMASS_fr))+40.0_pr*tauc*RHO*aRho0Tau0 &
             + 55.0_pr*tauc*RHO2*hRho0Tau0+27.0_pr*RHO**2*hRho1Rho1+45.0_pr*tauc*RHO2*hRho1Tau0 &
             + 45.0_pr*tauc*RHO2*hRho1Tau1+40.0_pr*tauc*RHO*Ctau(0) +5.0_pr*u*tauc*RHO*daRho0Tau0 &
             + 9.0_pr*u*RHO*daRho1Rho1+15.0_pr*u*tauc*RHO*daRho1Tau1 &
             + 5.0_pr*u*tauc*RHO2*dhRho0Tau0+9.0_pr*u*RHO2*dhRho1Rho1 &
             + 9.0_pr*u*tauc*RHO2*dhRho1Tau0 +15.0_pr*u*tauc*RHO2*dhRho1Tau1))/(27.0_pr*sigma)
    Ctau(1) = (hbzero-hbzero*(VMASS-VMASS_fr)+RHO*(aRho0Tau0-aRho1Tau1+RHO*hRho0Tau0-RHO*hRho1Tau1+Ctau(0)))/RHO

  End Subroutine calculate_C_from_NM
  !======================================================================
  !> This routine calculates the contribution to the nuclear matter properties comming from the
  !> finite range part of the functional. The formulas for E,P,K,ASS,LASS and KASS are from
  !> \cite sellahewa2014isovector. The formulas for SMASS and VMASS where derived from the energy
  !> density as linear combinations of the derivatives with respect of \f$\tau_n\f$ and \f$\tau_p\f$.
  !> The formulas are as follows
  !> \f{align*}{
  !>   E^{\rm NM}(\rho_0) & =
  !>     \frac{1}{2}\sum_{i=1,n_g}\left\{A_0^i\rho_0+B_0^ig_0(\mu_i k_F)
  !>     \right\}, \\
  !>   P^{\rm NM}(\rho_0) &=
  !>     \sum_{i=1,n_g} \left\{\frac{1}{2} A_0^i\rho_0^2 +
  !>     B_0^ip_0(\mu_ik_F)\rho_0\right\}, \\
  !>   K^{\rm NM}(\rho_0) &= -3\sum_{i=1,n_g} B_0^ik_0(\mu_ik_F), \\
  !>   M_s(\rho_0) &= \frac{2m_N}{\hbar^2}\frac{1}{2}\sum_{i=1,n_g}
  !>     B_0^im_0(\mu_i k_F)\frac{3}{5k_F^2}, \\
  !>   a^{\rm NM}_{\rm sym}(\rho_0) &= \sum_{i=1,n_g}\left\{\frac{1}{2}
  !>     A_1^i\rho_0+\frac{1}{6}\left[B_{nn}^i s_1(\mu_i k_F) + B_{np}^i
  !>     s_2(\mu_i k_F) \right] \right\}, \\
  !>   L^{\rm NM}_{\rm sym}(\rho_0) &= \sum_{i=1,n_g}\left\{\frac{3}{2}
  !>     A_1^i\rho_0+\frac{1}{6}\left[B_{nn}^i l_1(\mu_i k_F) + B_{np}^i
  !>     l_2(\mu_i k_F) \right] \right\}, \\
  !>   \Delta K^{\rm NM}(\rho_0) &= -\frac{1}{6}\sum_{i=1,n_g}\left\{
  !>     B_{nn}^i k_1(\mu_i k_F) + B_{np}^ik_2(\mu_i k_F)\right\}, \\
  !>   M_v(\rho_0) &= \frac{2m_N}{\hbar^2}\frac{1}{2}\sum_{i=1,n_g}
  !>    \left\{ B_0^im_0(\mu_i k_F) + B_{nn}^im_1(\mu_i k_F) +
  !>    B_{np}^im_2(\mu_i k_F) \right\}\frac{3}{5k_F^2},
  !> \f}
  !> where
  !> \f{align*}{
  !>  k_F &= \left(\frac{3\pi^2}{2} \rho_0 \right)^{1/3}, \\
  !>  A_0^i &= \frac{\pi^{3/2}\mu_i^3}{4}(4W_i+2B_i-2H_i-M_i), \\
  !>  A_1^i &= \frac{\pi^{3/2}\mu_i^3}{4}(-2H_i-M_i), \\
  !>  B_{nn}^i &= -\frac{1}{\sqrt{\pi}}(W_i+2B_i-H_i-2M_i), \\
  !>  B_{np}^i &= -\frac{1}{\sqrt{\pi}}(-H_i-2M_i), \\
  !>  B_0^i &= B_{nn}^i + B_{np}^i, \\
  !>  g_0(q) &= \frac{1}{q^3}\left[2-3q^2-\left(2-q^2\right)e^{-q^2}
  !>     \right] + \sqrt{\pi} {\rm erf}(q), \\
  !>  p_0(q) &= \frac{1}{2q^3}\left[-2+ q^2+\left(2+q^2\right)e^{-q^2}
  !>     \right], \\
  !>  k_0(q) &= \frac{1}{q^3}\left[-6+ 2q^2+\left(6+4q^2+q^4\right)
  !>     e^{-q^2} \right], \\
  !>  m_0(q) &= \frac{1}{2q^3}\left[-6+ 3q^2+\left(6+3q^2\right)e^{-q^2}
  !>     \right], \\
  !>  s_1(q) &=\frac{1}{q}\left[1-\left(1+q^2\right)e^{-q^2}\right], \\
  !>  s_2(q) &=\frac{1}{q}\left[1-q^2-e^{-q^2}\right], \\
  !>  l_1(q) &=\frac{1}{q}\left[-1+\left(1+q^2+2q^4\right)e^{-q^2}
  !>     \right], \\
  !>  l_2(q) &=\frac{1}{q}\left[-1-q^2+\left(1+2q^2\right)e^{-q^2}
  !>     \right], \\
  !>  k_1(q) &=\frac{1}{2q}\left[-2+\left(2+2q^2+q^4+2q^6\right)e^{-q^2}
  !>     \right], \\
  !>  k_2(q) &=\frac{1}{2q}\left[-2-q^2+\left(2+3q^2+2q^4\right)e^{-q^2}
  !>     \right], \\
  !>  m_1(q) &=\frac{1}{2q^3}\left[-10+3q^2+\left(10+7q^2+2q^4\right)
  !>     e^{-q^2} \right], \\
  !>  m_2(q) &=\frac{1}{q^3}\left[-2+q^4+\left(2+2q^2\right)e^{-q^2}
  !>     \right].
  !> \f}
  !>
  !>  @param[in] rho - Saturation density of nuclear matter
  !>  @param[in] E - Energy per nucleon of nuclear matter
  !>  @param[in] P - Pressure of nuclear matter
  !>  @param[in] K - Incompressibility parameter
  !>  @param[in] SMASS -Inverse of the scalar effective mass \f$ 1/m_{s}^{*} \f$
  !>  @param[in] ASS - Asymmetry coefficient \f$ a_{\rm sym} \f$ of nuclear matter
  !>  @param[in] LASS - Slope of the energy per nucleon curve \f$ L_{\rm sym} \f$
  !>  @param[in] VMASS - Vector effective mass
  !======================================================================
  Subroutine calculate_finite_range_NM(rho,E,P,K,SMASS,ASS,LASS,KASS,VMASS)
    Real(pr), Intent(In) :: rho
    Real(pr), intent(Out) :: E,P,K,SMASS,ASS,LASS,KASS,VMASS

    Integer(ipr) :: ig
    Real(pr) :: A0,A1,Bnn,Bnp,B0,W,B,H,M,mui,kF,q,tauc
    Real(pr) :: g0,p0,k0,m0,s1,s2,l1,l2,k1,k2,m1,m2

    kF = (1.5_pr*pi**2*rho)**(1.0_pr/3.0_pr)
    tauc = 3*kF**2/5.0_pr
    E=0.0_pr; P=0.0_pr; K=0.0_pr; SMASS=0.0_pr
    ASS=0.0_pr; LASS=0.0_pr; KASS=0.0_pr; VMASS=0.0_pr

    Do ig = 1,n_g
       mui = mu_g(ig)
       W = W_g(ig)
       B = B_g(ig)
       H = H_g(ig)
       M = M_g(ig)
       A0 =  0.25_pr*(sqrt(pi)*mui)**3*(4*W+2*B-2*H-M)
       A1 =  0.25_pr*(sqrt(pi)*mui)**3*(-2*H-M)
       If(force_is_dme) Then
          Bnn = 0._pr
          Bnp = 0._pr
       Else
          Bnn = -(W+2*B-H-2*M)/sqrt(pi)
          Bnp =  (H+2*M)/sqrt(pi)
       End If
       B0 = Bnn + Bnp
       q = mui*kF
       g0 = 2/q**3 - 3/q - (2/q**3 - 1/q)*exp(-q**2) + sqrt(pi)*erf(q)
       p0 = -1/q**3 + 1/(2*q) + (1/q**3 + 1/(2*q))*exp(-q**2)
       k0 = -6/q**3 + 2/q + (6/q**3 + 4/q + q)*exp(-q**2)
       m0 = 3*(-2+q**2+exp(-q**2)*(2+q**2))/(2*q**3)
       s1 = 1/q - (1/q + q)*exp(-q**2)
       s2 = 1/q - q -1/q*exp(-q**2)
       l1 = -1/q + (1/q + q + 2*q**3)*exp(-q**2)
       l2 = -1/q - q + (1/q + 2*q)*exp(-q**2)
       k1 = -1/q + (1/q + q + q**3/2._pr + q**5)*exp(-q**2)
       k2 = -1/q - q/2._pr + (1/q + 3/2._pr*q + q**3)*exp(-q**2)
       m1 = (-10+3*q**2+exp(-q**2)*(10+7*q**2+2*q**4))/(2*q**3)
       m2 = (-2+q**4+2*exp(-q**2)*(1+q**2))/(q**3)
       E = E + 0.5_pr*(A0*rho + B0*g0)
       P = P + 0.5_pr*A0*rho**2 + B0*p0*rho
       K = K - 3*B0*k0
       SMASS = SMASS + B0*m0/(2*tauc*hbzero)
       ASS = ASS + 0.5_pr*A1*rho + (Bnn*s1+Bnp*s2)/6._pr
       LASS = LASS + 1.5_pr*A1*rho + (Bnn*l1+Bnp*l2)/6._pr
       KASS = KASS - 2*(Bnn*k1+Bnp*k2)/3._pr
       VMASS = VMASS + (B0*m0+Bnn*m1+Bnp*m2)/(2*tauc*hbzero)
    End Do

  End Subroutine calculate_finite_range_NM
  !======================================================================
  !> This routine calculates the finite range contribution to the Nuclear Matter pressure.
  !> See \cite sellahewa2014isovector for details
  !> @result \f$ P^{\rm NM}(\rho_0) = \sum_{i=1,n_g} \left\{\frac{1}{2}
  !>          A_0^i\rho_0^2 + B_0^ip_0(\mu_ik_F)\rho_0\right\}, \f$
  !======================================================================
  Function P_SNM_FR(rho) result(P)
    Real(pr), Intent(In) :: rho

    Real(pr) :: P
    Integer(ipr) :: ig
    Real(pr) :: A0,B0,W,B,H,M,mui,kF,q,p0

    kF = (1.5_pr*pi**2*rho)**(1._pr/3._pr)
    P=0._pr
    Do ig = 1,n_g
       mui = mu_g(ig)
       W = W_g(ig)
       B = B_g(ig)
       H = H_g(ig)
       M = M_g(ig)
       A0 =  0.25_pr*(Sqrt(pi)*mui)**3*(4*W+2*B-2*H-M)
       If(force_is_dme) Then
          B0 = 0.0_pr
       Else
          B0 = -(W+2*B-2*H-4*M)/Sqrt(pi)
       End If
       q = mui*kF
       p0 = -1/q**3 + 1/(2*q) + (1/q**3 + 1/(2*q))*Exp(-q**2)
       P = P + 0.5_pr*A0*rho**2 + B0*p0*rho
    End Do

  End Function P_SNM_FR
  !=======================================================================
  !> This routine calculates the nuclear matter properties
  !=======================================================================
  Subroutine calculate_NM_properties()
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho_NM2
    Real(pr), Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    Real(pr) :: E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr

    RHO_NM=find_NM_RHOC()

    If(finite_range) Then
       Call calculate_finite_range_NM(rho_nm,E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr)
    Else
       E_fr=0_pr; P_fr=0_pr; K_fr=0_pr; SMASS_fr=0_pr
       ASS_fr=0_pr; LASS_fr=0_pr; KA_fr=0._pr; VMASS_fr=0_pr
    End If

    aRho0Rho0   = 0.50_pr*(aRhoRho   + bRhoRho)*mevfm
    aRho1Rho1   = 0.50_pr*(aRhoRho   - bRhoRho)*mevfm
    aRho0Tau0   = 0.50_pr*(aRhoTau   + bRhoTau)*mevfm
    aRho1Tau1   = 0.50_pr*(aRhoTau   - bRhoTau)*mevfm
    daRho0Rho0  = 0.50_pr*(daRhoRho  + dbRhoRho)*mevfm
    daRho1Rho1  = 0.50_pr*(daRhoRho  - dbRhoRho)*mevfm
    daRho0Tau0  = 0.50_pr*(daRhoTau  + dbRhoTau)*mevfm
    daRho1Tau1  = 0.50_pr*(daRhoTau  - dbRhoTau)*mevfm
    ddaRho0Rho0 = 0.50_pr*(ddaRhoRho + ddbRhoRho)*mevfm
    ddaRho1Rho1 = 0.50_pr*(ddaRhoRho - ddbRhoRho)*mevfm
    ddaRho0Tau0 = 0.50_pr*(ddaRhoTau + ddbRhoTau)*mevfm
    ddaRho1Tau1 = 0.50_pr*(ddaRhoTau - ddbRhoTau)*mevfm
    tauc=CK*RHO_NM**c23; u=(kfconst/mass_pion)*RHO_NM**c13; rho_NM2=rho_NM**2

    ! Symmetric Nuclear Matter
    E_NM = tauc*hbzero+RHO_NM*(aRho0Rho0+RHO_NM*hRho0Rho0+Crho(0)+RHO_NM**sigma*Cdrho(0)) &
         + tauc*RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0)) + E_fr
    P_NM = c13*RHO_NM**2*((2.0_pr*tauc*hbzero)/RHO_NM+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
         + 6.0_pr*RHO_NM*hRho0Rho0+8.0_pr*tauc*RHO_NM*hRho0Tau0+3.0_pr*Crho(0) &
         + 3.0_pr*(1+sigma)*RHO_NM**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
         + u*tauc*daRho0Tau0+u*RHO_NM*dhRho0Rho0+u*tauc*RHO_NM*dhRho0Tau0) + P_fr
    SMASS_NM = 1.0_pr+(RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0)))/hbzero + SMASS_fr
    K_NM = 9.0_pr*sigma*(1+sigma)*RHO_NM**(1+sigma)*Cdrho(0) &
         + (-2.0_pr*tauc*hbzero+10.0_pr*tauc*RHO_NM*aRho0Tau0+18.0_pr*RHO_NM2*hRho0Rho0 &
         + 40.0_pr*tauc*RHO_NM**2*hRho0Tau0+4.0_pr*u*RHO_NM*daRho0Rho0 &
         + RHO_NM*(10.0_pr*tauc*Ctau(0)+u*(8.0_pr*tauc*daRho0Tau0+u*ddaRho0Rho0 &
         + (10.0_pr*RHO_NM*dhRho0Rho0+14.0_pr*tauc*RHO_NM*dhRho0Tau0 &
         + (u*tauc*ddaRho0Tau0+u*RHO_NM*ddhRho0Rho0+u*tauc*RHO_NM*ddhRho0Tau0))))) + k_fr

    ! Asymmetric Nuclear Matter
    ASS_NM = RHO_NM2*hRho1Rho1+RHO_NM*(aRho1Rho1+Crho(1)+RHO_NM**sigma*Cdrho(1)) &
           + (tauc*(5.0_pr*hbzero+RHO_NM*(5.0_pr*aRho0Tau0+15.0_pr*aRho1Tau1+5.0_pr*RHO_NM*hRho0Tau0 &
           + 9.0_pr*RHO_NM*hRho1Tau0+5.0_pr*(3.0_pr*RHO_NM*hRho1Tau1+Ctau(0)+3.0_pr*Ctau(1)))))/9.0_pr + ASS_fr
    VMASS_NM = (hbzero+RHO_NM*(aRho0Tau0-aRho1Tau1+RHO_NM*hRho0Tau0-RHO_NM*hRho1Tau1+Ctau(0)-Ctau(1)))/hbzero + VMASS_fr
    LASS_NM = 6.0_pr*RHO_NM2*hRho1Rho1+3.0_pr*RHO_NM*(aRho1Rho1+Crho(1)+(1.0_pr+sigma)*RHO_NM**sigma*Cdrho(1)) &
            + u*RHO_NM*daRho1Rho1 +u*RHO_NM2*dhRho1Rho1 &
            + (tauc*(10.0_pr*hbzero+8.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
            + 25.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
            + 5.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
            + u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0+15.0_pr*dhRho1Tau1)))/9.0_pr + LASS_fr
    KA_NM = 18.0_pr*RHO_NM2*hRho1Rho1+9.0_pr*sigma*(1.0_pr+sigma)*RHO_NM**(1.0_pr+sigma)*Cdrho(1) &
          + 4.0_pr*u*RHO_NM*daRho1Rho1 +10.0_pr*u*RHO_NM2*dhRho1Rho1 &
          +  u**2*RHO_NM*ddaRho1Rho1+u**2*RHO_NM2*ddhRho1Rho1 &
          + (tauc*(-10.0_pr*hbzero+40.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
          + 50.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
          + 40.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
          + 14.0_pr*u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0 &
          + 15.0_pr*dhRho1Tau1)+5.0_pr*u**2*RHO_NM*(ddaRho0Tau0 &
          + 3.0_pr*ddaRho1Tau1)+u**2*RHO_NM2*(5.0_pr*ddhRho0Tau0+9*ddhRho1Tau0+15*ddhRho1Tau1)))/9._pr + KA_fr

  End Subroutine calculate_NM_properties
  !=======================================================================
  !> This routine finds the saturation density RHO_NM using the secant method
  !=======================================================================
  Real(pr) Function find_NM_RHOC()
    Integer(ipr) :: iter
    Real(pr) :: aRho0Rho0,daRho0Rho0
    Real(pr) :: aRho0Tau0,daRho0Tau0
    Real(pr) :: kfconstmpi,u,tauc
    Real(pr) :: rhom0,rhom,rhom2,w,w0,step,P_FR
    Real(pr), Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr

    kfconstmpi=kfconst/mass_pion; step=-0.0010_pr; iter=0
    ! initial point
    rhom=0.170_pr; tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2

    Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)

    If(finite_range) Then
       P_FR = P_SNM_FR(rhom)
    Else
       P_FR = 0
    End If
    aRho0Rho0 = 0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0 = 0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    aRho0Tau0 = 0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0 = 0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    w0 = c13*rhom2*((2.0_pr*tauc*hbzero)/rhom + 3.0_pr*aRho0Rho0 + 5.0_pr*tauc*aRho0Tau0 &
                   + 6.0_pr*rhom*hRho0Rho0 + 8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
                   + 3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0) + 5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
                   + u*tauc*daRho0Tau0 + u*rhom*dhRho0Rho0 + u*tauc*rhom*dhRho0Tau0) + P_FR
    rhom0 = rhom; rhom = rhom+step

    ! secant method
    Do While(Abs(step) >= eps*100.0_pr)
       iter = iter+1
       tauc = CK*rhom**c23; u = kfconstmpi*rhom**c13; rhom2 = rhom**2
       Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
       If(finite_range) Then
          P_FR = P_SNM_FR(rhom)
       Else
          P_FR = 0
       End If
       aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
       aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
       w=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom + 3.0_pr*aRho0Rho0 + 5.0_pr*tauc*aRho0Tau0 &
                    +6.0_pr*rhom*hRho0Rho0 + 8.0_pr*tauc*rhom*hRho0Tau0 + 3.0_pr*Crho(0) &
                    +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0) + 5.0_pr*tauc*Ctau(0) + u*daRho0Rho0 &
                    +u*tauc*daRho0Tau0 + u*rhom*dhRho0Rho0 + u*tauc*rhom*dhRho0Tau0) + P_FR
       step = -w*(rhom-rhom0)/(w-w0)
       rhom0 = rhom; w0 = w; rhom = rhom+step
       If(iter > 100) Stop 'STOP(In find_NM_RHOC)'
    End Do
    find_NM_RHOC=rhom

  End Function find_NM_RHOC
  !=======================================================================
  !> This routine defines the transformations from the \f$ (t,x) \f$ parametrization of the Skyrme
  !> pseudopotential to the Skyrme energy density functional. The latter is given by the following
  !> form
  !>   \f{align*}{
  !>     \mathcal{H}_{t}^{\text{even}} & =
  !>     C_{t}^{\rho\rho}\rho_{t}^{2}
  !>     +
  !>     C_{t}^{\rho\Delta\rho} \rho_{t}\Delta\rho_{t}
  !>     +
  !>     C_{t}^{\rho\tau} \rho_{t}\tau_{t}
  !>     +
  !>     C_{t}^{\rho\nabla J} \rho_{t}\gras{\nabla}\cdot\gras{J}_{t} \medskip\\
  !>     & +
  !>     C_{t}^{JJ} \sum_{\mu\nu} J^{t}_{\mu\nu}J_{t,\mu\nu}
  !>     +
  !>     C^{J\bar{J}}_t \,
  !>     \biggl[ \Bigl( \sum_\mu \tensor{J}_{t,\mu\mu}\Bigr)
  !>             \Bigl( \sum_\mu \tensor{J}_{t,\mu\mu}\Bigr)
  !>     +
  !>     \sum_{\mu\nu} \tensor{J}_{t,\mu\nu}\tensor{J}_{t,\nu\mu} \biggr]
  !>     \medskip\\
  !>     \mathcal{H}_{t}^{\text{odd}}  & =
  !>     C_{t}^{ss}\gras{s}_{t}^{2}
  !>     +
  !>     C_{t}^{s\Delta s} \gras{s}_{t}\cdot\Delta\gras{s}_{t}
  !>     +
  !>     C_{t}^{sj} \gras{j}^{2}_{t}
  !>     +
  !>     C_{t}^{sT} \gras{s}_{t}\cdot\gras{T}_{t}  \medskip\\
  !>     & +
  !>     C_{t}^{s\nabla j} \gras{s}_{t}\cdot\left( \gras{\nabla}\times\gras{j}_{t} \right)
  !>     +
  !>     C_{t}^{s\nabla s} \left( \nabla\cdot\gras{s}_{t} \right)^{2}
  !>     +
  !>     C_{t}^{sF} \gras{s}_{t}\cdot\gras{F}_{t}
  !>   \f}
  !> and the transformation coefficients can be found in Sec. 9.1.1 of \cite schunck2019energy.
  !> Note that for now, we use local gauge invariance to set the coupling constants \f$ C_{t}^{sj} \f$,
  !> \f$ C_{t}^{s\nabla j} \f$, \f$ C_{t}^{sT} \f$ and \f$ C_{t}^{sF} \f$; see also discussion about
  !> the tensor term in that same reference.
  !=======================================================================
  Subroutine C_from_t()

    ! Time-even coupling constants
    Crho(0)  =   3.0_pr/8.0_pr  * t0
    Crho(1)  = -(1.0_pr/4.0_pr) * t0*(0.50_pr+x0)
    Cdrho(0) =  (1.0_pr/16.0_pr)* t3
    Cdrho(1) = -(1.0_pr/24.0_pr)* t3*(0.50_pr+x3)
    Ctau(0)  =  (3.0_pr/16.0_pr)* t1 + (1.0_pr/4.0_pr)*t2*(5.0_pr/4.0_pr+x2)
    Ctau(1)  = -(1.0_pr/8.0_pr) * t1*(0.5+x1)+(1.0_pr/8.0_pr)*t2*(0.50_pr+x2)
    CrDr(0)  =  (1.0_pr/16.0_pr)* t2*(5.0_pr/4.0_pr+x2) - (9.0_pr/64.0_pr)*t1
    CrDr(1)  =  (3.0_pr/32.0_pr)* t1*(0.5+x1) + (1.0_pr/32.0_pr)*t2*(0.50_pr+x2)
    CJ(0)    = -(1.0_pr/16.0_pr)*(t1*(2.0_pr*x1-1.0_pr) + t2*(2.0_pr*x2+1)-5*te-15*to)
    CJ(1)    = -(1.0_pr/16.0_pr)*(t2 -t1 + 5.0_pr*te -5.0_pr*to )
    CJbar(0) = -(3.0_pr/8.0_pr) *(te+3.0_pr*to)
    CJbar(1) = -(3.0_pr/8.0_pr) *(te-to)
    CrdJ(0)  = -b4-(0.50_pr)*b4p
    CrdJ(1)  = -0.50_pr*b4p

    ! Time-odd coupling constants
    Cs(0)    = -(1.0_pr/4.0_pr) * t0*(0.5_pr - x0)
    Cs(1)    = -(1.0_pr/8.0_pr) * t0
    Cds(0)   = -(1.0_pr/24.0_pr)* t3*(0.5_pr - x3)
    Cds(1)   = -(1.0_pr/48.0_pr)* t3
    Csj(0)   = -Ctau(0)
    Csj(1)   = -Ctau(1)
    CsDs(0)  =  (3.0_pr/32.0_pr)* t1*(0.5_pr - x1) + (1.0_pr/32.0_pr)* t2*(0.5_pr + x2)
    CsDs(1)  =  (3.0_pr/64.0_pr)* t1 + (1.0_pr/64.0_pr)* t2
    CsdJ(0)  =  CrdJ(0)
    CsdJ(1)  =  CrdJ(1)
    cT(0)    = -CJ(0)
    cT(1)    = -CJ(1)
    CsNabs(0)= -(3.0_pr/32.0_pr)*(te-to)
    CsNabs(1)= -(1.0_pr/32.0_pr)*(3.0_pr*te+to)
    CF(0)    = -2.0_pr*CJbar(0)
    CF(1)    = -2.0_pr*CJbar(1)

  End Subroutine C_from_t
  !=======================================================================
  !> This routine computes the (t,x) parameterization of the Skyrme effective interaction from the coupling
  !> constants of the energy functional; see Chapter 9 of \cite schunck2019energy for additional discussion.
  !=======================================================================
  Subroutine t_from_C()

    t0  =  (8.0_pr/3)*Crho(0)
    t1  =  4.0_pr/3.0_pr*(Ctau(0)-4.0_pr*CrDr(0))
    t2  =  4.0_pr/3.0_pr*(3.0_pr*Ctau(0)-6.0_pr*Ctau(1)+4.0_pr*CrDr(0)-8.0_pr*CrDr(1))
    t3  =  16.0_pr*Cdrho(0)
    x0  = -0.50_pr*(3.0_pr*Crho(1)/Crho(0)+1.0_pr)
    x1  =  2.0_pr*(-Ctau(0)-3.0_pr*Ctau(1)+4.0_pr*CrDr(0)+12.0_pr*CrDr(1))/t1/3.0_pr
    x2  = -2.0_pr*(3.0_pr*Ctau(0)-15.0_pr*Ctau(1)+4.0_pr*CrDr(0)-20.0_pr*CrDr(1))/t2/3.0_pr
    x3  = -0.50_pr*(3.0_pr*Cdrho(1)/Cdrho(0)+1.0_pr)
    b4  =  CrdJ(1)-CrdJ(0)
    b4p = -2.0_pr*CrdJ(1)
    te  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)-9.0_pr*CJ(1)-4.0_pr*CrDr(0)+12.0_pr*CrDr(1)-2.0_pr*Ctau(0)+6.0_pr*Ctau(1))
    to  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)+3.0_pr*CJ(1)+4.0_pr*CrDr(0)+4.0_pr*CrDr(1))
    If(finite_range .And. .not.force_is_dme) Then
       x0 = zero
       x1 = zero
       x2 = zero
    End If
    ! Set public variables
    t0_pub = t0; t1_pub  = t1;  t2_pub = t2; t3_pub = t3
    x0_pub = x0; x1_pub  = x1;  x2_pub = x2; x3_pub = x3
    b4_pub = b4; b4p_pub = b4p; te_pub = te; to_pub = to

  End Subroutine t_from_C
  !=======================================================================
  !> This function defines an external field, here a one-body potential that has the form of an
  !> harmonic trap. In cylindrical coordinates \f$ (\rho,z) \f$, we have
  !>   \f[
  !>        V_{\rm ext.}(\boldsymbol{r}) = E_{0} + \frac{1}{2}\hbar\Omega_{z}z^2 + \hbar\Omega_{\perp}\rho^2
  !>   \f]
  !> where \f$ E_{0} \f$ gives the scale and \f$ (\Omega_{z},\Omega_{\perp}) \f$ give the slopes of the
  !> constraining potential.
  !> UNUSED FOR NOW
  !=======================================================================
  Elemental Function Vexternal(t,x,y,z)
    Integer(ipr), Intent(In) :: t  ! isospin index: 0=isoscalar, 1=isovector
    Real(pr), Intent(In) :: x,y,z  ! position in cartesian basis
    Real(pr) :: Vexternal

    Vexternal = Real(t,Kind=pr)*(x**2+y**2+z**2)*0.0_pr

  End Function Vexternal
  !=======================================================================
  !
  !=======================================================================
End Module UNEDF
