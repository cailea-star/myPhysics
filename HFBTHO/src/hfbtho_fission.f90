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
!                  FISSION FRAGMENT PROPERTIES PACKAGE                 !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------!
!>  This module contains a number of routines and functions dealing with the properties of
!>  the fission fragments. It contains the  determination of the position of the neck between
!>  the two prefragments, the calculation of the matrix elements and expectation value of the
!>  Gaussian neck operator, as well as the expectation value of the multipole moments in the
!>  intrinsic frame of each prefragment. For l = m = 0, the latter can be used to define the
!>  number of protons and total number of particles in each prefragment.
!>
!> @author
!> Nicolas Schunck
!----------------------------------------------------------------------
!  Subroutines: - neck_computeField(ib)
!               - QNFIND
!               - DEFMAS(NZMAXX,ZPOINT,SFACTO)
!               - center_of_mass(Z_POSI,CENLEF,CENRIG)
!               - QLMFRA(Z_POSI,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
!               - DEFSPH(LAMACT,COSTHE,SPHHAR)
!               - print_moments(iw)
!  Functions: - DERIVE(Z_POSI)
!             - QMOM_K(Z_POSI,KORDER)
!----------------------------------------------------------------------!
Module HFBTHO_fission_fragments

  Use HFBTHO_utilities

  Implicit None

  Public neck_computeField,QNFIND,test_density,neck,wave_localization,center_of_mass,QLMFRA,print_moments, &
         calc_qh,average_J_in_fragments

  Private DERIVE,QMOM_K,neck_expectation,DEFMAS,DEFSPH,calc_lx

  Logical, Public, Save :: qp_occupation = .True.

  ! Constraints on the neck
  Real(pr), Public, Save :: AN_VAL=1.0_pr !< Range in fermis in the Gaussian neck operator
  Real(pr), Public, Save :: Q_NECK=0.0_pr !< Value of the Gaussian neck operator
  Real(pr), Public, Save :: Z_NECK=0.0_pr !< Position of the neck (dimensionless)
  Real(pr), Public, Save :: CENLEF !< Position of the center of mass for the left fragment (dimensionless)
  Real(pr), Public, Save :: CENRIG !< Position of the center of mass for the right fragment (dimensionless)
  ! Fission fragments
  Real(pr), Allocatable, Public, Save :: QLMTOT(:,:) !< Arrays containing \f$ \langle \hat{Q}_{\ell} \rangle \f$ for
                                                     !< isoscalar pseudo-densities. First index is \f$ \ell \f$, second index is 0
                                                     !, for the left fragment, 1 for the right fragment
  Real(pr), Allocatable, Public, Save :: QLMPRO(:,:) !< Same as QLMTOT(:,:) but expectation values are computed with charge pseudo-densities

  ! Private variables
  Integer(ipr), Private, Save :: debug_fission=0 !< Internal flag to activate/deactivate debugging messages
  Integer(ipr), Private, Save :: NPOINT=500 !< Number of points used to interpolate the densities along the z-axis
  Real(pr), Private :: ZNMINI=0.0_pr !< z-lower bound to search for the neck position (will be overwritten)
  Real(pr), Private :: ZNMAXI=0.0_pr !< z-upper bound to search for the neck position (will be overwritten)

Contains
  !=======================================================================
  !> Calculates the matrix of the Gaussian neck operator in the current block.
  !=======================================================================
  Subroutine neck_computeField(ib,gaussian_neck)
    Use HFBTHO_Gauss, Only : nghl,fh
    Use HFBTHO, Only : n00max,bz,id,ia,NS,QHLA_opt

    Integer(ipr), Intent(in) :: ib !< - Index of the current block
    Real(pr), Allocatable, Intent(Inout) :: gaussian_neck(:) !< - Linearized matrix of the neck operator in the current block

    Integer(ipr) :: i,nd,nd2,ihli,ihil,im,n1,n2
    Integer(ipr) :: ja,nsa,nsb,nsab,ssu,ssd
    Real(pr) :: z,qhla,vh,fiun1,fiun2,fidn1,fidn2,vnhl
    Real(pr), Allocatable :: Vmom(:)
    Integer(ipr) :: ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax)

    ! Compute value of neck operator on integration mesh
    Allocate(Vmom(1:nghl))
    Do ihli = 1,nghl
       z = fh(ihli)
       Vmom(ihli) = Exp(-((z-Z_NECK*bz)/AN_VAL)**2)
    End Do !ihli

    ! Form matrix of the Gaussian neck operator in HO basis
    nd=id(ib); nd2=nd*nd; im=ia(ib)
    ! sum over gauss integration points
    Do ihil=1,nghl
       vnhl=Vmom(ihil)
       ! scan over basis states
       Do n1=1,nd
          ja=n1+im; nsa=NS(ja); ssu=Max(nsa,0); ssd=Max(-nsa,0)
          QHLA=QHLA_opt(ja,ihil)
          OMPFIU(N1)=QHLA*ssu
          OMPFID(N1)=QHLA*ssd
       End Do
       i=0
       Do n1=1,nd
          ja=n1+im; nsa=NS(ja)
          fiun1=OMPFIU(N1); fidn1=OMPFID(N1)
          Do n2=1,n1
             i=i+1; nsb=NS(n2+im); nsab=nsa+nsb; vh=0.0_pr
             If(nsab /= 0) Then
                If(nsb > 0) Then
                   fiun2 = OMPFIU(N2)
                   vh    = fiun1*fiun2
                Else
                   fidn2 = OMPFID(N2)
                   vh    = fidn1*fidn2
                End If
                gaussian_neck(i)=gaussian_neck(i)+vh*vnhl
             End If
          End Do !n2
       End Do !n1
    End Do !ihil

    Deallocate(Vmom)

  End Subroutine neck_computeField
  !====================================================================
  !>  The subroutine finds the position of the neck, where
  !>  \f$ \langle \hat{Q}_{N}\rangle \f$ is the lowest.  At the minimum,
  !>  all partial derivatives of \f$ \langle \hat{Q}_{N}\rangle \f$ with
  !>  respect to \f$ x_{N}, y_{N}, z_{N} \f$ (coordinates of the neck)
  !>  should be 0. We find these 0 by using a variant of the Newton method.
  !====================================================================
  Subroutine QNFIND()
    Use math, Only : ZBRENT
    Use HFBTHO, Only : iiter,neckValue,Parity,lambdaMax,neck_constraints,xmix,neckLag
    USe HFBTHO_multipole_moments, Only : qmoment

    Logical :: symmetrical,enforced_symmetry,any_odd_multipole
    Integer(ipr) :: i,l,IFOUND,N_SCAN,KORDER,found
    Real(pr) :: QN_VAL,QBEFOR,ZBEFOR,DQNPRE,DQNCUR,XLOWER,XUPPER
    Real(pr) :: dq,Z_INIT_back,Qtotal
    Real(pr) :: Z_MINI,Z_MAXI,Z_INIT,Q_INIT,Z_POSI,TOLERA

    ! Initialization of the neck coordinates from previously found solutions.
    ! At the very first iteration, (X_NECK,Y_NECK,Z_NECK)=(0,0,0)
    If(debug_fission >= 1) write(*,'("Iteration: ",i4," initial zN: ",f20.14," qN: ", 2F20.14)') iiter,Z_NECK,neckValue
    Z_INIT=Z_NECK; Q_INIT=0.0_pr

    symmetrical = .True.

    ! Test if symmetry flags enforce symmetric nucleus
       enforced_symmetry = .True.
    If(.Not.Parity) Then
       enforced_symmetry = .False.
    End If
    ! Test if all multipole moments are actually zero
    any_odd_multipole = .False.
    Do l=1,lambdaMax,2
       Qtotal = qmoment(l,1) + qmoment(l,2)
       If(Abs(Qtotal) > 1.e-4) Then
          any_odd_multipole =.True.
       End If
    End Do
    ! System is symmetric either if symmetry flags enforce it, or if self-consistency makes it so
    If(.Not.enforced_symmetry) Then
       symmetrical = .Not.any_odd_multipole
    Else
       symmetrical = enforced_symmetry
    End If

    KORDER=0

    ! Initial scanning of the values of QN: we need this to get a somewhat
    ! reliable initial estimate of the position of the neck (only if
    ! reflection symmetry is broken). We search around the last known
    ! position of the neck.
    If(.Not.symmetrical) Then

       DQNPRE=0.1_pr; N_SCAN=100
       Z_MINI=-5.0_pr; Z_MAXI=+5.0_pr
       dq = (Z_MAXI-Z_MINI)/Real(N_SCAN-1,Kind=pr)
       ZBEFOR=Z_MINI; QBEFOR=QMOM_K(ZBEFOR-dq,KORDER)

       found = 0
       Do i=1,N_SCAN
          Z_POSI=Z_MINI+dq*Real(i-1,Kind=pr)
          QN_VAL=QMOM_K(Z_POSI,KORDER)
          DQNCUR=DERIVE(Z_POSI)
          If(debug_fission >= 2) &
             Write(*,'("z=",f20.14," QN=",f20.14," dQN(analytic)=",f20.14," dQN(numerical)=",f20.14)') &
                        Z_POSI,QN_VAL,DQNCUR,(QN_VAL-QBEFOR)/dq
          ! Found the minimum
          If(DQNCUR > 0.0_pr .And. DQNPRE < 0.0_pr) Then
             Z_INIT=0.5_pr*(Z_POSI+ZBEFOR)
             ZNMINI=Z_INIT-0.5_pr*dq; ZNMAXI=Z_INIT+0.5_pr*dq
             Q_INIT=QN_VAL
             found =found + 1
          End If
          ZBEFOR=Z_POSI
          QBEFOR=QN_VAL
          DQNPRE=DQNCUR
       End Do

       If(debug_fission >= 1) Then
          Write(*,'("NUMITE=",i4," Z_NECK=",f20.14)') iiter,Z_NECK
          Write(*,'("ZNMINI=",f20.14," ZNMAXI=",f20.14)') ZNMINI,ZNMAXI
          Write(*,'("found=",i4," Z_POSI=",f20.14," QN=",f20.14)') found,Z_INIT,Q_INIT
       End If

       ! Search for the z-coordinate of the neck (most likely to be different from 0)
       If(found > 0) Then
          Z_INIT_back = Z_INIT ! backup of rough estimate
          XLOWER=ZNMINI; XUPPER=ZNMAXI; TOLERA=1.D-5; IFOUND=0
          Z_INIT=ZBRENT(DERIVE,XLOWER,XUPPER,TOLERA,IFOUND)
          If(debug_fission >= 1) Write(*,'("found (Newton)=",i4," Z_INIT=",f20.14)') ifound,Z_INIT
          ! If ZBRENT gives a value too different from the rough estimate, there is a problem...
          If(Abs(Z_INIT-Z_INIT_back) > 0.5_pr*Abs(dq)) Then
             If(debug_fission >= 1) Write(*,'("Initial guess: ",f20.14," ZBRENT found: ", F20.14)') Z_INIT_back,Z_INIT
             Z_INIT = Z_INIT_back
          End If
          If(IFOUND == 0) Z_INIT = Z_INIT_back
       Else
          ! If no candidate for the neck position was found, we default to 0
          Z_INIT=0.0_pr
       End If

       If(neck_constraints) Then
          If(iiter >= 1) Then
             Z_NECK = xmix*Z_INIT + (one-xmix)*Z_NECK
          Else
             If(Abs(neckLag) > 1.e-12_pr) Then
                Z_NECK = xmix*Z_INIT + (one-xmix)*Z_NECK
             Else
                Z_NECK = Z_INIT
             End If
          End If
       Else
          Z_NECK = Z_INIT
       End If
       Q_NECK = QMOM_K(Z_NECK,KORDER)
       If(debug_fission >= 1) Write(*,'("Iteration: ",i4," final zN: ",f20.14," qN: ", 2F20.14)') iiter,Z_NECK,Q_NECK

    Else
       ! Reflection symmetry is conserved, the neck is at z=0, we compute QN
       Z_NECK=0.0_pr
       Q_NECK=QMOM_K(Z_NECK,KORDER)
    End If

    neckValue = Q_NECK

  End Subroutine QNFIND
  !====================================================================
  !>  This function computes the first derivative of the Gaussian neck
  !> expectation value with respect to the z-coordinate of the neck.
  !====================================================================
  Real(pr) Function DERIVE(Z_POSI)
    Use HFBTHO, Only : bz

    Real(pr), Intent(In) :: Z_POSI !< - Current position of the neck (dimensionless)

    Integer(ipr) :: KORDER
    Real(pr) :: DERIV1,PREFAC

    DERIV1 = 0.0_pr
    PREFAC = 2.0_pr*bz*bz/AN_VAL**2
    KORDER = 1; DERIV1 = QMOM_K(Z_POSI,KORDER)
    KORDER = 0; DERIV1 = DERIV1 - Z_POSI*QMOM_K(Z_POSI,KORDER)
    DERIVE = PREFAC*DERIV1

  End Function DERIVE
  !====================================================================
  !>  This function computes the moments \f$ I_k \f$ of order k of the
  !>  Gaussian neck operator, i.e., the object
  !>    \f[
  !>        I_k = \int_{-\infty}^{+\infty} d\xi\, \xi^{k} \rho(\gras{r}) \hat{Q}_{\rm N}
  !>    \f]
  !====================================================================
  Real(pr) Function QMOM_K(Z_POSI,KORDER)
    Use HFBTHO_Gauss, Only : ngl,ngh,xh,wdcor
    Use HFBTHO, Only : bz,ro

    Integer(ipr), Intent(In) :: KORDER
    Real(pr), Intent(In) :: Z_POSI

    Integer(ipr) :: il,ih,ihli
    Real(pr) :: QN_VAL,ZNECKZ,COORDO,PREFAC

    QN_VAL=0.0_pr
    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          ZNECKZ=Exp(-((xh(ih)-Z_POSI)*bz/AN_VAL)**2)
          COORDO=xh(ih)
          If(KORDER == 0) Then
             PREFAC=1.0_pr
          Else
             PREFAC=COORDO**KORDER
          End If
          QN_VAL=QN_VAL+wdcor(ihli)*(ro(ihli,1)+ro(ihli,2))*ZNECKZ*PREFAC
       End Do
    End Do

    QMOM_K=QN_VAL

  End Function QMOM_K
  !====================================================================
  !> This function simply tests the normalization of the total density.
  !====================================================================
  Real(pr) Function test_density()
    Use HFBTHO_Gauss, Only : ngl,ngh,wdcor
    Use HFBTHO, Only : ro

    Integer(ipr) :: il,ih,ihli
    Real(pr) :: QN_VAL

    QN_VAL=0.0_pr
    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          QN_VAL=QN_VAL+wdcor(ihli)*(ro(ihli,1)+ro(ihli,2))
       End Do
    End Do

    test_density=QN_VAL

  End Function test_density
  !=======================================================================
  !> This routine computes the contribution of a given K-block to the expectation value of the Gaussian
  !> neck operator in configuration space, that is
  !>  \f[
  !>      \langle \hat{Q}_{\rm N} \rangle_{K} = \sum_{ij=1}^{n_{K}} \hat{Q}_{\rm N,ij}\rho^{K}_{ji}
  !>  \f]
  !=======================================================================
  Subroutine neck_expectation(ib,qval,rho,dd)
    Use HFBTHO, Only : id

    Integer(ipr), Intent(In) :: ib
    Real(pr), Allocatable, Intent(In) :: rho(:,:)
    Real(pr), Intent(Inout) :: qval,dd

    Integer(ipr) :: nd,nd2,nhfb,j,n1,n2
    Real(pr) :: hla
    Real(pr), Allocatable :: dblmul(:,:),gaussian_neck(:),qblock(:,:)

    nd=id(ib); nd2=nd*nd; nhfb=nd+nd;
    Do n1=1,nd
       dd=dd+rho(n1,n1)
    End Do

    Allocate(dblmul(nd,nd)); dblmul=zero
    If(Allocated(gaussian_neck)) Deallocate(gaussian_neck)
    Allocate(gaussian_neck(1:nd2))
    gaussian_neck=zero
    Call neck_computeField(ib,gaussian_neck)
    j=0
    Do n1=1,nd
       Do n2=1,n1
          j=j+1;hla=gaussian_neck(j)
          dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
       End Do
    End Do
    Allocate(qblock(nd,nd))
    Call dgemm('n','n',nd,nd,nd,one,rho,nd,dblmul,nd,zero,qblock,nd)
    Do n1=1,nd
       qval=qval+qblock(n1,n1)
    End Do

  End Subroutine neck_expectation
  !====================================================================
  !> This routine computes the expectation value of the Gaussian neck operator in configuration space,
  !> that is, as \f$ \langle \hat{Q}_{\rm N} \rangle = \text{Tr} \big( \hat{Q}_{\rm N}\rho \big) \f$.
  !=======================================================================
  Subroutine neck()
    Use HFBTHO, Only : nb,id,kd,ka,RVqpN,RUqpN,KpwiN,RVqpP,RUqpP,KpwiP

    Integer(ipr) :: it,ib,nd,nd2,n_qp,k,n2
    Real(pr) :: qvala,dd2
    Real(pr), Allocatable :: rho(:,:), Umatr(:,:),Vmatr(:,:)

    qvala=zero; dd2=zero

    it = 1

    Do ib = 1,nb
       ! Get characteristics of the current block: size of HO basis and number of qp
       nd=id(ib); nd2=nd*nd; n_qp=kd(ib,it)
       If(n_qp > 0) Then
          Allocate(Umatr(nd,n_qp)); Umatr=zero
          Allocate(Vmatr(nd,n_qp)); Vmatr=zero
          Do k=1,n_qp
             Do n2=1,nd
                Vmatr(n2,k)=RVqpN(KpwiN(ka(ib,it)+k)+n2)
                Umatr(n2,k)=RUqpN(KpwiN(ka(ib,it)+k)+n2)
             End Do
          End Do
          ! In configuration space, test that Tr(rho) = N (contained in dd1 and dd2), and that
          ! Tr(Q*rho) = Q (in qvala and qvalb)
          Allocate(rho(nd,nd)); rho=zero
          Call dgemm('n','t',nd,nd,n_qp,one,Vmatr,nd,Vmatr,nd,zero,rho,nd)
          Call neck_expectation(ib,qvala,rho,dd2)
          Deallocate(rho)
          Deallocate(Umatr,Vmatr)
       End If
    End Do

    it = 2

    Do ib = 1,nb
       ! Get characteristics of the current block: size of HO basis and number of qp
       nd=id(ib); nd2=nd*nd; n_qp=kd(ib,it)
       If(n_qp > 0) Then
          Allocate(Umatr(nd,n_qp)); Umatr=zero
          Allocate(Vmatr(nd,n_qp)); Vmatr=zero
          Do k=1,n_qp
             Do n2=1,nd
                Vmatr(n2,k)=RVqpP(KpwiP(ka(ib,it)+k)+n2)
                Umatr(n2,k)=RUqpP(KpwiP(ka(ib,it)+k)+n2)
             End Do
          End Do
          ! In configuration space, test that Tr(rho) = N (contained in dd1 and dd2), and that
          ! Tr(Q*rho) = Q (in qvala and qvalb)
          Allocate(rho(nd,nd)); rho=zero
          Call dgemm('n','t',nd,nd,n_qp,one,Vmatr,nd,Vmatr,nd,zero,rho,nd)
          Call neck_expectation(ib,qvala,rho,dd2)
          Deallocate(rho)
          Deallocate(Umatr,Vmatr)
       End If
    End Do

    Write(*,'("qvala=",f14.7," dd2=",f14.7)') 2*qvala,2*dd2

  End Subroutine neck
  !====================================================================
  !
  !=======================================================================
  Subroutine wave_localization(ZPOINT)
    Use HFBTHO, Only : nzx,nb,id,kd,ka,ia,nr,nl,ns,nz,RVqpN,KpwiN,SFACTO
    USe HFBTHO_gauss, Only : DEVHER

    Real(pr), Intent(In) :: ZPOINT !< - Position of the neck \f$ \xi \f$ (dimensionless)
    Integer(ipr) :: it,ib,nd,n_qp,k,n1,n2,im
    Integer(ipr) :: nr2,nl2,ns2,nz2,nr1,nl1,ns1,nz1
    Real(pr) :: cN

    Call DEVHER(nzx)

    If(Allocated(SFACTO)) Deallocate(SFACTO)
    Allocate(SFACTO(0:nzx,0:nzx)); SFACTO=0.0_pr
    Call DEFMAS(nzx,ZPOINT,SFACTO)

    it = 1

    Do ib = 1,nb
       ! Get characteristics of the current block: size of HO basis, number of qp and running index im
       nd=id(ib); n_qp=kd(ib,it); im=ia(ib)
       If(n_qp > 0) Then
          Do k=1,n_qp
             cN = zero
             Do n2=1,nd
                nr2=nr(im+n2); nl2=nl(im+n2); ns2=ns(im+n2); nz2=nz(im+n2)
                Do n1=1,nd
                   nr1=nr(im+n1); nl1=nl(im+n1); ns1=ns(im+n1); nz1=nz(im+n1)
                   If(nr1 == nr2 .And. nl1 == nl2 .And. ns1 == ns2) Then
                      cN = cN + RVqpN(KpwiN(ka(ib,it)+k)+n1)*RVqpN(KpwiN(ka(ib,it)+k)+n2)*SFACTO(nz1,nz2)
                   End If
                End Do
             End Do
          End Do
       End If
    End Do

  End Subroutine wave_localization
  !====================================================================
  !>  This subroutine computes
  !>   \f[
  !>     S_{m_{z}n_{z}} = \sum_{k}^{n_z+m_z} C_{n_{z}m_{z}}^{k}(00)
  !>              \int_{z_{N}}^{+\infty} d\xi H_{k}(\xi) e^{-\xi^2}
  !>   \f]
  !>  This is needed to evaluate the  expectation value of the
  !>  density in the fission fragments.
  !====================================================================
  Subroutine DEFMAS(NZMAXX,ZPOINT,SFACTO)
    Use math, Only : PROINT
    Use HFBTHO_Gauss, Only : D_HERM,COEF00,HERFAC

    Integer(ipr), Intent(In) :: NZMAXX !< - Maximum value for either \f$ n_{z} \f$ or \f$ m_{z} \f$
    Real(pr), Intent(In) :: ZPOINT !< - Position of the neck \f$ \xi \f$ (dimensionless)
    Real(pr), Allocatable, Intent(Inout) :: SFACTO(:,:) !< - Array containing \f$ S_{m_{z}n_{z}} \f$

    Integer(ipr) :: NORDER,NZ,MZ,K
    Real(pr) :: VALINT,XPOINT,SUMVAL,F_INTE
    Real(pr), Allocatable :: PHERMI(:),DHERMI(:),ZHERMI(:)

    ! Computing the value of the probability distribution function
    VALINT=0.0_pr; VALINT=PROINT(ZPOINT)

    ! Compute value of Hermite polynomials up to order Max(Nz+Mz) at z=0
    XPOINT=zero; NORDER=2*NZMAXX
    Allocate(PHERMI(1:NORDER+1)); Allocate(DHERMI(1:NORDER+1))
    Call D_HERM(XPOINT,NORDER,PHERMI,DHERMI,NORDER+1)

    ! Compute value of Hermite polynomials up to order Max(Nz+Mz) at z=Abs(ZPOINT)
    XPOINT=Abs(ZPOINT); NORDER=2*NZMAXX
    Allocate(ZHERMI(1:NORDER+1))
    Call D_HERM(XPOINT,NORDER,ZHERMI,DHERMI,NORDER+1)

    SUMVAL = 0.0_pr
    ! Loop over bra and ket z quantum number
    Do NZ=0,NZMAXX
       Do MZ=0,NZMAXX
          ! Initial value of the sum for K = 0
          SUMVAL = COEF00(0,MZ,NZ)*0.5_pr*Sqrt(Sqrt(pi))*(1.0_pr - VALINT)
          If(ZPOINT < 0.0_pr) SUMVAL = COEF00(0,MZ,NZ)*half*Sqrt(Sqrt(pi))*(1.0_pr + VALINT)
          Do K=1,NZ+MZ
             ! Hermite polynomial of order k is stored in PHERMI(k+1), but this shift does
             ! *not* apply to its normalization coefficient, which is stored in HERFAC(k)
             ! as expected. Below, we need Hermite polynomials of order k-1 and normalization
             ! coefficients of order k
             If(Mod(K,2) == 1) Then
                F_INTE = ZHERMI(K)*Exp(-ZPOINT*ZPOINT)/HERFAC(K)
             Else
                F_INTE =(ZHERMI(K)*Exp(-ZPOINT*ZPOINT) - PHERMI(K)) / HERFAC(K)
                If(ZPOINT < 0.0_pr) F_INTE = - F_INTE
             End If
             SUMVAL = SUMVAL + COEF00(K,MZ,NZ)*F_INTE
          End Do
          SFACTO(MZ,NZ)=SUMVAL
       End Do
    End Do

  End Subroutine DEFMAS
  !====================================================================
  !> This function computes position of the center of mass of the left
  !> and right fragments. The resulting quantity is in fermis.
  !====================================================================
  Subroutine center_of_mass(Z_POSI,CENLEF,CENRIG)
    Use math, Only : integrate_simpson,SPLINT,SPLINE
    Use HFBTHO_Gauss, Only : ngl,ngh,xh,wl
    Use HFBTHO, Only : bz,bp,ro

    Real(pr), Intent(In) :: Z_POSI !< - Position of the neck (dimensionless)
    Real(pr), Intent(Inout) :: CENLEF !< - Center of mass of the left fragment (dimensionless)
    Real(pr), Intent(Inout) :: CENRIG !< - Center of mass of the right fragment (dimensionless)

    Integer(ipr) :: IERROR,i,ih,il,ihli
    Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,RHOTMP,ZPOTMP,dh
    Real(pr) :: Z_ACTU,ZPOINT,DENLOC,W_HERM,RHONOR,Z_NORM
    Real(pr), Allocatable :: FUNCTI(:),GUNCTI(:),XVALUE(:),ZINTER(:),AUXSTO(:)
    Real(pr), Allocatable :: DINTER(:,:)

    Z_ACTU=Z_POSI

    ! Integrations over the transverse coordinates are done using the dimensionless units
    ! xi_x and xi_y. In the longitudinal direction, the density is defined as a function
    ! of the dimensionless variable xi_z = bz * z. The operator O(z) = z must thus be
    ! defined on the grid z = xi_z/bz with integration limits given below.
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=xh(1)*bz + 1.D-14
    XFINIS=Z_ACTU*bz
    ! Boundary test: make sure there is a valid integration interval
    If(Z_ACTU < XBEGIN) Then
       Write(*,'("ATTENTION! In center_of_mass(): changing neck position from ",f20.14," to ",f20.14)') &
                  Z_ACTU,XBEGIN + 1.D-12
       Z_ACTU = XBEGIN + 1.D-12
    End If

    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          XVALUE(ih)=xh(ih)*bz
          FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
       End Do
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))    /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR == 0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(*,'("In center_of_mass() - Error in interpolating the function (right fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)

    ! Computing the integrals \int \rho(z) and \int z\rho(z)
    !          -\infty < x < +\infty
    !          -\infty < y < +\infty
    !          -\infty < z < +z_{N}
    ! Integrations over x and y are performed "exactly" by Gauss-Hermite
    ! quadratures, integration over z is performed numerically using the
    ! Simpson 3/8 rule
    Allocate(FUNCTI(1:NPOINT),GUNCTI(1:NPOINT))
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; RHOTMP=0.0_pr
       GUNCTI(i)=0.0_pr; ZPOTMP=0.0_pr
       Do il=1,ngl
          ZPOINT=ZINTER(i)
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          RHOTMP=RHOTMP+W_HERM*DENLOC
          ZPOTMP=ZPOTMP+W_HERM*DENLOC*ZPOINT
       End Do
       FUNCTI(i)=RHOTMP; GUNCTI(i)=ZPOTMP
    End Do
    ! Integrating over z, from -infty to zN
    RHONOR=0.0_pr; Call integrate_simpson(FUNCTI,NPOINT,X_STEP,RHONOR)
    Z_NORM=0.0_pr; Call integrate_simpson(GUNCTI,NPOINT,X_STEP,Z_NORM)
    CENLEF=Z_NORM/RHONOR
    Deallocate(FUNCTI,GUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)

    ! Repeating this whole procedure for the right fragment
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=Z_ACTU*bz
    XFINIS=xh(ngh)*bz - 1.D-14
    If(Z_ACTU > XFINIS) Then
       Write(*,'("ATTENTION! In center_of_mass(): changing neck position from ",f20.14," to ",f20.14)') &
                  Z_ACTU,XFINIS - 1.D-12
       Z_ACTU = XFINIS - 1.D-12
    End If
    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          XVALUE(ih)=xh(ih)*bz
          FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
       End Do
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))       /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR == 0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(*,'("In center_of_mass() - Error in interpolating the function (left fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)

    Allocate(FUNCTI(1:NPOINT),GUNCTI(1:NPOINT))
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; RHOTMP=0.0_pr
       GUNCTI(i)=0.0_pr; ZPOTMP=0.0_pr
       Do il=1,ngl
          ZPOINT=ZINTER(i)
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          RHOTMP=RHOTMP+W_HERM*DENLOC
          ZPOTMP=ZPOTMP+W_HERM*DENLOC*ZPOINT
       End Do
       FUNCTI(i)=RHOTMP; GUNCTI(i)=ZPOTMP
    End Do
    RHONOR=0.0_pr; Call integrate_simpson(FUNCTI,NPOINT,X_STEP,RHONOR)
    Z_NORM=0.0_pr; Call integrate_simpson(GUNCTI,NPOINT,X_STEP,Z_NORM)
    CENRIG=Z_NORM/RHONOR

    Deallocate(FUNCTI,GUNCTI,AUXSTO,ZINTER,DINTER,XVALUE)

  End Subroutine center_of_mass
  !====================================================================
  !> This function computes the expectation value of the multipole moment
  !> operators \f$ Q_{\ell m} \f$ in the left and right fragments,
  !> assuming the position of the neck is \f$ z_{N} \f$.
  !====================================================================
  Subroutine QLMFRA(Z_POSI,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
    Use math, Only : integrate_simpson,SPLINT,SPLINE
    Use HFBTHO_Gauss, Only : ngl,ngh,xh,xl,wl
    Use HFBTHO, Only : bz,bp,ro

    Integer(ipr), Intent(In) :: LAMACT !< - \f$ \ell \f$ value
    Integer(ipr), Intent(In) :: I_TYPE !< - m value
    Real(pr), Intent(In) :: Z_POSI !< - Position of the neck (dimensionless)
    Real(pr), Intent(In) :: CENLEF !< - Center of mass of the left fragment (dimensionless)
    Real(pr), Intent(In) :: CENRIG !< - Center of mass of the right fragment (dimensionless)
    Real(pr), Intent(Inout) :: QLMLEF !< - \f$ \langle \hat{Q}_{\ell m} \rangle \f$ for the left fragment (dimensionless)
    Real(pr), Intent(Inout) :: QLMRIG !< - \f$ \langle \hat{Q}_{\ell m} \rangle \f$ for the right fragment (dimensionless)

    Integer(ipr) :: IERROR,i,ih,il,ihli
    Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,dh
    Real(pr) :: Z_ACTU,ZPOINT,DENLOC,W_HERM,QLMSUM,QLMVAL,RESULT
    Real(pr) :: RADIUS,COSTHE,epsilon,one
    Real(pr), Dimension(1:2) :: SPHHAR
    Real(pr), Allocatable :: FUNCTI(:),XVALUE(:),ZINTER(:),AUXSTO(:)
    Real(pr), Allocatable :: DINTER(:,:)

    Z_ACTU=Z_POSI

    ! Define boundaries for the integration over z, see comments in subroutine center_of_mass()
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=xh(1)*bz + 1.0e-14_pr
    XFINIS=Z_ACTU*bz
    ! Boundary test: make sure there is a valid integration interval
    If(Z_ACTU < XBEGIN) Then
       Write(*,'("ATTENTION! In QLMFRA(): changing neck position from ",f20.14," to ",f20.14)') Z_ACTU,XBEGIN+1.0e-12_pr
       Z_ACTU = XBEGIN + 1.0e-12_pr
    End If

    Do il=1,ngl
       ! See comments in subroutine center_of_mass()
       ! Total density
       If(I_TYPE == 1) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
          End Do
       End If
       ! Proton (=charge in HFBTHO) density
       If(I_TYPE == 2) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,2)
          End Do
       End If
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))    /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR == 0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(*,'("In QLMFRA() - Error in interpolating the function (right fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)

    ! Computing the expectation values of multipole moment for the left fragment
    Allocate(FUNCTI(1:NPOINT))
    epsilon=1.0e-14_pr; one=1.0_pr
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; QLMSUM=0.0_pr
       Do il=1,ngl
          ! Coordinates (x,y,z) for the multipole moments must be in fermis
          ZPOINT=ZINTER(i)-CENLEF ! shift with respect to the c.o.m. of left fragment
          RADIUS=Sqrt(xl(il)*bp**2+ZPOINT**2)
          ! Angle theta
          If(RADIUS <= epsilon) Then
             COSTHE=0.0_pr
          Else
             COSTHE=ZPOINT/RADIUS
          End If
          If(Abs(COSTHE-one) <= epsilon) Then
             COSTHE=one-2.0_pr*epsilon
          End If
          If(Abs(COSTHE+one) <= epsilon) Then
             COSTHE=2.0_pr*epsilon-one
          End If
          Call DEFSPH(LAMACT,COSTHE,SPHHAR)
          If(LAMACT == 0) Then
             QLMVAL=1.0_pr
          Else
             If(RADIUS >= epsilon) Then
                QLMVAL=RADIUS**(LAMACT)*SPHHAR(1)
             Else
                QLMVAL=0.0_pr
             End If
          End If
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          QLMSUM=QLMSUM+W_HERM*DENLOC*QLMVAL
       End Do
       FUNCTI(i)=QLMSUM
    End Do
    ! Integrating over z, from -infty to zN
    Call integrate_simpson(FUNCTI,NPOINT,X_STEP,RESULT)
    QLMLEF=RESULT
    Deallocate(FUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)

    ! Repeating the procedure for the right fragment
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=Z_ACTU*bz
    XFINIS=xh(ngh)*bz - 1.0e-14_pr
    ! Boundary test: make sure there is a valid integration interval
    If(Z_ACTU > XFINIS) Then
       Write(6,'("ATTENTION! In QLMFRA(): changing neck position from ",f20.14," to ",f20.14)') Z_ACTU,XFINIS-1.0e-12_pr
       Z_ACTU = XFINIS - 1.0e-12_pr
    End If
    Do il=1,ngl
       If(I_TYPE == 1) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
          End Do
       End If
       ! Proton (=charge in HFBTHO) density
       If(I_TYPE == 2) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,2)
          End Do
       End If
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))    /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR == 0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(*,'("In QLMFRA() - Error in interpolating the function (left fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)

    Allocate(FUNCTI(1:NPOINT))
    epsilon=1.0e-14_pr; one=1.0_pr
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; QLMSUM=0.0_pr
       Do il=1,ngl
          ! Coordinates (x,y,z) for the multipole moments must be in fermis
          ZPOINT=ZINTER(i)-CENRIG ! shift with respect to the c.o.m. of left fragment
          RADIUS=Sqrt(xl(il)*bp**2+ZPOINT**2)
          ! Angle theta
          If(RADIUS <= epsilon) Then
             COSTHE=0.0_pr
          Else
             COSTHE=ZPOINT/RADIUS
          End If
          If(Abs(COSTHE-one) <= epsilon) Then
             COSTHE=one-2.0_pr*epsilon
          End If
          If(Abs(COSTHE+one) <= epsilon) Then
             COSTHE=2.0_pr*epsilon-one
          End If
          Call DEFSPH(LAMACT,COSTHE,SPHHAR)
          If(LAMACT == 0) Then
             QLMVAL=1.0_pr
          Else
             If(RADIUS >= epsilon) Then
                QLMVAL=RADIUS**(LAMACT)*SPHHAR(1)
             Else
                QLMVAL=0.0_pr
             End If
          End If
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          QLMSUM=QLMSUM+W_HERM*DENLOC*QLMVAL
       End Do
       FUNCTI(i)=QLMSUM
    End Do
    Call integrate_simpson(FUNCTI,NPOINT,X_STEP,RESULT)
    QLMRIG=RESULT

    Deallocate(FUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)

  End Subroutine QLMFRA
  !====================================================================
  !> This routine computes the value of the multipole moment operators
  !> \f$ Q_{\ell 0}(\theta, \phi) = P_{\ell}(\cos\theta) \f$.
  !====================================================================
  Subroutine DEFSPH(LAMACT,COSTHE,SPHHAR)
    USe math, Only : DEFLEG
    Use HFBTHO_multipole_moments, Only : q_units

    Integer(ipr), Intent(In) :: LAMACT !< - \f$ \ell \f$ value
    Real(pr), Intent(In) :: COSTHE !< - \f$ \cos\theta \f$ value
    Real(pr), Dimension(1:2), Intent(Inout) :: SPHHAR !< - Value \f$ P_{\ell}(\cos\theta) \f$

    Integer(ipr) :: i
    Real(pr) :: ZLEGPO,PIARGU,FACMUL
    Real(pr), Allocatable :: FACTOR(:)

    ! Computing P_{l,m}(cos(theta))
    ZLEGPO=DEFLEG(LAMACT,0,COSTHE)

    ! Defining the factorials
    Allocate(FACTOR(0:LAMACT))
    FACTOR(0)=1.0_pr
    Do i=1,LAMACT
       FACTOR(I)=FACTOR(i-1)*Real(i,Kind=pr)
    End Do

    ! Computing the spherical harmonics
    PIARGU=4.0_pr*Atan(1.0_pr)

    FACMUL = Sqrt(0.25_pr*Real(2*LAMACT+1,Kind=pr)/PIARGU) * q_units(LAMACT)

    SPHHAR(1) = FACMUL*ZLEGPO
    SPHHAR(2) = 0.0_pr

    Deallocate(FACTOR)

  End Subroutine DEFSPH
  !====================================================================
  !> The routine prints the characteristics of the fission fragments.
  !====================================================================
  Subroutine print_moments(iw)
    Use HFBTHO, Only : neckValue,jxsq_rho_right,jxsq_rho_left,jxsq_kappa_left,jxsq_kappa_right

    Integer(ipr), Intent(in) :: iw !< - Fortran unit where data is printed

    Write(iw,'(2X,"Fission fragments characteristics")')
    Write(iw,'(2X,"=================================")')
    Write(iw,'(2X,"Gaussian neck .......",2f12.4)') neckValue,Z_NECK
    Write(iw,'(2X,"Centers of mass .....",2f12.4)') CENLEF,CENRIG
    Write(iw,'(2X,"     Observable         Left fragment  Right fragment")')
    Write(iw,'(2X,"Charge Z ............",f12.4,f15.4)') QLMPRO(0,0),QLMPRO(0,1)
    Write(iw,'(2X,"Mass A ..............",f12.4,f15.4)') QLMTOT(0,0),QLMTOT(0,1)
    Write(iw,'(2X,"q10 [b^1/2] .........",2f15.6)') QLMTOT(1,0),QLMTOT(1,1)
    Write(iw,'(2X,"q20 [b] .............",2f15.6)') QLMTOT(2,0),QLMTOT(2,1)
    Write(iw,'(2X,"q30 [b^3/2] .........",2f15.6)') QLMTOT(3,0),QLMTOT(3,1)
    Write(iw,'(2X,"q40 [b^2] ...........",2f15.6)') QLMTOT(4,0),QLMTOT(4,1)
    Write(iw,'(2X,"q50 [b^5/2] .........",2f15.6)') QLMTOT(5,0),QLMTOT(5,1)
    Write(iw,'(2X,"q60 [b^3] ...........",2f15.6)') QLMTOT(6,0),QLMTOT(6,1)
    Write(iw,'(2X,"q70 [b^7/2] .........",2f15.6)') QLMTOT(7,0),QLMTOT(7,1)
    Write(iw,'(2X,"q80 [b^4] ...........",2f15.6)') QLMTOT(8,0),QLMTOT(8,1)
    Write(iw,'(2X,"<Jx^2> from rho .....",2f15.6)') jxsq_rho_left,jxsq_rho_right
    Write(iw,'(2X,"<Jx^2> from kappa ...",2f15.6)') jxsq_kappa_left,jxsq_kappa_right
    Write(iw,'(2X,"<J^2>=2*<Jx^2> ......",2f15.6)') two*(jxsq_rho_left+jxsq_kappa_left),&
                                                    two*(jxsq_rho_right+jxsq_kappa_right)

  End Subroutine print_moments
  !==============================================================
  !> This subroutine calculates the z-component of the wave function at any given point.
  !==============================================================
  Subroutine calc_qh(zpoint)
    Use math, Only : sq,sqi
    Use HFBTHO, Only : qh_zneck,bz,nzm

    Real(pr), Intent(In) :: zpoint !< -  zpoint - z coordinate (dimensionless)

    Integer(ipr) :: n
    Real(pr) :: w4pii,z,w0

    If(Allocated(qh_zneck)) Deallocate(qh_zneck)
    Allocate(qh_zneck(0:nzm))

    w4pii=pi**(-0.250_pr)

    z=zpoint
    w0=pi**(-0.250_pr)*Exp(-half*z*z)/sqrt(bz)
    qh_zneck(0)=w0
    qh_zneck(1)=w0*sq(2)*z
    Do n=2,nzm
       qh_zneck(n)=sqi(n)*(sq(2)*z*qh_zneck(n-1)-sq(n-1)*qh_zneck(n-2))
    End Do

  End Subroutine calc_qh
  !==========================================================================
  !> This subroutine calculates the expectation value of the \f$ \hat{\ell}_x \f$ operator in the HO basis
  !> with respect to the center of mass zc of the fission fragment.
  !==========================================================================
  Subroutine calc_lx(nz1,nr1,nl1,nz2,nr2,nl2,zc_right,rlx_right,zc_left,rlx_left)
    Use math, Only : sq,sqi
    Use HFBTHO_Gauss, Only : ngl,sxl,ql,ql1
    Use HFBTHO, Only : bz,bp,qh_zneck,SFACTO,nzm

    Integer(ipr), Intent(In) :: nz1  !< - \f$ n_{z} \f$ quantum number of the left (bra) state
    Integer(ipr), Intent(In) :: nr1  !< - \f$ n_{r} \f$ quantum numbers of the left (bra) state
    Integer(ipr), Intent(In) :: nl1  !< - \f$ \Lambda \f$ quantum numbers of the left (bra) state
    Integer(ipr), Intent(In) :: nz2  !< - \f$ n_{z} \f$ quantum numbers of the right (ket) state
    Integer(ipr), Intent(In) :: nr2  !< - \f$ n_{r} \f$ quantum numbers of the right (ket) state
    Integer(ipr), Intent(In) :: nl2  !< - \f$ \Lambda \f$ quantum numbers of the right (ket) state
    Real(pr), Intent(In) :: zc_right !< - z coordinate of the CoM of the right fragment (in fm)
    Real(pr), Intent(In) :: zc_left  !< - z coordinate of the CoM of the left fragment (in fm)
    Real(pr), Intent(Out) :: rlx_right !< - matrix element of the \f$ \hat{\ell}_x \f$ operator in the right fragment
    Real(pr), Intent(Out) :: rlx_left  !< - matrix element of the \f$ \hat{\ell}_x \f$ operator in the left fragment

    Integer(ipr) :: il
    Real(pr) :: fz,fp,heavyside_at_zero,delnl1,delnl2,delnz1,delnz2,delnz3,fac1,fac2,fac3,facr1,facr2,facr3,&
                facz1_right,facz2_right,facz3_right,facz1_left,facz2_left,facz3_left

    rlx_right=zero; rlx_left=zero

    If(nl1 /= nl2+1 .And. nl1 /= nl2-1) Return

    fz=1.0_pr/bz; fp=1.0_pr/bp; heavyside_at_zero=0.5_pr
    delnl1=0.0_pr; delnl2=0.0_pr
    If(nl1 == nl2+1) delnl1=1.0_pr
    If(nl1 == nl2-1) delnl2=1.0_pr

    !----------------------------------------------------------
    ! integrals in rperp direction, the same for both fragments
    ! ---------------------------------------------------------
    facr1=zero; facr2=zero; facr3=zero
    Do il=1,NGL
       facr1 = facr1 + ql(nr1,nl1,il)*ql(nr2,nl2,il)*sxl(il)  ! linear in rperp
       facr2 = facr2 + ql(nr1,nl1,il)*ql(nr2,nl2,il)/sxl(il)  ! 1/rperp
       facr3 = facr3 + ql(nr1,nl1,il)*ql1(nr2,nl2,il)/sxl(il) ! derivation over rperp
    End Do

    facr1=facr1*bp; facr2=facr2*fp; facr3=facr3*fp        ! proper dimensions
    ! ----------------------------------------
    ! integrals in z direction, right fragment
    ! ----------------------------------------
    facz1_right = 0.0_pr
    If(nz2 < nzm) facz1_right = facz1_right - sq(nz2+1)*SFACTO(nz1,nz2+1)
    If(nz2>0)     facz1_right = facz1_right + sq(nz2)  *SFACTO(nz1,nz2-1)
    facz2_right = qh_zneck(nz1)*qh_zneck(nz2)*heavyside_at_zero
    facz3_right = -sq(2)*zc_right*fz*SFACTO(nz1,nz2)
    If(nz2 < nzm) facz3_right = facz3_right + sq(nz2+1)*SFACTO(nz1,nz2+1)
    If(nz2 > 0)   facz3_right = facz3_right + sq(nz2)  *SFACTO(nz1,nz2-1)

    facz1_right=facz1_right*sqi(2)*fz; facz3_right=facz3_right*sqi(2)*bz

    ! ---------------------------------------
    ! integrals in z direction, left fragment
    ! ---------------------------------------
    delnz1=0.0_pr; delnz2=0.0_pr; delnz3=0.0_pr
    If(nz1 == nz2+1) delnz1=1.0_pr
    If(nz1 == nz2-1) delnz2=1.0_pr
    If(nz1 == nz2) delnz3=1.0_pr

    facz1_left = 0.0_pr
    If(nz2 < nzm) facz1_left = facz1_left - sq(nz2+1)*(delnz1-SFACTO(nz1,nz2+1))
    If(nz2 > 0)   facz1_left = facz1_left + sq(nz2)  *(delnz2-SFACTO(nz1,nz2-1))
    facz2_left = -qh_zneck(nz1)*qh_zneck(nz2)*(1.0_pr-heavyside_at_zero)
    facz3_left = -sq(2)*zc_left*fz*(delnz3-SFACTO(nz1,nz2))
    If(nz2 < nzm) facz3_left = facz3_left + sq(nz2+1)*(delnz1-SFACTO(nz1,nz2+1))
    If(nz2 > 0  ) facz3_left = facz3_left + sq(nz2)  *(delnz2-SFACTO(nz1,nz2-1))

    facz1_left=facz1_left*sqi(2)*fz; facz3_left=facz3_left*sqi(2)*bz

    ! -------------------------------------
    ! matrix element for the right fragment
    ! -------------------------------------
    fac1=-0.5_pr*(delnl1-delnl2)*facr1*(facz1_right+facz2_right)
    fac2=-0.5_pr*(delnl1+delnl2)*nl2*facr2*facz3_right
    fac3=+0.5_pr*(delnl1-delnl2)*facr3*facz3_right
    rlx_right=fac1+fac2+fac3
    ! ------------------------------------
    ! matrix element for the left fragment
    ! ------------------------------------
    fac1=-0.5_pr*(delnl1-delnl2)*facr1*(facz1_left+facz2_left)
    fac2=-0.5_pr*(delnl1+delnl2)*nl2*facr2*facz3_left
    fac3=+0.5_pr*(delnl1-delnl2)*facr3*facz3_left
    rlx_left=fac1+fac2+fac3

  End Subroutine calc_lx
  !===========================================================================
  !> This subroutine calculates the mean square angular momenta of fission fragments.
  !===========================================================================
  Subroutine average_J_in_fragments()
    Use HFBTHO, Only : nt,nb,id,ia,nr,nl,nz,ns,ka,kd,RVqpN,RUqpN,KpwiN,RVqpP,RUqpP,KpwiP,SFACTO,qh_zneck, &
                       jxsq_rho_right,jxsq_rho_left,jxsq_kappa_left,jxsq_kappa_right

    Integer(ipr) :: n1,n2,ib1,nd1,im1,nz1,nr1,nl1,ns1,ib2,nd2,im2,nz2,nr2,nl2,ns2,iosc1,iosc2,it,ib,nd,n_qp,im,k
    Real(pr) :: rlx_left,rlx_right,rsx_left,rsx_right,jxsq_rhoN_left,jxsq_rhoN_right,jxsq_rhoP_left,jxsq_rhoP_right,&
                jxsq_kappaN_right,jxsq_kappaN_left,jxsq_kappaP_right,jxsq_kappaP_left,delnz
    Real(pr), Allocatable:: rhoN(:,:),rhoNm(:,:),rhoP(:,:),rhoPm(:,:),kappaN(:,:),kappaP(:,:),aux1(:,:),aux2(:,:),aux3(:,:)
    Real(pr), Allocatable:: jx_left(:,:),jx_right(:,:)

    Allocate(jx_left(nt,nt),jx_right(nt,nt))
    jx_left=zero; jx_right=zero

    Allocate(rhoN(nt,nt),rhoNm(nt,nt),rhoP(nt,nt),rhoPm(nt,nt),kappaN(nt,nt),kappaP(nt,nt),aux1(nt,nt),aux2(nt,nt),aux3(nt,nt))
    rhoN=zero; rhoNm=zero; rhoP=zero; rhoPm=zero; kappaN=zero; kappaP=zero

    ! ---------------------------------------------------------
    ! Calculate matrix elements of the jx operator in fragments
    ! ---------------------------------------------------------

    Call calc_qh(z_neck)

    Do ib1=1,nb
       nd1=id(ib1); im1=ia(ib1)
       Do n1=1,nd1
          iosc1=im1+n1; nr1=nr(iosc1); nl1=nl(iosc1); nz1=nz(iosc1); ns1=ns(iosc1)

          Do ib2=1,nb
             nd2=id(ib2); im2=ia(ib2)
             Do n2=1,nd2
                iosc2=im2+n2; nr2=nr(iosc2); nl2=nl(iosc2); nz2=nz(iosc2); ns2=ns(iosc2)

                ! Orbital angular momentum part
                rlx_right=zero; rlx_left=zero
                If(ns1 == ns2) Call calc_lx(nz1,nr1,nl1,nz2,nr2,nl2,cenrig,rlx_right,cenlef,rlx_left)

                ! Spin part
                rsx_right=zero; rsx_left=zero
                If(nr1 == nr2 .And. nl1 == nl2 .And. ns1 == -ns2) Then
                   rsx_right=0.5_pr*SFACTO(nz1,nz2)
                   delnz=0.0_pr
                   If(nz1 == nz2) delnz=1.0_pr
                   rsx_left=0.5_pr*(delnz-SFACTO(nz1,nz2))
                End If

                ! Matrix elements
                jx_right(iosc1,iosc2)=rlx_right+rsx_right
                jx_left(iosc1,iosc2)=rlx_left+rsx_left

             End Do
          End Do
          !
       End Do
    End Do
    !
    ! --------------------------------------------
    ! Calculate total density matrix of the system
    ! --------------------------------------------
    Do it=1,2
       Do ib=1,nb
          nd=id(ib); n_qp=kd(ib,it); im=ia(ib)
          Do n2=1,nd
             iosc2=im+n2
             Do n1=1,nd
                iosc1=im+n1
                ! Sum over qps in the block
                Do k=1,n_qp
                   If(it == 1) Then ! neutrons
                      rhoN(iosc1,iosc2)=rhoN(iosc1,iosc2)+RVqpN(KpwiN(ka(ib,it)+k)+n1)*RVqpN(KpwiN(ka(ib,it)+k)+n2)
                      kappaN(iosc1,iosc2)=kappaN(iosc1,iosc2)+RVqpN(KpwiN(ka(ib,it)+k)+n1)*RUqpN(KpwiN(ka(ib,it)+k)+n2)
                   Else If(it == 2) Then ! protons
                      rhoP(iosc1,iosc2)=rhoP(iosc1,iosc2)+RVqpP(KpwiP(ka(ib,it)+k)+n1)*RVqpP(KpwiP(ka(ib,it)+k)+n2)
                      kappaP(iosc1,iosc2)=kappaP(iosc1,iosc2)+RVqpP(KpwiP(ka(ib,it)+k)+n1)*RUqpP(KpwiP(ka(ib,it)+k)+n2)
                   End If
                End Do ! k
                !
             End Do ! n1
          End Do ! n2
       End Do ! ib
    End Do ! it

    rhoNm=-rhoN
    rhoPm=-rhoP
    Do iosc1=1,nt
       rhoNm(iosc1,iosc1)=rhoNm(iosc1,iosc1)+1.0_pr
       rhoPm(iosc1,iosc1)=rhoPm(iosc1,iosc1)+1.0_pr
    End Do

    ! ---------------------------------------------------
    ! Calculate angular momentum Jx of the right fragment
    ! ---------------------------------------------------
    ! rho contribution
    ! neutrons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_right,nt,rhoN,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,one,jx_right,nt,rhoNm,nt,zero,aux1,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux1,nt,aux2,nt,zero,aux3,nt)
    jxsq_rhoN_right=zero
    Do iosc1=1,nt
       jxsq_rhoN_right=jxsq_rhoN_right+aux3(iosc1,iosc1)
    End Do
    ! protons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_right,nt,rhoP,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,one,jx_right,nt,rhoPm,nt,zero,aux1,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux1,nt,aux2,nt,zero,aux3,nt)
    jxsq_rhoP_right=zero
    Do iosc1=1,nt
       jxsq_rhoP_right=jxsq_rhoP_right+aux3(iosc1,iosc1)
    End Do
    ! kappa contribution
    ! neutrons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_right,nt,kappaN,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux2,nt,aux2,nt,zero,aux3,nt)
    jxsq_kappaN_right=zero
    Do iosc1=1,nt
       jxsq_kappaN_right=jxsq_kappaN_right+aux3(iosc1,iosc1)
    End Do
    ! protons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_right,nt,kappaP,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux2,nt,aux2,nt,zero,aux3,nt)
    jxsq_kappaP_right=zero
    Do iosc1=1,nt
       jxsq_kappaP_right=jxsq_kappaP_right+aux3(iosc1,iosc1)
    End Do

    jxsq_rho_right=jxsq_rhoN_right+jxsq_rhoP_right
    jxsq_kappa_right=-jxsq_kappaN_right-jxsq_kappaP_right

    ! --------------------------------------------------
    ! Calculate angular momentum Jx of the left fragment
    ! --------------------------------------------------
    ! rho contribution
    ! neutrons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_left,nt,rhoN,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,one,jx_left,nt,rhoNm,nt,zero,aux1,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux1,nt,aux2,nt,zero,aux3,nt)
    jxsq_rhoN_left=zero
    Do iosc1=1,nt
       jxsq_rhoN_left=jxsq_rhoN_left+aux3(iosc1,iosc1)
    End Do
    ! protons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_left,nt,rhoP,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,one,jx_left,nt,rhoPm,nt,zero,aux1,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux1,nt,aux2,nt,zero,aux3,nt)
    jxsq_rhoP_left=zero
    Do iosc1=1,nt
       jxsq_rhoP_left=jxsq_rhoP_left+aux3(iosc1,iosc1)
    End Do
    ! kappa contribution
    ! neutrons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_left,nt,kappaN,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux2,nt,aux2,nt,zero,aux3,nt)
    jxsq_kappaN_left=zero
    Do iosc1=1,nt
       jxsq_kappaN_left=jxsq_kappaN_left+aux3(iosc1,iosc1)
    End Do
    ! protons
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_left,nt,kappaP,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux2,nt,aux2,nt,zero,aux3,nt)
    jxsq_kappaP_left=zero
    Do iosc1=1,nt
       jxsq_kappaP_left=jxsq_kappaP_left+aux3(iosc1,iosc1)
    End Do

    jxsq_rho_left=jxsq_rhoN_left+jxsq_rhoP_left
    jxsq_kappa_left=-jxsq_kappaN_left-jxsq_kappaP_left

    Deallocate(jx_left,jx_right)
    Deallocate(rhoN,rhoNm,rhoP,rhoPm,kappaN,kappaP,aux1,aux2,aux3)
    If(Allocated(qh_zneck)) Deallocate(qh_zneck)

  End Subroutine average_J_in_fragments
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_fission_fragments
