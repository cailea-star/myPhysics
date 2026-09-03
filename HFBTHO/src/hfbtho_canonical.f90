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
!                        CANONICAL BASIS                               !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!>  This module defines routines and functions to work with the canonical basis
!>
!>  @author
!>  Nicolas Schunck, Mario Stoitsov
!----------------------------------------------------------------------
!  Subroutines: - Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
!               - canonical_on_grid(derivatives)
!               - sample_occupations()
!----------------------------------------------------------------------
Module HFBTHO_canonical

  Use HFBTHO_utilities

  Implicit None

  Public Canonical,canonical_on_grid,sample_occupations

  Private

  Logical, Public, Save :: statistical_PNP = .False.

  ! Public variables coming from HFBTHO
  Integer(ipr), Allocatable, Public :: numax(:,:)
  Real(pr), Dimension(2), Public :: Dispersion,v2min,v2minv
  Real(pr), Allocatable, Target, Public :: ddc(:,:,:)  !< canonical wave functions \f$ |\varphi_{i}\rangle \f$ (in the HO basis)
  Real(pr), Allocatable, Target, Public :: vk(:,:)     !< canonical occupations \f$ v_{i}^{2} \f$
  Real(pr), Allocatable, Target, Public :: ek(:,:)     !< canonical single-particle energies \f$ e_{i} \f$
  Real(pr), Allocatable, Target, Public :: dk(:,:)     !< canonical pairing gaps \f$ \Delta_{i} \f$
  Real(pr), Allocatable, Target, Public :: hfbcan(:,:) !< canonical wave functions \f$ |\varphi_{i}\rangle \f$ (in a given Omega block)
  Real(pr), Allocatable, Target, Public :: evvkcan(:)  !< canonical single-particle energies \f$ e_{i} \f$ (in a given Omega block)
  Real(pr), Allocatable, Target, Public :: vkmax(:,:)
  ! Public variables read from file in this module
  Integer(ipr), Public :: ntot
  Integer(ipr), Dimension(2), Public :: n_states !< number of canonical states for protons (it=1) and for neutrons (it=2)

  ! Private variables
  Integer(ipr), Private :: lcc

Contains
  !=======================================================================
  !> This routines determines the canonical basis by diagonalization of the density matrix.
  !=======================================================================
  Subroutine Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
    Use HFBTHO, Only : RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP,KpwiN,KpwiP,KqpN,KqpP,ept,switch_on_temperature, &
                       blocanon,nhhdim,brin,pwiblo,eqpmin,blomax,bloall,bloblo,blo123,blok1k2,bloqpdif, &
                       alwork,lwork,ialwork,ilwork,ierror_flag,ierror_info,ndx,ala,iparenti

    Integer(ipr), Intent(In) :: it,icanon,k2,k1,nd,i0,ib,ibiblo,ibroib
    Integer(ipr), Intent(Inout) :: lc,m

    Integer(ipr) :: i,k,kk,n1,n2,nd1,ier
    Integer(ipr) :: nhhph,nhhpp,ibro,il,iu,NUMFOU
    Real(pr) :: s1,s2,vx,h1,d1,h2,d2,ddn1,ddn2,vl,vu,ABSTOL
    Real(pr), Allocatable :: hh(:,:),de(:,:)
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
#if(USE_ESSL==0)
    Integer(ipr), Allocatable :: ISUPPZ(:)
    Real(pr), Allocatable :: eigenv(:),eigenf(:,:)
#else
    Integer(ipr), Allocatable :: ifail(:)
    Real(pr), Allocatable :: ewavef(:,:)
#endif
    Real(pr), External :: DLAMCH

    If(it == 1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
    End If

#if(USE_ESSL==0)
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    !ialwork=26*nd; ilwork=10*nd
    ialwork=1+6*nd+2*nd**2; ilwork=3+5*nd ! DSYEVD
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 1
#else
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    ialwork=8*nd; ilwork=5*nd;
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 0
#endif

    ABSTOL=2.0_pr*DLAMCH('S')
    m=0

    !------------------------------------------------------
    ! Vanishing pairing => just taking the HF states
    !------------------------------------------------------
    If(Abs(ept(it)) < 0.0001_pr .And. (.Not.switch_on_temperature)) Then
       Do k=1,nd
          kk=k1+k-1; lc=lc+1                              ! total number of the canonical states
          ddc(1:nd,lc,it)=zero; vk(lc,it)=zero            ! zeros: nd could be larger then k2-k1+1
          If(kk > k2) Cycle
          vx=zero
          Do i=1,nd
             h1=VqpPo(KpwiPo(kk)+i)**2; vx=vx+h1
          End Do
          If(vx <= zero) vx=zero                         ! roundoff errors
          If(vx >= one ) vx=one
          Do i=1,nd
             If(vx >= half) Then
                ddc(i,lc,it)=VqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             Else
                ddc(i,lc,it)=UqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             End If
          End Do
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx)  ! internal P/N Dispersion
          If(Abs(vx-half) <= v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx        ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx                                    ! (ph) s.p. occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon /= 0) Then
             ek(lc,it)=EqpPo(KqpPo(kk))*(one-two*vx)+ala(it)      ! (ph) s.p. energies
             dk(lc,it)=zero                                       ! (ph) s.p. deltas
          End If
       End Do !k
    Else
       !------------------------------------------------------
       ! Pairing => calculate canonical basis
       !------------------------------------------------------
#if(USE_ESSL==0)
       VL=0.0_pr; VU=0.0_pr; IL=1; IU=nd; NUMFOU=0
       Allocate(ISUPPZ(2*nd))
       Allocate(eigenv(nd)); eigenv(1:nd)=0.0_pr
       Allocate(eigenf(nd,nd)); eigenf(1:nd,1:nd)=0.0_pr
       ier=0; Call DSYEVD('V','L',nd,hfbcan,ndx,eigenv,ALWORK,ialwork,LWORK,ilwork,ier)
       !ier=0; Call DSYEVR('V','A','L',nd,hfbcan,ndx,VL,VU,IL,IU,ABSTOL,NUMFOU, &
       !                   eigenv,eigenf,nd,ISUPPZ,ALWORK,ialwork,LWORK,ilwork,ier)
       !hfbcan(1:nd,1:nd) = eigenf(1:nd,1:nd)
       evvkcan(1:nd) = eigenv(1:nd)
       Deallocate(eigenv,eigenf,ISUPPZ)
#else
       vl=0.0_pr; vu=0.0_pr; il=1; iu=1; m=0; abstol=0.0_pr
       Allocate(ifail(nd),ewavef(ndx,ndx))
       ier=0; Call DSYEVX('V','A','L',nd,hfbcan,ndx,vl,vu,il,iu,abstol,m, &
                          evvkcan,ewavef,ndx,ALWORK,ialwork,LWORK,ifail,ier)
       hfbcan(1:ndx,1:ndx) = ewavef(1:ndx,1:ndx)
       Deallocate(ifail,ewavef)
#endif
       !------------------------------------------------------
       ! Eigenvalues and wavefunctions
       !------------------------------------------------------
       Do k=1,nd
          lc=lc+1 ! total number of canonical states
          Do i=1,nd
             ddc(i,lc,it)=hfbcan(i,k) ! (ph) canon orbitals in conf.space
          End Do
          vx=evvkcan(k)*half
          If(vx <= zero) vx=zero ! roundoff errors
          If(vx >= one ) vx=one
          ! blocking
          If(ibiblo == ib .And. vx > 0.49_pr .And. vx <= 0.51_pr) blocanon(it)=lc
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx) ! internal P/N Dispersion
          If(Abs(vx-half) <= v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx ! (ph) canon occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon /= 0) Then
             ! canonical energies and deltas (no physical meaning in PNP)
             nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
             Allocate(hh(nd,nd),de(nd,nd))
             ibro=ibroib
             Do n1=1,nd
                Do n2=1,n1
                   ibro=ibro+1
                   vx=brin(nhhph+ibro); hh(n2,n1)=vx; hh(n1,n2)=vx
                   vx=brin(nhhpp+ibro); de(n2,n1)=vx; de(n1,n2)=vx
                End Do
             End Do
             h1=zero; d1=zero
             Do n2=1,nd
                h2=zero; d2=zero
                Do n1=1,nd
                   ddn1=hfbcan(n1,k)
                   h2=h2+ddn1*hh(n1,n2)
                   d2=d2+ddn1*de(n1,n2)
                End Do
                ddn2=hfbcan(n2,k)
                h1=h1+h2*ddn2
                d1=d1+d2*ddn2
             End Do
             ek(lc,it)=h1 ! (ph) canon s.p. energies
             dk(lc,it)=d1 ! (ph) canon s.p. deltas
             Deallocate(hh,de)
          End If
       End Do !k
    End If

    !------------------------------------------------------
    ! RESU only
    !------------------------------------------------------
    If(icanon /= 0) Then
       !------------------------------------------------------
       ! Find maximal HO components of all qp states
       !------------------------------------------------------
       Do k=k1,k2
          s1=zero
          Do n1=1,nd
             nd1=nd+n1
             s2=Max(s1,Abs(VqpPo(KpwiPo(k)+n1)),Abs(UqpPo(KpwiPo(k)+n1)))
             If(s2 > s1) Then
                s1=s2
                vkmax(k,it)=s1    ! maximal overlap
                numax(k,it)=n1+i0 ! its number in k[k1,k2] numbering
             End If
          End Do
       End Do
       !------------------------------------------------------
       ! Searching for possible blocking candidates
       !------------------------------------------------------
       If(iparenti(1) == 0 .And. iparenti(2) == 0) Then
          n1=0
          Do k=k1,k2
             n1=n1+1
              ! Search within |(1-2*N)*Eqpe| lover than 'pwiblo'
              ! The levels number n1 is 1,2,3,... for the given block ([123] numbering)
             If(Abs(EqpPo(KqpPo(k))-eqpmin(it)) <= pwiblo) Then
                blomax(it)=blomax(it)+1                         ! blocked state #, maximel # of block candidates
                If(blomax(it) > bloall) Then
                   ierror_flag=ierror_flag+1
                   ierror_info(ierror_flag)='Too many blocking candidates! Increase bloall and run again'
                   Return
                End If
                bloblo(blomax(it),it)=ib ! Omega block where to block
                blo123(blomax(it),it)=n1 ! index of the blocking candidate within the block
                blok1k2(blomax(it),it)=k ! index of the blocking candidate amongst all quasiparticles
                bloqpdif(blomax(it),it)=Abs(EqpPo(KqpPo(k))-eqpmin(it)) ! difference in energy between current qp and one with lowest energy
             End If
          End Do
       End If
    End If

  End Subroutine Canonical
  !=======================================================================
  !> This routine calculates the value of canonical wave functions on the Gauss quadrature grid
  !> and stores the result in a file.
  !=======================================================================
  Subroutine canonical_on_grid(derivatives)
    Use HFBTHO_Gauss, Only : ngh,ngl,xh,xl,qh,ql,wl,wh,qh1,ql1
#if(USE_MPI==0)
    Use HFBTHO, Only : bz,bp,pwi,nb,id,ia,ala,kd,nz,nr,nl,ns
#else
    Use HFBTHO, Only : do_print,bz,bp,pwi,nb,id,ia,ala,kd,nz,nr,nl,ns
#endif

    Logical, Intent(In), Optional :: derivatives

    Logical :: with_threshold,write_derivatives
    Integer(ipr) :: it,k0,ib,nd,im,k,lcanwf,i,ih,il,iosc,nza,nra,nla,nsa
    Real(pr) :: ekk,xqhla_up,xqhla_dw,xfi1r_up,xfi1r_dw,xfi1z_up,xfi1z_dw,xfi2d_up,xfi2d_dw
    Real(pr) :: bpi,bpi2,bzi,bzi2,cnraa,cnzaa,sml2,v2,v4
    Character(50) :: canwfname

    ! include all the canonical states if set to .False.
    with_threshold=.True.
    If(Present(derivatives)) Then
       write_derivatives = derivatives
    Else
       write_derivatives = .False.
    End If

    lcanwf=20
     Write(canwfname, '("canonwf",a4)') '.dat'
#if(USE_MPI==0)
     Open(lcanwf,file=canwfname,status='unknown')
#else
     If(do_print == 1) Open(lcanwf,file=canwfname,status='unknown')
#endif
    !
    Write(lcanwf,'(A)') '#  -----------------------------------------------'
    Write(lcanwf,'(A)') '# | Canonical wave functions on a coordinate grid |'
    Write(lcanwf,'(A)') '#  -----------------------------------------------'
    Write(lcanwf,'(A)') '# tau......: id of the isospin'
    Write(lcanwf,'(A)') '# omega....: value of omega=lambda+ms'
    Write(lcanwf,'(A)') '# k........: id of the states in the canonical basis'
    Write(lcanwf,'(A)') '# Esp_c....: canonical s.p. energies'
    Write(lcanwf,'(A)') '# Eqp_c....: canonical q.p. energies'
    Write(lcanwf,'(A)') '# delta....: pairing gap'
    Write(lcanwf,'(A)') '# P_occ....: occupation factor v^2'
    Write(lcanwf,'(A)') '# z........: z-coordinate in [fm]'
    Write(lcanwf,'(A)') '# rperp....: distance to z-axis in [fm]'
    Write(lcanwf,'(A)') '# sigma....: spin coordinate'
    Write(lcanwf,'(A)') '# canon_wf.: dataset of the canonical wavefunctions'
    Write(lcanwf,'(A)') '#            layout: one row per state, same order as tau, omega, k, etc'
    Write(lcanwf,'(A)') '#                    one column per node in the grid, in the same order'

    ! Print global configuration
    Write(lcanwf,'(/,A)') '[global]'
    Write(lcanwf,'(A,ES12.5)') 'Eqp_thresh = ', pwi
    Write(lcanwf,'(A)') 'isospins   =      neutron       proton'
    ! Loop over isospins
    Do it=1,2
       k0 = 0
       n_states(it) = 0
       ! Loop over Omega blocks
       Do ib=1,nb
          nd = id(ib)
          ! Canonical states
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             n_states(it) = n_states(it) + 1
          End Do
       End Do
    End Do
    Write(lcanwf,'(A,I12,1X,I12)')       'n_states   = ', n_states(1), n_states(2)
    Write(lcanwf,'(A,ES12.5,1X,ES12.5)') 'lambda     = ', ala(1), ala(2)

    ! Print the spatial discretization grid
    Write(lcanwf,'(/,A)') '[grid]'
    Write(lcanwf,'(A)') 'regular = yes'
    Write(lcanwf,'(A)') 'order   = rp theta z sigma'
    Write(lcanwf,'(A)') 'type    = gauss-hermite custom gauss-laguerre-sq ones'
    Write(lcanwf,'(A)',advance='no') 'rp      ='
    Do il=1,NGL
       Write(lcanwf,'(1X,ES22.15)',advance='no') sqrt(xl(il))*bp
    End Do
    Write(lcanwf,'()')
    Write(lcanwf,'(A,ES22.15)') 'rp:b    = ', bp
    Write(lcanwf,'(A,I12)')     'rp:n    = ', NGL
    Write(lcanwf,'(A,ES22.15)') 'theta   = ', 0.0
    Write(lcanwf,'(A,ES22.15)') 'theta:w = ', 2.0*pi
    Write(lcanwf,'(A)',advance='no') 'z       ='
    Do ih=1,NGH
       Write(lcanwf,'(1X,ES22.15)',advance='no') xh(ih)*bz
    End Do
    Write(lcanwf,'()')
    Write(lcanwf,'(A,ES22.15)') 'z:b     = ', bz
    Write(lcanwf,'(A,I12)')     'z:n     = ', NGH
    Write(lcanwf,'(A)')         'sigma   = +1/2 -1/2'

    Write(lcanwf,'(/,A)') '[states]'

    ! Print isospin
    Write(lcanwf,'(A)',advance='no') 'tau   ='
    Do it=1,2
       k0 = 0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             Write(lcanwf,'(1X,I12)',advance='no') it-1
          End Do
       End Do
    End Do
    Write(lcanwf,'()')

    ! Print omega block
    Write(lcanwf,'(A)',advance='no') 'omega ='
    Do it=1,2
       k0 = 0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             Write(lcanwf,'(1X,I10,"/2")',advance='no') 2*ib-1
          End Do
       End Do
    End Do
    Write(lcanwf,'()')

    ! Print state id in omega block
    Write(lcanwf,'(A)',advance='no') 'k     ='
    Do it=1,2
       k0 = 0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             Write(lcanwf,'(1X,I12)',advance='no') k-1
          End Do
       End Do
    End Do
    Write(lcanwf,'()')

    ! Print canonical s.p. energies
    Write(lcanwf,'(A)',advance='no') 'Esp_c ='
    Do it=1,2
       k0 = 0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             Write(lcanwf,'(1X,ES12.5)',advance='no') ek(k0,it)
          End Do
       End Do
    End Do
    Write(lcanwf,'()')

    ! Print canonical q.p. energies
    Write(lcanwf,'(A)',advance='no') 'Eqp_c ='
    Do it=1,2
       k0 = 0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             ekk = Sqrt((ek(k0,it)-ala(it))**2+dk(k0,it)**2)+ala(it)
             Write(lcanwf,'(1X,ES12.5)',advance='no') ekk
          End Do
       End Do
    End Do
    Write(lcanwf,'()')

    ! Print pairing strength delta
    Write(lcanwf,'(A)',advance='no') 'delta ='
    Do it=1,2
       k0 = 0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             Write(lcanwf,'(1X,ES12.5)',advance='no') dk(k0,it)
          End Do
       End Do
    End Do
    Write(lcanwf,'()')

    ! Print occupation probabilities v^2
    Write(lcanwf,'(A)',advance='no') 'P_occ ='
    Do it=1,2
       k0 = 0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          Do k=1,nd
             k0 = k0 + 1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             Write(lcanwf,'(1X,ES12.5)',advance='no') vk(k0,it)
          End Do
       End Do
    End Do
    Write(lcanwf,'(/)')

    ! Print wave functions, derivatives
    Write(lcanwf,'(A)') '<canon_wf'
    Do it=1,2
       k0=0
       Do ib=1,nb
          If(kd(ib,it) <= 0) Cycle
          nd = id(ib)
          im = ia(ib)
          Do k=1,nd
             k0=k0+1
             If(with_threshold .And. ek(k0,it) >= pwi) Cycle
             ! Coordinate space
             Do il=1,NGL
                Do ih=1,NGH
                   xqhla_up = 0.0_pr
                   xqhla_dw = 0.0_pr
                   ! Oscillator states
                   Do i=1,nd
                      iosc=im+i; nza=nz(iosc); nra=nr(iosc); nla=nl(iosc); nsa=ns(iosc)
                      sml2 = nla*nla; cnzaa = nza+nza+1; cnraa = nra+nra+nla+1
                      !v2=half/xl(il); v4=v2*v2
                      If(nsa == +1) Then
                         xqhla_up = xqhla_up + ddc(i,k0,it)*qh(nza,ih)*ql(nra,nla,il)
                      End If
                      If(nsa == -1) Then
                         xqhla_dw = xqhla_dw + ddc(i,k0,it)*qh(nza,ih)*ql(nra,nla,il)
                      End If
                   End Do
                   xqhla_up = xqhla_up / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                   xqhla_dw = xqhla_dw / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                   Write(lcanwf,'(1X,ES19.12,1X,ES19.12)',advance='no') xqhla_up, xqhla_dw
                End Do
             End Do
             Write(lcanwf,'()')
          End Do
       End Do
    End Do
    Write(lcanwf,'(A)') '>'

    If(write_derivatives) Then
       ! Print derivatives and Laplacian
       Write(lcanwf,'(A)') '<deriv_wf'
       bpi=one/bp; bpi2=bpi*bpi; bzi=one/bz; bzi2=bzi*bzi
       Do it=1,2
          k0=0
          Do ib=1,nb
             If(kd(ib,it) <= 0) Cycle
             nd = id(ib)
             im = ia(ib)
             Do k=1,nd
                k0=k0+1
                If(with_threshold .And. ek(k0,it) >= pwi) Cycle
                ! Coordinate space
                Do il=1,ngl
                   Do ih=1,ngh
                      xfi1r_up = 0.0_pr; xfi1r_dw = 0.0_pr
                      xfi1z_up = 0.0_pr; xfi1z_dw = 0.0_pr
                      xfi2d_up = 0.0_pr; xfi2d_dw = 0.0_pr
                      ! Oscillator states
                      Do i=1,nd
                         iosc=im+i; nza=nz(iosc); nra=nr(iosc); nla=nl(iosc); nsa=ns(iosc)
                         sml2 = nla*nla; cnzaa = nza+nza+1; cnraa = nra+nra+nla+1
                         v2=half/xl(il); v4=v2*v2
                         If(nsa == +1) Then
                            xfi1r_up = xfi1r_up + ddc(i,k0,it)*qh(nza,ih)*ql1(nra,nla,il)*(two*sqrt(xl(il))*bpi)*half/xl(il)
                            xfi1z_up = xfi1z_up + ddc(i,k0,it)*qh1(nza,ih)*ql(nra,nla,il)/bz
                            xfi2d_up = xfi2d_up + ddc(i,k0,it)*qh(nza,ih)*ql(nra,nla,il)*((xh(ih)**2-cnzaa)*bzi2+four*(p14-cnraa*v2+sml2*v4)*xl(il)*bpi2 )
                         End If
                         If(nsa == -1) Then
                            xfi1r_dw = xfi1r_dw + ddc(i,k0,it)*qh(nza,ih)*ql1(nra,nla,il)*(two*sqrt(xl(il))*bpi)*half/xl(il)
                            xfi1z_dw = xfi1z_dw + ddc(i,k0,it)*qh1(nza,ih)*ql(nra,nla,il)/bz
                            xfi2d_dw = xfi2d_dw + ddc(i,k0,it)*qh(nza,ih)*ql(nra,nla,il)*((xh(ih)**2-cnzaa)*bzi2+four*(p14-cnraa*v2+sml2*v4)*xl(il)*bpi2 )
                         End If
                      End Do
                      xfi1r_up = xfi1r_up / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                      xfi1r_dw = xfi1r_dw / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                      xfi1z_up = xfi1z_up / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                      xfi1z_dw = xfi1z_dw / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                      xfi2d_up = xfi2d_up / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                      xfi2d_dw = xfi2d_dw / Sqrt(wh(ih)) / Sqrt(wl(il)) / Sqrt(2.0_pr*pi) / Sqrt(bz) / bp * Sqrt(2.0_pr)
                      Write(lcanwf,'(3(1X,ES19.12,1X,ES19.12))',advance='no') xfi1r_up,xfi1r_dw,xfi1z_up,xfi1z_dw,xfi2d_up,xfi2d_dw
                   End Do
                End Do
                Write(lcanwf,'()')
             End Do
          End Do
       End Do
       Write(lcanwf,'(A)') '>'
    End If

    Close(lcanwf)

  End Subroutine canonical_on_grid
  !=======================================================================
  !> This routine records on disk the occupations both of canonical states and of
  !> quasiparticles. This information is used to extract the number of particles in
  !> fission fragments using the technique discussed in \cite verriere2019numbera.
  !=======================================================================
  Subroutine sample_occupations()
    Use HFBTHO, Only : npr,itmin,itmax,nb,id,ka,kd,uk

    Character(Len=1) :: type_part
    Integer(ipr) :: it,ib,k,k1,k2,ierr,k0,nd,nstate
    Real(pr) :: pn,ek0,vk0

    Open(55,file='occup_canonical.dat',status='unknown',form='formatted',IOSTAT=ierr)
    Write(55,'(2i4)') npr(1),npr(2)
    Write(55,'("Aneck   0.0000000000000E+0.000")')
    Write(55,'("Sneck   0.0000000000000E+0.000")')
    Write(55,'("Qneck 0 0.0000000000000E+0.000")')
    Write(55,'("Qneck 1 0.0000000000000E+0.000")')
    Write(55,'("n_gh   001")')
    Do it=itmin,itmax
       If(it == 1) type_part='0'
       If(it == 2) type_part='1'
       k0=0; nstate=0
       Do ib=1,nb
          nd=id(ib)
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1 <= k2) Then
             Do k=1,nd                   ! print active states only
                k0=k0+1
                ek0=ek(k0,it)            ! canonical s.p. energies
                vk0=vk(k0,it)            ! canonical occupations v^2
                If(vk0 > -1.0e-4_pr) Then   ! print If signIficant v^2
                   nstate=nstate+1
                End If
             End Do
          End If
       End Do
       Write(55,'("basis_",a1,i5)') type_part,nstate
    End Do
    Do it=itmin,itmax
       If(it == 1) type_part='0'
       If(it == 2) type_part='1'
       k0=0
       Do ib=1,nb
          nd=id(ib)
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1 <= k2) Then
             Do k=1,nd                   ! print active states only
                k0=k0+1
                ek0=ek(k0,it)            ! canonical s.p. energies
                vk0=vk(k0,it)            ! canonical occupations v^2
                If(vk0 > -1.0e-4_pr) Then   ! print If signIficant v^2
                   pn=vk(k0,it)          ! qp probabilities
                   Write(55,'(a1,1x,e20.13)') type_part,pn
                End If
             End Do
          End If
       End Do !ib
    End Do !it
    Close(55)
    ! Write qp energies and occupations
    Open(55,file='occup_qp.dat',status='unknown',form='formatted',IOSTAT=ierr)
    Write(55,'(2i4)') npr(1),npr(2)
    Write(55,'("Aneck   0.0000000000000E+0.000")')
    Write(55,'("Sneck   0.0000000000000E+0.000")')
    Write(55,'("Qneck 0 0.0000000000000E+0.000")')
    Write(55,'("Qneck 1 0.0000000000000E+0.000")')
    Write(55,'("n_gh   001")')
    Do it=itmin,itmax
       If(it == 1) type_part='0'
       If(it == 2) type_part='1'
       nstate=0
       Do ib=1,nb
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1 <= k2) Then
             Do k=k1,k2      ! print active states only
                pn=uk(k,it)  ! qp probabilities
                nstate=nstate+1
             End Do
          End If
       End Do
       Write(55,'("basis_",a1,i5)') type_part,nstate
    End Do
    Do it=itmin,itmax
       If(it == 1) type_part='0'
       If(it == 2) type_part='1'
       Do ib=1,nb
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1 <= k2) Then
             Do k=k1,k2     ! print active states only
                pn=uk(k,it) ! qp probabilities
                Write(55,'(a1,1x,e20.13)') type_part,pn
             End Do
          End If
       End Do !ib
    End Do !it
    Close(55)

  End Subroutine sample_occupations
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_canonical