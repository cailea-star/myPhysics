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
!                      FINITE-RANGE GOGNY PACKAGE                      !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!>  This module provides all the routines needed to compute the matrix
!>  elements and expectation value of the finite-range Gogny force.
!>  The separability of the matrix elements in axial and radial
!>  components is widely used to signinificantly reduce computation time.
!>  The matrix elements are computed using the Gogny expansions outlined
!>  in \cite younes2009gaussian.
!>
!>  @author
!>  Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - allocate_fr()
!               - gamdel_gogny()
!               - gogny_matrix_elements()
!               - CalculateVzGogny()
!               - CalculateVrGogny()
!               - calculateTz()
!               - calculateTr()
!               - calculateCpolar2cartesian()
!               - calculateME1D()
!               - radial_matrix_elements(ni,li,nj,lj,nk,lk,nl,ll)
!               - calculate_Zblock()
!               - calculate_Nblock()
!               - trace_product_2(A,B,tr1,tr2)
!               - bdiag_trace(A,tr1,tr2)
!               - bdiag_print(A)
!               - test_HOWF_gauss
!               - LaguerreL(n,alpha,x,Ln,Lnp,Lnm1)
!               - GaussLaguerreWX(alfa,w,x)
!  Functions: - TrCoefficient(n1,k1,n2,k2,n,k)
!             - TrSumTerms(n1,k1,n2,k2,n,k)
!             - MatrixElement_z(ni,nj,nk,nl,mu,b)
!             - Ibarz(m,n,Gz)
!             - MatrixElement_r(ni,li,nj,lj,nk,lk,nl,ll,mu,b)
!             - Ibarr(n1,k1,n2,k2,Gp)
!             - factrl(n)
!             - HyperGeom2F1(a,b,c,x)
!             - binomialco(m,n)
!             - upperfactrl(x,i)
!             - l_block
!             - zindex(nzi,nzj,nzk,nzl)
!             - rindex(nri,nrj,nrk,nrl,li,lj,lk,ll,n)
!             - trace_product(A,B)
!             - MatrixElement_zZR(ni,nj,nk,nl)
!             - MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll)
!             - N_radial(nr,l)
!             - N_axial(nz)
!             - HermiteH(n,x)
!----------------------------------------------------------------------
Module HFBTHO_Gogny

  Use HFBTHO_utilities

  Implicit None

  Public allocate_fr,gamdel_gogny,gogny_matrix_elements,zindex,rindex,trace_product,trace_product_2

  Private CalculateVzGogny,CalculateVrGogny,calculateTz,calculateTr,TrCoefficient,TrSumTerms,calculateCpolar2cartesian, &
          calculateME1D,MatrixElement_z,Ibarz,radial_matrix_elements,MatrixElement_r,Ibarr,factrl,HyperGeom2F1,binomialco, &
          upperfactrl,calculate_Zblock,calculate_Nblock,l_block,bdiag_trace,bdiag_print,test_HOWF_gauss,MatrixElement_zZR, &
          MatrixElement_rZR,N_radial,N_axial,LaguerreL,GaussLaguerreWX,HermiteH

  ! Public variables
  Integer(ipr), Public, Save :: NumVz !< Number of matrix elements in the z-direction
  Integer(ipr), Public, Save :: NumVr !< Number of matrix elements in the perpendicular direction
  Integer(ipr), Public, Allocatable, Save :: ib_zrls(:,:,:,:) !<
  Integer(ipr), Public, Allocatable, Save :: i_zrls(:,:,:,:)
  Real(pr), Public, Allocatable, Save :: VzGogny(:,:) !< Matrix elements of the Gogny force in the z-direction
                                                      !< \f$ \langle n_z m_z | \hat{V} | n'_z m'_z \rangle \f$
  Real(pr), Public, Allocatable, Save :: VrGogny(:,:) !< Matrix elements of the Gogny force in the perpendicular direction
                                                      !< \f$ \langle n_r \Lambda_n m_r \Lambda_m | \hat{V} | n'_r \Lambda'_n m'_r \Lambda'_m \rangle \f$
  Real(pr), Public, Allocatable, Save :: gamma_g_dir(:,:) !< Matrix of the direct part of the mean field \f$ \Gamma^{(\mathrm{dir.})}_{nm} \f$ for the Gogny force
  Real(pr), Public, Allocatable, Save :: gamma_g_exc(:,:) !< Matrix of the exchange part of the mean field \f$ \Gamma^{(\mathrm{exc.})}_{nm} \f$ for the Gogny force
  Real(pr), Public, Allocatable, Save :: delta_g_dir(:,:) !< Matrix of the pairing field \f$ \Delta_{nm} \f$ for the Gogny force
  Real(pr), Public, Allocatable, Save :: coulf_g_dir(:,:) !< Matrix of the direct part of the mean field \f$ \Gamma^{(\mathrm{dir.})}_{nm} \f$ for the Coulomb force
  Real(pr), Public, Allocatable, Save :: coulf_g_exc(:,:) !< Matrix of the exchange part of the mean field \f$ \Gamma^{(\mathrm{exc.})}_{nm} \f$ for the Coulomb force
  Real(pr), Public, Allocatable, Save :: coulf_d_dir(:,:) !< Matrix of the pairing field \f$ \Delta_{nm} \f$ for the Coulomb force
#if(GOGNY_SYMMETRIES==0)
  Real(pr), Public, allocatable, Save  :: Vz_Gogny(:,:,:,:,:) !< Matrix elements of the Gogny force in the z-direction
                                                              !< \f$ \langle n_z m_z | \hat{V} | n'_z m'_z \rangle \f$
  Real(pr), Public, allocatable, Save  :: Vr_Gogny(:,:,:,:,:,:,:,:,:) !< Matrix elements of the Gogny force in the perpendicular direction
                                                                      !< \f$ \langle n_r \Lambda_n m_r \Lambda_m | \hat{V} | n'_r \Lambda'_n m'_r \Lambda'_m \rangle \f$
#endif
  logical, Public, Save :: matrix_elements_calculated !< Logical flag that informs if matrix elements have been calculated

  ! Private variables
  Integer(ipr), Private, Allocatable, Save :: Zblock(:,:,:,:) !< Indices corresponding to non-zero and different axial matrix elements
  Integer(ipr), Private, Allocatable, Save :: NBlock(:,:,:,:,:,:) !< Indices corresponding to non-zero and different radial matrix elements
  Real(pr), Private, Allocatable, Save :: T_z(:,:,:) !< \f$ T_{n_1 n_2}^n \f$= Tz(n1,n2,n) coefficients of the axial Gogny expansion
  Real(pr), Private, Allocatable, Save :: T_r(:,:,:,:,:) !< \f$ T_{n_1,k_1,n_2,k_2}^{n,k1+k2}\f$ = Tr(n1,k1,n2,k2,n) coefficients of the radial Gogny expansion
  Real(pr), Private, Allocatable, Save :: Cp2c(:,:,:) !< \f$ C^{n_r k}_{n_y}\f$= Cp2c(nr,k,ny) coefficients of the radial to cartesian
                                                           !< coordinates transformation of the Harmonic Oscillator wavefunction
  Real(pr), Private, Allocatable, Save :: ME1D(:,:,:,:,:) !< Two body potential matrix elements in one dimension
                                                               !< \f$ \langle n_1 n_2|\hat{V}|n_3 n_4\rangle \f$ = ME1D(n1,n2,n3,n4)
  Real(pr), Private, Allocatable, Save :: Vr_ig(:)

!$omp threadprivate(Vr_ig)
Contains
  !======================================================================
  !> This routine allocates arrays related to finite-range potentials.
  !======================================================================
  Subroutine allocate_fr()
    Use HFBTHO, Only : ndx,nbx

    If(Allocated(gamma_g_dir)) Deallocate(gamma_g_dir,gamma_g_exc,delta_g_dir,coulf_g_dir,coulf_g_exc,coulf_d_dir)
    Allocate(gamma_g_dir(ndx**2,2*nbx),gamma_g_exc(ndx**2,2*nbx),delta_g_dir(ndx**2,2*nbx),coulf_g_dir(ndx**2,2*nbx), &
             coulf_g_exc(ndx**2,2*nbx),coulf_d_dir(ndx**2,2*nbx))
    gamma_g_dir = zero
    gamma_g_exc = zero
    delta_g_dir = zero
    coulf_g_dir = zero
    coulf_g_exc = zero
    coulf_d_dir = zero

  End Subroutine allocate_fr
  !=======================================================================
  !> This routine calculates the HFB fields \f$ \Gamma \f$ and \f$ \Delta \f$ for a finite-range
  !> Gogny-type potential using the axial and radial components of the two-body potential matrix
  !> elements calculated with \ref gogny_matrix_elements().
  !>
  !> If the preprocessor variable GOGNY_SYMMETRIES is set to 1, the calculation takes advantage
  !> of the separability of the matrix elements and first contracts the axial matrix elements with
  !> the density and pairing matrices and creates the objects Zrho and Zkappa. Subsequently these
  !> objects are contracted with the radial matrix elements to obtain the \f$ \Gamma \f$ and
  !> \f$ \Delta \f$ fields. This significantly reduces the number of operations and speeds up the
  !> calculation. In this mode several other symmetries of the axial matrix element are used to
  !> reduce the computation time.
  !>
  !> For any other value of GOGNY_SYMMETRIES the loops to calculate the HFB fields go over every
  !> combination of quantum numbers. This mode of course requires a larger computation time and is
  !> only used when debugging the calculation of the fields.
  !=======================================================================
  Subroutine gamdel_gogny()
    Use UNEDF, Only: e2charg,force_is_dme,n_g,n_g_all,W_g,B_g,H_g,M_g,V_g_coul
    Use HFBTHO, Only : nrx,nlx,nbx,nzx,nttx,nzzx,nox,nzz,nrr,nll,nss,noo,ia,id,rk,ak,nz,nr,ns,nl,icou,allibro, &
                       nhhdim,nhhdim2,nhhdim3,brout,ntx,nrlx

    Integer(ipr) :: jb,i,ig
    Integer(ipr) :: i1,i2,i3
    Integer(ipr) :: nlb,nld,nzac
    Integer(ipr) :: iz_abcd,iz_abdc,iz_acbd
    Integer(ipr) :: ir_abcd,ir_abdc,ir_acbd
    Integer(ipr) :: jr_abcd,jr_abdc,jr_acbd
    Real(pr) :: rho_n,rho_p,kap_n,kap_p
    Real(pr) :: gamma_n_dir,gamma_n_exc,gamma_p_dir,gamma_p_exc
    Real(pr) :: delta_n_dir,delta_p_dir
    Real(pr) :: gamma_c_dir,gamma_c_exc,delta_c_dir
    Real(pr) :: V_abcd,V_abdc,V_acbd
    Real(pr) :: U_abcd,U_abdc,U_acbd
    Real(pr) :: Vdir,Udir,Vexc,Uexc,Vdel,Udel
    Real(pr) :: W,B,H,M
    Integer(ipr) :: ob,la,lb,lc,ld,iba,ibb,ibc,ibd
    Real(pr) :: Zrho_n_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    Real(pr) :: Zrho_p_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    Real(pr) :: Zrho_n_exc(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    Real(pr) :: Zrho_p_exc(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    Real(pr) :: Zkap_n_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    Real(pr) :: Zkap_p_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    Real(pr) :: zrho_nd,zrho_pd,zrho_ne,zrho_pe,zkap_n,zkap_p
    Integer(ipr) :: ra,za,sa,na,ida,ima,ibax
    Integer(ipr) :: rb,zb,sb,nb,idb,imb,ibbx
    Integer(ipr) :: rc,zc,sc,nc,imc
    Integer(ipr) :: rd,zd,sd,jd,nd,imd
    Integer(ipr) :: nac,ndb,sac,sdb,n_rlx
    Integer(ipr) :: ita,itb,itc,itd,npb,npd,izb,izd

    gamma_g_dir = zero
    gamma_g_exc = zero
    delta_g_dir = zero
    coulf_g_dir = zero
    coulf_g_exc = zero
    coulf_d_dir = zero

#if(GOGNY_SYMMETRIES==1)
    n_rlx = Max(2*nrx,nlx)
    If(.not.force_is_dme) Then
       !$OMP Parallel Default(None)&
       !$OMP& SHARED(nzx,nttx,nrr,nll,nss,noo,nzzx,nzz,ib_zrls,id,ia,nbx,&
       !$OMP&        i_zrls,rk,ak,VzGogny,ntx,nz,nr,nl,ns,nox,nrx,nrlx,n_rlx,&
       !$OMP&        W_g,B_g,H_g,M_g,VrGogny,allibro,nhhdim,nhhdim2,nhhdim3,&
       !$OMP&        gamma_g_dir,gamma_g_exc,delta_g_dir,brout,n_g_all,n_g, &
       !$OMP&        V_g_coul,coulf_g_dir,coulf_g_exc,coulf_d_dir,icou,e2charg) &
       !$OMP& PRIVATE(nzac,Zrho_n_dir,Zrho_p_dir,Zrho_n_exc,Zrho_p_exc,&
       !$OMP&         Zkap_n_dir,Zkap_p_dir,zc,za,itb,rb,lb,sb,npb,nlb,itd,npd,&
       !$OMP&         rd,ld,sd,nld,izb,zb,ibb,idb,imb,ibbx,jb,nb,izd,zd,ibd,jd,&
       !$OMP&         imd,nd,ndb,rho_n,rho_p,kap_n,kap_p,iz_abdc,iz_abcd,&
       !$OMP&         iz_acbd,ig,V_abcd,V_abdc,V_acbd,ita,ra,la,sa,iba,ida,ima,&
       !$OMP&         ibax,na,itc,rc,lc,sc,ibc,imc,nc,sac,gamma_n_dir,&
       !$OMP&         gamma_p_dir,gamma_n_exc,gamma_p_exc,delta_n_dir,&
       !$OMP&         delta_p_dir,ob,sdb,ir_abdc,ir_abcd,ir_acbd,jr_abdc,&
       !$OMP&         jr_abcd,jr_acbd,W,B,H,M,zrho_nd,zrho_pd,zrho_ne,zrho_pe,&
       !$OMP&         zkap_n,zkap_p,U_abcd,U_abdc,U_acbd,Vdir,Udir,Vexc,Uexc,&
       !$OMP&         Vdel,Udel,i,i1,i2,i3,nac,&
       !$OMP&         gamma_c_dir,gamma_c_exc,delta_c_dir)
       !$OMP DO SCHEDULE(DYNAMIC)
       Do nzac = 1,(nzx+1)**2
          Zrho_n_dir = 0
          Zrho_p_dir = 0
          Zrho_n_exc = 0
          Zrho_p_exc = 0
          Zkap_n_dir = 0
          Zkap_p_dir = 0
          zc = mod(nzac-1,nzx+1)
          za = (nzac-zc-1)/(nzx+1)
          Do itb = 1,nttx
             rb = nrr(itb); lb = nll(itb); sb = nss(itb); npb = noo(itb)
             nlb = mod(npb-lb+1,2)
             Do itd = 1,nttx
                npd = noo(itd)
                If(npb /= npd) Cycle
                rd = nrr(itd); ld = nll(itd); sd = nss(itd)
                nld = mod(npd-ld+1,2)
                Do izb = 1,nzzx(itb)
                   zb = nzz(itb,izb)
                   ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
                   idb = id(ibb)
                   imb=ia(ibb); ibbx=ibb+nbx
                   jb = i_zrls(zb,rb,lb,(sb+1)/2)
                   nb = jb-imb
                   Do izd = 1,nzzx(itd)
                      zd = nzz(itd,izd)
                      If(mod(za+zb+zc+zd,2) /= 0) Cycle
                      ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                      If(ibb /= ibd) Cycle
                      jd = i_zrls(zd,rd,ld,(sd+1)/2)
                      imd=ia(ibd)
                      nd = jd - imd
                      ndb = nd + (nb-1)*idb
                      rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                      rho_p = rk(ndb,ibbx)*0.5_pr
                      kap_n = ak(ndb,ibb )        !kappa_db
                      kap_p = ak(ndb,ibbx)
                      iz_abdc=zindex(za,zb,zd,zc)
                      If(zc == zd) Then
                         iz_abcd=iz_abdc
                      Else
                         iz_abcd=zindex(za,zb,zc,zd)
                      End If
                      If(zd == zb) Then
                         iz_acbd=iz_abdc
                      Else If(zc == zb) Then
                         iz_acbd=iz_abcd
                      Else
                         iz_acbd=zindex(za,zc,zb,zd)
                      End If
                      Do ig = 1,n_g_all
                         V_abcd = VzGogny(ig,iz_abcd)
                         V_abdc = VzGogny(ig,iz_abdc)
                         V_acbd = VzGogny(ig,iz_acbd)
                         Zrho_n_dir(ig,rd,rb,nld,nlb,npb) = Zrho_n_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_n
                         Zrho_p_dir(ig,rd,rb,nld,nlb,npb) = Zrho_p_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_p
                         Zrho_n_exc(ig,rd,rb,nld,nlb,npb) = Zrho_n_exc(ig,rd,rb,nld,nlb,npb)+V_abdc*rho_n
                         Zrho_p_exc(ig,rd,rb,nld,nlb,npb) = Zrho_p_exc(ig,rd,rb,nld,nlb,npb)+V_abdc*rho_p
                         Zkap_n_dir(ig,rd,rb,nld,nlb,npb) = Zkap_n_dir(ig,rd,rb,nld,nlb,npb)+V_acbd*kap_n
                         Zkap_p_dir(ig,rd,rb,nld,nlb,npb) = Zkap_p_dir(ig,rd,rb,nld,nlb,npb)+V_acbd*kap_p
                      End Do !ig
                   End Do ! izd
                End Do !izb
             End Do !itd
          End Do !itb
          Do ita = 1,ntx
             If(za /= nz(ita)) Cycle
             ra = nr(ita); la = nl(ita); sa = ns(ita)
             iba = ib_zrls(za,ra,la,(sa+1)/2)
             ida=id(iba); ima=ia(iba); ibax=iba+nbx
             na = ita-ima
             Do itc = 1, ita
                If(zc /= nz(itc)) Cycle
                rc = nr(itc); lc = nl(itc); sc = ns(itc)
                ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
                If(ibc /= iba) Cycle
                imc=ia(ibc)
                nc = itc-imc
                sac = sa + sc
                gamma_n_dir = zero
                gamma_p_dir = zero
                gamma_n_exc = zero
                gamma_p_exc = zero
                delta_n_dir = zero
                delta_p_dir = zero
                gamma_c_dir = zero
                gamma_c_exc = zero
                delta_c_dir = zero
                Do ob = 0,nox
                   ibb = ob + 1
                   Do nlb = 0,1
                      lb = ob+nlb
                      sb = -2*nlb+1
                      Do nld = 0,1
                         ld = ob+nld
                         If(la+lb /= lc+ld .And. la-lb /= lc-ld) Cycle
                         sd = -2*nld+1
                         sdb = sd + sb
                         Do rb = 0,nrx
                            If(2*rb+lb > nrlx) Cycle
                            Do rd = 0,nrx
                               If(2*rd+ld > nrlx) Cycle
                               If(sac /= 0) Then
                                  ir_abdc=rindex(ra,rb,rd,rc,la, lb, ld, lc,n_rlx)
                                  If(rb == rd .And. lb == ld) Then
                                     ir_acbd=ir_abdc
                                  Else
                                     ir_acbd=rindex(ra,rc,rb,rd,la,-lc,lb,-ld,n_rlx)
                                  End If
                                  If(rc == rd .And. lc == ld) Then
                                     ir_abcd = ir_abdc
                                  Else If(rc == rb .And. lc == lb) Then
                                     ir_abcd = ir_acbd
                                  Else
                                     ir_abcd=rindex(ra,rb,rc,rd,la,lb,lc,ld,n_rlx)
                                  End If
                                  If(lb == 0 .And. ld == 0) Then
                                     jr_abdc = ir_abdc
                                     jr_acbd = ir_acbd
                                     jr_abcd = ir_abcd
                                  Else
                                     jr_abdc=rindex(ra,rb,rd,rc,la,-lb,-ld,lc,n_rlx)
                                     If(rb == rd .And. lb == ld) Then
                                        jr_acbd = jr_abdc
                                     Else
                                        jr_acbd=rindex(ra,rc,rb,rd,la,-lc,-lb,ld,n_rlx)
                                     End If
                                     If(rc == rb .And. lb == lc) Then
                                        jr_abcd = ir_acbd
                                     Else If(rc == rd .And. lc == ld) Then
                                        jr_abcd = ir_abdc
                                     Else
                                        jr_abcd=rindex(ra,rb,rc,rd,la,-lb,lc,-ld,n_rlx)
                                     End If
                                  End If
                               Else
                                  If(sa == sb) Then
                                     jr_abdc=rindex(ra,rb,rd,rc,la,-lb,-ld, lc,n_rlx)
                                     ir_acbd=rindex(ra,rc,rb,rd,la,-lc, lb,-ld,n_rlx)
                                  Else
                                     ir_abdc=rindex(ra,rb,rd,rc,la, lb, ld,lc,n_rlx)
                                     jr_acbd=rindex(ra,rc,rb,rd,la,-lc,-lb,ld,n_rlx)
                                  End If
                               End If
                               Do ig = 1,n_g
                                  W = W_g(ig)
                                  B = B_g(ig)
                                  H = H_g(ig)
                                  M = M_g(ig)
                                  zrho_nd = Zrho_n_dir(ig,rd,rb,nld,nlb,ibb)
                                  zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                                  zrho_ne = Zrho_n_exc(ig,rd,rb,nld,nlb,ibb)
                                  zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                  zkap_n  = Zkap_n_dir(ig,rd,rb,nld,nlb,ibb)
                                  zkap_p  = Zkap_p_dir(ig,rd,rb,nld,nlb,ibb)
                                  If(sac /= 0 .And. sdb /= 0) Then
                                     V_abdc=VrGogny(ig,ir_abdc)
                                     U_abdc=VrGogny(ig,jr_abdc)
                                     V_abcd=VrGogny(ig,ir_abcd)
                                     U_abcd=VrGogny(ig,jr_abcd)
                                     V_acbd=VrGogny(ig,ir_acbd)
                                     U_acbd=VrGogny(ig,jr_acbd)
                                     If(sa == sb) Then
                                        !up up up up|or|down down down down
                                        Vdir = V_abcd
                                        Udir = U_abcd
                                        Vexc = V_abdc
                                        Uexc = U_abdc
                                        Vdel = V_acbd
                                        Udel = U_acbd
                                     Else
                                        !up up down down|or|down down up up
                                        Vdir = U_abcd
                                        Udir = V_abcd
                                        Vexc = U_abdc
                                        Uexc = V_abdc
                                        Vdel = U_acbd
                                        Udel = V_acbd
                                     End If
                                     gamma_n_dir = gamma_n_dir + Vdir*((W+B-M-H)*zrho_nd + (W+B)*zrho_pd) &
                                                               + Udir*(    (W-H)*zrho_nd + W*zrho_pd)
                                     gamma_p_dir = gamma_p_dir + Vdir*((W+B-M-H)*zrho_pd + (W+B)*zrho_nd) &
                                                               + Udir*(    (W-H)*zrho_pd + W*zrho_nd)
                                     gamma_n_exc = gamma_n_exc + Vexc*((M+H-W-B)*zrho_ne + (M+H)*zrho_pe) &
                                                               + Uexc*(    (M-B)*zrho_ne + M*zrho_pe)
                                     gamma_p_exc = gamma_p_exc + Vexc*((M+H-W-B)*zrho_pe + (M+H)*zrho_ne) &
                                                               + Uexc*(    (M-B)*zrho_pe + M*zrho_ne)
                                     delta_n_dir = delta_n_dir + (Vdel*(W-H)+Udel*(M-B))*zkap_n
                                     delta_p_dir = delta_p_dir + (Vdel*(W-H)+Udel*(M-B))*zkap_p
                                  Else If(sac == 0 .And. sdb == 0) Then
                                     If(sa == sb) Then
                                        !up down down up|or|down up up down
                                        Vexc = -VrGogny(ig,jr_abdc)
                                        Vdel =  VrGogny(ig,ir_acbd)
                                     Else
                                        !up down up down|or|down up down up
                                        Vexc =  VrGogny(ig,ir_abdc)
                                        Vdel = -VrGogny(ig,jr_acbd)
                                     End If
                                     gamma_n_exc = gamma_n_exc + Vexc*((H-W)*zrho_ne + H*zrho_pe)
                                     gamma_p_exc = gamma_p_exc + Vexc*((H-W)*zrho_pe + H*zrho_ne)
                                     delta_n_dir = delta_n_dir + Vdel*zkap_n*(W-H+B-M)
                                     delta_p_dir = delta_p_dir + Vdel*zkap_p*(W-H+B-M)
                                  End If
                               End Do !ig
                               Do ig = n_g+1,n_g_all
                                  W = V_g_coul(ig-n_g) * e2charg
                                  If(icou == -5) Then
                                     zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                                     zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                     zkap_p  = Zkap_p_dir(ig,rd,rb,nld,nlb,ibb)
                                  Else If(icou == -4) Then
                                     zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                                     zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                     zkap_p  = 0.0_pr
                                  Else If(icou == -3) Then
                                     zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                                     zrho_pe = 0.0_pr
                                     zkap_p  = 0.0_pr
                                  Else If(icou == -2) Then
                                     zrho_pd = 0.0_pr
                                     zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                     zkap_p  = 0.0_pr
                                  Else
                                     zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                                     zrho_pe = 0.0_pr
                                     zkap_p  = 0.0_pr
                                  End If
                                  If(sac /= 0 .And. sdb /= 0) Then
                                     V_acbd=VrGogny(ig,ir_acbd)
                                     U_acbd=VrGogny(ig,jr_acbd)
                                     If(sa == sb) Then
                                        !up up up up|or|down down down down
                                        Vdir = VrGogny(ig,ir_abcd)
                                        Udir = VrGogny(ig,jr_abcd)
                                        If(icou == -5) Then
                                           Vexc = VrGogny(ig,ir_abdc)
                                           Vdel = V_acbd
                                        Else If(icou == -4) Then
                                           Vexc = VrGogny(ig,ir_abdc)
                                           Vdel = 0.0_pr
                                        Else If(icou == -3) Then
                                           Vexc = 0.0_pr
                                           Vdel = 0.0_pr
                                        Else If(icou == -2) Then
                                           Vexc = VrGogny(ig,ir_abdc)
                                           Vdel = 0.0_pr
                                        Else
                                           Vexc = 0.0_pr
                                           Vdel = 0.0_pr
                                        End If
                                     Else
                                        !up up down down|or|down down up up
                                        Vdir = VrGogny(ig,jr_abcd)
                                        Udir = VrGogny(ig,ir_abcd)
                                        If(icou == -5) Then
                                           Vexc = VrGogny(ig,jr_abdc)
                                           Vdel = U_acbd
                                        Else If(icou == -4) Then
                                           Vexc = VrGogny(ig,jr_abdc)
                                           Vdel = 0.0_pr
                                        Else If(icou == -3) Then
                                           Vexc = 0.0_pr
                                           Vdel = 0.0_pr
                                        Else If(icou == -2) Then
                                           Vexc = VrGogny(ig,jr_abdc)
                                           Vdel = 0.0_pr
                                        Else
                                           Vexc = 0.0_pr
                                           Vdel = 0.0_pr
                                        End If
                                     End If
                                     gamma_c_dir = gamma_c_dir + (Vdir+Udir)*W*zrho_pd
                                     gamma_c_exc = gamma_c_exc - Vexc*W*zrho_pe
                                     delta_c_dir = delta_c_dir + Vdel*W*zkap_p
                                  Else If(sac == 0 .And. sdb == 0) Then
                                     If(sa == sb) Then
                                        !up down down up|or|down up up down
                                        If(icou == -5) Then
                                           Vexc = -VrGogny(ig,jr_abdc)
                                           Vdel =  VrGogny(ig,ir_acbd)
                                        Else If(icou == -4) Then
                                           Vexc = -VrGogny(ig,jr_abdc)
                                           Vdel =  0.0_pr
                                        Else If(icou == -3) Then
                                           Vexc =  0.0_pr
                                           Vdel =  0.0_pr
                                        Else If(icou == -2) Then
                                           Vexc = -VrGogny(ig,jr_abdc)
                                           Vdel =  0.0_pr
                                        Else
                                           Vexc =  0.0_pr
                                           Vdel =  0.0_pr
                                        End If
                                     Else
                                        !up down up down|or|down up down up
                                        If(icou == -5) Then
                                           Vexc =  VrGogny(ig,ir_abdc)
                                           Vdel = -VrGogny(ig,jr_acbd)
                                        Else If(icou == -4) Then
                                           Vexc =  VrGogny(ig,ir_abdc)
                                           Vdel =  0.0_pr
                                        Else If(icou == -3) Then
                                           Vexc =  0.0_pr
                                           Vdel =  0.0_pr
                                        Else If(icou == -2) Then
                                           Vexc =  VrGogny(ig,ir_abdc)
                                           Vdel =  0.0_pr
                                        Else
                                           Vexc =  0.0_pr
                                           Vdel =  0.0_pr
                                        End If
                                     End If
                                     gamma_c_exc = gamma_c_exc - Vexc*W*zrho_pe
                                     delta_c_dir = delta_c_dir + Vdel*W*zkap_p
                                  End If
                               End Do !ig
                            End Do !rd
                         End Do !rb
                      End Do !nld
                   End Do !nlb
                End Do !ob
                i=allibro(iba)+na*(na-1)/2+nc
                i1=i+nhhdim
                i2=i+nhhdim2
                i3=i+nhhdim3
                nac = na+(nc-1)*ida;
                gamma_g_dir(nac,iba )=gamma_n_dir
                gamma_g_dir(nac,ibax)=gamma_p_dir
                gamma_g_exc(nac,iba )=gamma_n_exc
                gamma_g_exc(nac,ibax)=gamma_p_exc
                delta_g_dir(nac,iba )=delta_n_dir
                delta_g_dir(nac,ibax)=delta_p_dir
                coulf_g_dir(nac,ibax)=gamma_c_dir
                coulf_g_exc(nac,ibax)=gamma_c_exc
                coulf_d_dir(nac,ibax)=delta_c_dir
                brout(i )=brout(i )+gamma_n_dir+gamma_n_exc
                brout(i1)=brout(i1)+gamma_p_dir+gamma_p_exc+gamma_c_dir+gamma_c_exc
                brout(i2)=brout(i2)+delta_n_dir
                brout(i3)=brout(i3)+delta_p_dir+delta_c_dir
             End Do !itc
          End Do !ita
       End Do !nzac
       !$OMP End Do
       !$OMP End Parallel
    Else
       !$OMP Parallel Default(None)&
       !$OMP& SHARED(nzx,nttx,nrr,nll,nss,noo,nzzx,nzz,ib_zrls,id,ia,nbx,&
       !$OMP&        i_zrls,rk,ak,VzGogny,ntx,nz,nr,nl,ns,nox,nrx,nrlx,n_rlx,&
       !$OMP&        W_g,B_g,H_g,M_g,VrGogny,allibro,nhhdim,nhhdim2,nhhdim3,&
       !$OMP&        gamma_g_dir,gamma_g_exc,delta_g_dir,brout,n_g,force_is_dme) &
       !$OMP& PRIVATE(nzac,Zrho_n_dir,Zrho_p_dir,Zrho_n_exc,Zrho_p_exc,&
       !$OMP&         Zkap_n_dir,Zkap_p_dir,zc,za,itb,rb,lb,sb,npb,nlb,itd,npd,&
       !$OMP&         rd,ld,sd,nld,izb,zb,ibb,idb,imb,ibbx,jb,nb,izd,zd,ibd,jd,&
       !$OMP&         imd,nd,ndb,rho_n,rho_p,kap_n,kap_p,iz_abdc,iz_abcd,&
       !$OMP&         iz_acbd,ig,V_abcd,V_abdc,V_acbd,ita,ra,la,sa,iba,ida,ima,&
       !$OMP&         ibax,na,itc,rc,lc,sc,ibc,imc,nc,sac,gamma_n_dir,&
       !$OMP&         gamma_p_dir,gamma_n_exc,gamma_p_exc,delta_n_dir,&
       !$OMP&         delta_p_dir,ob,sdb,ir_abdc,ir_abcd,ir_acbd,jr_abdc,&
       !$OMP&         jr_abcd,jr_acbd,W,B,H,M,zrho_nd,zrho_pd,zrho_ne,zrho_pe,&
       !$OMP&         zkap_n,zkap_p,U_abcd,U_abdc,U_acbd,Vdir,Udir,Vexc,Uexc,&
       !$OMP&         Vdel,Udel,i,i1,i2,i3,nac)
       !$OMP DO SCHEDULE(DYNAMIC)
       Do nzac = 1,(nzx+1)**2
          Zrho_n_dir = 0
          Zrho_p_dir = 0
          zc = mod(nzac-1,nzx+1)
          za = (nzac-zc-1)/(nzx+1)
          Do itb = 1,nttx
             rb = nrr(itb); lb = nll(itb); sb = nss(itb); npb = noo(itb)
             nlb = mod(npb-lb+1,2)
             Do itd = 1,nttx
                npd = noo(itd)
                If(npb /= npd) Cycle
                rd = nrr(itd); ld = nll(itd); sd = nss(itd)
                nld = mod(npd-ld+1,2)
                Do izb = 1,nzzx(itb)
                   zb = nzz(itb,izb)
                   ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
                   idb = id(ibb)
                   imb=ia(ibb); ibbx=ibb+nbx
                   jb = i_zrls(zb,rb,lb,(sb+1)/2)
                   nb = jb-imb
                   Do izd = 1,nzzx(itd)
                      zd = nzz(itd,izd)
                      If(mod(za+zb+zc+zd,2) /= 0) Cycle
                      ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                      If(ibb /= ibd) Cycle
                      jd = i_zrls(zd,rd,ld,(sd+1)/2)
                      imd=ia(ibd)
                      nd = jd - imd
                      ndb = nd + (nb-1)*idb
                      rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                      rho_p = rk(ndb,ibbx)*0.5_pr
                      iz_abcd=zindex(za,zb,zc,zd)
                      Do ig = 1,n_g
                         V_abcd = VzGogny(ig,iz_abcd)
                         Zrho_n_dir(ig,rd,rb,nld,nlb,npb) = Zrho_n_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_n
                         Zrho_p_dir(ig,rd,rb,nld,nlb,npb) = Zrho_p_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_p
                      End Do !ig
                   End Do ! izd
                End Do !izb
             End Do !itd
          End Do !itb
          Do ita = 1,ntx
             If(za /= nz(ita)) Cycle
             ra = nr(ita); la = nl(ita); sa = ns(ita)
             iba = ib_zrls(za,ra,la,(sa+1)/2)
             ida=id(iba); ima=ia(iba); ibax=iba+nbx
             na = ita-ima
             Do itc = 1, ita
                If(zc /= nz(itc)) Cycle
                rc = nr(itc); lc = nl(itc); sc = ns(itc)
                ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
                If(ibc /= iba) Cycle
                imc=ia(ibc)
                nc = itc-imc
                sac = sa + sc
                If(sac == 0) Cycle
                gamma_n_dir = zero
                gamma_p_dir = zero
                Do ob = 0,nox
                   ibb = ob + 1
                   Do nlb = 0,1
                      lb = ob+nlb
                      sb = -2*nlb+1
                      Do nld = 0,1
                         ld = ob+nld
                         If(la+lb /= lc+ld .And. la-lb /= lc-ld) Cycle
                         sd = -2*nld+1
                         sdb = sd + sb
                         If(sdb == 0) Cycle
                         Do rb = 0,nrx
                            If(2*rb+lb > nrlx) Cycle
                            Do rd = 0,nrx
                               If(2*rd+ld > nrlx) Cycle
                               ir_abcd=rindex(ra,rb,rc,rd,la,lb,lc,ld,n_rlx)
                               If(lb == 0 .And. ld == 0) Then
                                  jr_abcd = ir_abcd
                               Else
                                  jr_abcd=rindex(ra,rb,rc,rd,la,-lb,lc,-ld,n_rlx)
                               End If
                               Do ig = 1,n_g
                                  W = W_g(ig)
                                  B = B_g(ig)
                                  H = H_g(ig)
                                  M = M_g(ig)
                                  zrho_nd = Zrho_n_dir(ig,rd,rb,nld,nlb,ibb)
                                  zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                                  V_abcd=VrGogny(ig,ir_abcd)
                                  U_abcd=VrGogny(ig,jr_abcd)
                                  If(sa == sb) Then
                                     !up up up up|or|down down down down
                                     Vdir = V_abcd
                                     Udir = U_abcd
                                  Else
                                     !up up down down|or|down down up up
                                     Vdir = U_abcd
                                     Udir = V_abcd
                                  End If
                                  gamma_n_dir = gamma_n_dir + Vdir*((W+B-M-H)*zrho_nd + (W+B)*zrho_pd)&
                                                            + Udir*(    (W-H)*zrho_nd + W*zrho_pd)
                                  gamma_p_dir = gamma_p_dir + Vdir*((W+B-M-H)*zrho_pd + (W+B)*zrho_nd)&
                                                            + Udir*(    (W-H)*zrho_pd + W*zrho_nd)
                               End Do !ig
                            End Do !rd
                         End Do !rb
                      End Do !nld
                   End Do !nlb
                End Do !ob
                i=allibro(iba)+na*(na-1)/2+nc
                i1=i+nhhdim
                i2=i+nhhdim2
                i3=i+nhhdim3
                nac = na+(nc-1)*ida;
                gamma_g_dir(nac,iba )=gamma_n_dir
                gamma_g_dir(nac,ibax)=gamma_p_dir
                brout(i )=brout(i )+gamma_n_dir
                brout(i1)=brout(i1)+gamma_p_dir
             End Do !itc
          End Do !ita
       End Do !nzac
       !$OMP End Do
       !$OMP End Parallel
    End If
#else
    Do ita = 1,ntx
       ra = nr(ita); za = nz(ita); la = nl(ita); sa = ns(ita)
       iba = ib_zrls(za,ra,la,(sa+1)/2)
       ida=id(iba); ima=ia(iba); ibax=iba+nbx
       na = ita-ima
       Do itc = 1, ita
          rc = nr(itc); zc = nz(itc); lc = nl(itc); sc = ns(itc)
          ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
          idc=id(ibc); imc=ia(ibc); ibcx=ibc+nbx
          nc = itc-imc
          sac = sa + sc
          If(ibc /= iba) Cycle
          gamma_n_dir = zero
          gamma_p_dir = zero
          gamma_n_exc = zero
          gamma_p_exc = zero
          delta_n_dir = zero
          delta_p_dir = zero
          Do itb = 1,ntx
             rb = nr(itb); zb = nz(itb); lb = nl(itb); sb = ns(itb)
             ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
             idb=id(ibb); imb=ia(ibb); ibbx=ibb+nbx
             nb = itb-imb
             Do itd = 1,ntx
                rd = nr(itd); zd = nz(itd); ld = nl(itd); sd = ns(itd)
                ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                idd=id(ibd); imd=ia(ibd); ibdx=ibd+nbx
                nd = itd-imd
                sdb = sd + sb
                If(ibb /= ibd) Cycle
                ndb = nd + (nb-1)*idd
                rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                rho_p = rk(ndb,ibbx)*0.5_pr
                kap_n = ak(ndb,ibb )        !kappa_db
                kap_p = ak(ndb,ibbx)
                Do ig = 1,n_g
                   W = W_g(ig)
                   B = B_g(ig)
                   H = H_g(ig)
                   M = M_g(ig)
                   If(sac /= 0 .And. sdb /= 0) Then
                      V_abdc=Vr_Gogny(ig,ra,la,rb, lb,rd, ld,rc, lc)
                      U_abdc=Vr_Gogny(ig,ra,la,rb,-lb,rd,-ld,rc, lc)
                      V_abcd=Vr_Gogny(ig,ra,la,rb, lb,rc, lc,rd, ld)
                      U_abcd=Vr_Gogny(ig,ra,la,rb,-lb,rc, lc,rd,-ld)
                      V_acbd=Vr_Gogny(ig,ra,la,rc,-lc,rb, lb,rd,-ld)
                      U_acbd=Vr_Gogny(ig,ra,la,rc,-lc,rb,-lb,rd, ld)
                      If(sa == sb) Then
                         ! up up   up up   |or| down down down down
                         Vdir = V_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Udir = U_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Vexc = V_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Uexc = U_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel = V_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                         Udel = U_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                      Else
                         !up up  down down |or| down down up up
                         Vdir = U_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Udir = V_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Vexc = U_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Uexc = V_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel = U_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                         Udel = V_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                      End If
                      gamma_n_dir = gamma_n_dir + Vdir*((W+B-M-H)*rho_n + (W+B)*rho_p)&
                                                + Udir*(    (W-H)*rho_n +     W*rho_p)
                      gamma_p_dir = gamma_p_dir + Vdir*((W+B-M-H)*rho_p + (W+B)*rho_n)&
                                                + Udir*(    (W-H)*rho_p +     W*rho_n)
                      gamma_n_exc = gamma_n_exc + Vexc*((M+H-W-B)*rho_n + (M+H)*rho_p)&
                                                + Uexc*(    (M-B)*rho_n +     M*rho_p)
                      gamma_p_exc = gamma_p_exc + Vexc*((M+H-W-B)*rho_p + (M+H)*rho_n)&
                                                + Uexc*(    (M-B)*rho_p +     M*rho_n)
                      delta_n_dir = delta_n_dir + (Vdel*(W-H)+Udel*(M-B))*kap_n
                      delta_p_dir = delta_p_dir + (Vdel*(W-H)+Udel*(M-B))*kap_p
                   Else If(sac == 0 .And. sdb == 0) Then
                      If(sa == sb) Then
                         !up down down up|or| down up up down
                         Vexc = -Vr_Gogny(ig,ra,la,rb,-lb,rd,-ld,rc,lc)*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel =  Vr_Gogny(ig,ra,la,rc,-lc,rb,lb,rd,-ld)*Vz_Gogny(ig,za,zc,zb,zd)
                      Else
                         !up down up down |or| down up down up
                         Vexc = Vr_Gogny(ig,ra,la,rb,lb,rd,ld,rc,lc)*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel =-Vr_Gogny(ig,ra,la,rc,-lc,rb,-lb,rd,ld)*Vz_Gogny(ig,za,zc,zb,zd)
                      End If
                      gamma_n_exc = gamma_n_exc + Vexc*((H-W)*rho_n + H*rho_p)
                      gamma_p_exc = gamma_p_exc + Vexc*((H-W)*rho_p + H*rho_n)
                      delta_n_dir = delta_n_dir + Vdel*kap_n*(W-H+B-M)
                      delta_p_dir = delta_p_dir + Vdel*kap_p*(W-H+B-M)
                   End If
                End Do !ig
             End Do !itd
          End Do !itb
          If(force_is_dme) Then
             !only the Hartree Part
             gamma_n_exc = zero
             gamma_p_exc = zero
             delta_n_dir = zero
             delta_p_dir = zero
          End If
          i=allibro(iba)+na*(na-1)/2+nc
          i1=i+nhhdim
          i2=i+nhhdim2
          i3=i+nhhdim3
          nac = na+(nc-1)*ida;
          gamma_g_dir(nac,iba )=gamma_n_dir
          gamma_g_dir(nac,ibax)=gamma_p_dir
          gamma_g_exc(nac,iba )=gamma_n_exc
          gamma_g_exc(nac,ibax)=gamma_p_exc
          delta_g_dir(nac,iba )=delta_n_dir
          delta_g_dir(nac,ibax)=delta_p_dir
          brout(i )=brout(i )+gamma_n_dir+gamma_n_exc
          brout(i1)=brout(i1)+gamma_p_dir+gamma_p_exc
          brout(i2)=brout(i2)+delta_n_dir
          brout(i3)=brout(i3)+delta_p_dir
       End Do !itc
    End Do !ita
#endif

  End Subroutine gamdel_gogny
  !======================================================================
  !> This routine calculates the two-body matrix elements \f$ \langle n_1 n_2 |  V_i(r) | n_3 n_4 \rangle \f$
  !> for a Gaussian potential \f$ V_i(r) = e^{(\mathbf{r}_1-\mathbf{r}_2)^2/\mu_i^2}\f$. The matrix elements
  !> are separated as the product of axial \f$ V^z_{n_1n_2n_3n_4}\f$ and radial \f$V^r_{n_1n_2n_3n_4}\f$
  !> components.
  !>
  !> If the preprocessor directive GOGNY_HYPER is set to 1, \f$ V^z_{n_1n_2n_3n_4} \f$ is calculated with
  !> an expansion based on the hypergeometric function \f$_2F_1\f$. This expansion maintains accuracy with
  !> increasing basis size; see \ref matrixelement_z() for details. \f$V^r_{n_1n_2n_3n_4} \f$ is calculated
  !> by expanding each two-dimensional harmonic oscillator function as a sum of products of two one-dimensional
  !> harmonic oscillator functions; see \ref radial_matrix_elements() for more details. This expansion allows
  !> to calculate \f$V^r\f$ with the same subroutine used to calculate \f$V^z\f$ in order to preserve accuracy
  !> as the basis size increases.
  !>
  !> For any other value of GOGNY_HYPER a set of more direct expansions is used to calculate \f$V^z\f$ and
  !> \f$V^r\f$; see \ref matrixelement_z() and \ref matrixelement_r() for more details. These expansions
  !> start to loose accuracy as the basis size increases.
  !>
  !> The calculation of the radial components is timed and the required wall clock time is printed as output.
  !======================================================================
  Subroutine gogny_matrix_elements()
    Use HFBTHO, Only : do_print

    Integer(ipr) ::  t1,t2,countrate,countmax,iw

    If(matrix_elements_calculated) Return

    If(do_print == 1) Then
       Do iw = lout,lfile
          Write(iw,*)
          Write(iw,*) ' Calculating the finite range matrix elements'
       End Do
    End If

    Call calculateTz()
#if(GOGNY_HYPER==1)
    Call calculateCpolar2cartesian()
    Call calculateME1D()
#else
    Call calculateTr()
#endif
    Call calculateVzGogny()
    Call system_clock(t1,countrate,countmax)
    Call calculateVrGogny()
    Call system_clock(t2,countrate,countmax)

    If(do_print == 1) Then
       Do iw = lout,lfile
          Write(iw,'(a33,f15.4)') '  Matrix elements wall clock time', (t2-t1)/real(countrate,kind=pr)
       End Do
    End If

    matrix_elements_calculated = .True.

  End Subroutine gogny_matrix_elements
  !======================================================================
  !> This routine calculates the the necessary two-body potential matrix elements
  !> \f$\langle n_{z_i} n_{z_j}|V_z|n_{z_k} n_{z_l} \rangle \f$.
  !>
  !> If the preprocessor variable GOGNY_SYMMETRIES is set to 1 only states with \f$ n_{z_k} \geq n_{z_i} \f$,
  !> \f$ n_{z_l} \geq n_{z_j} \f$ are calculated and stored. In any case, states where
  !> \f$ n_{z_i}+n_{z_j}+n_{z_k}+n_{z_l} \f$ is not an even number are not calculated since those matrix
  !> elements are zero.
  !======================================================================
  Subroutine CalculateVzGogny()
    Use UNEDF, Only : n_g_all,mu_g_all
    Use HFBTHO, Only : nzx,bz

    Real(pr) :: Vz
    Integer(ipr) :: ii,N,nzi,nzj,nzk,nzl,ig

    n = nzx
#if(GOGNY_SYMMETRIES==1)
    Call calculate_Zblock()
    ii = Zblock(n,n,n,n); NumVz = ii
    If(Allocated(VzGogny)) Deallocate(VzGogny)
    Allocate(VzGogny(1:n_g_all,0:ii))
    VzGogny(1:n_g_all,0) = zero
    !$OMP Parallel Default(None) &
    !$OMP& SHARED(n,Zblock,VzGogny,mu_g_all,bz,n_g_all) &
    !$OMP& PRIVATE(nzi,nzj,nzk,nzl,ii,ig,Vz)
    !$OMP DO SCHEDULE(DYNAMIC)
    Do nzi = 0,n
       Do nzj = 0,n
           Do nzk = 0,n
             Do nzl = 0,n
                If(nzk < nzi) Cycle
                If(nzl < nzj) Cycle
                If(mod(nzi+nzj+nzk+nzl,2) /= 0) Cycle
                ii = Zblock(nzi,nzj,nzk,nzl)
                Do ig = 1,n_g_all
                   Vz = MatrixElement_z(nzi,nzj,nzk,nzl,mu_g_all(ig),bz)
                   VzGogny(ig,ii) = Vz
                End Do
             End Do
          End Do
       End Do
    End Do
    !$OMP End Do
    !$OMP End Parallel
#else
    If(Allocated(Vz_Gogny)) Deallocate(Vz_Gogny)
    Allocate(Vz_Gogny(1:n_g_all,0:n,0:n,0:n,0:n))
    Vz_gogny = zero
    Do nzi = 0,n
       Do nzj = 0,n
           Do nzk = 0,n
             Do nzl = 0,n
                If(mod(nzi+nzj+nzk+nzl,2) /= 0) Cycle
                Do ig = 1,n_g_all
                   Vz = MatrixElement_z(nzi,nzj,nzk,nzl,mu_g_all(ig),bz)
                   Vz_Gogny(ig,nzi,nzj,nzk,nzl) = Vz
                End Do
             End Do
          End Do
       End Do
    End Do
#endif

  End Subroutine CalculateVzGogny
  !======================================================================
  !> This routine calculates the necessary two-body potential matrix elements
  !> \f$ \langle n_{r_i} \Lambda_i, n_{r_j} \Lambda_j |V_p | n_{r_k} \Lambda_k, n_{r_l} \Lambda_l \rangle\f$.
  !>
  !> If the preprocessor variable GOGNY_SYMMETRIES is set to 1 only the matrix elements that will be
  !> used in the calculation of the HFB fields are calculated. In any case, states where \f$ -\Lambda_i
  !> - \Lambda_j + \Lambda_k + \Lambda_l \neq 0 \f$ are not calculated since the matrix element is zero.
  !======================================================================
  Subroutine CalculateVrGogny()
    Use UNEDF, Only : n_g_all,force_is_DME
    Use HFBTHO, Only : nr,nz,nl,ns,ntx,nrx,nlx,nttx,nrr,nll,nss,noo

#if(GOGNY_SYMMETRIES==1)
    Integer(ipr) :: ii,ig,im,n
    Integer(ipr) :: nza,nzc,nra,nrb,nrc,nrd
    Integer(ipr) :: nla,nlb,nlc,nld,nsa,nsc,nsac,nsb,nsd,nsdb,ir_abcd
    Integer(ipr) :: ita,itb,itc,itd,iba,ibc,ibd,ibb,ir
    Integer(ipr) , Allocatable :: index_flag(:), nr_flag(:,:), nl_flag(:,:),ir_flag(:)
    Integer(ipr) :: ilauf
#else
    Integer(ipr) :: n,nra,nrb,nrc,nrd,nla,nlb,nlc,nld,ig
#endif

    N = Max(2*nrx,nlx)
    im = mod(N,2)
#if(GOGNY_SYMMETRIES==1)
    Call calculate_Nblock()
    ii = nblock(n/2,n/2,n/2,n/2,im,-im)+l_block(n/2,n/2,im,-im,-im,n); NumVr = ii
    If(Allocated(VrGogny)) Deallocate(VrGogny)
    Allocate(VrGogny(1:n_g_all,0:ii)); VrGogny=0.0_pr
    Allocate(index_flag(0:ii)); index_flag=0
    Allocate(nr_flag(1:4,0:ii))
    Allocate(nl_flag(1:4,0:ii))
    Allocate(ir_flag(0:ii))

    ! We go over the loop that the field calculation will go through in
    ! order to identify the matrix elements that will be necesary and
    ! calculate only those on a subsequent loop.
    ilauf = 0
    If(.not.force_is_DME) Then
    Do ita = 1,ntx
       nra = nr(ita); nza = nz(ita); nla = nl(ita); nsa = ns(ita)
       If(nza > 1) Cycle
       iba = ib_zrls(nza,nra,nla,(nsa+1)/2)
       Do itc = 1, ita
          nrc = nr(itc); nzc = nz(itc); nlc = nl(itc); nsc = ns(itc)
          If(nzc > 1) Cycle
          ibc = ib_zrls(nzc,nrc,nlc,(nsc+1)/2)
          nsac = nsa + nsc
          If(ibc /= iba) Cycle
          Do itb = 1,nttx
             nrb=nrr(itb); nlb=nll(itb); nsb=nss(itb); ibb=noo(itb)
             Do itd = 1,nttx
                nrd=nrr(itd); nld=nll(itd); nsd=nss(itd); ibd=noo(itd)
                If(ibb /= ibd) Cycle
                If(nla+nlb /= nlc+nld .And. nla-nlb /= nlc-nld) Cycle
                nsdb = nsd + nsb
                If(nsac /= 0) Then
                   ir=rindex(nra,nrb,nrd,nrc,nla,nlb,nld,nlc,n)
                   If(index_flag(ir) == 0) Then
                      nr_flag(1,ilauf) = nra
                      nr_flag(2,ilauf) = nrb
                      nr_flag(3,ilauf) = nrd
                      nr_flag(4,ilauf) = nrc
                      nl_flag(1,ilauf) = nla
                      nl_flag(2,ilauf) = nlb
                      nl_flag(3,ilauf) = nld
                      nl_flag(4,ilauf) = nlc
                      ir_flag(ilauf) = ir
                      index_flag(ir) = 1
                      ilauf = ilauf + 1
                   End If
                   If(nrb /= nrd.or.nlb /= nld) Then
                      ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,nlb,-nld,n)
                      If(index_flag(ir) == 0) Then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrc
                         nr_flag(3,ilauf) = nrb
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlc
                         nl_flag(3,ilauf) = nlb
                         nl_flag(4,ilauf) =-nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      End If
                   End If
                   If((nrc /= nrd .Or. nlc /= nld) .And. (nrc /= nrb .Or. nlc /= nlb)) Then
                      ir=rindex(nra,nrb,nrc,nrd,nla,nlb,nlc,nld,n)
                      If(index_flag(ir) == 0) Then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrc
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) = nlb
                         nl_flag(3,ilauf) = nlc
                         nl_flag(4,ilauf) = nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      End If
                   End If
                   If(nlb /= 0 .Or. nld /= 0) Then
                      ir=rindex(nra,nrb,nrd,nrc,nla,-nlb,-nld,nlc,n)
                      If(index_flag(ir) == 0)then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrd
                         nr_flag(4,ilauf) = nrc
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlb
                         nl_flag(3,ilauf) =-nld
                         nl_flag(4,ilauf) = nlc
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf+1
                      End If
                      If(nrb /= nrd .Or. nlb /= nld) Then
                         ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,-nlb,nld,n)
                         If(index_flag(ir) == 0) Then
                            nr_flag(1,ilauf) = nra
                            nr_flag(2,ilauf) = nrc
                            nr_flag(3,ilauf) = nrb
                            nr_flag(4,ilauf) = nrd
                            nl_flag(1,ilauf) = nla
                            nl_flag(2,ilauf) =-nlc
                            nl_flag(3,ilauf) =-nlb
                            nl_flag(4,ilauf) = nld
                            ir_flag(ilauf) = ir
                            index_flag(ir) = 1
                            ilauf = ilauf + 1
                         End If
                      End If
                      If((nrc /= nrb .Or. nlb /= nlc) .And. (nrc /= nrd .Or. nlc /= nld)) Then
                         ir=rindex(nra,nrb,nrc,nrd,nla,-nlb,nlc,-nld,n)
                         If(index_flag(ir) == 0) Then
                            nr_flag(1,ilauf)= nra
                            nr_flag(2,ilauf)= nrb
                            nr_flag(3,ilauf)= nrc
                            nr_flag(4,ilauf)= nrd
                            nl_flag(1,ilauf)= nla
                            nl_flag(2,ilauf)=-nlb
                            nl_flag(3,ilauf)= nlc
                            nl_flag(4,ilauf)=-nld
                            ir_flag(ilauf)=ir
                            index_flag(ir)=1
                            ilauf = ilauf + 1
                         End If
                      End If
                   End If
                Else
                   If(nsa == nsb) Then
                      ir=rindex(nra,nrb,nrd,nrc,nla,-nlb,-nld,nlc,n)
                      If(index_flag(ir) == 0)then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrd
                         nr_flag(4,ilauf) = nrc
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlb
                         nl_flag(3,ilauf) =-nld
                         nl_flag(4,ilauf) = nlc
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf+1
                      End If
                      ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,nlb,-nld,n)
                      If(index_flag(ir) == 0) Then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrc
                         nr_flag(3,ilauf) = nrb
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlc
                         nl_flag(3,ilauf) = nlb
                         nl_flag(4,ilauf) =-nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      End If
                   Else
                      ir=rindex(nra,nrb,nrd,nrc,nla,nlb,nld,nlc,n)
                      If(index_flag(ir) == 0) Then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrd
                         nr_flag(4,ilauf) = nrc
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) = nlb
                         nl_flag(3,ilauf) = nld
                         nl_flag(4,ilauf) = nlc
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      End If
                      ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,-nlb,nld,n)
                      If(index_flag(ir) == 0) Then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrc
                         nr_flag(3,ilauf) = nrb
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlc
                         nl_flag(3,ilauf) =-nlb
                         nl_flag(4,ilauf) = nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      End If
                   End If
                End If
             End Do !itd
          End Do !itb
       End Do !itc
    End Do !ita
    Else
    Do ita = 1,ntx
       nra = nr(ita); nza = nz(ita); nla = nl(ita); nsa = ns(ita)
       If(nza > 1) Cycle
       iba = ib_zrls(nza,nra,nla,(nsa+1)/2)
       Do itc = 1, ita
          nrc = nr(itc); nzc = nz(itc); nlc = nl(itc); nsc = ns(itc)
          If(nzc > 1) Cycle
          ibc = ib_zrls(nzc,nrc,nlc,(nsc+1)/2)
          nsac = nsa + nsc
          If(ibc /= iba) Cycle
          Do itb = 1,nttx
             nrb=nrr(itb); nlb=nll(itb); nsb=nss(itb); ibb=noo(itb)
             Do itd = 1,nttx
                nrd=nrr(itd); nld=nll(itd); nsd=nss(itd); ibd=noo(itd)
                If(ibb /= ibd) Cycle
                If(nla+nlb /= nlc+nld .And. nla-nlb /= nlc-nld) Cycle
                nsdb = nsd + nsb
                If(nsac /= 0) Then
                   ir=rindex(nra,nrb,nrc,nrd,nla,nlb,nlc,nld,n)
                   If(index_flag(ir) == 0) Then
                      nr_flag(1,ilauf) = nra
                      nr_flag(2,ilauf) = nrb
                      nr_flag(3,ilauf) = nrc
                      nr_flag(4,ilauf) = nrd
                      nl_flag(1,ilauf) = nla
                      nl_flag(2,ilauf) = nlb
                      nl_flag(3,ilauf) = nlc
                      nl_flag(4,ilauf) = nld
                      ir_flag(ilauf) = ir
                      index_flag(ir) = 1
                      ilauf = ilauf + 1
                   End If
                   If(nlb /= 0 .Or. nld /= 0) Then
                      ir=rindex(nra,nrb,nrc,nrd,nla,-nlb,nlc,-nld,n)
                      If(index_flag(ir) == 0) Then
                         nr_flag(1,ilauf)= nra
                         nr_flag(2,ilauf)= nrb
                         nr_flag(3,ilauf)= nrc
                         nr_flag(4,ilauf)= nrd
                         nl_flag(1,ilauf)= nla
                         nl_flag(2,ilauf)=-nlb
                         nl_flag(3,ilauf)= nlc
                         nl_flag(4,ilauf)=-nld
                         ir_flag(ilauf)=ir
                         index_flag(ir)=1
                         ilauf = ilauf + 1
                      End If
                   End If
                End If
             End Do !itd
          End Do !itb
       End Do !itc
    End Do !ita
    End If

    ! The actual calculation of the matrix elements is done below using
    ! the flags that were set on the previous loop

    !$OMP Parallel Default(None) &
    !$OMP& SHARED(ilauf,nr_flag,nl_flag,ir_flag,VrGogny,n_g_all) &
    !$OMP& PRIVATE(ii,nra,nrb,nrc,nrd,nla,nlb,nlc,nld,ir_abcd,ig)
    !$OMP DO SCHEDULE(DYNAMIC)
    Do ii = 0,ilauf-1
       nra = nr_flag(1,ii)
       nrb = nr_flag(2,ii)
       nrc = nr_flag(3,ii)
       nrd = nr_flag(4,ii)
       nla = nl_flag(1,ii)
       nlb = nl_flag(2,ii)
       nlc = nl_flag(3,ii)
       nld = nl_flag(4,ii)
       ir_abcd = ir_flag(ii)
       Call radial_matrix_elements(nra,nla,nrb,nlb,nrc,nlc,nrd,nld)
       Do ig = 1,n_g_all
          VrGogny(ig,ir_abcd) = Vr_ig(ig)
       End Do
    End Do
    !$OMP End Do
    !$OMP End Parallel
#else
    If(Allocated(Vr_Gogny)) Deallocate(Vr_Gogny)
    Allocate(Vr_Gogny(1:n_g_all,0:n/2,-n:n,0:n/2,-n:n,0:n/2,-n:n,0:n/2,-n:n))
    Vr_gogny = zero
    Do nra = 0,n/2
       Do nrb = 0,n/2
          Do nrc = 0,n/2
             Do nrd = 0,n/2
                Do nla = -n+2*nra,n-2*nra
                   Do nlb = -n+2*nrb,n-2*nrb
                      Do nlc = -n+2*nrc,n-2*nrc
                         Do nld = -n+2*nrd,n-2*nrd
                            Call radial_matrix_elements(nra,nla,nrb,nlb,nrc,nlc,nrd,nld)
                            Do ig = 1,n_g_all
                               Vr_gogny(ig,nra,nla,nrb,nlb,nrc,nlc,nrd,nld) = Vr_ig(ig)
                            End Do
                         End Do
                      End Do
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
#endif

  End Subroutine CalculateVrGogny
  !======================================================================
  !> This routine calculates all the non-zero \f$ T_{n_1 n_2}^n \f$ coefficients. These coefficients
  !> expand a product of two one-dimensional harmonic oscillator wave functions (HOWF) into a linear
  !> combination of single one-dimensional HOWF's.
  !>
  !> In particular, these coefficients are calculated as
  !>    \f[
  !>        T_{n_1 n_2}^n \equiv \frac{\sqrt{n_1!n_2!n!}}
  !>          {\displaystyle\left(\frac{n-n+n}{2}\right)! \left(\frac{n-n+n}{2}\right)! \left(\frac{n-n+n}{2}\right)!}.
  !>    \f]
  !> The coefficients are stored in the private array Tz(n1,n2,n); see \cite younes2009gaussian for a
  !> derivation of these coefficients.
  !======================================================================
  Subroutine calculateTz()
    Use HFBTHO, Only : nrx,nlx,nzx

    Integer(ipr) :: n1,n2,n,nz
    Real(pr) :: n1f,n2f,nf,m1f,m2f,m3f

#if(GOGNY_HYPER==1)
    nz = Max(nzx,Max(2*nrx,nlx))
#else
    nz = nzx
#endif
    If(Allocated(T_z)) Deallocate(T_z)
    Allocate(T_z(0:nz,0:nz,0:2*nz))
    T_z = zero
    !$OMP Parallel Default(None) &
    !$OMP& SHARED(nz,T_z) &
    !$OMP& PRIVATE(n1,n2,n,n1f,n2f,nf,m1f,m2f,m3f)
    !$OMP DO SCHEDULE(DYNAMIC)
    Do n1 = 0, nz
       Do n2 = 0, nz
          Do n = 0, 2*nz
             If(n2 < n1) Cycle
             If(n < Abs(n2-n1) .Or. n > n2+n1) Cycle
             If(Mod(n,2) /= Mod(n2+n1,2))      Cycle
             n1f = factrl(n1)
             n2f = factrl(n2)
             nf = factrl(n)
             m1f = factrl((-n1+n2+n)/2)
             m2f = factrl((+n1-n2+n)/2)
             m3f = factrl((+n1+n2-n)/2)
             T_z(n1,n2,n)=(Sqrt(n2f)/m1f)*(Sqrt(nf)/m2f)*(Sqrt(n1f)/m3f)
             T_z(n2,n1,n)=T_z(n1,n2,n)
          End Do
       End Do
    End Do
    !$OMP End Do
    !$OMP End Parallel

  End Subroutine calculateTz
  !======================================================================
  !> This routine calculates all the nonzero \f$ T_{n_1,k_1,n_2,k_2}^{n,k1+k2}\f$ coefficients to transform
  !> a product of two two-dimensional harmonic oscillator wave functions in radial coordinates into a linear
  !> combination of single two dimensional HO wave functions in radial coordinates.
  !>
  !> In particular, these coefficients are calculated as
  !>   \f[
  !>      T_{n_1,k_1;n_2,k_2}^{n,k_1+k_2} = (-1)^{n_1+n_2-n}
  !>      \sqrt{\frac{n!(n_1+|k_1|)!(n_2+|k_2|)!}{n_1!n_2!(n+|k_1+k_2)!}}
  !>      \sum_{m_1=0}^{n_1} \sum_{m_2=0}^{n_2} (-1)^{m_1+m_2}
  !>      {{n_1}\choose{m_1}} {{n_2}\choose{m_2}}
  !>      {{n_{1,2}-m_1-m_2}\choose{n}}\frac{(n_{1,2}+|k_1+k_2|-m_1-m_2)!}
  !>      {(n_1+|k_1|-m_1)!(n_2+|k_2|-m_2)!},
  !>   \f]
  !> where \f$ n_{1,2} = n_1+n_2+\frac{|k_1|+|k_2|-|k_1+k_2|}{2}\f$. The coefficients are stored in the
  !> private array Tr(n1,k1,n2,k2,n); see \cite younes2009gaussian for a derivation of these coefficients.
  !======================================================================
  Subroutine calculateTr()
    Use HFBTHO, Only : nrx,nlx

    Integer(ipr) :: n1,n2,n,k1,k2,nm12,nmax,nm2
    Real(pr) :: Tpp,Tmp

    nmax = Max(2*nrx,nlx)
    nm2 = nmax/2
    If(Allocated(T_r)) Deallocate(T_r)
    Allocate(T_r(0:nm2,-nmax:nmax,0:nm2,-nmax:nmax,0:nmax))
    Do n1 = 0, nm2
       Do k1 = 0,nmax-2*n1
          Do n2 = n1, nm2
             Do k2 = 0,nmax-2*n2
                nm12 = n1 + n2 + (k1+k2-Abs(-k1+k2))/2
                Do n = 0, nm12
                   Tpp = TrCoefficient(n1, k1,n2,k2,n, k1+k2)
                   Tmp = TrCoefficient(n1,-k1,n2,k2,n,-k1+k2)
                   T_r(n1, k1,n2, k2,n) = Tpp
                   T_r(n2, k2,n1, k1,n) = Tpp
                   T_r(n1,-k1,n2,-k2,n) = Tpp
                   T_r(n2,-k2,n1,-k1,n) = Tpp
                   T_r(n1,-k1,n2, k2,n) = Tmp
                   T_r(n2, k2,n1,-k1,n) = Tmp
                   T_r(n1, k1,n2,-k2,n) = Tmp
                   T_r(n2,-k2,n1, k1,n) = Tmp
                End Do
             End Do
          End Do
       End Do
    End Do

  End Subroutine calculateTr
  !======================================================================
  !> This function calculates the \f$ T_{n_1,k_1,n_2,k_2}^{n,k} \f$ coefficient. Since
  !> \f$ T_{n_1,k_1,n_2,k_2}^{n,k} = T_{n,k,n_2,k_2}^{n_1,k_1} = T_{n_1,k_1,n,k}^{n_2,k_2} \f$, the
  !> combination with the least number of operations is calculated in order to reduce loss of
  !> accuracy. The routine \ref trsumterms() is used to calculate the number of operations for
  !> each combination; see Eqs. (C.9) and (C.10) of \cite younes2009gaussian for details.
  !>
  !> @result \f$ \displaystyle T_{n_1,k_1,n_2,k_2}^{n,k} = (-1)^{n_1+n_2+n}
  !>          \sqrt{\frac{n!(n_1+|k_1|)!(n_2+|k_2|)!}{n_1!n_2!(n+|k|)!}}
  !>          \sum_{m_1=0}^{n_1} \sum_{m_2=0}^{n_2} (-1)^{m_1+m_2}
  !>         {{n_1}\choose{m_1}} {{n_2}\choose{m_2}} {{n_{1,2}-m_1-m_2}\choose{n}}
  !>         \frac{(n_{1,2}+|k|-m_1-m_2)!}{(n_1+|k_1|-m_1)!(n_2+|k_2|-m_2)!}\f$, where
  !>         \f$ n_{1,2} = n_1+n_2+\frac{|k_1|+|k_2|-|k|}{2} \f$.
  !======================================================================
  Function TrCoefficient(n1,k1,n2,k2,n,k) result(Tr)
    Integer(ipr), Intent(In) :: n1 !< - first  radial principal quantum number
    Integer(ipr), Intent(In) :: n2 !< - second radial principal quantum number
    Integer(ipr), Intent(In) :: n  !< - third  radial principal quantum number
    Integer(ipr), Intent(In) :: k1 !< - first  radial orbital quantum number
    Integer(ipr), Intent(In) :: k2 !< - second radial orbital quantum number
    Integer(ipr), Intent(In) :: k  !< - third  radial orbital quantum number
    Real(pr) :: Tr

    Integer(ipr), Dimension(3) :: Nsums
    Integer(ipr) :: minp,in1,in2,in,ik1,ik2,ik
    Integer(ipr) :: n12,m1,m2
    Real(pr) :: n1f,nk1f,n2f,nk2f,nf,nkf,fac,Ti

    If(Mod(k1+k2+k,2) /= 0) Stop 'wrong k values in TrCoefficient'

    Nsums(1) = TrSumTerms(n ,k ,n2,k2,n1,k1)
    Nsums(2) = TrSumTerms(n1,k1,n ,k ,n2,k2)
    NSums(3) = TrSumTerms(n1,k1,n2,k2,n ,k )
    minp = minloc(Nsums,1)

    Select Case(minp)
    Case(1)
       in1 = n
       ik1 = k
       in2 = n2
       ik2 = k2
       in  = n1
       ik  = k1
    Case(2)
       in1 = n1
       ik1 = k1
       in2 = n
       ik2 = k
       in  = n2
       ik  = k2
    Case(3)
       in1 = n1
       ik1 = k1
       in2 = n2
       ik2 = k2
       in  = n
       ik  = k
    End Select

    n12  = in1 + in2 + (Abs(ik1)+Abs(ik2)-Abs(ik))/2
    n1f  = factrl(in1)
    nk1f = factrl(in1 + Abs(ik1))
    n2f  = factrl(in2)
    nk2f = factrl(in2 + Abs(ik2))
    nf   = factrl(in)
    nkf  = factrl(in+Abs(ik))

    Tr = 0.0_pr
    Do m1 = 0,Min(in1,n12-in)
       Do m2 = 0,Min(in2,n12-in-m1)
          Ti = (-1)**(m1+m2)*binomialco(in1,m1)*binomialco(in2,m2)*binomialco(n12-m1-m2,in)&
                            *(factrl(n12+Abs(ik)-m1-m2)/factrl(in1+Abs(ik1)-m1))/factrl(in2+Abs(ik2)-m2)
          Tr = Tr + Ti
       End Do
    End Do
    Fac = (-1)**(in1+in2-in)*Sqrt(nf/nkf)*Sqrt(nk1f/n1f)*Sqrt(nk2f/n2f)
    Tr = Fac*Tr

  End Function TrCoefficient
  !======================================================================
  !> This function calculates the number of summations required to calculate the
  !> \f$ T_{n_1,k_1,n_2,k_2}^{n,k} \f$  coefficient.
  !======================================================================
  Function TrSumTerms(n1,k1,n2,k2,n,k) result(NT)
    Integer(ipr), Intent(In) :: n1 !< - first  radial principal quantum number
    Integer(ipr), Intent(In) :: n2 !< - second radial principal quantum number
    Integer(ipr), Intent(In) :: n  !< - third  radial principal quantum number
    Integer(ipr), Intent(In) :: k1 !< - first  radial orbital quantum number
    Integer(ipr), Intent(In) :: k2 !< - second radial orbital quantum number
    Integer(ipr), Intent(In) :: k  !< - third  radial orbital quantum number
    Integer(ipr) :: NT

    Integer(ipr) :: n12,mu,ab_min,ab_max

    n12 = n1 + n2 + (Abs(k1)+Abs(k2)-Abs(k))/2
    mu = n12-n
    If(mu < 0) Then
       NT = 0
       Return
    End If
    ab_min = Min(n1,n2)
    ab_max = Max(n1,n2)
    If(mu <= ab_min) Then
       NT = (mu+1)*(mu+2)/2
    Else If(mu <= ab_max) Then
       NT = (mu+1)*(mu+2)/2 - (mu-ab_min)*(mu-ab_min+1)/2
    Else
       If(mu <= n1+n2) Then
          NT = (mu+1)*(mu+2)/2 - (mu-ab_min)*(mu-ab_min+1)/2 - (mu-ab_max)*(mu-ab_max+1)/2
       Else
          NT = (n1+1)*(n2+1)
       End If
    End If

  End Function TrSumTerms
  !======================================================================
  !> This routine calculates all the necessary \f$C_{n_x n_y}^{n k}i^{n_y}\f$ coefficients with
  !> \f$n_x+n_y=2n+|k| \f$ to transform a two-dimensional harmonic oscillator wave function (HOWF)
  !> in radial coordinates into a linear combination of products of one-dimensional HOWF in
  !> cartesian coordinates. The coefficients are calculated with
  !>    \f[
  !>        C_{n_x n_y}^{n k}i^{n_y} = \frac{2^{-n-|k|/2}(-1)^n
  !>             \sqrt{n!(n+|k|)!}}{\sqrt{(2n+|k|-n_y)!n_y!}}
  !>             \sum_{q=0}^{\min(n_y,n+(|k|-k)/2)}
  !>             {{2n+|k|-n_y}\choose{n-q+\frac{|k|-k}{2}}}
  !>             {{n_y}\choose{q}} (-1)^{n_y-q}
  !>    \f]
  !> The \f$ i^{n_y}\f$ factor makes the result a real number. These transformation coefficients
  !> are only calculated when the preprocessor directive GOGNY_HYPER is set to 1. The coefficients
  !> are stored in the private array Cp2c(nr,k,ny); see \cite younes2009gaussian for a derivation
  !> and additional details.
  !======================================================================
  Subroutine calculateCpolar2cartesian()
    Use HFBTHO, Only : nrx,nlx

    Integer(ipr) :: n,k,ny,nx,q,qmax,nshells
    Real(pr) :: A,xsum

    nshells = Max(2*nrx,nlx)
    If(Allocated(Cp2c)) Deallocate(Cp2c)
    Allocate(cp2c(0:Nshells,-2*Nshells:2*Nshells,0:2*Nshells))
    cp2c = zero
    !$OMP Parallel Default(None) &
    !$OMP& SHARED(Nshells,Cp2c) &
    !$OMP& PRIVATE(n,k,ny,nx,A,xsum,qmax,q)
    !$OMP DO SCHEDULE(DYNAMIC)
    Do n = 0,Nshells
       Do k = -2*nshells,2*nshells
          Do ny = 0,2*nshells
             If(Abs(k) > 2*nshells-2*n) Cycle
             If(ny > 2*n+Abs(k))        Cycle
             nx = 2*n+Abs(k)-ny
             A = (-1)**n*(2**(-n-Abs(k)*0.5_pr))*Sqrt(factrl(n+Abs(k)))*(Sqrt(factrl(n))/(Sqrt(factrl(nx))*Sqrt(factrl(ny))))
             xsum = zero
             qmax = Min(ny,n+(Abs(k)-k)/2)
             Do q = 0,qmax
                xsum = xsum + BinomialCo(nx,n-q+(Abs(k)-k)/2)*BinomialCo(ny,q)*(-1)**(ny-q)
             End Do
             Cp2c(n,k,ny) = A*xsum
          End Do
       End Do
    End Do
    !$OMP End Do
    !$OMP End Parallel

  End Subroutine calculateCpolar2cartesian
  !======================================================================
  !> This routine calculates all the non-zero axial two-body matrix elements
  !> \f$ \langle n_1 n_2|V|n_3 n_4 \rangle \f$ using the radial oscillator parameter \f$b_{\perp}\f$.
  !>
  !> These matrix elements are only calculated when the preprocessor directive GOGNY_HYPER is set to 1.
  !> The matrix elements are stored in the private array ME1D(n1,n2,n3,n4).
  !======================================================================
  Subroutine calculateME1D()
    Use UNEDF, Only : n_g_all,mu_g_all
    Use HFBTHO, Only : nrx,nlx,bp

    Integer(ipr) :: ni,nj,nk,nl,n,ig
    Real(pr) :: ME

    n = Max(2*nrx,nlx)
    If(Allocated(ME1D)) Deallocate(ME1D)
    Allocate(ME1D(1:n_g_all,0:n,0:n,0:n,0:n))
    ME1D = zero
    !$OMP Parallel Default(None) &
    !$OMP& SHARED(n,ME1D,mu_g_all,bp,n_g_all) &
    !$OMP& PRIVATE(ni,nj,nk,nl,ig,ME)
    !$OMP DO SCHEDULE(DYNAMIC)
    Do ni = 0,n
       Do nj = 0,n
          Do nk = 0,n
             Do nl = 0,n
                If(nk < ni) Cycle
                If(nl < nj) Cycle
                If(Mod(ni+nj+nk+nl,2) /= 0) Cycle
                Do ig = 1,n_g_all
                   ME = MatrixElement_z(ni,nj,nk,nl,mu_g_all(ig),bp)
                   ME1D(ig,ni,nj,nk,nl) = ME
                   ME1D(ig,nk,nj,ni,nl) = ME
                   ME1D(ig,ni,nl,nk,nj) = ME
                   ME1D(ig,nk,nl,ni,nj) = ME
                End Do
             End Do
          End Do
       End Do
    End Do
    !$OMP End Do
    !$OMP End Parallel

  End Subroutine calculateME1D
  !======================================================================
  !> This function calculates the axial component of the two-body potential matrix element for a
  !> Gaussian potential with axial symmetry. If the preproccesor directive GOGNY_HYPER is set
  !> to 1 the matrix element is calculated using the modified Gogny expansion that uses a
  !> hypergeometric function and preserves accuracy with a large basis size; see Eq. (D.12) in
  !> \cite younes2009gaussian.
  !>    \f[
  !>       \langle n_i n_j| \hat{V}_z|n_k n_l \rangle
  !>       = \frac{\mu}{\sqrt{2\pi^3}b_z}\sum_{n=|n_j-n_l|,2}^{n_j+n_l}
  !>         T_{n_j n_l}^{n}\bar{F}_{n_i n_k}^{n},
  !>    \f]
  !> where
  !>    \f[
  !>       \bar{F}_{n_i n_k}^{n} = \frac{\Gamma(\xi-n_i)\Gamma(\xi-n_k)
  !>         \Gamma(\xi-n)}{(1+\mu^2/(2b_z^2))^\xi\sqrt{n!n_i!n_k!}}
  !>         \ _2F_1(-n_i,n_k;-\xi+n+1;-\mu^2/(2b_z^2))
  !>    \f]
  !> with \f$ \xi = \frac{n_i+n_k+n+1}{2} \f$.
  !>
  !> For any other value of GOGNY_HYPER the matrix element is calculated using the direct Gogny
  !> transformation which looses numerical accuracy with a large basis size; see Eq. (10) in
  !> \cite younes2009gaussian.
  !>    \f[
  !>       \langle n_i,n_j|V_z|n_k,n_l \rangle = \sqrt{\frac{G_z-1}{G_z+1}}
  !>       \sum_{m=|n_i-n_k|,2}^{n_i+n_k}\sum_{n=|n_j-n_l|,2}^{n_j+n_l}
  !>       T_{n_i,n_k}^{m} T_{n_j,n_l}^{n} \bar{I}(m,n)
  !>    \f]
  !======================================================================
  Function MatrixElement_z(ni,nj,nk,nl,mu,b) result(Vz)
    Integer(ipr), Intent(In) :: ni !< - first  z-component quantum number
    Integer(ipr), Intent(In) :: nj !< - second z-component quantum number
    Integer(ipr), Intent(In) :: nk !< - third  z-component quantum number
    Integer(ipr), Intent(In) :: nl !< - fourth z-component quantum number
    Real(pr), Intent(In) :: mu,b
    Real(pr) :: Vz

#if(GOGNY_HYPER==1)
    Integer(ipr) :: nz
    Integer(ipr), Dimension(2) :: Nsums
    Integer(ipr) :: minp,ini,inj,ink,inl
    Real(pr) :: Fbar,xi,z
#else
    Integer(ipr) :: mz,nz
    Real(pr) :: Gz
#endif

    Vz = 0.0_pr

#if(GOGNY_HYPER==1)
    Nsums(1) = (nj+nl-Abs(nj-nl))/2 + 1
    Nsums(2) = (ni+nk-Abs(ni-nk))/2 + 1
    If(nsums(1) == nsums(2)) Then
       nsums(1) = Min(ni,nk)
       nsums(2) = Min(nj,nl)
    End If
    minp = minloc(Nsums,1)

    Select Case(minp)
    Case(1)
       ini = ni
       inj = nj
       ink = nk
       inl = nl
    Case(2)
       ini = nj
       inj = ni
       ink = nl
       inl = nk
    End Select

    z = one+mu**2/(two*b**2)
    Vz = zero
    Do nz = Abs(inj-inl),inj+inl,2
       If(Mod(ini+ink+nz,2) /= 0) exit
       xi = (ini+ink+nz+1)*0.5_pr
       Fbar = gamma(xi-ini)*gamma(xi-ink)*gamma(xi-nz)/(z**xi*Sqrt(factrl(ini)*factrl(ink)*factrl(nz))) &
             *HyperGeom2F1(-ini,-ink,nz+1-xi,1-z)
       Vz = Vz + T_z(inj,inl,nz)*Fbar
    End Do
    Vz = mu/(Sqrt(two*pi**3)*b)*Vz
#else
    Gz = one + (mu/b)**2
    Do mz = Abs(ni-nk),ni+nk,2
       Do nz = Abs(nj-nl),nj+nl,2
          Vz = Vz + T_z(ni,nk,mz)*T_z(nj,nl,nz)*Ibarz(mz,nz,Gz)
       End Do
    End Do
    Vz = Sqrt((Gz-1)/(Gz+1))*Vz
#endif

  End Function MatrixElement_z
  !======================================================================
  !> This function calculates the \f$ \bar{I}(m,n) \f$ coefficient necessary to calculate the axial
  !> component of the two-body potential matrix element for a Gaussian potential with axial symmetry;
  !> see \cite younes2009gaussian.
  !>
  !> @result \f$ \displaystyle \bar{I}(m,n) = \sqrt{\frac{m!n!}{2^{m+n}}}
  !>                            \frac{(-1)^{(m-n)/2}}{\left(\frac{m+n}{2}\right)!(1+G_z)^{(m+n)/2}}
  !>                            {{m+n}\choose{n}}\f$
  !======================================================================
  Function Ibarz(m,n,Gz) result(Iz)
    Integer(ipr), Intent(In) :: m
    Integer(ipr), Intent(In) :: n
    Real(pr), Intent(In) :: Gz
    Real(pr) :: Iz

    If(Mod(m+n,2) /= 0) Then
       Iz = 0.0_pr
       Return
    End If
    Iz = (-1)**((m-n)/2)*(Sqrt(factrl(m))/factrl((m+n)/2))*Sqrt(factrl(n)/(2+2*Gz)**(m+n))*binomialco(m+n,n)

  End Function Ibarz
  !======================================================================
  !> This routine calculates the radial component of the two-body potential matrix element for a given
  !> set of radial quantum numbers. The matrix elements are calculated for all the Gaussians in the
  !> finite-range interaction since each Gaussian has its own range parameter in the mu_g_all(:) array.
  !>
  !> If the preprocessor variable GOGNY_HYPER is set to 1 the matrix elements are obtained by
  !> transforming the HO wavefunctions from radial into cartesian coordinates and separating the
  !> integration into a product of two one-dimensional two body matrix elements. The transformation
  !> is given by
  !>   \f{multline*}{ \langle n_{r_i} \Lambda_i,n_{r_j}\Lambda_j|\hat{V}_p|
  !>       n_{r_k} \Lambda_k,n_{r_l}\Lambda_l \rangle =
  !>       \sum_{n_{y_i}=0}^{2n_{r_i}+|\Lambda_i|}
  !>       \sum_{n_{y_j}=0}^{2n_{r_j}+|\Lambda_j|}
  !>       \sum_{n_{y_k}=0}^{2n_{r_k}+|\Lambda_k|}
  !>       \sum_{n_{y_l}=0}^{2n_{r_l}+|\Lambda_l|}
  !>       C_{n_{x_i} n_{y_i}}^{n_{r_i}\Lambda_i *}
  !>       C_{n_{x_j} n_{y_j}}^{n_{r_j}\Lambda_j *}
  !>       C_{n_{x_k}n_{y_k}}^{n_{r_k}\Lambda_k}
  !>       C_{n_{x_l}n_{y_l}}^{n_{r_l}\Lambda_l} \\
  !>       \times\langle n_{x_i}n_{x_j}|\hat{V}_{\rm 1D}|n_{x_k}n_{x_l}\rangle
  !>       \langle n_{y_i}n_{y_j}|\hat{V}_{\rm 1D}|n_{y_k}n_{y_l}\rangle
  !>   \f}
  !> The one-dimensional matrix elements have already been calculated by the \ref calculateme1d()
  !> subroutine and stored in the ME1D array. Since \ref calculateme1d() uses the modified Gogny
  !> transformation appropriate for large basis size, this transformation is also used when preserving
  !> accuracy with large basis.
  !>
  !> For any other value of GOGNY_HYPER the function \ref matrixelement_r(), which uses the direct
  !> Gogny transformation, is used.
  !======================================================================
  Subroutine radial_matrix_elements(ni,li,nj,lj,nk,lk,nl,ll)
    Use UNEDF, Only : n_g_all

    Integer(ipr), Intent(In) :: ni !< - first  radial principal quantum number
    Integer(ipr), Intent(In) :: nj !< - second radial principal quantum number
    Integer(ipr), Intent(In) :: nk !< - third  radial principal quantum number
    Integer(ipr), Intent(In) :: nl !< - fourth radial principal quantum number
    Integer(ipr), Intent(In) :: li !< - first  radial orbital quantum number
    Integer(ipr), Intent(In) :: lj !< - second radial orbital quantum number
    Integer(ipr), Intent(In) :: lk !< - third  radial orbital quantum number
    Integer(ipr), Intent(In) :: ll !< - fourth radial orbital quantum number

    Real(pr) :: Ci,Cj,Ck,Cl
    Integer(ipr) :: nyi,nyj,nyk,nyl,nxi,nxj,nxk,nxl,ig

    If(Allocated(Vr_ig)) Deallocate(Vr_ig)
    Allocate(Vr_ig(1:n_g_all))
    Vr_ig = zero
    If(-li-lj+lk+ll /= 0) Return
#if(GOGNY_HYPER==1)
    Do nyi = 0,2*ni+Abs(li)
       nxi = 2*ni+Abs(li) - nyi
       Ci = Cp2c(ni,li,nyi)
       Do nyj = 0,2*nj+Abs(lj)
          nxj = 2*nj+Abs(lj)-nyj
          Cj = Cp2c(nj,lj,nyj)
          Do nyk = 0,2*nk+Abs(lk)
             nxk = 2*nk+Abs(lk) - nyk
             Ck = Cp2c(nk,lk,nyk)
             Do nyl = Mod(nyi+nyj+nyk,2),2*nl+Abs(ll),2
                nxl = 2*nl+Abs(ll)-nyl
                Cl = Cp2c(nl,ll,nyl)
                Do ig = 1,n_g_all
                   Vr_ig(ig) = Vr_ig(ig) + (-1)**(nyi+nyj+(nyi+nyj+nyk+nyl)/2)*Ci*Cj*Ck*Cl &
                                         *ME1D(ig,nxi,nxj,nxk,nxl)*ME1D(ig,nyi,nyj,nyk,nyl)
                End Do
             End Do
          End Do
       End Do
    End Do
#else
    Do ig = 1,n_g_all
       Vr_ig(ig) = MatrixElement_r(ni,li,nj,lj,nk,lk,nl,ll,mu_g_all(ig),bp)
    End Do
#endif

  End Subroutine radial_matrix_elements
  !======================================================================
  !> This function calculates the radial component of the twp-body potential matrix element for a
  !> Gaussian potential with axial symmetry using the direct Gogny transformation. This transformation
  !> looses accuracy as the basis size increases and should not be used in production runs; see
  !> \cite younes2009gaussian
  !>
  !> @result \f[ \displaystyle \langle n_{r_i},l_i,n_{r_j},l_j|V_p|n_{r_k},l_k,n_{r_l},l_l \rangle =
  !>         \frac{G_p-1}{G_p+1} \sum_{n_r=0}^{n_{\bar{j},l}}
  !>         \sum_{n=0}^{n_{\bar{i},k}}
  !>         T_{n_{r_i},-l_i,n_{r_k},l_k}^{n,-l_i+l_k}
  !>         T_{n_{r_j},-l_j,n_{r_l},l_l}^{n_r,-l_j+l_l}
  !>         \bar{I}(n_r,-l_j+l_l,n,-l_i+l_k)
  !>         \delta_{-l_j+l_l,l_i+l_k} \f]
  !======================================================================
  Function MatrixElement_r(ni,li,nj,lj,nk,lk,nl,ll,mu,b) result(Vr)
    Integer(ipr), Intent(In) :: ni!< - first  radial principal quantum number
    Integer(ipr), Intent(In) :: nj!< - second radial principal quantum number
    Integer(ipr), Intent(In) :: nk!< - third  radial principal quantum number
    Integer(ipr), Intent(In) :: nl!< - fourth radial principal quantum number
    Integer(ipr), Intent(In) :: li!< - first  radial orbital quantum number
    Integer(ipr), Intent(In) :: lj!< - second radial orbital quantum number
    Integer(ipr), Intent(In) :: lk!< - third  radial orbital quantum number
    Integer(ipr), Intent(In) :: ll!< - fourth radial orbital quantum number
    Real(pr), Intent(In) :: mu!< - range of the Gaussian potential
    Real(pr), Intent(In) :: b!< - Harmonic Oscillator length
    Real(pr) :: Vr

    Real(pr) :: Gp
    Integer(ipr) :: mr,nr,nmjl,nmik

    Vr = 0.0_pr
    If(-li-lj+lk+ll /= 0) Return
    Gp = one + (mu/b)**2
    nmjl = nj+nl+(Abs(-lj)+Abs(ll)-Abs(-lj+ll))/2
    nmik = ni+nk+(Abs(-li)+Abs(lk)-Abs(-li+lk))/2
    Do mr = 0,nmik
       Do nr = 0,nmjl
          Vr = Vr + T_r(ni,-li,nk,lk,mr)*T_r(nj,-lj,nl,ll,nr)*Ibarr(nr,-lj+ll,mr,-li+lk,Gp)
       End Do
    End Do
    Vr = (Gp-1)/(Gp+1)*Vr

  End Function MatrixElement_r
  !======================================================================
  !> This function calculates the \f$ \bar{I}(n_1,k_1,n_2,k_2) \f$ coefficient necessary to calculate
  !> the radial component of the two-body potential matrix element component for a Gaussian potential
  !> with axial symmetry; see \cite younes2009gaussian.
  !>
  !> @result \f$ \bar{I}(n_1,k_1,n_2,k_2) = \sqrt{\frac{(n_1+|k_1|)!n_2!}
  !>         {n_1!(n_2+|k_2|)!}} {{n_1+n_2+|k_1|}\choose{n_2}}
  !>         \frac{\delta_{k_1,-k_2}}{(G_p+1)^{n_1+n_2+|k_1|}}\f$
  !======================================================================
  Function Ibarr(n1,k1,n2,k2,Gp) result(Ir)
    Integer(ipr), Intent(In):: n1 !< - first  radial principal quantum number
    Integer(ipr), Intent(In):: n2 !< - second radial principal quantum number
    Integer(ipr), Intent(In):: k1 !< - first  radial orbital quantum number
    Integer(ipr), Intent(In):: k2 !< - second radial orbital quantum number
    Real(pr), Intent(In) :: Gp
    Real(pr) :: Ir
    Integer(ipr) :: k
    Real(pr) :: n1f,n1kf,n2f,n2kf

    If(k1 /= -k2) Then
       Ir = 0.0_pr
       Return
    End If
    k = Abs(k1)
    n1f = factrl(n1)
    n2f = factrl(n2)
    n1kf = factrl(n1+k)
    n2kf = factrl(n2+k)
    Ir = Sqrt(n1kf/n1f)*Sqrt(n2f/n2kf)*(binomialco(n1+n2+k,n2)/((Gp+1)**(n1+n2+k)))

  End Function Ibarr
  !======================================================================
  !> This function recursively calculates and stores in an array (to avoid recalculating in future
  !> calls) \f$ 0!, 1!, 2!, \ldots, n! \f$. If \f$n > 170\f$, \f$ n! = \Gamma(n+1) \f$ is used
  !> although it will overflow on most machines. Note that this routine duplicates \ref hfbtho_math.gfv()
  !>
  !> @result  \f$ n! \f$
  !======================================================================
  Function factrl(n) result(fact)
    Integer(ipr), Intent(In) :: n !< - An integer
    Real(pr) :: fact

    Integer(ipr), save :: ntop = 0
    Integer(ipr), parameter :: nmax = 170
    Integer(ipr) :: i
    Real(pr), Dimension(0:nmax), save :: a=0._pr

    If(n < 0) Then
       Write(*,*) 'negative integer in factrl'
       Stop
    End If
    If(ntop == 0) a(0)=one
    If(n <= ntop) Then
       fact = a(n)
    Else If(n <= nmax) Then
       Do i = ntop+1,n
          a(i) = Real(i,Kind=pr)*a(i-1)
       End Do
       ntop = n
       fact = a(n)
    Else
       fact = Exp(log_gamma(n+1._pr))
    End If

  End Function factrl
  !======================================================================
  !> This function calculates the Hypergeometric function \f$_2F_1(a,b,c;x) \f$ in the very particular
  !> case where \f$ a \leq 0\f$, \f$ b \leq 0\f$ and c is NOT a negative integer such that
  !> \f$ |c| \leq \min(|a|,|b|)\f$
  !>
  !> @result
  !>     \f[
  !>         _2F_1(a,b,c;x) = \sum_{i=0}^d {{d}\choose{i}}{{e}
  !>          \choose{i}} \frac{i!}{(c)_i}x^i
  !>     \f]
  !> where \f$ d = \min(|a|,|b|)\f$, \f$ e = \max(|a|,|b|)\f$ and \f$ (x)_i \f$ is the rising
  !> Pochhammer Symbol (also known as upper factorial)
  !======================================================================
  Function HyperGeom2F1(a,b,c,x) result(HG2F1)
    Integer(ipr), Intent(In) :: a !< - First parameter, has to be negative
    Integer(ipr), Intent(In) :: b !< - Second parameter, has to be negative
    Real(pr), Intent(In) :: c !< - Third parameter, cannot be a negative integer such that \f$ |c| \leq \min(|a|,|b|)\f$
    Real(pr), Intent(In) :: x !< - Value where the function is evaluated
    Real(pr) :: HG2F1

    Integer(ipr) :: d,e,i
    Real(pr) :: HGi

    If(a > 0 .Or. b > 0) Then
       Write(*,*) 'a or b is not a negative integer in HyperGeom2F1'
       Stop
    End If
    d = Abs(Max(a,b))
    e = Abs(Min(a,b))
    If(c < 0.0_pr .And. Abs(c-nint(c)) < 1.0e-14_pr .And. Abs(c) <= d) Then
       Write(*,*) 'c is a negative integer greater or equal than Max(a,b)in HyperGeom2F1'
       Stop
    End If
    HG2F1 = zero
    Do i = 0,d
       HGi = + binomialco(d,i)*binomialco(e,i)*x**i*factrl(i)/upperfactrl(c,i)
       HG2F1 = HG2F1 + HGi
    End Do

  End Function HyperGeom2F1
  !======================================================================
  !> This function calculates the binomial coefficient \f$ {{m}\choose{n}} \f$. Overflows are avoided
  !> by changing the factorials for logarithms of the \f$ \Gamma \f$ function and taking the exponential
  !>
  !> @result  \f$ \displaystyle {{m}\choose{n}}=\frac{m!}{n!(m-n)!)}\f$
  !======================================================================
  Function binomialco(m,n) result(bc)
    Integer(ipr), Intent(In):: m !< - A positive integer
    Integer(ipr), Intent(In):: n !< - An integer
    Real(pr) :: bc

    If(m < 0) Then
       Write(*,*) 'negative integer in binomialco'
       Stop
    End If
    If(n < 0.or.n > m) Then
       bc = zero
       Return
    End If
    If(n == 0 .Or. n == m) Then
       bc = one
       Return
    End If
    If(n == 1 .Or. n == m-1) Then
       bc = Real(m,Kind=pr)
       Return
    End If
    If(m <= 170) Then
       If(n <= m/2) Then
          bc = (factrl(m)/factrl(m-n))/factrl(n)
       Else
          bc = (factrl(m)/factrl(n))/factrl(m-n)
       End If
    Else
       bc = exp(log_gamma(m+1.0_pr) - log_gamma(n+1.0_pr) - log_gamma(m-n+1.0_pr))
    End If

  End Function binomialco
  !======================================================================
  !> This function calculates the upper factorial, also known as the rising Pochhamer symbol
  !> \f$ (x)_i \f$
  !> @result  \f$ (x)_i =  x(x+1)\ldots (x+i-1)\f$
  !======================================================================
  Function upperfactrl(x,i) result(upf)
    Real(pr), Intent(In) :: x !< - Value where the upper factorial is evaluated
    Integer(ipr), Intent(In) :: i !< - An integer
    Real(pr) :: upf

    Integer(ipr) :: j

    If(i < 0) Then
       Write(*,*) 'negative integer in upperfactrl'
       Stop
    End If
    upf = one
    If(i == 0) return
    Do j = 0,i-1
       upf = upf*(x+Real(j,Kind=pr))
    End Do

  End Function upperfactrl
  !======================================================================
  !> This routine calculates and stores in memory what is called the ZBlock. ZBlock is a (N x N x N x N)
  !> array that is used to obtain the  transformation from the axial component quantum numbers (\f$ n_z \f$'s)
  !> into the index of an array that contains only the the non-zero axial two body potential matrix elements.
  !> IF \f$N\f$ is the number of shells, the ZBlock is
  !>    \f{eqnarray*}{
  !>     NB & = &\sum_{i=0}^{n_{z_i}-1} \sum_{j=0}^{N} \sum_{k=i}^{N}
  !>        \sum_{\substack{l=j\\ i+j+k+l\ {\rm is\ even}}}^{N} 1
  !>       +\sum_{j=0}^{n_{z_j}-1} \sum_{k=n_{z_i}}^{N}
  !>        \sum_{\substack{l=j\\ n_{z_i}+j+k+l\ {\rm is\ even}}}^{N} 1 \\
  !>        & + & \sum_{k=n_{z_i}}^{n_{z_k}-1}
  !>              \sum_{\substack{l=n_{z_j}\\ n_{z_i}+n_{z_j}+k+
  !>                    l\ {\rm is\ even}}}^{N} 1
  !>             +\sum_{\substack{l=n_{z_j}\\ n_{z_i}+n_{z_j}+n_{z_k}
  !>                    +l\ {\rm is\ even}}}^{n_{z_l}-1} 1  + 1
  !>    \f}
  !======================================================================
  Subroutine calculate_Zblock()
    Use HFBTHO, Only : nzx

    Integer(ipr) :: sni,snj,snk,snl
    Integer(ipr) :: xni,xnj,xnk,xnl
    Integer(ipr) :: ii,jj,kk,ll,n

    n = nzx
    If(Allocated(ZBlock)) Deallocate(ZBlock)
    Allocate(ZBlock(0:n,0:n,0:n,0:n))
    Zblock = 0
    xni = 0
    Do ii = 0,n
       sni = xni
       xnj = 0
       Do jj = 0,n
          snj = xnj
          xnk = 0
          Do kk = ii,n
             snk = xnk
             xnl = 0
             Do ll = jj+Mod(ii+kk,2),n,2
                snl = xnl
                ZBlock(ii,jj,kk,ll)=sni+snj+snk+snl+1
                xni = xni + 1
                xnj = xnj + 1
                xnk = xnk + 1
                xnl = xnl + 1
             End Do
          End Do
       End Do
    End Do

  End Subroutine calculate_Zblock
  !======================================================================
  !> This routine calculates and stores in memory what is called the NBlock. NBlock is an
  !> (N/2 x N/2 x N/2 x N/2 x 2N x 2N) array that is used to obtain the transformation from the
  !> radial component quantum numbers (\f$ n \f$'s and \f$ \Lambda \f$'s) into the index of an
  !> array that contains only the non-zero two-body potential matrix elements. The Nblock is the
  !> part of that index that only depends on the \f$n\f$ quantum numbers and the first two
  !> \f$\Lambda\f$ quantum number. This part takes the most time to calculate (and is therefore
  !> stored in memory so that it is only calculated once). The other part of that index is calculated
  !> with \ref l_block(). If \f$N\f$ is the number of shells, the NBlock is
  !>   \f{eqnarray*}{
  !>    NB &=& \sum_{i=0}^{n_{ri} -1} \sum_{j=0}^{N/2} \sum_{k=0}^{N/2}
  !>           \sum_{l=0}^{N/2} \sum_{p=-N+2i}^{N-2i}\;
  !>           \sum_{q=-N+2j}^{N-2j}\;
  !>           \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !>               ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>       &+& \sum_{j=0}^{n_{rj}-1} \sum_{k=0}^{N/2}
  !>           \sum_{l=0}^{N/2} \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>           \sum_{q=-N+2j}^{N-2j}\;
  !>           \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !>               ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>       &+& \sum_{k=0}^{n_{rk}-1} \sum_{l=0}^{N/2}
  !>           \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>           \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>           \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !>               ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>       &+& \sum_{l=0}^{n_{rl}-1} \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>           \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>           \sum_{r=\max(-N+2n_{rk},-\min(N-l,-q)+p+q)}
  !>               ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>       &+& \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>           \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>           \sum_{r=\max(-N+2n_{rk},-\min(N-n_{rl},-q)+p+q)}
  !>               ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>       &+& \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>           \sum_{r=\max(-N+2n_{rk},-\min(N-n_{rl},-q)+\Lambda_i+q)}
  !>               ^{\min(\min(N-2k,-\Lambda_i),N-2l+\Lambda_i+q)} 1
  !>   \f}
  !======================================================================
  Subroutine calculate_Nblock()
    Use HFBTHO, Only : nrx,nlx

    Integer(ipr) :: sni,snj,snk,snl,sli,slj
    Integer(ipr) :: xni,xnj,xnk,xnl,xli,xlj
    Integer(ipr) :: ii,jj,kk,ll,i,j,k,n2,n

    n = Max(2*nrx,nlx)
    n2 = n/2
    If(Allocated(NBlock)) Deallocate(NBlock)
    Allocate(NBlock(0:n2,0:n2,0:n2,0:n2,-n:n,-n:n))
    NBlock = 0
    xni = 0
    Do ii = 0,n2
       sni = xni
       xnj = 0
       Do jj = 0,n2
          snj = xnj
          xnk = 0
          Do kk = 0,n2
             snk = xnk
             xnl = 0
             Do ll = 0,n2
                snl = xnl
                xli = 0
                Do i = -N+2*ii,N-2*ii
                   sli = xli
                   xlj = 0
                   Do j = -N+2*jj,N-2*jj
                      slj = xlj
                      Nblock(ii,jj,kk,ll,i,j)=sni+snj+snk+snl+sli+slj
                      Do k = Max(-n+2*kk,-Min(n-2*ll,-j)+i+j),Min(Min(n-2*kk,-i),n-2*ll+i+j)
                         xni = xni + 1
                         xnj = xnj + 1
                         xnk = xnk + 1
                         xnl = xnl + 1
                         xli = xli + 1
                         xlj = xlj + 1
                      End Do
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do

  End Subroutine calculate_Nblock
  !======================================================================
  !> This routine calculates what is called the LBlock. LBlock is a number that is used to obtain the
  !> transformation from the radial component quantum numbers (\f$n\f$'s and \f$\Lambda\f$'s) into the
  !> index of an array that contains only the non-zero Two body potential matrix elements. The Lblock
  !> is the part of that index that depends on five of the radial quantum  numbers (which would require
  !> an \f$2 N^5\f$ array to store in memory) but can be calculated on the fly rather quickly. The other
  !> part of that index is called the NBlock and is calculated and stored in memory by \ref calculate_nblock().
  !>
  !> @result \f$
  !>  LB = 1 + \Lambda_k - \max(-N+2n_{rk},-\min(N-2n_{rl},-\Lambda_j)
  !>                              +\Lambda_i+\Lambda_j) \f$,
  !> where \f$N\f$ is the number of shells.
  !======================================================================
  Function l_block(kn,ln,il,jl,kl,n) result(lb)
    Integer(ipr), Intent(In) :: kn !< - third radial quantum number
    Integer(ipr), Intent(In) :: ln !< - fourth radial quantum number
    Integer(ipr), Intent(In) :: il !< - first orbital quantum number
    Integer(ipr), Intent(In) :: jl !< - second orbital quantum number
    Integer(ipr), Intent(In) :: kl !< - third orbital quantum number
    Integer(ipr), Intent(In) :: n !< - total number of shells
    Integer(ipr) :: lb

    lb = 1 + kl - Max(-N+2*kn,-Min(n-2*ln,-jl)+il+jl)

    !The code comented bellow is from when the Lblock contained also
    !the sumations depending on in and kn. It is no longer needed
    !but is kept in case it becomes necessary later.

    ! mi  = Min(N-2*kn,-il)
    ! mi1 = Min(jl-1,-n+2*ln)
    ! ma1 = Max(-N+2*jn,mi1+1)
    ! mi2 = Min(mi1,2*(kn-ln)-il)
    ! ma2 = Max(-N+2*jn,mi2+1)
    ! mi3 = Min(jl-1,(2*kn-n-il-Mod(Abs(2*kn-n-il),2))/2)
    ! ma3 = Max(ma1,mi3+1)
    ! mi4 = Min(mi2,mi-n+2*ln-il)
    ! ma4 = Max(-n+2*jn,mi4+1)
    ! mi5 = Min(mi1,mi-n+2*ln-il)
    ! ma5 = Max(ma2,mi5+1)
    ! mi6 = Min(mi3,mi-n+2*ln-il)
    ! ma6 = Max(ma1,mi6+1)
    ! mi7 = Min(jl-1,mi-n+2*ln-il)
    ! ma7 = Max(ma3,mi7+1)
    ! ma8 = Max(-N+2*jn,2*(kn+ln-n)-il)
    ! mi8 = Min(mi1,mi+n-2*ln-il)
    ! ma9 = Max(ma1,2*(kn+ln-n)-il)
    ! mi9 = Min(mi7,n-2*ln)
    ! mi10= Min(jl-1,(mi-il-Mod(Abs(mi-il),2))/2)
    ! if(ma8 <= mi4) Then
    !    lb = lb -((-1+ma8-mi4)*(2+2*il-4*kn-4*ln+ma8+mi4+4*n))/2
    ! endif
    ! if(ma4 <= mi2.and.-N+2*kn <= mi) Then
    !    lb = lb + (1 - ma4 + mi2)*(1 - 2*kn + mi + n)
    ! endif
    ! if(ma2 <= mi5.and.-n+2*ln <= n-2*ln) Then
    !    lb = lb + (1 - ma2 + mi5)*(1 - 4*ln + 2*n)
    ! endif
    ! if(ma5 <= mi8) Then
    !    lb = lb + ((-1+ma5-mi8)*(-2+2*il+4*ln+ma5-2*mi+mi8-2*n))/2
    ! endif
    ! if(ma9 <= mi6) Then
    !    lb = lb -((-1+ma9-mi6)*(2+2*il-4*kn-4*ln+ma9+mi6+4*n))/2
    ! endif
    ! if(ma6 <= mi3.and.-N+2*kn <= mi) Then
    !    lb = lb + (1 - ma6 + mi3)*(1 - 2*kn + mi + n)
    ! endif
    ! if(ma3 <= mi9) Then
    !    lb = lb + ((-1+ma3-mi9)*(-2+4*ln+ma3+mi9-2*n))/2
    ! endif
    ! if(ma7 <= mi10) Then
    !    lb = lb + (-1 + ma7 - mi10)*(-1 + il + ma7 - mi + mi10)
    ! endif
    ! do i = -n+2*in,il-1
    !    Do j = -n+2*jn,n-2*jn
    !       Do k = Max(-n+2*kn,-min(n-2*ln,-j)+i+j),min(min(n-2*kn,-i),n-2*ln+i+j)
    !          lb = lb + 1
    !       End Do
    !    End Do
    ! enddo
    ! ma = Max(-N+2*jn,-N+2*ln+1)
    ! mi1 = Min(il-1,-n+2*kn)
    ! ma1 = Max(-n+2*in,mi1+1)
    ! mi2 = Min(mi1,-n+2*(kn+ln-jn))
    ! ma2 = Max(-N+2*in,-n+2*kn)
    ! if(ma2 <= mi2) Then
    !    lb = lb -((-1+ma2-mi2)*(6+12*kn**2+4*ma2+ma2**2+5*mi2+&
    !         ma2*mi2+mi2**2+3*(3+ma2+mi2)*n+3*n**2-&
    !         6*kn*(3+ma2+mi2+2*n)))/6 !si1
    ! endif
    ! mi2 = Min(N-2*jn,(-n+2*kn-mi1-Mod(-n+2*kn-mi1,2))/2)
    ! ma2 = mi2+1
    ! mi3 = Min(mi2,2*(ln-kn)-mi1)
    ! ma3 = Max(ma,mi3+1)
    ! mi4 = Min(mi3,-n+2*(kn+ln-in))
    ! ma4 = Max(ma,mi4+1)
    ! ma5 = Max(ma,-2*(n-kn-ln)-mi1)
    ! if(ma5 <= mi4) Then
    !    lb = lb -((-1+ma5-mi4)*(6+12*kn**2+12*ln**2+4*ma5+ma5**2+&
    !         9*mi1+3*ma5*mi1+3*mi1**2+5*mi4+ma5*mi4+3*mi1*mi4+&
    !         mi4**2+6*kn*(-3+4*ln-ma5-2*mi1-mi4-4*n)+&
    !         6*(3+ma5+2*mi1+mi4)*n+12*n**2-&
    !         6*ln*(3+ma5+2*mi1+mi4+4*n)))/6 !si2
    ! endif
    ! if(ma4 <= mi3.and.-n+2*in <= mi1) Then
    !    lb = lb+((-1+ma4-mi3)*(-1+2*in-mi1-n)*&
    !         (2+2*in-4*kn-4*ln+ma4+mi1+mi3+3*n))/2 !si3
    ! endif
    ! mi4 = Min(mi2,N+2*(ln-kn-in))
    ! ma4 = Max(ma3,mi4+1)
    ! mi5 = Min(mi4,-n+2*(kn+ln-in))
    ! ma5 = Max(ma3,mi5+1)
    ! if(ma3 <= mi5) Then
    !    lb = lb + (1-ma3+mi5)*(-1+4*kn-2*n)*(-1+2*kn-n) !si4
    ! endif
    ! if(ma5 <= mi4) Then
    !    lb = lb +((-1+ma5-mi4)*(-6+12*in**2-36*kn**2-6*ln+12*ln**2+&
    !         ma5-6*ln*ma5+ma5**2+2*mi4-6*ln*mi4+ma5*mi4+mi4**2+&
    !         3*(-5-4*ln+ma5+mi4)*n-9*n**2+&
    !         6*in*(1-4*kn-4*ln+ma5+mi4+2*n)+&
    !         6*kn*(5+4*ln-ma5-mi4+6*n)))/6 !si5
    ! endif
    ! if(ma3 <= mi4) Then
    !    lb = lb +((-1+ma3-mi4)*(4*kn-4*ln+ma3+2*mi1+mi4)*&
    !         (-1+4*kn-2*n))/2 !si6
    ! endif
    ! if(ma4 <= mi2.and.-N+2*in <= mi1) Then
    !    lb = lb + (1-ma4+mi2)*(1-2*in+mi1+n)*(1-4*kn+2*n)!si7
    ! endif
    ! mi3 = Min(N-2*jn,kn-in)
    ! ma3 = Max(ma2,mi3+1)
    ! mi4 = Min(mi3,-n+4*kn-2*ln)
    ! ma4 = Max(ma2,mi4+1)
    ! mi5 = Min(mi4,n+2*(ln-kn-in))
    ! ma5 = Max(ma2,mi5+1)
    ! mi6 = Min(mi5,-n+2*(kn+ln-in))
    ! ma6 = Max(ma2,mi6+1)
    ! if(ma2 <= mi6) Then
    !    lb = lb + (1-ma2+mi6)*(-1+4*kn-2*n)*(-1+2*kn-n)!si8
    ! endif
    ! if(ma6 <= mi5) Then
    !    lb = lb+((-1+ma6-mi5)*(-6+12*in**2-36*kn**2-6*ln+12*ln**2+&
    !         ma6-6*ln*ma6+ma6**2+2*mi5-6*ln*mi5+ma6*mi5+mi5**2+&
    !         3*(-5-4*ln+ma6+mi5)*n-9*n**2+&
    !         6*in*(1-4*kn-4*ln+ma6+mi5+2*n)+&
    !         6*kn*(5+4*ln-ma6-mi5+6*n)))/6!si9
    ! endif
    ! mi6 = Min(mi5,-n+4*kn-2*ln-1)
    ! if(ma2 <= mi6) Then
    !    lb = lb -((-1+ma2-mi6)*(-1+4*kn-2*n)*&
    !         (-8*kn+4*ln+ma2+mi6+2*n))/2 !si10
    ! endif
    ! if(ma5 <= mi4) Then
    !    lb = lb -((-1+ma5-mi4)*(-1+2*in-2*kn+ma5+mi4)*&
    !         (-1+4*kn-2*n))!si11
    ! endif
    ! mi5 = Min(mi3,-n+2*(kn+ln-in))
    ! ma5 = Max(ma4,mi5+1)
    ! mi6 = Min(mi5,n-2*ln)
    ! if(ma4 <= mi6) Then
    !    lb = lb -((-1+ma4-mi6)*(6+12*ln**2-5*ma4+ma4**2-4*mi6+&
    !       ma4*mi6+mi6**2+6*ln*(-3+ma4+mi6-2*n)-&
    !       3*(-3+ma4+mi6)*n+3*n**2))/6 !si12
    ! endif
    ! if(ma5 <= mi3) Then
    !    lb = lb +(-1+ma5-mi3)*(-1+2*in-2*kn+ma5+mi3)*&
    !         (1+in-kn-2*ln+n)!si13
    ! endif
    ! mi4 = Min(mi3,2*(ln-kn)-mi1)
    ! ma4 = Max(ma2,mi4+1)
    ! mi5 = Min(mi4,n-2*ln)
    ! if(ma2 <= mi5) Then
    !    lb = lb +((-1+ma2-mi5)*(4*ma2**2+3*ma2*mi1+4*ma2*mi5+3*mi1*mi5+&
    !      4*mi5**2-6*kn*(-2+4*ln+ma2+mi5-2*n)-&
    !      3*(ma2+2*mi1+mi5)*n-6*n**2+12*ln*(ma2+mi1+mi5+n)-&
    !      2*(4*ma2+3*mi1+2*mi5+3*n)))/6 !si14
    ! endif
    ! mi5 = Min(mi3,-n+4*kn-2*ln)
    ! ma5 = Max(ma4,mi5+1)
    ! mi6 = Min(mi5,(n-2*kn-mi1-Mod(n-2*kn-mi1,2))/2)
    ! ma6 = Max(ma4,mi6+1)
    ! if(ma4 <= mi6) Then
    !    lb = lb +((-1+ma4-mi6)*(-36*kn**2-5*ma4+4*ma4**2-3*mi1+&
    !         6*ma4*mi1+3*mi1**2-mi6+4*ma4*mi6+6*mi1*mi6+4*mi6**2-&
    !         3*n-6*(ma4+mi1+mi6)*n-9*n**2+&
    !         6*kn*(1+2*ma4+2*mi1+2*mi6+6*n)))/6!si15
    ! endif
    ! if(ma6 <= mi5.and.1 <= 2*(n-2*kn)) Then
    !    lb = lb +(1-ma6+mi5)*(-1+4*kn-2*n)*(2*kn-n)!si16
    ! endif
    ! mi4 = Min(mi3,n-2*ln)
    ! if(ma5 <= mi4) Then
    !    lb = lb +((-1+ma5-mi4)*(12*ln**2+(-2+ma5)*ma5-mi4+ma5*mi4+&
    !         mi4**2+6*ln*(-1+ma5+mi4)-&
    !         6*kn*(-2+4*ln+ma5+mi4-2*n)-3*n*(1+n)))/3 !si17
    ! endif
    ! mi6 = Min(mi3,(n-2*kn-mi1-Mod(n-2*kn-mi1,2))/2)
    ! ma6 = Max(ma5,mi6+1)
    ! if(ma5 <= mi6) Then
    !    lb = lb +((-1+ma5-mi6)*(4*kn**2-4*ln**2-ma5+ma5**2-mi1+&
    !         2*ma5*mi1+mi1**2+ma5*mi6+2*mi1*mi6+mi6**2+&
    !         kn*(-2+4*ma5+4*mi1+4*mi6-4*n)-&
    !         2*ln*(-1+ma5+mi6-2*n)-(ma5+2*mi1+mi6)*n))/2!si18
    ! endif
    ! mi4 = Min(mi3,n-2*ln-1)
    ! if(ma6 <= mi4) Then
    !    lb = lb -((-1+ma6-mi4)*(12*ln**2-2*ma6+ma6**2-mi4+ma6*mi4+&
    !         mi4**2+6*ln*(-1+ma6+mi4-2*n)-3*(-1+ma6+mi4)*n+&
    !         3*n**2))/6!si19
    ! endif
    ! mi4 = Min(n-2*jn,2*(ln-kn)-mi1)
    ! ma4 = Max(ma3,mi4+1)
    ! mi5 = Min(mi4,n-2*ln)
    ! if(ma3 <= mi5.and.-N+2*in <= mi1) Then
    !    lb = lb +((-1+ma3-mi5)*(-2+4*ln+ma3+mi5-2*n)*&
    !         (1-2*in+mi1+n))/2!si20
    ! endif
    ! mi5 = Min(N-2*jn,N+2*(ln-kn-in))
    ! ma5 = Max(ma4,mi5+1)
    ! mi6 = Min(mi5,n-2*ln)
    ! if(ma4 <= mi6) Then
    !    lb = lb -((-1+ma4-mi6)*(6-24*ln**2-7*ma4+2*ma4**2-5*mi6+&
    !         2*ma4*mi6+2*mi6**2+6*in*(-2+4*ln+ma4+mi6-2*n)+&
    !         6*kn*(-2+4*ln+ma4+mi6-2*n)-6*(-2+ma4+mi6)*n+&
    !         6*n**2))/6 !si21
    ! endif
    ! mi6 = Min(mi5,(n-2*kn-mi1-Mod(n-2*kn-mi1,2))/2)
    ! ma6 = Max(ma4,mi6+1)
    ! if(ma4 <= mi6) Then
    !    lb = lb +((-1+ma4-mi6)*(4*kn**2-4*ln**2-ma4+ma4**2-mi1+&
    !         2*ma4*mi1+mi1**2+ma4*mi6+2*mi1*mi6+mi6**2+&
    !         kn*(-2+4*ma4+4*mi1+4*mi6-4*n)-&
    !         2*ln*(-1+ma4+mi6-2*n)-(ma4+2*mi1+mi6)*n))/2!si22
    ! endif
    ! mi6 = Min(mi5,n-2*ln-1)
    ! if(ma6 <= mi6) Then
    !    lb = lb -((-1+ma6-mi6)*(12*ln**2-2*ma6+ma6**2-mi6+ma6*mi6+&
    !         mi6**2+6*ln*(-1+ma6+mi6-2*n)-3*(-1+ma6+mi6)*n+&
    !         3*n**2))/6!si23
    ! endif
    ! mi6 = Min(n-2*jn,(n-2*kn-mi1-Mod(n-2*kn-mi1,2))/2)
    ! ma6 = Max(ma5,mi6+1)
    ! if(ma5 <= mi6.and.-N+2*in <= mi1) Then
    !    lb = lb -((-1+ma5-mi6)*(-2+2*in+4*kn+2*ma5+mi1+2*mi6-3*n)*&
    !     (-1+2*in-mi1-n))/2!si24
    ! endif
    ! mi6 = N+min(-2*jn,-kn-in)
    ! if(ma6 <= mi6) Then
    !    lb = lb -((-1+ma6-mi6)*(6-18*kn-11*ma6-7*mi6+18*n+&
    !       2*(6*in**2+3*in*(-3+4*kn+2*ma6+2*mi6-4*n)+&
    !          2*(3*kn**2+ma6**2+ma6*mi6+mi6**2+&
    !            3*kn*(ma6+mi6-2*n)-3*(ma6+mi6)*n+3*n**2))))/6!si25
    ! endif
    ! mi2 = Min(il-1,0)
    ! mi3 = Min(mi2,n-4*ln+2*kn)
    ! ma3 = Max(ma1,mi3+1)
    ! mi4 = Min(mi3,-n+2*(kn+ln-jn))
    ! ma4 = Max(ma1,mi4+1)
    ! if(ma1 <= mi4) Then
    !    lb = lb -((-1+ma1-mi4)*(6+12*kn**2+4*ma1+ma1**2+5*mi4+&
    !         ma1*mi4+mi4**2+3*(3+ma1+mi4)*n+3*n**2-&
    !         6*kn*(3+ma1+mi4+2*n)))/6!si26
    ! endif
    ! if(ma4 <= mi3.and.jn <= ln) Then
    !    lb = lb +((-1+2*jn-2*ln)*(-1+ma4-mi3)* &
    !         (2+2*jn-4*kn-2*ln+ma4+mi3+2*n))/2!si27
    ! endif
    ! mi4 = Min(mi2,n+2*(kn-ln-jn))
    ! ma4 = Max(ma3,mi4+1)
    ! mi5 = Min(mi4,-n+2*(kn+ln-jn))
    ! ma5 = Max(ma3,mi5+1)
    ! if(ma3 <= mi5) Then
    !    lb = lb +(1-ma3+mi5)*(-1+4*ln-2*n)*(-1+2*ln-n)!si28
    ! endif
    ! if(ma5 <= mi4) Then
    !    lb = lb +((-1+ma5-mi4)*(-6+12*jn**2+12*kn**2+30*ln-36*ln**2+&
    !         ma5-6*ln*ma5+ma5**2+2*mi4-6*ln*mi4+ma5*mi4+mi4**2+&
    !         6*kn*(-1+4*ln-ma5-mi4-2*n)+&
    !         3*(-5+12*ln+ma5+mi4)*n-9*n**2+&
    !         6*jn*(1-4*kn-4*ln+ma5+mi4+2*n)))/6!si29
    ! endif
    ! if(ma3 <= mi4) Then
    !    lb = lb +((-1+ma3-mi4)*(-1+4*ln-2*n)*&
    !         (-4*kn+8*ln+ma3+mi4-2*n))/2!si30
    ! endif
    ! if(ma4 <= mi2.and.jn <= ln) Then
    !    lb = lb +(1-2*jn+2*ln)*(1-ma4+mi2)*(1-4*ln+2*n)!si31
    ! endif
    ! mi3 = Min(il-1,ln-jn)
    ! ma3 = Max(mi2,mi3)+1
    ! mi4 = Min(mi3,-n+4*ln-2*kn)
    ! ma4 = Max(mi2,mi4)+1
    ! mi5 = Min(mi4,n+2*(kn-ln-jn))
    ! ma5 = Max(mi2,mi5)+1
    ! mi6 = Min(mi5,-n+2*(kn+ln-jn))
    ! ma6 = Max(mi2,mi6)+1
    ! if(mi2+1 <= mi6) Then
    !    lb = lb +(-mi2+mi6)*(-1+4*ln-2*n)*(-1+2*ln-n)!si32
    ! endif
    ! if(ma6 <= mi5) Then
    !    lb = lb +((-1+ma6-mi5)*(-6+12*jn**2+12*kn**2+30*ln-36*ln**2+&
    !         ma6-6*ln*ma6+ma6**2+2*mi5-6*ln*mi5+ma6*mi5+mi5**2+&
    !         6*kn*(-1+4*ln-ma6-mi5-2*n)+&
    !         3*(-5+12*ln+ma6+mi5)*n-9*n**2+&
    !         6*jn*(1-4*kn-4*ln+ma6+mi5+2*n)))/6!si33
    ! endif
    ! mi6 = Min(mi5,-n+4*ln-2*kn-1)
    ! if(mi2+1 <= mi6) Then
    !    lb = lb -((mi2-mi6)*(-1+4*ln-2*n)*&
    !         (1+4*kn-8*ln+mi2+mi6+2*n))/2!si34
    ! endif
    ! if(ma5 <= mi4) Then
    !    lb = lb -((-1+ma5-mi4)*(-1+2*jn-2*ln+ma5+mi4)*&
    !         (-1+4*ln-2*n))!si35
    ! endif
    ! mi5 = Min(mi3,-n+2*(kn+ln-jn))
    ! ma5 = Max(ma4,mi5+1)
    ! if(ma4 <= mi5) Then
    !    lb = lb -((-1+ma4-mi5)*(6+12*kn**2-5*ma4+ma4**2-4*mi5+&
    !         ma4*mi5+mi5**2+6*kn*(-3+ma4+mi5-2*n)-&
    !         3*(-3+ma4+mi5)*n+3*n**2))/6!si36
    ! endif
    ! if(ma5 <= mi3) Then
    !    lb = lb +(-1+ma5-mi3)*(-1+2*jn-2*ln+ma5+mi3)*&
    !         (1+jn-2*kn-ln+n)!si37
    ! endif
    ! mi4 = Min(mi3,n-4*ln+2*kn)
    ! ma4 = Max(mi2,mi4)+1
    ! if(mi2+1 <= mi4) Then
    !    lb = lb+((mi2-mi4)*(6*kn*(1+mi2+mi4)+&
    !         2*(-1+mi2**2+mi2*mi4+mi4**2)-3*(1+mi2+mi4)*n))/3!si38
    ! endif
    ! mi5 = Min(mi3,-n+4*ln-2*kn)
    ! ma5 = Max(ma4,mi5+1)
    ! mi6 = Min(mi5,n-2*ln)
    ! ma6 = Max(ma4,mi6+1)
    ! if(ma4 <= mi6) Then
    !    lb = lb +((-1+ma4-mi6)*(4*ma4**2-mi6+&
    !         ma4*(-5+24*ln+4*mi6-12*n)+4*mi6*(6*ln+mi6-3*n)))/6!si39
    ! endif
    ! if(ma6 <= mi5.and.1 <= 2*(n-2*ln)) Then
    !    lb = lb + (1 - ma6 + mi5)*(-1 + 4*ln - 2*n)*(2*ln - n)!si40
    ! endif
    ! if(ma5 <= mi3) Then
    !    lb = lb +((-1+ma5-mi3)*(12*kn**2+(-2+ma5)*ma5-mi3+ma5*mi3+&
    !      mi3**2+6*kn*(-1-4*ln+ma5+mi3)-&
    !      6*ln*(-2+ma5+mi3-2*n)-3*n*(1+n)))/3!si41
    ! endif
    ! mi6 = Min(mi3,n-2*ln)
    ! ma6 = Max(ma5,mi6+1)
    ! if(ma5 <= mi6) Then
    !    lb = lb +((-1+ma5-mi6)*(-4*kn**2+4*ln*(-1+4*ln)-ma5+ma5**2+&
    !         ma5*mi6+mi6**2-2*kn*(-1+ma5+mi6-2*n)+&
    !         8*ln*(ma5+mi6-2*n)+n-3*(ma5+mi6)*n+3*n**2))/2!si42
    ! endif
    ! if(ma6 <= mi3) Then
    !    lb = lb -((-1+ma6-mi3)*(12*kn**2-2*ma6+ma6**2-mi3+ma6*mi3+&
    !         mi3**2+6*kn*(-1+ma6+mi3-2*n)-3*(-1+ma6+mi3)*n+&
    !         3*n**2))/6!si43
    ! endif
    ! mi4 = Min(il-1,n-4*ln+2*kn)
    ! ma4 = Max(ma3,mi4+1)
    ! if(ma3 <= mi4.and.jn <= ln) Then
    !    lb = lb -((-1+2*jn-2*ln)*(-1+ma3-mi4)*&
    !         (-2+4*kn+ma3+mi4-2*n))/2!si44
    ! endif
    ! mi5 = Min(il-1,n+2*(kn-ln-jn))
    ! ma5 = Max(ma4,mi5+1)
    ! if(ma4 <= mi5) Then
    !    lb = lb -((-1+ma4-mi5)*(6-24*kn**2+24*kn*ln-7*ma4-5*mi5+&
    !         6*jn*(-2+4*kn+ma4+mi5-2*n)+12*n+&
    !         2*(ma4**2+ma4*mi5+mi5**2+3*ln*(-2+ma4+mi5-2*n)-&
    !         3*(ma4+mi5)*n+3*n**2)))/6!si45
    ! endif
    ! mi6 = Min(mi5,n-2*ln)
    ! ma6 = Max(ma4,mi6+1)
    ! if(ma4 <= mi6) Then
    !    lb = lb+((-1+ma4-mi6)*(-4*kn**2+4*ln*(-1+4*ln)-ma4+ma4**2+&
    !         ma4*mi6+mi6**2-2*kn*(-1+ma4+mi6-2*n)+&
    !         8*ln*(ma4+mi6-2*n)+n-3*(ma4+mi6)*n+3*n**2))/2!si46
    ! endif
    ! if(ma6 <= mi5) Then
    !    lb = lb -((-1+ma6-mi5)*(12*kn**2-2*ma6+ma6**2-mi5+ma6*mi5+&
    !         mi5**2+6*kn*(-1+ma6+mi5-2*n)-3*(-1+ma6+mi5)*n+&
    !         3*n**2))/6!si47
    ! endif
    ! mi6 = Min(il-1,n-2*ln)
    ! ma6 = Max(ma5,mi6+1)
    ! if(ma5 <= mi6.and.jn <= ln) Then
    !    lb = lb -((-1+2*jn-2*ln)*(-1+ma5-mi6)*&
    !         (-1+jn+3*ln+ma5+mi6-2*n))!si48
    ! endif
    ! if(ma6 <= il-1) Then
    !    lb = lb+((il-ma6)*(17+4*il**2-30*jn-30*ln-15*ma6+&
    !         il*(-15+12*jn+12*ln+4*ma6-12*n)+30*n+&
    !         4*(3*jn**2+6*jn*ln+3*ln**2+3*jn*ma6+3*ln*ma6+ma6**2-&
    !         3*(2*(jn+ln)+ma6)*n+3*n**2)))/6!si49
    ! endif
    ! mi2 = Min(il-1,-n+ln+jn)
    ! ma2 = Max(ma1,mi2+1)
    ! mi3 = Min(N-2*jn,(-n-mi2-Mod(Abs(-n-mi2),2))/2+kn)
    ! ma3 = Max(ma,mi3+1)
    ! if(ma <= mi3.and.ma1 <= mi2) Then
    !    lb = lb +((-1+ma1-mi2)*(-1+ma-mi3)*&
    !         (2-4*kn-4*ln+ma+ma1+mi2+mi3+4*n))/2!si50
    ! endif
    ! mi4 = Min(N-2*jn,(-n-ma1-Mod(Abs(-n-ma1),2))/2+kn)
    ! ma4 = Max(ma3,mi4+1)
    ! if(ma3 <= mi4) Then
    !    lb = lb -((-1+ma3-mi4)*(-2+2*kn+4*ln-ma1-3*n)*&
    !         (-1-2*kn+ma1+ma3+mi4+n))/2!si51
    ! endif
    ! if(ma3 <= mi4) Then
    !    lb = lb +((-1+ma3-mi4)*(4*ma3**2+3*ma3*mi2+4*ma3*mi4+3*mi2*mi4+&
    !         4*mi4**2-6*kn*(-2+4*ln+ma3+mi4-2*n)-&
    !         3*(ma3+2*mi2+mi4)*n-6*n**2+12*ln*(ma3+mi2+mi4+n)-&
    !         2*(4*ma3+3*mi2+2*mi4+3*n)))/6!si52
    ! endif
    ! mi5 = Min(N-2*jn,n-2*ln)
    ! if(ma4 <= mi5.and.ma1 <= mi2) Then
    !    lb = lb -((-1+ma1-mi2)*(-1+ma4-mi5)*&
    !         (-2+4*ln+ma4+mi5-2*n))/2!si53
    ! endif
    ! mi3 = Min(il-1,(-n+2*ln-ma-Mod(Abs(-n+2*ln-ma),2))/2)
    ! ma3 = Max(ma2,mi3+1)
    ! if(ma2 == mi3) Then
    !    mi4 = Min(-n+2*ln-2*ma2,(-n-ma2-Mod(Abs(-n-ma2),2))/2+kn)
    !    ma4 = Max(ma,mi4+1)
    !    If(ma <= mi4) Then
    !       lb = lb -((-1+ma-mi4)*(2-4*kn-4*ln+ma+2*ma2+mi4+4*n))/2!si54
    !    End If
    !    mi5 = Min(-n+2*ln-2*ma2,n-2*ln)
    !    If(ma4 <= mi5) Then
    !       lb = lb +((-1+ma4-mi5)*(-2+4*ln+ma4+mi5-2*n))/2!si55
    !    End If
    !    mi4 = Min(N-2*jn,(-n-ma2-Mod(Abs(-n-ma2),2))/2+kn)
    !    ma4 = Max(-n+2*ln-2*ma2,mi4)+1
    !    If(-n+2*ln-2*ma2+1 <= mi4) Then
    !       lb = lb +(1-2*kn-ma2+n)*(-2*ln+2*ma2+mi4+n)!si56
    !    End If
    !    mi5 = Min(N-2*jn,-ma2)
    !    If(ma4 <= mi5) Then
    !       lb = lb + (-1 + ma4 - mi5)*(-1 + 2*ma2 + ma4 + mi5)!si57
    !    End If
    ! Else If(ma2 < mi3) Then
    !    If(Mod(-mi3,2) == 1) Then
    !       mi4 = Min(-n+2*ln-2*mi3,(-n-mi3-Mod(Abs(-n-mi3),2))/2+kn)
    !       ma4 = Max(ma,mi4+1)
    !       If(ma <= mi4) Then
    !          lb = lb -((-1+ma-mi4)*&
    !               (2-4*kn-4*ln+ma+2*mi3+mi4+4*n))/2!si58
    !       End If
    !       mi5 = Min(-n+2*ln-2*mi3,n-2*ln)
    !       If(ma4 <= mi5) Then
    !          lb = lb +((-1+ma4-mi5)*(-2+4*ln+ma4+mi5-2*n))/2!si59
    !       End If
    !       mi4 = Min(n-2*jn,(-n-mi3-Mod(Abs(-n-mi3),2))/2+kn)
    !       ma4 = Max(-n+2*ln-2*mi3+1,mi4+1)
    !       If(-n+2*ln-2*mi3+1 <= mi4) Then
    !          lb = lb +(1-2*kn-mi3+n)*(-2*ln+2*mi3+mi4+n)!si60
    !       End If
    !       mi5 = Min(N-2*jn,-mi3)
    !       If(ma4 <= mi5) Then
    !          lb = lb +(-1 + ma4 - mi5)*(-1 + ma4 + 2*mi3 + mi5)!si61
    !       End If
    !    End If
    !    mi4 = (-ma2-Mod(-ma2+1,2))/2
    !    ma4 = (-mi3+Mod(-mi3,2))/2
    !    mn2 = (-n-Mod(N,2))/2
    !    mi5 = Min(mi4,(N-2*ln+kn+mn2-Mod(Abs(N-2*ln+kn+mn2),3))/3)
    !    ma5 = Max(ma4,mi5+1)
    !    If(ma4 <= mi5) Then
    !       lb = lb +((-1 + ma4 - mi5)*(-2 + 4*kn + 2*ln - ma - 3*n)*&
    !            (1 + 2*ln - ma + 2*ma4 + 2*mi5 - n))/2 !sx1
    !    End If
    !    mi6 = Min(mi4,ma-kn-mn2-1)
    !    ma6 = Max(ma5,mi6+1)
    !    mi7 = Min(mi6,(n-Mod(n,2))/2-ln)
    !    ma7 = Max(ma5,mi7+1)
    !    If(ma5 <= mi7) Then
    !       lb = lb+((-1+ma5-mi7)*(36*ln**2-3*(-3+ma)*ma+&
    !            2*(-3+ma5*(-7+8*ma5)+mi7+8*ma5*mi7+8*mi7**2)+&
    !            6*ln*(1-2*ma+8*ma5+8*mi7-6*n)-3*n+&
    !            6*(ma - 4*(ma5+mi7))*n+9*n**2))/6 !sx2
    !    End If
    !    If(ma7 <= mi6.and.ma <= n-2*ln) Then
    !       lb = lb+((1-ma7+mi6)*(-2+2*ln+ma-n)*(-1+2*ln+ma-n))/2! sx3
    !    End If
    !    If(ma6 <= mi4) Then
    !       lb = lb+((-1+ma6-mi4)*(-2+3*kn**2+ma+ma**2-2*ma*ma6+ma6**2+&
    !            mi4-2*ma*mi4+ma6*mi4+mi4**2-3*mn2+ma6*mn2+mi4*mn2-&
    !            mn2**2+2*ln*(2-2*ma+ma6+mi4+2*mn2)+&
    !            kn*(1+4*ln-4*ma+3*ma6+3*mi4+2*mn2-4*n)+&
    !            2*(-2+2*ma-ma6-mi4-2*mn2)*n))/2!sx4
    !    End If
    !    mi7 = Min(mi4,(n-Mod(n,2))/2-ln)
    !    ma7 = Max(ma6,mi7+1)
    !    If(ma6 <= mi7) Then
    !       lb = lb +((-1+ma6-mi7)*(-kn**2+12*ln**2-4*ma6+5*ma6**2+mi7+&
    !            5*ma6*mi7+5*mi7**2+mn2-ma6*mn2-mi7*mn2-mn2**2+&
    !            2*ln*(-1+7*ma6+7*mi7-2*mn2-6*n)-&
    !            kn*(-1+4*ln+ma6+mi7+2*mn2-2*n)+n-7*ma6*n-&
    !            7*mi7*n+2*mn2*n+3*n**2))/2 !sx5
    !    End If
    !    If(ma7 <= mi4) Then
    !       lb = lb -((-1+ma7-mi4)*(3*kn**2+12*ln**2-2*ma7+ma7**2-mi4+&
    !            ma7*mi4+mi4**2-3*mn2+3*ma7*mn2+3*mi4*mn2+3*mn2**2+&
    !            6*ln*(-1+ma7+mi4+2*mn2-2*n)+&
    !            3*kn*(-1+4*ln+ma7+mi4+2*mn2-2*n)-&
    !            3*(-1+ma7+mi4+2*mn2)*n+3*n**2))/6!sx6
    !    End If
    !    mi5 = Min(mi4,(mn2+kn+n-2*ln-1-Mod(Abs(mn2+kn+n-2*ln-1),3))/3)
    !    ma5 = Max(ma4,mi5+1)
    !    If(ma4 <= mi5) Then
    !       lb = lb +((-1+ma4-mi5)*(4*kn**2+4*ln+ma4+5*mi5+&
    !            4*(ma4*(ln+ma4)+(ln+ma4)*mi5+mi5**2)-&
    !            2*kn*(1+4*ln+4*ma4+4*mi5-2*mn2)-2*(ma4+mi5)*mn2+&
    !            (2*kn+4*ln+ma4+mi5-2*mn2)*n-2*n**2-2*(mn2+n)))/2!sx7
    !    End If
    !    mi6 = Min(mi5,(n-Mod(n,2))/2-jn)
    !    ma6 = Max(ma4,mi6+1)
    !    If(ma4 <= mi6) Then
    !       lb = lb -((-1+ma4-mi6)*(-ma4+mi6+&
    !            2*(3*kn**2+ma4**2+ma4*mi6+mi6**2-&
    !            3*kn*(ma4+mi6-2*mn2)-3*(ma4+mi6)*mn2+3*mn2**2)))/6!sx8
    !    End If
    !    If(ma6 <= mi5) Then
    !       lb = lb +((-1+ma6-mi5)*(8*jn**2-ma6+mi5+8*jn*(ma6+mi5-n)+&
    !            2*(-kn**2+ma6**2+ma6*mi5+mi5**2+kn*(ma6+mi5-2*mn2)+&
    !            ma6*mn2+mi5*mn2-mn2**2-2*(ma6+mi5)*n+n**2)))/2!sx9
    !    End If
    !    mi6 = Min(mi4,(n-Mod(n,2))/2-jn)
    !    ma6 = Max(ma5,mi6+1)
    !    mi7 = Min(mi6,(n-1-Mod(n-1,2))/2-ln)
    !    If(ma5 <= mi7) Then
    !       lb = lb -((-1+ma5-mi7)*(-2*ma5+2*mi7+&
    !            4*(3*ln**2+ma5**2+ma5*mi7+mi7**2+3*ln*(ma5+mi7))-&
    !            6*(2*ln+ma5+mi7)*n+3*n**2))/3!sx10
    !    End If
    !    If(ma6 <= mi4) Then
    !       lb = lb +4*(jn-ln)*(-1+ma6-mi4)*(jn+ln+ma6+mi4-n)!sx11
    !    End If
    !    mn2 = (-n+1-Mod(N+1,2))/2
    !    mi5 = Min(mi4,(N-2*ln+kn+mn2-2-Mod(Abs(N-2*ln+kn+mn2-2),3))/3)
    !    ma5 = Max(ma4,mi5+1)
    !    If(ma4 <= mi5) Then
    !       lb = lb +((-1 + ma4 - mi5)*(-2 + 4*kn + 2*ln - ma - 3*n)*&
    !            (3 + 2*ln - ma + 2*ma4 + 2*mi5 - n))/2!sy1
    !    End If
    !    mi6 = Min(mi4,ma-kn-mn2-1)
    !    ma6 = Max(ma5,mi6+1)
    !    mi7 = Min(mi6,(n-1-Mod(n-1,2))/2-ln)
    !    ma7 = Max(ma5,mi7+1)
    !    If(ma5 <= mi7) Then
    !       lb = lb +((-1+ma5-mi7)*(36*ln**2+9*ma-3*ma**2+10*ma5+26*mi7+&
    !            16*(ma5**2+ma5*mi7+mi7**2)-27*n+&
    !            6*(ma-4*(ma5+mi7))*n+9*n**2-&
    !            6*ln*(-9+2*ma-8*ma5-8*mi7+6*n)))/6 !sy2
    !    End If
    !    If(ma7 <= mi6.and.ma <= n-2*ln) Then
    !       lb = lb +((1-ma7+mi6)*(-2+2*ln+ma-n)*(-1+2*ln+ma-n))/2!sy3
    !    End If
    !    If(ma6 <= mi4) Then
    !       lb = lb +((-1+ma6-mi4)*(3*kn**2-ma+ma**2+&
    !            ma6-2*ma*ma6+ma6**2+&
    !            2*mi4-2*ma*mi4+ma6*mi4+mi4**2-mn2+ma6*mn2+mi4*mn2-&
    !            mn2**2+2*ln*(2-2*ma+ma6+mi4+2*mn2)+&
    !            kn*(3+4*ln-4*ma+3*ma6+3*mi4+2*mn2-4*n)+&
    !            2*(-2+2*ma-ma6-mi4-2*mn2)*n))/2 !sy4
    !    End If
    !    mi7 = Min(mi4,(n-1-Mod(n-1,2))/2-ln)
    !    ma7 = Max(ma6,mi7+1)
    !    If(ma6 <= mi7) Then
    !       lb = lb +((-1+ma6-mi7)*(2-kn**2+12*ln**2+&
    !            4*ma6+5*ma6**2+9*mi7+&
    !            5*ma6*mi7+5*mi7**2+mn2-ma6*mn2-mi7*mn2-mn2**2+&
    !            2*ln*(7+7*ma6+7*mi7-2*mn2-6*n)-&
    !            kn*(-1+4*ln+ma6+mi7+2*mn2-2*n)-&
    !            (7*(1+ma6+mi7)-2*mn2)*n+3*n**2))/2 !sy5
    !    End If
    !    If(ma7 <= mi4) Then
    !       lb = lb -((-1+ma7-mi4)*(3*kn**2+12*ln**2-2*ma7+ma7**2-mi4+&
    !            ma7*mi4+mi4**2-3*mn2+3*ma7*mn2+3*mi4*mn2+3*mn2**2+&
    !            6*ln*(-1+ma7+mi4+2*mn2-2*n)+&
    !            3*kn*(-1+4*ln+ma7+mi4+2*mn2-2*n)-&
    !            3*(-1+ma7+mi4+2*mn2)*n+3*n**2))/6!sy6
    !    End If
    !    mi5 = Min(mi4,(mn2+kn+n-2*ln-Mod(Abs(mn2+kn+n-2*ln),3))/3-1)
    !    ma5 = Max(ma4,mi5+1)
    !    If(ma4 <= mi5) Then
    !       lb = lb +((-1+ma4-mi5)*(4*kn**2+&
    !            4*(2+2*ma4+ma4**2+(3+ma4)*mi5+mi5**2+&
    !            ln*(2+ma4+mi5))-2*(2+ma4+mi5)*mn2+&
    !            (4*ln+ma4+mi5-2*mn2)*n-2*n**2+&
    !            2*kn*(-6-4*ln-4*ma4-4*mi5+2*mn2+n)))/2!sy7
    !    End If
    !    mi6 = Min(mi5,(n-1-Mod(n-1,2))/2-jn)
    !    ma6 = Max(ma4,mi6+1)
    !    If(ma4 <= mi6) Then
    !       lb = lb -((-1+ma4-mi6)*(6+6*kn**2+5*ma4+2*ma4**2+7*mi6+&
    !            2*ma4*mi6+2*mi6**2-6*kn*(2+ma4+mi6-2*mn2)-&
    !            6*(2+ma4+mi6)*mn2+6*mn2**2))/6!sy8
    !    End If
    !    If(ma6 <= mi5) Then
    !       lb = lb+((-1+ma6-mi5)*(8*jn**2+ma6+3*mi5+4*mn2+&
    !            8*jn*(1+ma6+mi5-n)-4*n+&
    !            2*(-kn**2+ma6**2+ma6*mi5+mi5**2+&
    !            kn*(2+ma6+mi5-2*mn2)+ma6*mn2+mi5*mn2-mn2**2-&
    !            2*(ma6+mi5)*n+n**2)))/2!sy9
    !    End If
    !    mi6 = Min(mi4,(n-1-Mod(n-1,2))/2-jn)
    !    ma6 = Max(ma5,mi6+1)
    !    mi7 = Min(mi6,(n-Mod(n,2))/2-ln-1)
    !    If(ma5 <= mi7) Then
    !       lb = lb -((-1+ma5-mi7)*(3+4*ma5+8*mi7+&
    !            4*(3*ln**2+ma5**2+ma5*mi7+mi7**2+&
    !            3*ln*(1+ma5+mi7))-6*n-6*(2*ln+ma5+mi7)*n+3*n**2))/3!sy10
    !    End If
    !    If(ma6 <= mi4) Then
    !       lb = lb +4*(jn-ln)*(-1+ma6-mi4)*(1+jn+ln+ma6+mi4-n)!sy11
    !    End If
    !    If(Mod(-ma2,2) == 0) Then
    !       mi4 = Min(-n+2*ln-2*ma2,(-n-ma2-Mod(Abs(-n-ma2),2))/2+kn)
    !       ma4 = Max(ma,mi4+1)
    !       If(ma <= mi4) Then
    !          lb = lb -((-1+ma-mi4)*&
    !               (2-4*kn-4*ln+ma+2*ma2+mi4+4*n))/2!si62
    !       End If
    !       mi5 = Min(-n+2*ln-2*ma2,n-2*ln)
    !       If(ma4 <= mi5) Then
    !          lb = lb +((-1+ma4-mi5)*(-2+4*ln+ma4+mi5-2*n))/2!si63
    !       End If
    !       mi4 = Min(n-2*jn,(-n-ma2-Mod(Abs(-n-ma2),2))/2+kn)
    !       ma4 = Max(-n+2*ln-2*ma2,mi4)+1
    !       If(-n+2*ln-2*ma2+1 <= mi4) Then
    !          lb = lb +(1-2*kn-ma2+n)*(-2*ln+2*ma2+mi4+n)!si64
    !       End If
    !       mi5 = Min(N-2*jn,-ma2)
    !       If(ma4 <= mi5) Then
    !          lb = lb +(-1+ma4-mi5)*(-1+2*ma2+ma4+mi5)!si65
    !       End If
    !    End If
    ! endif
    ! mi4 = Min(N-2*jn,(-n-il+1-Mod(Abs(-n-il+1),2))/2+kn)
    ! ma4 = Max(ma,mi4+1)
    ! if(ma <= mi4.and.ma3 <= il-1) Then
    !    lb = lb +((il-ma3)*(-1+ma-mi4)*(-3+il+4*kn+ma3-2*n))/2!si66
    ! endif
    ! mi5 = Min(N-2*jn,(-n-ma3-Mod(Abs(-n-ma3),2))/2+kn)
    ! ma5 = Max(ma4,mi5+1)
    ! if(ma4 <= mi5) Then
    !    lb = lb +((-1+ma4-mi5)*(36*kn**2-3*(-3+ma3)*ma3+ma4+4*ma4**2+&
    !         4*ma4*mi5+mi5*(5+4*mi5)+6*(ma3+2*(ma4+mi5))*n+&
    !         9*n**2-3*(2+n)-6*kn*(-1+2*ma3+4*ma4+4*mi5+6*n)))/6!si67
    ! endif
    ! mi6 = Min(mi5,-il+1)
    ! ma6 = Max(ma4,mi6+1)
    ! if(ma4 <= mi6) Then
    !    lb = lb +(-1+ma4-mi6)*(-1+il+2*kn-n)*(-1+il-2*kn+ma4+mi6+n)!si68
    ! endif
    ! if(ma6 <= mi5) Then
    !    lb = lb -((-1+ma6-mi5)*(-ma6+mi5+&
    !         2*(12*kn**2+ma6**2+ma6*mi5+mi5**2+3*(ma6+mi5)*n+&
    !         3*n**2-6*kn*(ma6+mi5+2*n))))/6!si69
    ! endif
    ! mi6 = Min(N-2*jn,-il+1)
    ! ma6 = Max(ma5,mi6+1)
    ! if(ma5 <= mi6.and.ma3 <= il-1) Then
    !    lb = lb +(il-ma3)*(-1+ma5-mi6)*(-2+il+ma3+ma5+mi6)!si70
    ! endif
    ! mi6 = Min(N-2*jn,-ma3)
    ! if(ma6 <= mi6) Then
    !    lb = lb -((-1+ma6-mi6)*(6-7*ma6-5*mi6+&
    !         2*(3*ma3**2+ma6**2+ma6*mi6+mi6**2+&
    !         3*ma3*(-2+ma6+mi6))))/6!si71
    ! endif
  End Function l_block
  !======================================================================
  !> Given four axial quantum numbers, this function returns the corresponding index in the one-dimensional
  !> array that contains all non-zero and different axial components. If the sum of the four quantum numbers
  !> is not an even number, in which case \f$ V^z_{ijkl} = 0\f$, the function returns 0 as the index.
  !======================================================================
  Function zindex(nzi,nzj,nzk,nzl) result(iz)
    Integer(ipr) :: nzi !< - First  axial quantum number
    Integer(ipr) :: nzj !< - Second axial quantum number
    Integer(ipr) :: nzk !< - Third  axial quantum number
    Integer(ipr) :: nzl !< - Fourth axial quantum number
    Integer(ipr) :: iz

    If(Mod(nzi+nzj+nzk+nzl,2) == 0) Then
       iz = ZBlock(Min(nzi,nzk),Min(nzj,nzl),Max(nzi,nzk),Max(nzj,nzl))
    Else
       iz = 0
    End If

  End Function zindex
  !======================================================================
  !> Given four radial and four angular quantum numbers, this function returns the corresponding index
  !> in the one-dimensional array that contains all non-zero and different radial components of two-body
  !> potential matrix elements for a gaussian potential. If the quantum numbers correspond to a matrix
  !> element such that \f$ V^r_{ijkl} = 0\f$, the function returns 0 as the index.
  !=====================================================================
  Function rindex(nri,nrj,nrk,nrl,li,lj,lk,ll,n) result(ir)
    Integer(ipr) :: nri !< - First  radial quantum number
    Integer(ipr) :: nrj !< - Second radial quantum number
    Integer(ipr) :: nrk !< - Third  radial quantum number
    Integer(ipr) :: nrl !< - Fourth radial quantum number
    Integer(ipr) :: li !< - First  angular quantum number
    Integer(ipr) :: lj !< - Second angular quantum number
    Integer(ipr) :: lk !< - Third  angular quantum number
    Integer(ipr) :: ll !< - Fourth angular quantum number
    Integer(ipr) :: n !< - Total number of shells
    Integer(ipr) :: ir

    If(li+lj == ll+lk) Then
       If(lk <= -li .And. ll <= -lj) Then
         ir=nblock(nri,nrj,nrk,nrl,li,lj)+l_block(nrk,nrl,li,lj,lk,n)
       Else If(ll <= -lj) Then
         ir=nblock(nrk,nrj,nri,nrl,-lk,lj)+l_block(nri,nrl,-lk,lj,-li,n)
       Else If(lk <= -li) Then
         ir=nblock(nri,nrl,nrk,nrj,li,-ll)+l_block(nrk,nrj,li,-ll,lk,n)
       Else
         ir=nblock(nrk,nrl,nri,nrj,-lk,-ll)+l_block(nri,nrj,-lk,-ll,-li,n)
       End If
    Else
       ir = 0
    End If

  End Function rindex
  !======================================================================
  !> Given two equaly-sized, square, symmetric and block-diagonal matrices, this function returns the
  !> trace of the product of such matrices. The matrices are given as two-dimensional arrays. The
  !> second index indicates the block within the matrix. The first index gives the matrix element
  !> within each block. In particular the \f$ (i,j) \f$ pair is transformed into an index given by
  !> \f$ i + (j-1)*d \f$, where \f$ d \f$ is the size of the block.
  !======================================================================
  Function trace_product(A,B) result(tr)
    Use HFBTHO, Only : nb,nbx,id

    Real(pr), Allocatable :: A(:,:) !< - First  Matrix
    Real(pr), Allocatable :: B(:,:) !< - Second Matrix
    Real(pr) :: tr,fac

    Integer(ipr) :: ib,ibx,n1,n2,nd,n12

    tr = 0.0_pr
    If(.not. Allocated(A) .Or. .not. Allocated(B)) Return

    tr = 0.0_pr
    Do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       Do n1 = 1,nd
          Do n2 = 1,n1
             If(n1 == n2) Then
                fac = one
             Else
                fac = two
             End If
             n12 = n1+(n2-1)*nd
             tr = tr + fac*A(n12,ib )*B(n12,ib )
             tr = tr + fac*A(n12,ibx)*B(n12,ibx)
          End Do!n2
       End Do!n1
    End Do!ib

  End Function trace_product
  !======================================================================
  !> Given two equaly-sized, square, symmetric and block-diagonal matrices, this function returns the
  !> trace of the product of such matrices. The matrices used with this subroutine are composed of two
  !> major blocks (one for neutrons and one for protons) which in turn are formed by smaller blocks.
  !> The subroutine returns a trace for each major block and the total trace is the sum of both traces.
  !> The matrices are given as two-dimensional arrays. The second index indicates the block within the
  !> matrix. The first index gives the matrix element within each block. In particular the \f$ (i,j) \f$
  !> pair is transformed into an index given by \f$ i + (j-1)*d \f$, where \f$ d \f$ is the size of the
  !> block.
  !======================================================================
  Subroutine trace_product_2(A,B,tr1,tr2)
    Use HFBTHO, Only : nb,nbx,id

    Real(pr), Allocatable , Dimension(:,:) :: A !< - First  Matrix
    Real(pr), Allocatable , Dimension(:,:) :: B !< - Second Matrix
    Real(pr), Intent(Out) :: tr1 !< - First half of the trace
    Real(pr), Intent(Out) :: tr2 !< - Second half of the trace

    Real(pr) :: fac
    Integer(ipr) :: ib,ibx,n1,n2,nd,n12

    tr1 = zero
    tr2 = zero
    If(.Not. Allocated(A) .Or. .Not. Allocated(B)) Return

    tr1 = zero
    tr2 = zero
    Do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       Do n1 = 1,nd
          Do n2 = 1,n1
             If(n1 == n2) Then
                fac = one
             Else
                fac = two
             End If
             n12 = n1+(n2-1)*nd
             tr1 = tr1 + fac*A(n12,ib )*B(n12,ib )
             tr2 = tr2 + fac*A(n12,ibx)*B(n12,ibx)
          End Do!n2
       End Do!n1
    End Do!i

  End Subroutine trace_product_2
  !======================================================================
  !> This subroutine calculates the the trace of a block diagonal matrix. The matrices used with this
  !> subroutine are composed of two major blocks (one for neutrons and one for protons) which in turn
  !> are formed by smaller blocks. The subroutine returns a trace for each major block and the total
  !> trace is the sum of both traces. The matrix is given as a two-dimensional array. The second index
  !> indicates the block within the matrix. The first index gives the matrix element within each block.
  !> In particular the \f$ (i,j) \f$ pair is transformed into an index given by \f$ i + (j-1)*d \f$,
  !> where \f$ d \f$ is the size of the block.
  !======================================================================
  Subroutine bdiag_trace(A,tr1,tr2)
    Use HFBTHO, Only : nb,nbx,id

    Real(pr), Allocatable , Dimension(:,:), Intent(In) :: A !< - Block diagonal matrix
    Real(pr) , intent(out) :: tr1 !< - First  half of the trace
    Real(pr) , intent(out) :: tr2 !< - Second half of the trace

    Integer(ipr) :: ib,ibx,n1,nd,n11

    tr1 = zero
    tr2 = zero
    If(.not. Allocated(A)) Return

    Do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       Do n1 = 1,nd
          n11 = n1+(n1-1)*nd
          tr1 = tr1 + A(n11,ib )
          tr2 = tr2 + A(n11,ibx)
       End Do!n1
    End Do!ib

  End Subroutine bdiag_trace
  !======================================================================
  !> This subroutine simply prints a given block diagonal matrix. The matrix is given as a two-dimensional
  !> array. The second index indicates the block within the matrix. The first index gives the matrix
  !> element within each block. In particular the \f$ (i,j) \f$ pair is transformed into an index given
  !> by \f$ i + (j-1)*d \f$, where \f$ d \f$ is the size of the block.
  !======================================================================
  Subroutine bdiag_print(A)
    Use HFBTHO, Only : nb,nbx,id

    Real(pr), Allocatable , Dimension(:,:) :: A !< - Block diagonal matrix
    Integer(ipr) :: ib,ibx,n1,n2,nd,n12

    If(.not. Allocated(A)) Return

    Do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       Do n1 = 1,nd
          Do n2 = 1,nd
             n12 = n2+(n1-1)*nd
             Write(*,*) ib ,n1,n2,n12,A(n12,ib )
          End Do
       End Do!n1
    End Do!ib

  End Subroutine bdiag_print
  !======================================================================
  !> This routine tests the orthonormality of the Harmonic oscillator wave functions in axial symmetry by
  !> printing the integral of all products of two wave functions. The integrals are calculated using the
  !> appropriate Gaussian quadrature technique. The two-body matrix elements of a contact interaction are
  !> also printed. This subroutine is used for debugging only and is not used in the actual DFT calculations.
  !======================================================================
  Subroutine test_HOWF_gauss()
    Use HFBTHO_Gauss, Only : ngl,ngh,xh,wh
    Use HFBTHO, Only : bz,bp,n00

    Integer(ipr) :: i,ni,nj,nk,nl,li,lj,lk,ll
    Real(pr) :: s,Lni,Lnj,Nr_i,Nr_j
    Real(pr) :: Hni,Hnj,Nz_i,Nz_j,wher
    Real(pr), Dimension(1:ngl) :: wlag,xlag

    Do ni = 0,n00/2
       Do nj = 0,n00/2
          Do li = 0,n00-2*ni
             lj = li
             Nr_i = N_radial(ni,li)
             Nr_j = N_radial(nj,lj)
             s = 0._pr
             Call GaussLaguerreWX((li+lj)*0.5_pr,wlag,xlag)
             Do i = 1,ngl
                Call LaguerreL(ni,Real(li,Kind=pr),xlag(i),Lni)
                Call LaguerreL(nj,Real(lj,Kind=pr),xlag(i),Lnj)
                s = s + wlag(i)*lni*lnj
             End Do
             Write(*,'(3i3,f15.8)') ni,nj,li,s*Nr_i*Nr_j*bp**2/2._pr
          End Do
       End Do
    End Do
    Write(*,*)
    Do ni = 0,n00/2
       Do nj = 0,n00/2
          Do nk = 0,n00/2
             Do nl = 0,n00/2
                Do li = -n00+2*ni,n00-2*ni
                   Do lj = -n00+2*nj,n00-2*nj
                      Do lk = -n00+2*nk,n00-2*nk
                         ll = li+lj-lk
                         If(Abs(ll) > n00-2*nl) Cycle
                         s = MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll)
                         Write(*,'(8i3,f15.8)') ni,nj,nk,nl,li,lj,lk,ll,s
                      End Do
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    Write(*,*)
    Do ni = 0,n00
       Do nj = 0,n00
          s = 0._pr
          Nz_i = N_axial(ni)
          Nz_j = N_axial(nj)
          Do i = 1,ngh
             wher = wh(i)*Exp(-xh(i)**2)
             Hni = HermiteH(ni,xh(i))
             Hnj = HermiteH(nj,xh(i))
             s = s + wher*Hni*Hnj*bz
          End Do
          Write(*,'(2i3,f15.8)') ni,nj,s*Nz_i*Nz_j
       End Do
    End Do
    Write(*,*)
    Do ni = 0,n00
       Do nj = 0,n00
          Do nk = 0,n00
             Do nl = 0,n00
                If(Mod(ni+nj+nk+nl,2) /= 0) Cycle
                Write(*,'(4i3,f15.8)') ni,nj,nk,nl,MatrixElement_zZR(ni,nj,nk,nl)
             End Do
          End Do
       End Do
    End Do
    Write(*,*)

  End Subroutine test_HOWF_gauss
  !======================================================================
  !> Given four principal quantum numbers, this function calculates the two-body matrix element of a
  !> contact interaction using one-dimensional harmonic oscillator wave functions. The matrix element
  !> integral is calculated using Gauss-Hermite quadrature. This function is used for debugging only
  !> and is not used in the actual DFT calculations.
  !>
  !> @result \f$ V_{n_i,n_j,n_k,n_l} = \int \phi_{n_i}(z_1;b_z)
  !>         \phi_{n_j}(z_2;b_z) \delta(z_1-z_2) \phi_{n_k}(z_1;b_z)
  !>         \phi_{n_l}(z_2;b_z) dz_1 dz_2 \f$
  !======================================================================
  Function MatrixElement_zZR(ni,nj,nk,nl) result(ME)
    Use HFBTHO_Gauss, Only : ngh,xh,wh
    Use HFBTHO, Only : bz

    Integer(ipr), Intent(In) :: ni !< - first  axial quantum number
    Integer(ipr), Intent(In) :: nj !< - second axial quantum number
    Integer(ipr), Intent(In) :: nk !< - third  axial quantum number
    Integer(ipr), Intent(In) :: nl !< - fourth axial quantum number
    Real(pr) :: ME

    Integer(ipr) :: i
    Real(pr) :: Nz_i,Nz_j,Nz_k,Nz_l,wher,Hni,Hnj,Hnk,Hnl

    ME = 0._pr
    If(Mod(ni+nj+nk+nl,2) /= 0) Return
    Nz_i = N_axial(ni)
    Nz_j = N_axial(nj)
    Nz_k = N_axial(nk)
    Nz_l = N_axial(nl)
    Do i = 1,ngh
       wher = wh(i)*Exp(-xh(i)**2)
       Hni = HermiteH(ni,xh(i)/Sqrt(2.0_pr))
       Hnj = HermiteH(nj,xh(i)/Sqrt(2.0_pr))
       Hnk = HermiteH(nk,xh(i)/Sqrt(2.0_pr))
       Hnl = HermiteH(nl,xh(i)/Sqrt(2.0_pr))
       ME = ME + wher*Hni*Hnj*Hnk*Hnl*bz/Sqrt(2.0_pr)
    End Do
    ME = ME*Nz_i*Nz_j*Nz_k*Nz_l

  End Function MatrixElement_zZR
  !======================================================================
  !> Given four radial and four angular quantum numbers,this function calculates the two-body matrix
  !> element of a contact interaction using two-dimensional harmonic oscillator wave functions in radial
  !> coordinates. The matrix element integral is calculated using Gauss-Laguerre quadrature. This
  !> function is used for debugging only and is not used in the actual DFT calculations.
  !>
  !> @result \f$ V_{n_i,\Lambda_i,n_j,\Lambda_j,n_k,\Lambda_k,n_l,
  !>             \Lambda_l} = \int \phi_{n_i,\Lambda_i}(r_1,\varphi_1;
  !>             b_\perp) \phi_{n_j,\Lambda_j}(r_2,\varphi_2;b_\perp)
  !>             \delta(r_1-r_2) \delta(\varphi_1-\varphi_2)
  !>             \phi_{n_k,\Lambda_k}(r_1,\varphi_1;b_\perp)
  !>             \phi_{n_l,\Lambda_l}(r_2,\varphi_2;b_\perp) dr_1
  !>             d \varphi_2 dr_2 d \varphi_2 \f$
  !======================================================================
  Function MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll) result(ME)
    Use HFBTHO_Gauss, Only : ngl
    Use HFBTHO, Only : bp

    Integer(ipr), Intent(In) :: ni !< - first  radial principal quantum number
    Integer(ipr), Intent(In) :: nj !< - second radial principal quantum number
    Integer(ipr), Intent(In) :: nk !< - third  radial principal quantum number
    Integer(ipr), Intent(In) :: nl !< - fourth radial principal quantum number
    Integer(ipr), Intent(In) :: li !< - first  radial orbital quantum number
    Integer(ipr), Intent(In) :: lj !< - second radial orbital quantum number
    Integer(ipr), Intent(In) :: lk !< - third  radial orbital quantum number
    Integer(ipr), Intent(In) :: ll !< - fourth radial orbital quantum number
    Real(pr) :: ME

    Integer(ipr) :: i
    Real(pr) :: Nr_i,Nr_j,Nr_k,Nr_l,Lni,Lnj,Lnk,Lnl
    Real(pr), Dimension(1:ngl) :: wlag,xlag

    ME = 0._pr
    If(li+lj /= lk+ll) Return
    Nr_i = N_radial(ni,li)
    Nr_j = N_radial(nj,lj)
    Nr_k = N_radial(nk,lk)
    Nr_l = N_radial(nl,ll)
    Call GaussLaguerreWX((Abs(li)+Abs(lj)+Abs(lk)+Abs(ll))*0.5_pr,wlag,xlag)
    Do i = 1,ngl
       Call LaguerreL(ni,Real(Abs(li),Kind=pr),xlag(i)/2.0_pr,Lni)
       Call LaguerreL(nj,Real(Abs(lj),Kind=pr),xlag(i)/2.0_pr,Lnj)
       Call LaguerreL(nk,Real(Abs(lk),Kind=pr),xlag(i)/2.0_pr,Lnk)
       Call LaguerreL(nl,Real(Abs(ll),Kind=pr),xlag(i)/2.0_pr,Lnl)
       ME = ME + wlag(i)*Lni*Lnj*Lnk*Lnl
    End Do
    ME = ME*Nr_i*Nr_j*Nr_k*Nr_l*bp**2*2._pr**(-(Abs(li)+Abs(lj)+Abs(lk)+Abs(ll)+4)*0.5_pr)/(2*pi)

  End Function MatrixElement_rZR
  !======================================================================
  !> Given a pair of radial and angular quantum numbers, this function calculates the normalization
  !> constant of a two-dimensional harmonic oscillator wave function in radial coordinates. This
  !> function is used for debugging only and is not used in the actual DFT calculations.
  !>
  !> @result \f$ N_{n_r,\Lambda} = \frac{1}{b_\perp} \left[
  !>             \frac{2n_r!}{(n_r+|\Lambda|)!} \right]^{1/2} \f$
  !======================================================================
  Function N_radial(nr,l) result(N_r)
    Use HFBTHO, Only : bp

    Integer(ipr), Intent(In) :: nr !< - radial main quantum number
    Integer(ipr), Intent(In) :: l  !< - radial orbital quantum number
    Real(pr) :: N_r

    Real(pr) :: nrf, nrlf

    nrf  = factrl(nr)
    nrlf = factrl(nr+Abs(l))
    N_r = Sqrt(2*nrf/nrlf)/bp

  End Function N_radial
  !======================================================================
  !> Given a quantum number, this function calculates the normalization constant of a one-dimensional
  !> harmonic oscillator wave function. This function is used for debugging only and is not used in the
  !> actual DFT calculations.
  !>
  !> @result \f$ N_{n_z} = \frac{1}{(b_z\sqrt{\pi}2^{n_z}n_z!)^{1/2}}\f$
  !======================================================================
  Function N_axial(nz) result(N_z)
    Use HFBTHO, Only : bz

    Integer(ipr), Intent(In) :: nz !< - axial quanutum number
    Real(pr) :: N_z

    Real(pr) :: nzf

    nzf  = factrl(nz)
    N_z = one/Sqrt(bz*Sqrt(pi)*(2**nz)*nzf)

  End Function N_axial
  !======================================================================
  !> This routine calculates the generalized Laguerre polynomial \f$L_n^\alpha(x)\f$ along with
  !> \f$ \frac{\partial}{\partial x} L_n^\alpha(x) \f$ and \f$L_{n-1}^\alpha(x)\f$. The last two
  !> are optional. Note that this routine duplicates hfbtho_gauss.gaupol().
  !======================================================================
  Subroutine LaguerreL(n,alpha,x,Ln,Lnp,Lnm1)
    Integer(ipr), Intent(In) :: n !< - Order of the polynomial
    Real(pr), Intent(In) :: alpha !< - Parameter \f$ \alpha \f$ of the generalized polynomial \f$ L_{n}^{\alpha} \f$
    Real(pr), Intent(In) :: x !< - Value where the polynomials are evaluated
    Real(pr), intent(out) :: Ln !< - Generalized Laguerre Polynomial \f$ L_n^\alpha(x) \f$
    Real(pr), optional, intent(out) :: Lnp  !< - Derivative of Ln with respect of x
    Real(pr), optional, intent(out) :: Lnm1 !< - Generalized Laguerre Polynomial \f$ L_{n-1}^\alpha(x) \f$

    Real(pr) :: Ljm2, Ljm1, Lj
    Integer(ipr) :: j

    Lj = one
    Ljm1 = zero
    Do j = 1,n
       Ljm2 = Ljm1
       Ljm1 = Lj
       Lj = ((-x+2*j-1+alpha)*Ljm1-(j-1+alpha)*Ljm2)/real(j,kind=pr)
    End Do
    Ln = Lj
    If(present(Lnp))  Lnp  = (n*Ln-(n+alpha)*Ljm1)/x
    If(present(Lnm1)) Lnm1 = Ljm1

  End Subroutine LaguerreL
  !======================================================================
  !> This routine calculates the weights w and nodes x for the Gauss-Laguerre quadrature integration
  !> with a given parameter \f$ \alpha \f$. Note that this routine duplicates hfbtho_gauss.gausspoints().
  !======================================================================
  Subroutine GaussLaguerreWX(alfa,w,x)
    Real(pr), Intent(In) :: alfa !< - \f$\alpha\f$ parameter on the weight function
    Real(pr), intent(out), Dimension(:) :: x !< - Integration nodes
    Real(pr), intent(out), Dimension(:) :: w !< - Integration weights

    Integer(ipr), Parameter :: maxit = 20
    Integer(ipr) :: i,its,n,ai
    Real(pr), Parameter :: eps = 2.0e-14_pr
    Real(pr) :: z,z1,Ln,Lnp,Lnm1

    n = size(w)
    If(size(w) /= size(x)) Then !Check that both arrays have the same size
       Write(*,*) 'Arrays w and x must have the same size in GaussLaguerreWX'
       Stop
    End If
    Do i = 1,n
       !intial guess for every root of the Laguerre Polynomial
       If(i == 1) Then
          z=(1+alfa)*(3+0.92_pr*alfa)/(1+2.4_pr*n+1.8_pr*alfa)
       Else If(i == 2) Then
          z=z+(15+6.25_PR*alfa)/(1+0.9_pr*alfa+2.5_pr*n)
       Else
          ai=i-2
          z=z+((1+2.55_pr*ai)/(1.9_pr*ai)+1.26_pr*ai*alfa/(1+3.5_pr*ai))&
               *(z-x(i-2))/(1+0.3_pr*alfa)
       End If
       ! Newton method to refine the roots
       Do its = 1, maxit
          Call LaguerreL(n,alfa,z,Ln,Lnp,Lnm1)
          z1 = z
          z = z1 - Ln/Lnp
          If(Abs(z-z1) <= eps*z) Exit
       End Do
       If(its==maxit+1) write(*,*) 'maxit exceeded in GaussLaguerreWX',its,Abs(z-z1)/z,z
       ! Save root and calculate weight
       x(i) = z
       If(alfa == 0._pr) Then
          w(i) = -1/(n*Lnm1*Lnp)
       Else
          w(i) = -gamma(n+alfa)/(factrl(n)*Lnm1*Lnp)
       End If
    End Do

  End Subroutine GaussLaguerreWX
  !======================================================================
  !> This routine calculates the Hermite polynomial \f$H_n(x)\f$. It duplicates hfbtho_gauss.gaupol().
  !>
  !> @result \f$ H_n(x) = \left(2x-\frac{d}{dx} \right)^n \cdot 1 \f$
  !======================================================================
  Function HermiteH(n,x) result(Hn)
    Integer(ipr), Intent(In) :: n !< - Order of the polynomial
    Real(pr), Intent(In) :: x !< - Value where the polynomial is evaluated
    Real(pr) :: Hn  !< - Hermite Polynomial \f$ H_n(x) \f$

    Real(pr) :: Hjm2, Hjm1, Hj
    Integer(ipr) :: j

    Hj = one
    Hjm1 = zero
    Do j = 1,n
       Hjm2 = Hjm1
       Hjm1 = Hj
       Hj = two*x*Hjm1-two*Real(j-1,Kind=pr)*Hjm2
    End Do
    Hn = Hj

  End Function HermiteH
  !======================================================================
  !>
  !======================================================================
End Module HFBTHO_Gogny
