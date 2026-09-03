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
!                      GAUSS QUADRATURE PACKAGE                        !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module provides a set of routines to compute the weights and
!> nodes of several standard quadrature schemes. It also contains
!> several routines specifically related to Hermite polynomials.
!>
!> @authors
!> Mario Stoitsov, Nicolas Schunck
!-------------------------------------------------------------------
!   Subroutines: - gaupol(lpr)
!                - coordinateLST_HO()
!                - optHFBTHO_HO()
!                - gausspoints()
!                - Gaussq(kindi,n,alpha,beta,kpts,endpts,b,t,w)
!                - Class(kindi,N,ALPHA,BETA,B,A,MUZERO)
!                - GBTQL2(N,D,E,Z,IERR)
!                - D_HERM(X,N,HER,DHER,NDIM)
!                - DEVHER(NOSACT)
!   Functions: - GBSLVE(SHIFT,N,A,B)
!              - pr_gamma(x)
!----------------------------------------------------------------------!

Module HFBTHO_Gauss

  Use HFBTHO_utilities

  Implicit None

  Public gaupol,coordinateLST_HO,gausspoints,gauleg,D_HERM,DEVHER

  Private optHFBTHO_HO,Gaussq,GBSLVE,Class,GBTQL2,pr_gamma

  ! Public variables
  Integer(ipr), Public, Save :: ngh  !< Number of Gauss-Hermite integration points
  Integer(ipr), Public, Save :: ngl  !< Number of Gauss-Laguerre integration points
  Integer(ipr), Public, Save :: nghl !< Number of Gauss integration points, nghl = ngh * ngl
  Integer(ipr), Public, Save :: nleg !< Number of Gauss-Legendre integration points
  Real(pr), Allocatable, Public, Save :: xh(:)     !< Dimensionless nodes of the Gauss-Hermite quadrature \f$ \xi_h \f$
  Real(pr), Allocatable, Public, Save :: wh(:)     !< Weighted weights of the Gauss-Hermite quadrature \f$ w_h = \omega_h e^{-\xi_h^2} \f$
  Real(pr), Allocatable, Public, Save :: xl(:)     !< Dimensionless nodes of the Gauss-Laguerre quadrature \f$ \eta_l \f$
  Real(pr), Allocatable, Public, Save :: sxl(:)    !< Square roots of the nodes of the Gauss-Laguerre quadrature \f$ \sqrt{\eta_l} \f$
  Real(pr), Allocatable, Public, Save :: wl(:)     !< Weighted weights of the Gauss-Laguerre quadrature \f$ w_l = \omega_{l} e^{-\eta_l} \f$
  Real(pr), Allocatable, Public, Save :: xleg(:)   !< Dimensionless nodes of the Gauss-Legendre quadrature \f$ \mu_{leg} \f$
  Real(pr), Allocatable, Public, Save :: wleg(:)   !< Weights of the Gauss-Legendre quadrature \f$ w_{leg} \f$
  Real(pr), Allocatable, Public, Save :: fh(:)     !< Nodes (with dimension) of the Gauss-Hermite quadrature \f$ z_h = \beta_z \xi_h \f$
  Real(pr), Allocatable, Public, Save :: fl(:)     !< Nodes (with dimension of the Gauss-Laguerre quadrature \f$ z_l = \beta_{\perp} \sqrt{\eta_l} \f$
  Real(pr), Allocatable, Public, Save :: wdcor(:)  !< Normalization factor for quadrature rules \f$ \pi b_z b_{\perp}^{2} w_h w_l \f$
  Real(pr), Allocatable, Public, Save :: wdcori(:) !< Inverse of the normalization factor wdcori = 1/wdcor
  Real(pr), Allocatable, Public, Save :: COEF00(:,:,:),HERFAC(:)
  Real(pr), Allocatable, Public, Save, Target :: qh(:,:),qh1(:,:),ql(:,:,:),ql1(:,:,:)

  ! Private variables
  Integer(ipr), PRIVATE :: debug_gauss = 0

Contains
  !=======================================================================
  !> Subroutine gaupol() computes the harmonic oscillator wave functions
  !> and their first derivatives in cylindrical coordinates \f$ (\rho,\theta,z )\f$. The full HO
  !> wave function reads
  !>  \f[
  !>     \psi_{n_{r}\Lambda n_{z}}(\rho, \theta, z) =
  !>     \psi_{n_{r}}^{|\Lambda|}(\eta) \psi_{n_{z}}(\xi) \frac{e^{i\Lambda\theta}}{\sqrt{2\pi}}
  !>                  \chi_{\Sigma}(\sigma),
  !>  \f]
  !>  where the dimensionless coordinates are \f$ \eta = \beta_{\perp}^{2}\rho^2 \f$,
  !>  \f$ \xi = \beta_{z}z \f$, \f$ \Lambda \f$ is the eigenvalue of \f$ \hat{\ell}_{z} \f$,
  !>  \f$ \Sigma \f$ the eigenvalue of \f$ \hat{s}_{z} \f$ and \f$ \sigma \f$ is a number
  !>  that labels the spin eigenstate in coordinate space, by opposition to \f$ \chi_{\Sigma} \f$,
  !>  which refers to a vector in spin space. Under the action of time-reversal, the HO
  !>  basis states turn into \f$ \psi_{n_{r}-\Lambda n_{z}}(\rho, \theta, z) \f$, where
  !>  the spatial part proportional to \f$ \rho, z \f$ is unchanged.
  !>
  !>  For the z-coordinate, the HO wave function reads
  !>     \f[
  !>        \psi_{n_{z}}(z) = \mathcal{N}_{n_{z}}\beta_{z}^{1/2} e^{-\xi^{2}/2} H_{n_{z}}(\xi),
  !>     \f]
  !>  with the normalization factor
  !>     \f[
  !>        \mathcal{N}_{n_{z}}
  !>      = \left( \frac{1}{\sqrt{\pi}2^{n_{z}}n_{z}!} \right)^{1/2}.
  !>     \f]
  !>  The first derivative of the HO wave function is
  !>     \f[
  !>        \frac{\partial\psi_{n_{z}}}{\partial z} = \beta_{z}
  !>            \left[ \sqrt{2n_{z}} \psi_{n_{z}-1}(\xi) - \xi \psi_{n_{z}}(\xi) \right].
  !>     \f]
  !>  In the code, variable \f$ \mathtt{qh} \f$ contains \f$ (1/\beta^{1/2}_{z})\psi_{n_{z}} \f$
  !>  on the nodes of the Gauss-Hermite integration mesh; variable \f$ \mathtt{qh1} \f$
  !>  contains the value of \f$ (1/\beta^{3/2}_{z})\partial\psi_{n_{z}}/\partial z \f$ on that same
  !>  mesh. Note that this is consistent with a dimension analysis: \f$ \psi_{n_{z}}(z) \f$
  !>  has dimension \f$ [L]^{-1/2} \f$ so that the integral over space of the square
  !>  \f$ |\psi_{n_{z}}(z)|^2 \f$ is dimensionless and can be interpreted as a probability.
  !>  Similarly, the dimension of the derivative \f$ \partial\psi_{n_{z}}/\partial z \f$ is
  !>  \f$ [L]^{-3/2} \f$. In the code, objects \f$ \mathtt{qh} \f$ and \f$ \mathtt{qh1} \f$
  !>  are thus dimensionless.
  !>
  !>  For the \f$\rho\f$-coordinate, the HO wave function reads
  !>     \f[
  !>        \psi_{n_{r}}^{\Lambda}(\eta) =
  !>        \mathcal{N}_{n_{r}}^{\Lambda} \beta_{\perp}\sqrt{2}
  !>              \eta^{|\Lambda|/2} e^{-\eta/2} L_{n_{r}}^{|\Lambda|}(\eta),
  !>     \f]
  !>  with the normalization factor
  !>     \f[
  !>        \mathcal{N}_{n_{r}}^{\Lambda}
  !>      = \left( \frac{n_{r}}{ (n_{r} + |\Lambda|)!} \right)^{1/2}.
  !>     \f]
  !>  The first derivative of the HO wave function is
  !>     \f[
  !>        \frac{\partial\psi_{n_{r}}^{\Lambda}}{\partial\rho} = \frac{\beta_{\perp}}{\sqrt{\eta}}
  !>            \left[ (2n_{r} + |\Lambda| - \eta)\psi_{n_{r}}^{\Lambda}(\eta)
  !>                   -2\sqrt{n_{r}(n_{r}+|\Lambda|)}\psi_{n_{r}-1}^{\Lambda}(\eta) \right].
  !>     \f]
  !>  In the code, variable \f$ \mathtt{ql} \f$ contains the value of
  !>  \f$ (1/\beta_{\perp}\sqrt{2})\psi_{n_{r}}^{\Lambda} \f$ on the nodes of the Gauss-Laguerre
  !>  integration mesh; variable \f$ \mathtt{ql1} \f$ contains the value of
  !>  \f$ (\sqrt{\eta}/\beta_{\perp}^{2}\sqrt{2}) \partial\psi_{n_{r}}^{\Lambda}/\partial \rho \f$
  !>  on that same mesh.
  !=======================================================================
  Subroutine gaupol(lpr)
    Use math, Only : sq,sqi,wfi
    Use HFBTHO, Only : do_print,nzm,nlm,nrm

    Logical, Intent(In) :: lpr

    Real(pr) :: w0,z,x,s,s0,s1,w00,w4pii,dsq,d1,d2,d3,d4,hs0,hs1
    Integer(ipr) :: ih,il,iw,ix,n,l,n1,n2

    !-----------------------------------------------
    ! z-coordinate
    !-----------------------------------------------
    w4pii=pi**(-0.250_pr)
    Do ih=1,ngh
       z=xh(ih); w0=w4pii*Exp(-half*z*z)
       w0 = w0*Sqrt(wh(ih))
       qh(0,ih)=w0;       qh(1,ih)=sq(2)*w0*z
       qh1(0,ih)=-w0*z;   qh1(1,ih)=sq(2)*w0*(one-z*z)
       Do n=2,nzm
          qh(n,ih)=sqi(n)*(sq(2)*z*qh(n-1,ih)-sq(n-1)*qh(n-2,ih))
          qh1(n,ih)=sq(n+n)*qh(n-1,ih)-z*qh(n,ih)
       End Do
    End Do
    !-----------------------------------------------
    ! perpendicular coordinate
    !-----------------------------------------------
    Do il=1,ngl
       x=xl(il); w00=sq(2)*Exp(-half*x)
       Do l=0,nlm
          w0=w00*Sqrt(half*wl(il)*x**l)
          ql(0,l,il)=wfi(l)*w0;         ql(1,l,il)=(l+1-x)*wfi(l+1)*w0
          ql1(0,l,il)=(l-x)*wfi(l)*w0;  ql1(1,l,il)=(Real(l*l+l,Kind=pr) &
                                                  -x*Real(l+l+3,Kind=pr)+x*x)*wfi(l+1)*w0
          Do n=2,nrm
             dsq=sq(n)*sq(n+l); d1=Real(n+n+l-1,Kind=pr)-x
             d2=sq(n-1)*sq(n-1+l); d3=n+n+l-x; d4=two*dsq
             ql(n,l,il)=(d1*ql(n-1,l,il)-d2*ql(n-2,l,il))/dsq
             ql1(n,l,il)=d3*ql(n,l,il)-d4*ql(n-1,l,il)
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Hermite orthonormalization
    !-----------------------------------------------
    hs0=zero; hs1=two
    Do n1=0,nzm
       Do n2=0,n1
          If(Mod(n1-n2,2) == 0) Then
             s=zero
             Do ih=1,ngh
                s=s+qh(n1,ih)*qh(n2,ih)
             End Do
             If(n1 /= n2) Then
                hs0=Max(s,hs0)
             Else
                hs1=Min(s,hs1)
             End If
          End If
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Laguerre orthonormalization
    !-----------------------------------------------
    s0=zero; s1=two
    Do l=0,nlm
       Do n1=0,nrm
          Do n2=0,n1
             s=zero
             Do il=1,ngl
                s=s+ql(n1,l,il)*ql(n2,l,il)
             End Do
             If(n1 /= n2) Then
                s0=Max(s,s0)
             Else
                s1=Min(s,s1)
             End If
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! print accuracy
    !-----------------------------------------------
    If(lpr .And. do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a)')  '            Integration Meshes           '
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a,i3)') '  Number of Gauss-Hermite mesh points ngh ....: ',ngh
          Write(iw,'(a,i3)') '  Number of Gauss-Laguerre mesh points ngl ...: ',ngl
          Write(iw,'(a,i3)') '  Number of Gauss-Legendre mesh points nleg ..: ',nleg
          Write(iw,'(a)') '  Integration boundaries'
          Write(iw,'(2(a,f12.8))')  &
               '    Hermite  - from xh(1)  =',xh(1),  ' to xh(ngh)   =',xh(ngh)
          Write(iw,'(2(a,f12.8))')  &
               '    Laguerre - From xl(1)  =',xl(1),  ' to xl(ngl)   =',xl(ngl)
          If(nleg > 0) Then
             Write(iw,'(2(a,f12.8))')  &
                  '    Legendre - From xleg(1)=',xleg(1),' to xleg(nleg)=',xleg(nleg)
          End If
          Write(iw,*) ' Max.dev.in:     Orthogonality            Normalization'
          Write(iw,*) ' Hermite  ',hs0,Abs(one-hs1)
          Write(iw,*) ' Laguerre ',s0,Abs(one-s1)
       End Do
    End If
    !-----------------------------------------------
    ! debug
    !-----------------------------------------------
    If (lpr .And. debug_gauss >= 4 .And. do_print == 1) Then
       ix=3
       Do iw=lout,lfile
          Write(iw,*) ' nz    qh(nz,ih=1,...)'
          Do n=0,nzm
             Write(iw,'(i4,3f15.8)') n,(qh(n,ih),ih=1,ix)
             Write(iw,'(i4,3f15.8)') n,(qh1(n,ih),ih=1,ix)
             Write(iw,*) ' '
          End Do
          Do l=0,nlm
             Write(iw,*) ' nr ml    ql(nr,l,il=1,...)'
             Do n=0,nrm
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql(n,l,il),il=1,ix)
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql1(n,l,il),il=1,ix)
                Write(iw,*) ' '
             End Do
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Hermite polynomials normalization
       !-----------------------------------------------
       Do n1=0,nzm
          Do n2=0,n1
             If(Mod(n1-n2,2) == 0) Then
                s=zero
                Do ih=1,ngh
                   s=s+qh(n1,ih)*qh(n2,ih)
                End Do
                Do iw=lout,lfile
                   Write(iw,'(" Gauss-Hermite: n1=",i3,"  n2=",i3,f20.8)') n1,n2,s
                End Do
             End If
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Laguerre polynomials normalization
       !-----------------------------------------------
       Do l=0,nlm
          Do n1=0,nrm
             Do n2=0,n1
                s=zero
                Do il=1,ngl
                   s=s+ql(n1,l,il)*ql(n2,l,il)
                End Do
                Do iw=lout,lfile
                   Write(iw,'(" Gauss Laguerre: l=",i2," n1=",i3,"  n2=",i3,f20.8)') l,n1,n2,s
                End Do
             End Do
          End Do
       End Do
    End If

  End Subroutine gaupol
  !====================================================================
  !> Defines and stores the weights and nodes of integration in the
  !> dimensionless coordinates \f$ \xi\f$ and \f$ \eta \f$.
  !====================================================================
  Subroutine coordinateLST_HO()
    Use HFBTHO, Only : ierror_flag,bp,bz

    Integer(ipr) :: il,ih,i

    Do il=1,ngl
       Do ih=1,ngh
          i=ih+(il-1)*ngh
          fh(i)=bz*xh(ih)
          fl(i)=bp*Sqrt(xl(il))
          wdcor(i)=pi*wh(ih)*wl(il)*bz*bp*bp
          wdcori(i)=one/wdcor(i)
       End Do
    End Do

    Call optHFBTHO_HO() ! optimal HO/THO combinations
    If(ierror_flag /= 0) Return

  End Subroutine coordinateLST_HO
  !=======================================================================
  !> This routine defines and stores the derivatives and Laplacian of basis functions
  !>   - \f$ \mathtt{FI2D_opt(JA,ihil)} \f$ contains the Laplacian(r,z) HOwf
  !>   - \f$ \mathtt{FID2D-xlamy2*FID } \f$ contains the Laplacian(r,z,phy) FID
  !>   - \f$ \mathtt{FIU2D-xlapy2*FIU } \f$ contains the Laplacian(r,z,phy) FIU
  !=======================================================================
  Subroutine optHFBTHO_HO()
    Use HFBTHO, Only : bp,bz,ntx,NB,ID,ia,nl,nr,nz,ns,QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt

    Integer(ipr) :: ih,il,ihil,ib,nd,nza,nra,nla,nsa,im,ja,n1
    Real(pr) :: qla,v2,v4,yi,y,y2,fi1r,fi1z,fi2d,qha,qhla,qhl1a,qh1la,cnraa,cnzaa,sml2
    Real(pr) :: bpi,bpi2,bzi,bzi2,xh2

    bpi=one/bp; bpi2=bpi*bpi; bzi=one/bz; bzi2=bzi*bzi

    ! Allocate the optimization arrays
    If(Allocated(QHLA_opt)) Deallocate(QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt)
    Allocate(QHLA_opt(ntx,nghl),FI1R_opt(ntx,nghl),FI1Z_opt(ntx,nghl),FI2D_opt(ntx,nghl),y_opt(nghl))

    ! Start blocks
    Do ib=1,NB
       ND=ID(ib); IM=ia(ib)
       ! Sum over Gauss integration points
       Do IL=1,ngl
          v2=half/xl(il); v4=v2*v2
          Do IH=1,ngh
             ihil=ih+(il-1)*ngh; xh2=xh(ih)**2
             ! HO-basis
             yi=Sqrt(xl(il))*bp; y=one/yi; y2=y*y
             y_opt(ihil)=y
             ! Scan over basis states
             Do N1=1,ND
                ja=n1+im; nla=nl(ja); nra=nr(ja); nza=nz(ja); nsa=ns(ja)
                sml2 = nla*nla; cnzaa = nza+nza+1; cnraa = nra+nra+nla+1
                qha = qh(nza,ih); qla = ql(nra,nla,il)
                qhl1a = qha*ql1(nra,nla,il)*v2; qh1la = qh1(nza,ih)*qla
                ! ho-basis
                qhla = qha*qla ! wave function f_i
                fi1r = (two*sqrt(xl(il))*bpi)*qhl1a ! derivative with respect to , \nabla_r f_i
                fi1z = bzi*qh1la ! derivative with respect to z, \nabla_z f_i
                fi2d = ((xh2-cnzaa)*bzi2+four*(p14-cnraa*v2+sml2*v4)*xl(il)*bpi2 )*qhla ! laplacian, \delta f_i
                qhla_opt(ja,ihil)=qhla; fi2d_opt(ja,ihil)=fi2d; fi1r_opt(ja,ihil)=fi1r; fi1z_opt(ja,ihil)=fi1z
             End Do !N1
          End Do !IH
       End Do !IL
    End Do !IB

  End Subroutine optHFBTHO_HO
  !=======================================================================
  !> The routine determines the points and weights for Gauss quadratures
  !> in the cases of Gauss-Legendre, -Laguerre and -Hermite formulas.
  !> Note that for Gauss-Hermite and Gauss-Laguerre integrations, the
  !> weights are defined in such a way that they contain the weight function
  !> associated with the orthohonal polynomials. For example, in the case
  !> of Gauss-Hermite quadrature, we may write
  !>
  !>  \f{align}{
  !>      \int dx f(x)g(x) & = \int dx e^{-x^2} f(x)e^{x^2} g(x)e^{x^2} \\
  !>                       & = \sum_{i} w_i f(x_i)e^{x_i^2} g(x_i)e^{x_i^2} \\
  !>                       & = \sum_{i} w_i e^{x_i^2} f(x_i)g(x_i)
  !>  \f}
  !>
  !> The weights are thus multiplied by \f$ e^{x_i^2} \f$.
  !=======================================================================
  Subroutine gausspoints()
    Use HFBTHO, Only : Parity,ierror_flag

    Real(pr):: al,be,sparity
    Real(pr), Allocatable :: endpts(:),b(:),t(:),w(:)
    Integer(ipr) :: N,i,j,KINDI,kpts,nparity

    al=0.0_pr; be=0.0_pr; kpts=0

    !--------------------------------------------------------------------
    !------------------>> Gauss-Hermite (positive nodes) <<--------------
    !--------------------------------------------------------------------
    If(Parity) Then
       KINDI=4; N=2*ngh ! Parity conserved
       nparity=ngh; sparity=two
    Else
       KINDI=4; N=ngh   ! Parity not conserved
       nparity=0;   sparity=one
    End If
    Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
    Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
    If(ierror_flag /= 0) Return
    Do i=nparity+1,N
       j=i-nparity
       xh(j)=t(i)
       ! Build in the Gaussian weight function into the weights wh
       wh(j)=sparity*Exp(xh(j)*xh(j)+Log(w(i)))
    End Do
    Deallocate(endpts,b,t,w)
    !--------------------------------------------------------------------
    !---------------------------->> Gauss-Laguerre <<--------------------|
    !--------------------------------------------------------------------
    KINDI=6; N=ngl
    Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
    Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
    If(ierror_flag /= 0) Return
    Do j=1,ngl
       xl(j)=t(j)
       ! Build in the exponential weight function into the weights wl
       wl(j)=Exp(xl(j)+Log(w(j)))
       sxl(j)=Sqrt(xl(j))
    End Do
    Deallocate(endpts,b,t,w)
    !--------------------------------------------------------------------
    !----------------->> Gauss-Legendre (positive nodes) <<--------------
    !--------------------------------------------------------------------
    If(nleg > 0) Then
       KINDI=1; N=2*nleg
       Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
       Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
       If(ierror_flag /= 0) Return
       Do j=1,nleg
          i=nleg+j
          xleg(j)=t(i); wleg(j)=w(i)
       End Do
       Deallocate(endpts,b,t,w)
    End If

  End Subroutine gausspoints
  !=======================================================================
  !> This set of routines computes the nodes t(j) and weights w(j) for
  !> Gaussian-type quadrature rules with pre-assigned nodes. These are
  !> used when one wishes to approximate
  !>
  !>   \f[
  !>       \int_{a}^{b}  f(x) w(x) dx \approx \sum_{j=1}^{n} w_j f(t_j)
  !>   \f]
  !>
  !> (note w(x) and w(j) have no connection with each other). Here w(x)
  !> is one of six possible non-negative weight functions (listed below),
  !> and f(x) is the function to be integrated. Gaussian quadrature is
  !> particularly useful on infinite intervals (with appropriate weight
  !> functions), since then other techniques often fail. Associated with
  !> each weight function w(x) is a set of orthogonal polynomials. The
  !> nodes t(j) are just the zeroes of the proper n-th degree polynomial.
  !>
  !> Inputs (all real numbers are in double precision)
  !>  - kindi: an integer between 1 and 6 giving the type of quadrature
  !>       * 1:  Legendre quadrature, w(x) = 1 on \f$ [-1, 1] \f$
  !>       * 2:  Chebyshev quadrature of the first kind
  !>             \f$ w(x) = 1/\sqrt(1 - x^2) \f$ on \f$ [-1, +1] \f$
  !>       * 3:  Chebyshev quadrature of the second kind
  !>             \f$ w(x) = \sqrt(1 - x^2) \f$ on \f$ [-1, 1] \f$
  !>       * 4:  Hermite quadrature, \f$ w(x) = \exp(-x^2) \f$ on
  !>             \f$ ]-\infty, +\infty[ \f$
  !>       * 5:  Jacobi quadrature, \f$ w(x) = (1-x)^{\alpha}(1+x)^{\beta}\f$
  !>             on \f$ [-1, 1] \f$, \f$ \alpha, \beta > -1 \f$.
  !>             Note: kind=2 and 3 are a special case of this.
  !>       * 6:  generalized Laguerre quadrature,
  !>             \f$ w(x) = \exp(-x) x^{\alpha} \f$ on
  !>             \f$ [0, +\infty[ \f$, with \f$ alpha > -1 \f$.
  !>  - n: the number of points used for the quadrature rule
  !>  - alpha: real parameter used only for Gauss-Jacobi and Gauss-Laguerre quadrature (otherwise use 0.d0).
  !>  - beta: real parameter used only for Gauss-Jacobi quadrature(otherwise use 0.d0)
  !>  - kpts: (integer) normally 0, unless the left or right end-
  !>          point (or both) of the interval is required to be a
  !>          node (this is called Gauss-Radau or Gauss-Lobatto
  !>          quadrature). Then kpts is the number of fixed
  !>          endpoints (1 or 2).
  !>  - endpts: real array of length 2. Contains the values of any fixed endpoints, if kpts = 1 or 2.
  !>  - b: real scratch array of length n
  !>
  !> Outputs (both double precision arrays of length n)
  !>  - t: the desired nodes.
  !>  - w: the desired weights w(j).
  !>
  !> NOTE: Underflow may sometimes occur, but is harmless.
  !>
  !> Adapted from GO library at www.netlib.org
  !=======================================================================
  Subroutine Gaussq(kindi,N,alpha,beta,kpts,endpts,b,t,w)
    Integer(ipr), Intent(In) :: N,kindi,kpts
    Real(pr), Intent(In) :: alpha,beta
    Real(pr), Dimension(2), Intent(In) :: endpts
    Real(pr), Dimension(n), Intent(Inout) :: b,t,w

    Integer(ipr) :: ierr
    Real(pr) :: MUZERO,GAM,T1

    Call Class(kindi,n,alpha,beta,b,t,muzero)

    If(KPTS == 0) Then
       w = 0.0_pr; w(1) = 1.0_pr
       Call GBTQL2(n,t,b,w,ierr)
       w = MUZERO*w*w
       Return
    End If

    If(KPTS == 2) Then
       GAM=GBSLVE(endpts(1),n,t,b)
       T1=((endpts(1)-endpts(2))/(GBSLVE(endpts(2),N,T,B)-GAM))
       B(n-1)=Sqrt(T1)
       T(n) = endpts(1)+GAM*T1
       w = 0.0_pr; w(1) = 1.0_pr
       Call GBTQL2(n,t,b,w,ierr)
       w = MUZERO*w*w
       Return
    End If

    T(n) = GBSLVE(endpts(1),n,t,b)*B(N-1)**2+endpts(1)
    w = 0.0_pr; w(1) = 1.0_pr
    Call GBTQL2(n,t,b,w,ierr)
    w = MUZERO*w*w

  End Subroutine Gaussq
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function GBSLVE(SHIFT,N,A,B)
    Integer(ipr), Intent(In) :: N
    Real(pr), Intent(In) :: SHIFT
    Real(pr), Dimension(N), Intent(In) :: A,B

    Integer(ipr) :: NM1,i
    Real(pr) :: ALPHA

    ALPHA=A(1)-SHIFT
    NM1=N-1
    Do I=2,NM1
       ALPHA=A(I)-SHIFT-B(I-1)**2/ALPHA
    End Do
    GBSLVE=1.0_pr/ALPHA

  End Function GBSLVE
  !=======================================================================
  !>  This procedure supplies the coefficients a(j), b(j) of the
  !>  recurrence relation
  !>
  !>   \f[
  !>       b_j p_j(x) = (x-a_j) p_{j-1}(x) - b_{j-1} p_{j-2}(x)
  !>   \f]
  !>
  !>  for the various classical (normalized) orthogonal polynomials,
  !>  and the zero-th moment
  !>
  !>   \f[
  !>       \mu_{0} = \int w(x) dx
  !>   \f]
  !>
  !>  of the given polynomial's weight function w(x). Since the
  !>  polynomials are orthonormalized, the tridiagonal matrix is
  !>  guaranteed to be symmetric.
  !>
  !>  The input parameter alpha is used only for Laguerre and Jacobi
  !>  polynomials, and the parameter beta is used only for Jacobi
  !>  polynomials. The Laguerre and Jacobi polynomials require the Gamma
  !>  function.
  !>
  !>  Adapted from GO library at www.netlib.org
  !=======================================================================
  Subroutine Class(kindi,N,ALPHA,BETA,B,A,MUZERO)
    Integer(ipr), Intent(In) :: kindi,N
    Real(pr), Intent(In) :: ALPHA,BETA
    Real(pr), Intent(Inout) :: MUZERO
    Real(pr), Dimension(N), Intent(Inout) :: A,B

    Integer(ipr) :: i,NM1
    Real(pr) :: ABI,DI20,AB,A2B2,FI

    NM1=N-1
    Select Case (kindi)
    Case (1) ! Legendre polynomials
       MUZERO=2.0_pr
       Do I=1,NM1
          A(I)=0.0_pr
          ABI=Real(I,Kind=pr)
          B(I)=ABI/Sqrt(4.0_pr*ABI*ABI-1.0_pr)
       End Do
       A(N)=0.0_pr
    Case (2) ! Chebyshev polynomials of the first kind
       MUZERO=PI
       Do I=1,NM1
          A(I)=0.0_pr
          B(I)=0.50_pr
       End Do
       B(1)=Sqrt(0.50_pr)
       A(N)=0.0_pr
    Case (3) ! Chebyshev polynomials of the second kind
       MUZERO=PI/2.0_pr
       Do I=1,NM1
          A(I)=0.0_pr
          B(I)=0.50_pr
       End Do
       A(N)=0.0_pr
    Case (4) ! Hermite polynomials
       MUZERO=Sqrt(PI)
       Do I=1,NM1
          A(I)=0.0_pr
          DI20=I/2.0_pr
          B(I)=Sqrt(DI20)
       End Do
       A(N)=0.0_pr
    Case (5) ! Jacobi polynomials
       AB=ALPHA+BETA
       ABI=2.0_pr+AB
       !MUZERO=2.0_pr**(AB+1.0_pr)*pr_gamma(ALPHA+1.0_pr)*pr_gamma(BETA+1.0_pr)/pr_gamma(ABI)
       MUZERO=Exp((AB+one)*Log(two)+Log(pr_gamma(ALPHA+one))+Log(pr_gamma(BETA+one))-Log(pr_gamma(ABI)))
       A(1)=(BETA-ALPHA)/ABI
       B(1)=Sqrt(4.0_pr*(1.0_pr+ALPHA)*(1.0_pr+BETA)/((ABI+1.0_pr)*ABI*ABI))
       A2B2=BETA*BETA-ALPHA*ALPHA
       Do I=2,NM1
          ABI=2.0_pr*I+AB
          A(I)=A2B2/((ABI-2.0_pr)*ABI)
          FI=I
          B(I)=Sqrt(4.0_pr*FI*(FI+ALPHA)*(FI+BETA)*(FI+AB)/((ABI*ABI-1.0_pr)*ABI*ABI))
       End Do
       ABI=2.0_pr*N+AB
       A(N)=A2B2/((ABI-2.0_pr)*ABI)
    Case (6) ! Laguerre polynomials
       MUZERO=pr_gamma(ALPHA+1.0_pr)
       Do I=1,NM1
          FI=I
          A(I)=2.0_pr*FI-1.0_pr+ALPHA
          B(I)=Sqrt(FI*(FI+ALPHA))
       End Do
       A(N)=2.0_pr*N-1.0_pr+ALPHA
    Case default
    End Select

  End Subroutine Class
  !=======================================================================
  !
  !=======================================================================
  Subroutine GBTQL2(N,D,E,Z,IERR)
    Integer(ipr), Intent(In) :: N
    Integer(ipr), Intent(Inout) :: IERR
    Real(pr), Dimension(N), Intent(Inout) :: D,E,Z

    Integer(ipr) :: I,J,K,L,M,II,MML
    Real(pr) :: MACHEP,P,G,R,S,C,F,B

    !MACHEP=16.0_pr**(-14)
    MACHEP=epsilon(1.0_pr)
    IERR=0
    If(N == 1) Return
    E(N)=0.0_pr
    Do L= 1,N
       J=0
       Do
          Do M=L,N
             If(M  ==  N) Exit
             If(Abs(E(M))  <=  MACHEP*(Abs(D(M))+Abs(D(M+1)))) Exit
             Continue
          End Do
          P=D(L)
          If(M  ==  L) Exit
          If(J  ==  30) Then
             IERR=L
             Return
          End If
          J=J+1
          G=(D(L+1)-P) / (2.0_pr*E(L))
          R=Sqrt(G*G+1.0_pr)
          G=D(M) - P + E(L)/(G+Sign(R,G))
          S=1.0_pr
          C=1.0_pr
          P=0.0_pr
          MML=M-L
          Do II=1, MML
             I=M-II
             F=S*E(I)
             B=C*E(I)
             If(Abs(F) >= Abs(G)) Then
                C=G/F
                R=Sqrt(C*C+1.0_pr)
                E(I+1)=F*R
                S=1.0_pr/R
                C=C*S
             Else
                S=F/G
                R=Sqrt(S*S+1.0_pr)
                E(I+1)=G*R
                C=1.0_pr/R
                S=S*C
             End If
             G=D(I+1)-P
             R=(D(I)-G)*S + 2.0_pr*C*B
             P=S*R
             D(I+1)=G+P
             G=C*R - B
             F=Z(I+1)
             Z(I+1)=S*Z(I) + C*F
             Z(I)=C*Z(I) - S*F
          End Do
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.0_pr
       End Do
    End Do
    Do II=2, N
       I=II-1
       K=I
       P=D(I)
       Do J=II,N
          If(D(J)  >=  P) Cycle
          K=J
          P=D(J)
       End Do
       If(K  ==  I) Cycle
       D(K)=D(I)
       D(I)=P
       P=Z(I)
       Z(I)=Z(K)
       Z(K)=P
    End Do

  End Subroutine GBTQL2
  !=======================================================================
  !  This routine  calculates the weights and knots of the
  !  Gauss-Legendre  integration formula of the N-th order
  !  and stores them in the  arrays X and W. X1 and X2 are
  !  integration limits. Normally, they are -1 and +1, but
  !  may be set to any values (X1 < X2), so you will n o t
  !  have to change variables in your program.
  !
  !                   You always have:
  !
  !   Integral(X1,X2)[f(x)dx] = Sum(i=1,N)[W(i)*f(X(i))]
  !
  !  Arrays  are stored in  ascending  order of knots X(i).
  !=======================================================================
  Subroutine gauleg(x1,x2,x,w,n)
    Integer(ipr), Intent(In) :: n
    Real(pr), Intent(In) :: x1,x2
    Real(pr), Dimension(n), Intent(Inout) :: x,w

    Integer(ipr) :: KX,m,k,j,i
    Real(pr) :: EPS,xm,xl,z,zx,z1,z2,zb,pp,p1,p2,p3

    EPS=3.0e-17_pr; KX=8
    m=(n+1)/2
    xm=half*(x2+x1)
    xl=half*(x2-x1)
    Do i=1,m
       z=Cos(pi*(Real(i,Kind=pr)-0.25_pr)/(Real(n,Kind=pr)+0.5_pr))
       zx=1.e+30_pr
       z2=one; k=0
       Do While (z2 > EPS .And. k < KX)
          p1=one
          p2=zero
          k=k+1
          Do j=1,n
             p3=p2
             p2=p1
             p1=((two*j-one)*z*p2-(j-one)*p3)/Real(j,Kind=pr)
          End Do
          pp=Real(n,Kind=pr)*(z*p1-p2)/(z*z-one)
          z1=z
          z=z1-p1/pp
          z2=Abs(z-z1)
          if (z2 <= zx) then
             zx=z2
             zb=z
          End if
       End Do
       x(i)=xm-xl*zb
       x(n+1-i)=xm+xl*zb
       w(i)=two*xl/((one-zb*zb)*pp*pp)
       w(n+1-i)=w(i)
    End Do

  End Subroutine gauleg
  !=======================================================================
  !> D_HERM calculates Hermite polynomials at point X up to N-th degree. It also calculates
  !> derivatives DHER(i) of these polynomials, DHER(I)=2*I*HER(I-1)
  !=======================================================================
  Subroutine D_HERM(X,N,HER,DHER,NDIM)
    Integer(ipr), Intent(In) :: N, NDIM
    Real(pr), Intent(In) :: X
    Real(pr), Dimension(1:NDIM), Intent(Inout) :: HER,DHER

    Integer(ipr) :: i
    Real(pr) :: F,DF

    HER (1)=one
    DHER(1)=zero

    If(N <= 0) Return

    HER (2)=X + X
    DHER(2)=two

    If((N-1) <= 0) Return

    Do I=2,N-1
       F=X*HER(I)-Real(I-1,Kind=pr)*HER(I-1)
       HER(I+1)=F+F
       DF=Real(I,Kind=pr)*HER(I)
       DHER(I+1)=DF+DF
    End Do

  End Subroutine D_HERM
  !=======================================================================
  !> The routine computes the expansion of a product of two Hermite polynomials up to order NOSACT
  !> as a linear combination of Hermite> polynomials. The expansion goes from 0 t0 2*NOSACT.
  !> Coefficients are stored in the global variable COEF00(0:NOSACT,0:NOSACT,0:ngh).The routine
  !> also computes the normalization factors of Hermite polynomials up to order 2*NOSACT, which
  !> are stored in HERFAC (also a global variable.
  !=======================================================================
  Subroutine DEVHER(NOSACT)
    Use HFBTHO, Only : nzm

    Integer(ipr), Intent(In) :: NOSACT

    Integer(ipr) :: I,IZEROS,NOSCIL,N,M,K,IGAUSS,NORDER
    Real(pr) :: XZER,RESULT,toto
    Real(pr), Dimension(1:2*NOSACT+1) :: PHERMI,DHERMI
    Real(pr), Dimension(0:2*NOSACT,1:ngh) :: HERPLN

    NORDER=2*NOSACT
    HERPLN = 0.0_pr; PHERMI = 0.0_pr

    If(Allocated(COEF00)) Deallocate(COEF00)
    Allocate(COEF00(0:NORDER,0:NOSACT,0:NOSACT)); COEF00 = 0.0_pr
    If(Allocated(HERFAC)) Deallocate(HERFAC)
    Allocate(HERFAC(0:NORDER)); HERFAC = 1.0_pr

    HERFAC(0)=Sqrt(Sqrt(pi))
    Do I=1,NORDER
       HERFAC(I)=HERFAC(I-1)*Sqrt(Real(2*I,Kind=pr))
    End Do

    ! Defining the values of the Hermite polynomials at Gauss zeros
    Do IZEROS = 1,ngh
       XZER = xh(IZEROS)
       Call D_HERM(XZER,NORDER,PHERMI,DHERMI,NORDER+1)
       Do NOSCIL = 0,NORDER
          HERPLN(NOSCIL,IZEROS) = PHERMI(NOSCIL+1)/HERFAC(NOSCIL)
       End Do
    End Do
    ! Testing orthonormalization of Hermite polynomials
    If(debug_gauss >= 1) Then
       Do IZEROS = 1,ngh
          Do NOSCIL = 0,nzm
             toto = Exp(half*xh(IZEROS)*xh(IZEROS))*qh(NOSCIL,IZEROS)/Sqrt(wh(IZEROS))
             Write(*,'("NOSCIL=",i4," IZEROS=",i4," xh=",e24.12," HERPLN-qh=",e24.12)') &
                        NOSCIL,IZEROS,xh(IZEROS),HERPLN(NOSCIL,IZEROS)-toto
          End Do
       End Do
    End If
    ! Calculating the expansion coefficients for the polynomials
    Do N=0,NOSACT
       Do M=0,NOSACT
          Do K=M+N,0,-2
             RESULT=zero
             Do IGAUSS=1,ngh
                RESULT=RESULT+wh(IGAUSS)*Exp(-xh(IGAUSS)*xh(IGAUSS))*HERPLN(N,IGAUSS)*HERPLN(M,IGAUSS)*HERPLN(K,IGAUSS)
             End Do
             COEF00(K,M,N)=RESULT
          End Do
       End Do
    End Do
    ! Testing orthonormalization of Hermite polynomials
    If(debug_gauss >= 1) Then
       Do N=0,NOSACT
          Do M=0,NOSACT
             RESULT=zero
             Do IGAUSS=1,ngh
                RESULT=RESULT+wh(IGAUSS)*Exp(-xh(IGAUSS)*xh(IGAUSS))*HERPLN(N,IGAUSS)*HERPLN(M,IGAUSS)
             End Do
             Write(*,'("N=",i3," M=",i3," \int e^(-x^2) H_n(x) H_m(x) = ",e24.12)') N,M,RESULT
          End Do
       End Do
    End If

  End Subroutine DEVHER
  !=======================================================================
  !>  pr_gamma evaluates \f$ \Gamma(x) \f$ for a real argument.
  !>
  !>  Discussion:
  !>    This function was originally named DGAMMA. However, a number of
  !>    Fortran compilers now include a library function of this name. To
  !>    avoid conflicts, this function was renamed pr_gamma. This routine
  !>    calculates the GAMMA function for a real argument X. Computation
  !>    is based on an algorithm outlined in reference 1 below. The
  !>    program uses rational functions that approximate the GAMMA
  !>    function to at least 20 significant decimal digits. Coefficients
  !>    for the approximation over the interval (1,2) are unpublished.
  !>    Those for the approximation for 12 <= X are from reference 2.
  !>
  !>  Licensing:
  !>    This code is distributed under the GNU LGPL license.
  !>
  !>  Modified:
  !>    18 January 2008
  !>
  !>  Author:
  !>    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !>    FORTRAN90 version by John Burkardt.
  !>
  !>  Reference:
  !>    - 1. William Cody, "An Overview of Software Development for Special
  !>         Functions," in Numerical Analysis Dundee, 1975, Edited by GA
  !>         Watson, Lecture Notes in Mathematics 506, Springer, 1976.
  !>    - 2. John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles
  !>        Mesztenyi, John Rice, Henry Thatcher, Christoph Witzgall,
  !>        "Computer Approximations," Wiley, 1968, LC: QA297.C64.
  !>
  !>  Parameters:
  !>    - Input, real ( kind = 8 ) X, the argument of the function.
  !>    - Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
  !=======================================================================
  Real(pr) Function pr_gamma(x)
    Real(pr), Intent(In) :: x

    ! Coefficients for minimax approximation over (12, INF).
    Logical :: parity
    Integer(ipr) :: i,n
    Real(pr) :: eps,fact,res,sum,xbig,xden,xinf,xminin,xnum,y,y1,ysq,z
    Real(pr), Parameter :: sqrtpi = 0.9189385332046727417803297_pr
    Real(pr), Dimension(7), Parameter :: c = [ -1.910444077728000000000e-03_pr, 8.417138778129500000000000e-04_pr, &
                                               -5.952379913043012000000e-04_pr, 7.936507935003502480000000e-04_pr, &
                                               -2.777777777777681622553e-03_pr, 8.333333333333333331554247e-02_pr, &
                                                                                5.708383526100000000000e-03_pr ]
    ! Numerator and denominator coefficients for rational minimax approximation over (1,2).
    Real(pr), Dimension(1:8), Parameter :: p = [ -1.71618513886549492533811e+00_pr,  2.47656508055759199108314e+01_pr, &
                                                 -3.79804256470945635097577e+02_pr,  6.29331155312818442661052e+02_pr, &
                                                  8.66966202790413211295064e+02_pr, -3.14512729688483675254357e+04_pr, &
                                                 -3.61444134186911729807069e+04_pr,  6.64561438202405440627855e+04_pr ]
    Real(pr), Dimension(1:8), Parameter :: q = [ -3.08402300119738975254353e+01_pr,  3.15350626979604161529144e+02_pr, &
                                                 -1.01515636749021914166146e+03_pr, -3.10777167157231109440444e+03_pr, &
                                                  2.25381184209801510330112e+04_pr,  4.75584627752788110767815e+03_pr, &
                                                 -1.34659959864969306392456e+05_pr, -1.15132259675553483497211e+05_pr ]

    parity = .False.; fact = one; n = 0; y = x
    xbig = 171.624e+00_pr
    xminin = Tiny(1.0_pr); eps = Epsilon(1.0_pr) ; xinf = Huge(1.0_pr)

    ! Argument is negative.
    If(y <= zero) Then

       y = - x
       y1 = Aint ( y )
       res = y - y1

       If(res/=zero) Then
          If(y1/=Aint(y1*half)*two) Then
             parity = .True.
          End If
          fact = - pi / Sin(pi*res)
          y = y + one
       Else
          res = xinf
          pr_gamma = res
          Return
       End If

    End If

    ! Argument is positive.
    if (y < eps) Then

       ! Argument < EPS.
       If(xminin <= y) Then
          res = one / y
       Else
          res = xinf
          pr_gamma = res
          Return
       End If

    Else If(y<pp12) Then

       y1 = y

       ! 0.0 < argument < 1.0.
       If(y < one) Then
          z = y
          y = y + one
       Else
          ! 1.0 < argument < 12.0.
          ! Reduce argument if necessary.
          n = Int(y) - 1
          y = y - Real(n,Kind=pr)
          z = y - one
       End If
       !
       ! Evaluate approximation for 1.0 < argument < 2.0.
       xnum = zero
       xden = one
       Do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
       End Do

       res = xnum / xden + one

       ! Adjust result for case  0.0 < argument < 1.0.
       If(y1<y) Then
          res = res / y1
       Else If(y<y1) Then
          ! Adjust result for case 2.0 < argument < 12.0.
          Do i = 1, n
             res = res * y
             y = y + one
          End Do
       End If

    Else

       !  Evaluate for 12.0 <= argument.
       If(y <= xbig) Then
          ysq = y * y
          sum = c(7)
          Do i = 1, 6
             sum = sum / ysq + c(i)
          End Do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - half ) * Log(y)
          res = Exp(sum)
       Else
          res = huge ( res )
          pr_gamma = res
          Return
       End If

    End If

    !  Final adjustments and Return.
    If(parity) Then
       res = - res
    End If

    If(fact/=one) Then
       res = fact / res
    End If

    pr_gamma = res

  End Function pr_gamma  ! Gamma function in double precision
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_Gauss
