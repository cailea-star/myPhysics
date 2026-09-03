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
!                      MULITPOLE MOMENTS PACKAGE                       !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module defines and computes the matrix elements and expectation
!> values of the axial multipole moments \f$ \hat{Q}_{\lambda 0} \f$
!> defined by
!>   \f[
!>        \hat{Q}_{\lambda 0}(r,\theta,\varphi)
!>        = \mathcal{N}_{\lambda}r^{\lambda} Y_{\lambda 0} (\theta,\varphi).
!>   \f]
!> Assuming axial symmetry, the multipole moments can be written
!>   \f[
!>        \hat{Q}_{\lambda 0}(r,\theta)
!>        = \mathcal{N}_{\lambda} \sqrt{\frac{2\lambda+1}{4\pi}}
!>                              r^{\lambda} P_{\lambda} (\cos\theta).
!>   \f]
!> In cylindrical coordinates, \f$ r^2 = \rho^{2} + z^2 \f$ and \f$ r\cos\theta = z \f$.
!> The multipole moments can be written as a function of \f$ \rho \f$ and \f$ z \f$ only.
!> The normalization constants \f$ \mathcal{N}_{\lambda} \f$ are the same as in the code
!> HFODD: all moments are expressed in powers of barns and the coefficients are thus
!>
!>   \f{align}{
!>        \mathcal{N}_{0} &  =\sqrt{4\pi} \\
!>        \mathcal{N}_{1} &  =\sqrt{\frac{4\pi}{3}}\frac{1}{10} \\
!>        \mathcal{N}_{2} &  =\sqrt{\frac{16\pi}{5}}\frac{1}{100} \\
!>        \mathcal{N}_{\lambda} &  =\frac{1}{10^{\lambda}},\ \ \lambda>2
!>   \f}
!>
!> @author
!> Nicolas Schunck
!-------------------------------------------------------------------
!  Subroutines: - moments_setUnits
!               - moments_computeValue
!               - moments_valueMesh(z,rrr,Qval)
!               - moments_computeField(lambda,ib,debug)
!               - moments_expectation(lambda,it,ib,qval,rho,dd)
!               - J_computeValue()
!               - calc_lx_tot(nz1,nr1,nl1,nz2,nr2,nl2,rlx_tot)
!----------------------------------------------------------------------!
Module HFBTHO_multipole_moments

  Use HFBTHO_utilities
  Use HFBTHO, Only : lambdaMax,ncons_max

  Implicit None

  Public moments_setUnits,moments_computeValue,moments_valueMesh,moments_computeField, &
         moments_expectation,J_computeValue

  Private calc_lx_tot

  Integer(ipr), Public, Save :: numberCons
  Integer(ipr), Dimension(1:ncons_max), Public, Save :: multLambda
  Real(pr), Dimension(0:lambdaMax), Public, Save :: q_units,multLag,multRequested
  Real(pr), Dimension(0:lambdaMax,1:3), Public, Save :: qmoment

Contains
  !=======================================================================
  !> This routine defines standard units for the multipole moments.
  !=======================================================================
  Subroutine moments_setUnits()
    Use HFBTHO, Only : lambdaMax

    Integer(ipr) :: lambda
    Real(pr) :: sqr4pi

    q_units=one
    sqr4pi=Sqrt(4.0_pr*pi)

    q_units(0)=+sqr4pi
    q_units(1)=+sqr4pi/Sqrt(three)
    q_units(2)=+sqr4pi/Sqrt(five)*two

    Do lambda=0,lambdaMax
       q_units(lambda)=q_units(lambda) / ten**lambda
    End Do

  End Subroutine moments_setUnits
  !=======================================================================
  !> This routine computes the expectation value of multipole moments in coordinate space.
  !=======================================================================
  Subroutine moments_computeValue()
    Use HFBTHO_Gauss, Only : nghl,wdcor,fh,fl
    Use HFBTHO, Only : lambdaMax,ro

    Integer(ipr) :: lambda,ihli
    Real(pr), Dimension(0:lambdaMax) :: Qval
    Real(pr) :: z,rrr,rown,rowp,whl,rn,rp

    qmoment=zero; Qval=zero

    Do ihli=1,nghl

       whl=wdcor(ihli)
       rn=ro(ihli,1); rp=ro(ihli,2)
       rown=whl*rn; rowp=whl*rp;
       z=fh(ihli); rrr=fl(ihli)**2

       Call moments_valueMesh(z,rrr,Qval)

       Do lambda=0,lambdaMax
          qmoment(lambda,1)=qmoment(lambda,1)+rown*Qval(lambda)
          qmoment(lambda,2)=qmoment(lambda,2)+rowp*Qval(lambda)
       End Do

    End Do

    Do lambda=0,lambdaMax
       qmoment(lambda,3)=qmoment(lambda,1)+qmoment(lambda,2)
    End Do

  End Subroutine moments_computeValue
  !=======================================================================
  !> This routine computes the value of multipole moments at point \f$ (\rho,z)\f$ of the quadrature
  !> grid.
  !>  @param[in] z - Value of z
  !>  @param[in] rrr - Value of radial coordinate \f$ \rho \f$
  !>  @param[inout] Qval - Value of \f$ \hat{Q}_{\lambda}(\rho,z) \f$
  !=======================================================================
  Subroutine moments_valueMesh(z,rrr,Qval)
    Use HFBTHO, Only : lambdaMax,Parity

    Real(pr), Intent(In) :: z,rrr
    Real(pr), Dimension(0:lambdaMax), Intent(Inout) :: Qval

    Integer(ipr) :: lambda
    Real(pr) :: sqr4pi,z2,z3,z4,z5,z6,z7,z8,rrr4,rrr6,fq,betac,rc,ac

    sqr4pi=one/Sqrt(4.0_pr*pi)

    z2=z*z; z3=z2*z; z4=z3*z; z5=z4*z; z6=z5*z; z7=z6*z; z8=z7*z
    rrr4=rrr*rrr; rrr6=rrr4*rrr

    Qval(0) =               sqr4pi
    Qval(1) = Sqrt(three)  *sqr4pi          * z
    Qval(2) = Sqrt(five)   *sqr4pi*half     * (two*z2-        rrr)
    Qval(3) = Sqrt(seven)  *sqr4pi*half     * (two*z3-three*z*rrr)
    Qval(4) = Sqrt(nine)   *sqr4pi*p18      * (eight*z4-24.0_pr*z2*rrr    + three *rrr4)
    Qval(5) = Sqrt(11.0_pr)*sqr4pi*p18      * (eight*z5-   pp40*z3*rrr    +pp15*z *rrr4)
    Qval(6) = Sqrt(13.0_pr)*sqr4pi/pp16     * (pp16*z6-120.0_pr*z4*rrr+ 90.0_pr*z2*rrr4-five     *rrr6)
    Qval(7) = Sqrt(15.0_pr)*sqr4pi/pp16     * (pp16*z7-168.0_pr*z5*rrr+210.0_pr*z3*rrr4-35.0_pr*z*rrr6)
    Qval(8) = Sqrt(17.0_pr)*sqr4pi/128.0_pr * (128.0_pr*z8-1792.0_pr*z6*rrr +3360.0_pr*z4*rrr4 &
                                                          -1120.0_pr*z2*rrr6+  35.0_pr*rrr4*rrr4)

    If(Parity) Then
       Qval(1)=zero; Qval(3)=zero;Qval(5)=zero; Qval(7)=zero
    End If

    fq = one

    Do lambda=0,lambdaMax
       Qval(lambda)=Qval(lambda)*q_units(lambda)*fq
    End Do

  End Subroutine moments_valueMesh
  !=======================================================================
  !> This routine computes the matrix elements of the multipole moment in the HO basis.
  !>  @param[in] lambda - Multipolarity \f$ \lambda \f$
  !>  @param[in] ib - Index of the current K block
  !>  @param[inout] multMatElems - Vector containing the matrix elements of \f$ \hat{Q}_{\lambda} \f$
  !>  @param[in] debug - Optional argument doing some additional tests
  !=======================================================================
  Subroutine moments_computeField(lambda,ib,multMatElems,debug)
    Use HFBTHO_Gauss, Only : nghl,fh,fl
    Use HFBTHO, Only : n00max,lambdaMax,id,ia,ns,QHLA_opt

    Integer(ipr), Intent(In) :: lambda,ib
    Real(pr), Allocatable, Intent(Inout) :: multMatElems(:)
    Logical, Optional, Intent(In) :: debug

    Integer(ipr) :: i,nd,nd2,ihli,ihil,im,n1,n2
    Integer(ipr) :: ja,nsa,nsb,nsab,ssu,ssd
    Real(pr) :: qhla,vh,fiun1,fiun2,fidn1,fidn2,vnhl
    Real(pr), Dimension(1:nghl) :: Vmom
    Real(pr), Dimension(0:lambdaMax) :: Qval
    Real(pr) :: z,rrr
    Integer(ipr), Parameter :: ndxmax=(n00max+2)*(n00max+2)/4
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax)

    Qval=zero

    ! Compute moment lambda on integration mesh
    Do ihli=1,nghl
       z=fh(ihli); rrr=fl(ihli)**2
       Call moments_valueMesh(z,rrr,Qval)
       Vmom(ihli)=Qval(lambda)
    End Do !ihli

    ! Form matrix of the multipole constraint lambda in HO basis
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
                If(.Not.present(debug)) Then
                   multMatElems(i)=multMatElems(i)+vh*vnhl
                End If
             End If
             If(present(debug)) Then
                If(n2 == n1) multMatElems(i)=multMatElems(i)+vh
             End If
          End Do !n2
       End Do !n1
    End Do !ihil

  End Subroutine moments_computeField
  !=======================================================================
  !> This routine computes the expectation value \f$ \langle \hat{Q}_{\lambda} \rangle \f$ of multipole
  !> moments in configuration space.
  !=======================================================================
  Subroutine moments_expectation(lambda,it,ib,qval,rho,dd,multMatElems)
    Use HFBTHO, Only : id,ia,nbx

    Integer(ipr), Intent(In) :: lambda,it,ib
    Real(pr), Allocatable, Intent(In) :: rho(:,:)
    Real(pr), Intent(Inout) :: qval,dd
    Real(pr), Allocatable, Intent(Inout) :: multMatElems(:)

    Integer(ipr) :: nd,nd2,nhfb,i0,m,j,n1,n2
    Real(pr) :: hla
    Real(pr), Allocatable :: dblmul(:,:),qblock(:,:)

    nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
    multMatElems=zero
    Call moments_computeField(lambda,ib,multMatElems)
    ! Trace of the density should give particle number
    Do n1=1,nd
       dd=dd+rho(n1,n1)
    End Do
    ! Expectation value of multipole moment in configuration space
    Allocate(dblmul(nd,nd));dblmul=zero
    j=0
    Do n1=1,nd
       Do n2=1,n1
          j=j+1;hla=multMatElems(j)
          dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
       End Do
    End Do
    Allocate(qblock(nd,nd))
    Call dgemm('n','n',nd,nd,nd,one,rho,nd,dblmul,nd,zero,qblock,nd)
    Do n1=1,nd
       qval=qval+qblock(n1,n1)
    End Do

  End Subroutine moments_expectation
  !=====================================================
  !> This routine computes the expectation value \f$ \langle \hat{J}^{2} \rangle \f$ of the
  !> square angular momenta.
  !=====================================================
  Subroutine J_computeValue()
    Use HFBTHO, Only : nb,id,ia,nr,nl,nz,ns,kd,nt,ka,RVqpN,RUqpN,KpwiN,RVqpP,RUqpP,KpwiP, &
                       jxsq_rhoN,jxsq_rhoP,jxsq_kappaN,jxsq_kappaP

    Integer(ipr) :: n1,n2,ib1,nd1,im1,nz1,nr1,nl1,ns1,ib2,nd2,im2,nz2,nr2,nl2,ns2,iosc1,iosc2,it,ib,nd,n_qp,im,k
    Real(pr) :: rlx_tot,rsx_tot
    Real(pr), Allocatable:: rhoN(:,:),rhoNm(:,:),rhoP(:,:),rhoPm(:,:),kappaN(:,:),kappaP(:,:),aux1(:,:),aux2(:,:),aux3(:,:)
    Real(pr), Allocatable:: jx_tot(:,:)

    Allocate(jx_tot(nt,nt)); jx_tot=zero
    Allocate(rhoN(nt,nt),rhoNm(nt,nt),rhoP(nt,nt),rhoPm(nt,nt),kappaN(nt,nt),kappaP(nt,nt),aux1(nt,nt),aux2(nt,nt),aux3(nt,nt))
    rhoN=zero; rhoNm=zero; rhoP=zero; rhoPm=zero; kappaN=zero; kappaP=zero

    ! ---------------------------------------------
    ! Calculate matrix elements of the jx operator
    ! ---------------------------------------------
    Do ib1=1,nb
       nd1=id(ib1); im1=ia(ib1)
       Do n1=1,nd1
          iosc1=im1+n1; nr1=nr(iosc1); nl1=nl(iosc1); nz1=nz(iosc1); ns1=ns(iosc1)

          Do ib2=1,nb
             nd2=id(ib2); im2=ia(ib2)
             Do n2=1,nd2
                iosc2=im2+n2; nr2=nr(iosc2); nl2=nl(iosc2); nz2=nz(iosc2); ns2=ns(iosc2)

                ! Orbital angular momentum part
                rlx_tot=zero
                If(ns1 == ns2) Then
                   Call calc_lx_tot(nz1,nr1,nl1,nz2,nr2,nl2,rlx_tot)
                End If

                ! Spin part
                rsx_tot=zero
                If(nr1 == nr2 .And. nz1 == nz2 .And. nl1 == nl2 .And. ns1 == -ns2) rsx_tot=0.5_pr

                ! Matrix elements
                jx_tot(iosc1,iosc2)=rlx_tot+rsx_tot

             End Do
          End Do

       End Do
    End Do

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
                      rhoN(iosc1,iosc2)  = rhoN(iosc1,iosc2)   + RVqpN(KpwiN(ka(ib,it)+k)+n1)*RVqpN(KpwiN(ka(ib,it)+k)+n2)
                      kappaN(iosc1,iosc2)= kappaN(iosc1,iosc2) + RVqpN(KpwiN(ka(ib,it)+k)+n1)*RUqpN(KpwiN(ka(ib,it)+k)+n2)
                   Else If(it == 2) Then ! protons
                      rhoP(iosc1,iosc2)  = rhoP(iosc1,iosc2)   + RVqpP(KpwiP(ka(ib,it)+k)+n1)*RVqpP(KpwiP(ka(ib,it)+k)+n2)
                      kappaP(iosc1,iosc2)= kappaP(iosc1,iosc2) + RVqpP(KpwiP(ka(ib,it)+k)+n1)*RUqpP(KpwiP(ka(ib,it)+k)+n2)
                   End If
                End Do ! k
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

    ! ----------------------------------------------------
    ! Calculate angular momentum Jx with respect to center
    ! ----------------------------------------------------
    ! rho contribution
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_tot,nt,rhoN,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,one,jx_tot,nt,rhoNm,nt,zero,aux1,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux1,nt,aux2,nt,zero,aux3,nt)
    jxsq_rhoN=zero
    Do iosc1=1,nt
       jxsq_rhoN=jxsq_rhoN+aux3(iosc1,iosc1)
    End Do

    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_tot,nt,rhoP,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,one,jx_tot,nt,rhoPm,nt,zero,aux1,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux1,nt,aux2,nt,zero,aux3,nt)
    jxsq_rhoP=zero
    Do iosc1=1,nt
       jxsq_rhoP=jxsq_rhoP+aux3(iosc1,iosc1)
    End Do

    ! kappa contribution
    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_tot,nt,kappaN,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux2,nt,aux2,nt,zero,aux3,nt)
    jxsq_kappaN=zero
    Do iosc1=1,nt
       jxsq_kappaN=jxsq_kappaN+aux3(iosc1,iosc1)
    End Do

    aux1=zero; aux2=zero; aux3=zero
    Call dgemm('N','N',nt,nt,nt,one,jx_tot,nt,kappaP,nt,zero,aux2,nt)
    Call dgemm('N','N',nt,nt,nt,two,aux2,nt,aux2,nt,zero,aux3,nt)
    jxsq_kappaP=zero
    Do iosc1=1,nt
       jxsq_kappaP=jxsq_kappaP+aux3(iosc1,iosc1)
    End Do

    Deallocate(jx_tot)
    Deallocate(rhoN,rhoNm,rhoP,rhoPm,kappaN,kappaP,aux1,aux2,aux3)

  End Subroutine J_computeValue
  !============================================================================
  !> This routine calculates the expectation value of the \f$ \hat{\ell}_x \f$ operator in the HO basis
  !> with respect to the center of mass of the coordinate system
  !==========================================================================
  Subroutine calc_lx_tot(nz1,nr1,nl1,nz2,nr2,nl2,rlx_tot)
    Use math, Only : sq,sqi
    Use HFBTHO_Gauss, Only : ngl,sxl,ql,ql1
    Use HFBTHO, Only : bz,bp

    Integer(ipr), Intent(In) :: nz1,nr1,nl1 !< - Quantum numbers of the left (bra) state
    Integer(ipr), Intent(In) :: nz2,nr2,nl2 !< - Quantum numbers of the right (ket) state
    Real(pr), Intent(Out) :: rlx_tot !< - Matrix element of \f$ \hat{\ell}_x \f$

    Integer(ipr) :: il
    Real(pr) :: fz,fp,fac1,fac2,fac3,facr1,facr2,facr3,facz1,facz2,delnl1,delnl2

    rlx_tot=zero

    If(nl1 /= nl2+1 .And. nl1 /= nl2-1) Return
    If(nz1 /= nz2+1 .And. nz1 /= nz2-1) Return

    fz=1.0_pr/bz; fp=1.0_pr/bp
    delnl1=0.0_pr; delnl2=0.0_pr
    If(nl1 == nl2+1) delnl1=1.0_pr
    If(nl1 == nl2-1) delnl2=1.0_pr

    !-----------------------------
    ! integrals in rperp direction
    ! ----------------------------
    facr1=zero; facr2=zero; facr3=zero
    Do il=1,ngl
       facr1 = facr1+ql(nr1,nl1,il)*ql(nr2,nl2,il)*sxl(il)  ! linear in rperp
       facr2 = facr2+ql(nr1,nl1,il)*ql(nr2,nl2,il)/sxl(il)  ! 1/rperp
       facr3 = facr3+ql(nr1,nl1,il)*ql1(nr2,nl2,il)/sxl(il) ! derivation over rperp
    End Do

    facr1=facr1*bp; facr2=facr2*fp; facr3=facr3*fp        ! proper dimensions
    ! ------------------------
    ! integrals in z direction
    ! ------------------------
    facz1=zero; facz2=zero
    If(nz1 == nz2+1) Then
       facz1 = -fz*sq(nz2+1)*sqi(2)
       facz2 = sq(nz2+1)*sqi(2)*bz
    Else If (nz1 == nz2-1) Then
       facz1 = +fz*sq(nz2)*sqi(2)
       facz2 = sq(nz2)*sqi(2)*bz
    End If
    ! ---------------
    ! matrix elements
    ! ---------------
    fac1 = -0.5_pr*(delnl1-delnl2)*facr1*facz1
    fac2 = -0.5_pr*(delnl1+delnl2)*nl2*facr2*facz2
    fac3 = +0.5_pr*(delnl1-delnl2)*facr3*facz2
    rlx_tot = fac1+fac2+fac3

  End Subroutine calc_lx_tot
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_multipole_moments
