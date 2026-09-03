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
!                PARTICLE NUMBER PROJECTION PACKAGE                    !
!                      (CANONICAL BASIS ONLY)                          !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module implements particle number projection of the HFB solution
!> in the canonical basis. It contains 3 main routines to compute the projected
!> densities, the expectation values of various contributions to the energy
!> on the projected HFB state, and the Coulomb energy.
!>
!>  @author
!>  Mario Stoitsov, Nicolas Schunck
!----------------------------------------------------------------------
!  Subroutines: - allocate_arrays_PNP()
!               - deallocate_arrays_PNP()
!               - expectpj(lpr)
!               - densitpj()
!               - coulompj()
!----------------------------------------------------------------------
Module HFBTHO_PNP

  Use HFBTHO_utilities

  Implicit None

  Public expectpj

  Private allocate_arrays_PNP,deallocate_arrays_PNP,densitpj,coulompj

  ! Public variables
  Integer(ipr), Public, Save :: keypjn,keypjp,ilpjn,ilpjp,ilpj,ilpjmax,ilpj2,ilnqx,ilnghl
  Integer(ipr), Public, Save :: iproj,npr1pj,npr2pj
  Real(pr), Public, Save :: rehfbcan,retotpj,depnp
  Complex(pr), Allocatable, Target, Public, Save :: phypj(:),sinphy(:),exp1iphy(:),exp1iphym(:),exp2iphy(:),exp2iphym(:)
  Complex(pr), Allocatable, Target, Public, Save :: coupj(:,:)
  Complex(pr), Allocatable, Target, Public, Save :: ropj(:,:,:),taupj(:,:,:),dropj(:,:,:),djpj(:,:,:),akapj(:,:,:),pjk(:,:)
  Complex(pr), Allocatable, Target, Public, Save :: SZFIpj(:,:,:),SFIZpj(:,:,:),SRFIpj(:,:,:),SFIRpj(:,:,:)
  Complex(pr), Allocatable, Target, Public, Save :: epj(:,:),ddepj(:,:,:),cpj(:,:,:),ypj(:,:,:),rpj(:,:,:)
  Real(pr), Allocatable, Target, Public, Save :: xg_norm(:,:),xg_phase(:,:)
  Real(pr), Private, Save :: number_parity = 1.0_pr

Contains
  !=======================================================================
  !> This routine allocates arrays needed for particle-number projection in the canonical basis.
  !=======================================================================
  Subroutine allocate_arrays_PNP()
    Use HFBTHO_Gauss, Only : nghl
    USe HFBTHO, Only : nqx,nb2x,nkx

    ilnqx=ilpjmax*nqx; ilnghl=ilpjmax*nghl
    If(Allocated(exp1iphy)) Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk,SZFIpj,SFIZpj,SRFIpj,SFIRpj, &
                                       epj,cpj,ypj,rpj,ddepj,phypj,sinphy,exp1iphy,exp2iphy,exp1iphym,exp2iphym, &
                                       xg_norm,xg_phase)
    Allocate(ropj(nghl,ilpjmax,2),taupj(nghl,ilpjmax,2),dropj(nghl,ilpjmax,2),djpj(nghl,ilpjmax,2),akapj(nghl,ilpjmax,2),  &
             SZFIpj(nghl,ilpjmax,2),SFIZpj(nghl,ilpjmax,2),SRFIpj(nghl,ilpjmax,2),SFIRpj(nghl,ilpjmax,2), &
             coupj(nghl,ilpjmax))
    Allocate(epj(ilpjmax,2),pjk(ilpjmax,2),xg_norm(ilpjmax,2),xg_phase(ilpjmax,2))
    Allocate(cpj(nkx,ilpjmax,2),ypj(nkx,ilpjmax,2),rpj(nkx,ilpjmax,2),ddepj(nqx,ilpjmax,nb2x))
    Allocate(phypj(ilpjmax),sinphy(ilpjmax),exp1iphy(ilpjmax),exp2iphy(ilpjmax),exp1iphym(ilpjmax),exp2iphym(ilpjmax))

  End Subroutine allocate_arrays_PNP
  !=======================================================================
  !> This routine deallocate arrays needed for particle-number projection in the canonical basis.
  !=======================================================================
  Subroutine deallocate_arrays_PNP()

    If(Allocated(exp1iphy)) Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk,SZFIpj,SFIZpj,SRFIpj,SFIRpj, &
                                       epj,cpj,ypj,rpj,ddepj,phypj,sinphy,exp1iphy,exp2iphy,exp1iphym,exp2iphym, &
                                       xg_norm,xg_phase)

  End Subroutine deallocate_arrays_PNP
  !=======================================================================
  !> This routine calculates the expectation values (tz is the particle number) in the case
  !> of particle-number projection in the canonical basis.
  !=======================================================================
  Subroutine expectpj(lpr)
    Use UNEDF, Only : CrdJ,CpV0,CpV1,sigma,use_j2terms,ta7,ta8,CExPar,hbzeron,hbzerop,ts1,ts2,t4o3,tv1,tv2,tv3,tv4,tv5,tv6
    Use HFBTHO_Gauss, Only : nghl,wdcor,fh,fl
    Use HFBTHO_canonical, Only : Dispersion,v2min
    Use HFBTHO_multipole_moments, Only : moments_computeValue,numberCons,multLambda,multLag,qmoment,multRequested
    Use HFBTHO, Only : rho_c,icou,tz,icstr,itmin,itmax,rms,q2,ass,bet,cex,amas,ffdef4,ffdef6,ffdef7, &
                       neckLag,neckValue,neckRequested,iLST1,pwi,keyblo,bloblo,ala,alast,iiter,iLST, &
                       b0,bz,bp,r00,hb0n,hb0p,ECMPAV,ehfb,ieresj,eresj,iasswrong,npr,inin,kindhfb,r02,r04, &
                       si,skyrme,do_print

    Logical, Intent(In) :: lpr !< - Prints (or not) results on screen

    Integer(ipr) :: i,j,it,ihli,iw,icons,lambda,ign,igp,k
    Integer(ipr), Allocatable :: in_vec(:),ip_vec(:)
    Complex(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Complex(pr) :: cekt(3),cdel(2),cept(3),cetot,etens,cq2pj(ilpjmax,ilpjmax)
    Complex(pr) :: cxn(2),crms(3),cq2(3),cq4(3),xnpj(2),rmspj(3),q2pj(3),q4pj(3)
    Complex(pr) :: evolpj,esurpj,ecdipj,ecexpj,ecoupj,ept1pj,ept2pj,epotpj,eki1pj, &
                   eki2pj,ekinpj,etotpj,espopj,epa1pj,epa2pj,epirpj,ede1pj,ede2pj, &
                   etenspj
    Complex(pr) :: eva,ev3,ev5,es5,eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,rt2,tnt,tpt,tt, &
                   dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul,pijk,row,cx,   &
                   dd1n,dd1p,ceptn,ceptp,cdeln,cdelp,cektn,cektp,e_rhorho,e_rhotau
    Complex(pr) :: rsa0
    Real(pr) :: whl,xn(3),q4(3),def(3),bet2(3),het4(3),r212,r222,rc,z,zz,rrr,p2,p3,p4
    Real(pr) :: rdelta(2),repair(3),rekin(3),revolpj,resurpj,respopj,recdipj,recexpj,retenspj
    Real(pr) :: RpRMSsq,RnRMSsq,DarwinFoldy

    ! Allocate necessary arrays
    Call allocate_arrays_PNP()

    ! Calculate complex densities and the direct Coulomb field
    Call densitpj()

    evolpj = c_zero; esurpj = c_zero; ecdipj = c_zero; ecexpj = c_zero; espopj  = c_zero;
    e_rhorho = c_zero; e_rhotau = c_zero
    ept1pj = c_zero; ept2pj = c_zero; eki1pj = c_zero; eki2pj = c_zero; epj     = c_zero;
    etotpj = c_zero; epa1pj = c_zero; epa2pj = c_zero; ede1pj = c_zero; ede2pj  = c_zero;
    xnpj   = c_zero; rmspj  = c_zero; q2pj   = c_zero; q4pj   = c_zero; etenspj = c_zero; cq2pj = c_zero

    Allocate(in_vec(1:ilpjn*ilpjp),ip_vec(1:ilpjn*ilpjp))
    i=0
    Do ign=1,ilpjn
       Do igp=1,ilpjp
          i=i+1
          in_vec(i) = ign
          ip_vec(i) = igp
       End Do
    End Do

    ! Loop over both neutron and proton gauge angles. Could be multithreaded in the future
    Do k=1,ilpjn*ilpjp
       i = in_vec(k)
       j = ip_vec(k)

       pijk = pjk(i,1)*pjk(j,2)

       cekt(:) = c_zero; cept(:) = c_zero; cdel(:) = c_zero
       cxn(:)  = c_zero; crms(:) = c_zero; cq2(:)  = c_zero; cq4(:) = c_zero

       eva   = c_zero; ev3   = c_zero; ev5   = c_zero; es5   = c_zero
       eso   = c_zero; ecodi = c_zero; ecoex = c_zero; etens = c_zero
       ceptn = c_zero; ceptp = c_zero; cdeln = c_zero; cdelp = c_zero
       cektn = c_zero; cektp = c_zero

       Do ihli = 1,nghl
          ! real
          whl = wdcor(ihli)
          z   = fh(ihli); zz = z*z; rrr = zz + fl(ihli)**2
          p2  = p32*zz   - half*rrr    !3/2 z*z-1/2 (z*z+r*r)=1/2(2 z*z-r*2)=1/2 Q
          p3  = p53*z*p2 - p23*rrr*z
          p4  = p74*z*p3 - p34*rrr*p2
          ! complex
          rn  = ropj(ihli,i,1); rp  = ropj(ihli,j,2); rnp2 = rn**2 + rp**2; rnp1=rn - rp
          ! ig - particle number, rms and deformations
          row = whl*rn; cxn(1)=cxn(1)+row; crms(1)=crms(1)+row*rrr; cq2(1)=cq2(1)+row*p2; cq4(1)=cq4(1)+row*p4
          row = whl*rp; cxn(2)=cxn(2)+row; crms(2)=crms(2)+row*rrr; cq2(2)=cq2(2)+row*p2; cq4(2)=cq4(2)+row*p4
          ! ig - energy contributions
          rt   = rn + rp;     rt2  = rt*rt
          tnt  = taupj(ihli,i,1); tpt  = taupj(ihli,j,2); tt = tnt + tpt
          dn   = dropj(ihli,i,1); dp   = dropj(ihli,j,2); dt = dn + dp
          akn  = akapj(ihli,i,1); akp  = akapj(ihli,j,2)
          akn2 = akn*akn;         akp2 = akp*akp
          adn  = akn*rn;          adp  = akp*rp
          ! ig-Pairing energy and delta
          rsa0 = (rt/rho_c)
          dd1n = CpV0(0)*(ONE-rsa0*CpV1(0))*whl
          dd1p = CpV0(1)*(ONE-rsa0*CpV1(1))*whl
          ceptn = ceptn + dd1n*akn2; ceptp = ceptp + dd1p*akp2
          cdeln = cdeln - dd1n*adn;  cdelp = cdelp - dd1p*adp
          ! Breakdown of energy contributions
          cektn = cektn + hb0n*whl*tnt; cektp = cektp + hb0p*whl*tpt   !kinetic (n and p)
          ev3 = ev3 + (tv1*rt2 - tv2*rnp2)*whl                         !volume
          eva = eva + (tv3*rt2-tv4*rnp2)*rt**sigma*whl
          ev5 = ev5 + (tv5*rt*tt + tv6*(rn*tnt + rp*tpt))*whl
          es5 = es5 + (ts1*rt*dt + ts2*(rn*dn + rp*dp))*whl            !surface
          eso = eso + (CrdJ(0)*rt+CrdJ(1)*rnp1)*djpj(ihli,i,1)*whl     !spin-orbit n
          eso = eso + (CrdJ(0)*rt-CrdJ(1)*rnp1)*djpj(ihli,j,2)*whl     !spin-orbit p
          If(icou >= 1 .Or. icou == -2) ecodi = ecodi + half*coupj(ihli,j)*rp*whl !Coulomb direct, substitution method
          If(icou == 2 .Or. icou == -3) ecoex = ecoex - cex*rp**t4o3*whl          !Coulomb exchange, Slater approximation
          If(use_j2terms) Then
             SFIZN=SFIZpj(IHLI,i,1); SFIRN=SFIRpj(IHLI,i,1); SZFIN=SZFIpj(IHLI,i,1); SRFIN=SRFIpj(IHLI,i,1)
             SFIZP=SFIZpj(IHLI,j,2); SFIRP=SFIRpj(IHLI,j,2); SZFIP=SZFIpj(IHLI,j,2); SRFIP=SRFIpj(IHLI,j,2)
             ETENS=ETENS+whl*(ta7*(SZFIN**2+SFIZN**2+SRFIN**2+SFIRN**2+SZFIP**2+SFIZP**2+SRFIP**2+SFIRP**2)&
                             +ta8*(SZFIN*SZFIP+SFIZN*SFIZP+SRFIN*SRFIP+SFIRN*SFIRP))
          End If
       End Do !ihli
       ! Intel compiler with -O2, -O3 does not like updating array element directly
       ! (like: cept(1)=cept(1)+...), hence the need for this construct
       cept(1) = ceptn; cept(2) = ceptp
       cdel(1) = cdeln; cdel(2) = cdelp
       cekt(1) = cektn; cekt(2) = cektp

       evol     = ev3 + eva + ev5; esurf = es5; ecoul = ecodi  + ecoex*CExPar
       cekt(3)  = cekt(1) + cekt(2); cept(3)  = cept(1) + cept(2)
       cetot    = cekt(3) + evol + esurf + eso + ecoul + cept(3)+ ETENS
       cdel(1)  = cdel(1)/tz(1); cdel(2)  = cdel(2)/tz(2)
       !------------------------------------------------
       ! half-projected energies required for the matrix elements
       !------------------------------------------------
       epj(i,1) = epj(i,1) + cetot*pjk(j,2)
       epj(j,2) = epj(j,2) + cetot*pjk(i,1)
       !------------------------------------------------
       ! for constraint contributions to half-projected energies
       !------------------------------------------------
       If(icstr /= 0) cq2pj(i,j)=two*(cq2(1)+cq2(2))
       !------------------------------------------------
       ! projected energies
       !------------------------------------------------
       evolpj = evolpj + pijk*evol;    esurpj = esurpj + pijk*esurf;   espopj = espopj + pijk*eso
       e_rhorho = e_rhorho + pijk*(ev3+eva); e_rhotau = e_rhotau + pijk*ev5
       epa1pj = epa1pj + pijk*cept(1); epa2pj = epa2pj + pijk*cept(2); epirpj = epa1pj + epa2pj
       ede1pj = ede1pj + pijk*cdel(1); ede2pj = ede2pj + pijk*cdel(2)
       ecdipj = ecdipj + pijk*ecodi;   ecexpj = ecexpj + pijk*ecoex;   ecoupj = ecdipj + ecexpj
       ept1pj = ept1pj + pijk*cept(1); ept2pj = ept2pj + pijk*cept(2); epotpj = ept1pj + ept2pj
       eki1pj = eki1pj + pijk*cekt(1); eki2pj = eki2pj + pijk*cekt(2); ekinpj = eki1pj + eki2pj
       etotpj = etotpj + pijk*cetot
       etenspj= etenspj+ pijk*etens
       !------------------------------------------------
       ! unprojected hfb total energy and constraint
       !------------------------------------------------
       If(i == 1 .And. j == 1) Then
          rehfbcan=Real(cetot,Kind=pr)
       End If

       ! projected particle numbers, rms, deformations
       If(j == 1) Then
          xnpj(1)  = xnpj(1)  + pjk(i,1)*cxn(1)
          rmspj(1) = rmspj(1) + pjk(i,1)*crms(1)
          q2pj(1)  = q2pj(1)  + pjk(i,1)*cq2(1)
          q4pj(1)  = q4pj(1)  + pjk(i,1)*cq4(1)
       End If
       If(i == 1) Then
          xnpj(2)  = xnpj(2)  + pjk(j,2)*cxn(2)
          rmspj(2) = rmspj(2) + pjk(j,2)*crms(2)
          q2pj(2)  = q2pj(2)  + pjk(j,2)*cq2(2)
          q4pj(2)  = q4pj(2)  + pjk(j,2)*cq4(2)
       End If

    End Do !k

    ! Real quantities at the end

    !------------------------------------------------
    ! Energies
    !------------------------------------------------
    rdelta(1) = Real(ede1pj,Kind=pr); rdelta(2) = Real(ede2pj,Kind=pr); retotpj   = Real(etotpj,Kind=pr);
    repair(1) = Real(epa1pj,Kind=pr); repair(2) = Real(epa2pj,Kind=pr); repair(3) = Real(epirpj,Kind=pr)
    rekin(1)  = Real(eki1pj,Kind=pr); rekin(2)  = Real(eki2pj,Kind=pr); rekin(3)  = Real(ekinpj,Kind=pr)
    revolpj   = Real(evolpj,Kind=pr); resurpj   = Real(esurpj,Kind=pr); respopj   = Real(espopj,Kind=pr);
    recdipj   = Real(ecdipj,Kind=pr); recexpj   = Real(ecexpj,Kind=pr); retenspj  = Real(etenspj,Kind=pr);
    depnp = retotpj - rehfbcan  !correlation energy due to projection
    !------------------------------------------------
    ! expectation values of multipole moments
    !------------------------------------------------
    Call moments_computeValue()
    !------------------------------------------------
    ! rms and deformations
    !------------------------------------------------
    Do it=itmin,itmax
       xn(it) = Real(xnpj(it),Kind=pr)
       rms(it)= Sqrt(Real(rmspj(it),Kind=pr)/xn(it))
       q2(it) = two*Real(q2pj(it),Kind=pr)    !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
       q4(it) = ffdef4*Real(q4pj(it),Kind=pr) !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
       def(it)= Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
    End Do
    r212    = rms(1)**2; r222 = rms(2)**2
    rms(3)  = Sqrt((xn(1)*r212+xn(2)*r222)/amas)
    q2(3)   = q2(1) + q2(2)  ! quadrupole moment
    q4(3)   = q4(1) + q4(2)  ! hexadecapole moment
    def(3)  = Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
    !------------------------------------------------
    ! other definitions of the same quantities
    !------------------------------------------------
    bet2(1) = ffdef6*q2(1)/(xn(1)*r02) !beta_n=Qn*Sqrt(5Pi)/(3N x^2)
    bet2(2) = ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
    bet2(3) = ffdef6*q2(3)/(amas*r02)
    het4(1) = ffdef7*q4(1)/(xn(1)*r04)
    het4(2) = ffdef7*q4(2)/(xn(2)*r04)
    het4(3) = ffdef7*q4(3)/(amas*r04)
    xn(3)   = xn(1) + xn(2)
    bet = def(3)
    !------------------------------------------------
    !  constraint constants and contributions to half-projected energies
    !------------------------------------------------
    If(icstr /= 0) Then
       cx=0.0_pr
       If(numberCons > 0) Then
          Do icons=1,numberCons
             lambda=multLambda(icons)
             If(lambda >= 1) Then
                cx = cx - multLag(lambda)*(qmoment(lambda,3)-multRequested(lambda))
             End If
             If(lambda == 0) Then
                cx = cx - neckLag*(neckValue-neckRequested)
             End If
          End Do
       End If
       Do i=1,ilpjn
          Do j=1,ilpjp
             epj(i,1) = epj(i,1) + cx*cq2pj(i,j)*pjk(j,2)
             epj(j,2) = epj(j,2) + cx*cq2pj(i,j)*pjk(i,1)
          End Do
       End Do
    End If

    If(lpr) Then
      !rc=Sqrt(r222+0.640_pr)
      ! Charge radius, from Adv. Nucl. Phys. 8, 219 (1975)
      RpRMSsq=0.769_pr
      RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
      DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
      rc = Sqrt(r222 + RpRMSsq + (xn(1)/xn(2))*RnRMSsq + DarwinFoldy)
      ! transitions to barn,barn^2,barn^4
      Do i=1,3
         q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
      End Do

      ! STORE to projected buffer 'eresj'
      ! ieresj=50 from module definitions
      ! ' si ','JININ'
      eresj(1)=si; eresj(2)=inin;
      ! ' A','   N ','   Z '
      eresj(3)=npr(1)+npr(2); eresj(4)=npr(1); eresj(5)=npr(2);
      ! ' Jln ',' Jlp '
      eresj(6)=alast(1); eresj(7)=alast(2);
      ! ,'JEtot','Jbett','Jbetn','Jbetp',' JQt ',' JQn ',' JQp '  &
      eresj(8)=retotpj; eresj(9)=def(3); eresj(10)=def(1); eresj(11)=def(2);
      eresj(12)=q2(3); eresj(13)=q2(1); eresj(14)=q2(2);
      ! ' JpEn',' JpEp',' JpDn',' JpDp',' JAsn',' JAsp'  &
      eresj(15)=repair(1); eresj(16)=repair(2);
      eresj(17)=rdelta(1); eresj(18)=rdelta(2); eresj(19)=ass(1); eresj(20)=ass(2);
      ! ,' Jrt ',' Jrn ',' Jrp ',' Jrc ',' Jht ',' Jhn ',' Jhp '  &
      eresj(21)=rms(3); eresj(22)=rms(1); eresj(23)=rms(2); eresj(24)=rc;
      eresj(25)=het4(3); eresj(26)=het4(1); eresj(27)=het4(2);
      ! ,' Jqht',' Jqhn',' Jqhp'  &
      eresj(28)=q4(3); eresj(29)=q4(1); eresj(30)=q4(2);
      ! ,' JKINt',' JKINn','JKINp',' JSO ','JCDIR',' JCEX','JDisn','JDisp'  &
      eresj(31)=rekin(3); eresj(32)=rekin(1); eresj(33)=rekin(2); eresj(34)=respopj;
      eresj(35)=recdipj; eresj(36)=recexpj; eresj(37)=Dispersion(1); eresj(38)=Dispersion(2);
      ! ,'JV2Mn','JV2Mp','JILST','JKIND','  JL '  &
      eresj(39)=v2min(1); eresj(40)=v2min(2)
      eresj(41)=iLST; eresj(42)=kindhfb; eresj(43)=iLpj;
      !  ,'JECMPAV','JECMPAV','JECMPAV'
      eresj(44)=ECMPAV(3); eresj(45)=ECMPAV(1); eresj(46)=ECMPAV(2);
      ! 'JA','JN',JZ'
      eresj(47)=Nint(xn(3)); eresj(48)=Nint(xn(1)); eresj(49)=Nint(xn(2));
      ! 'iter'
      eresj(50)=iiter
      ! nucleus with wrong asymptotic
      If(iasswrong(3) /= 0) eresj(21)=-eresj(21)

      ! WRITE to screen 'lout' and tape akzout.dat 'lfile'
      If(do_print == 1) Then
         Do iw=lout,lfile
            Write(iw,*)
            Write(iw,'(a,9x,a,/)')            '  NB! From expectpj (PNP PAV RESULTS)'
            Write(iw,*)
            If(iLST1 /= 0) Write(iw,'(a,6f15.6)') '  hfb decay const. ass ',ass
            Write(iw,'(a,8f15.6)') '  pairing: CpV0,CpV1,pwi... ',CpV0,CpV1,pwi
            Write(iw,'(a,a,a,2i5)') '  forces:   ',skyrme,',  Gauge points:',ilpjn,ilpjp
            If(keyblo(1) /= 0) Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ', bloblo(keyblo(1),1)
            If(keyblo(2) /= 0) Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ', bloblo(keyblo(2),2)
            Write(iw,*)
            Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
            Write(iw,'(a,6f15.6)') '  Requested part.numbs.',tz,Sum(tz)
            Write(iw,'(a,6f15.6)') '  Projected part.numbs.',xn
            Write(iw,'(a,3f15.6)') '  Dispersion dN2 ......',Dispersion
            Write(iw,'(a,6f15.6)') '  b0, bz, bp ..........',b0,bz,bp
            Write(iw,*)
            Write(iw,'(a,6f15.6)') '  lambda (ala) ........',ala
            Write(iw,'(a,6f15.6)') '  Lambda (alast) ......',alast
            Write(iw,'(a,6f15.6)') '  delta(n,p) ..........',rdelta
            Write(iw,'(a,6f15.6)') '  pairing energy ......',repair
            Write(iw,*)
            Write(iw,'(a,6f15.6)') '  rms-radius ..........',rms
            Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
            Write(iw,'(a,6f15.6)') '  deformation beta2 ...',def
            Write(iw,'(a,6f15.6)') '  quadrupole moment[b] ',q2
            Write(iw,'(a,6f15.6)') '  hexadecapole moment .',q4
            Write(iw,*)
            Write(iw,'(a,6f15.6)')     '  kinetic energy ......',rekin
            Write(iw,'(a,6f15.6)')     '  cmc-diagonal part ...',rekin(1)/hb0n*hbzeron-rekin(1),rekin(2)/hb0p*hbzerop-rekin(2), &
                                                                 rekin(1)/hb0n*hbzeron-rekin(1)+rekin(2)/hb0p*hbzerop-rekin(2)
            Write(iw,'(a,6f15.6)')     '  cmc-PAV .............',ECMPAV
            Write(iw,*)
            Write(iw,'(a,30x,6f15.6)') '  volume energy .......',revolpj
            Write(iw,'(a,30x,6f15.6)') '        rho_tau .......',Real(e_rhotau)
            Write(iw,'(a,30x,6f15.6)') '        rho_rho .......',Real(e_rhorho)
            Write(iw,'(a,30x,6f15.6)') '  volume energy .......',revolpj
            Write(iw,'(a,30x,6f15.6)') '  surface energy ......',resurpj
            Write(iw,'(a,30x,6f15.6)') '  spin-orbit energy ...',respopj
            Write(iw,'(a,30x,6f15.6)') '  coulomb direct ......',recdipj
            Write(iw,'(a,30x,6f15.6)') '  coulomb exchange ....',recexpj
            Write(iw,'(a,30x,6f15.6)') '  tensor energy .......',retenspj
            Write(iw,*)
            Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(qp) ....',ehfb
            Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(can,pj).',rehfbcan
            Write(iw,'(a,30x,f15.6)')  '  ehfb(qp)-ehfb(can,pj)',ehfb-rehfbcan
            Write(iw,'(a,30x,f15.6)')  '  Epj-ehfb(can,pj) ....',depnp
            Write(iw,'(a,30x,6f15.6)') '  Energy: Epj=E(PAV) ..',retotpj
            Write(iw,*)
         End Do
      End If

    End If

    ! Allocate necessary arrays
    Call deallocate_arrays_PNP()

  End Subroutine expectpj
  !=======================================================================
  !> This routine computes the gauge-dependent densities.
  !=======================================================================
  Subroutine densitpj()
    Use UNEDF, Only : CpV0,CpV1,sigma,use_j2terms,CExPar,hbzeron
    Use HFBTHO_Gauss, Only : ngl,xl,ngh,xh,nghl,wdcori
    Use HFBTHO_canonical, Only : vk,ddc
    Use HFBTHO, Only : varmas,itmin,itmax,ndx,nkx,npr,nb,id,ia,Parity,nbx,ka,kd,blocanon, &
                       y_opt,fiu,fiuz,fiur,fiud2n,fid,fidz,fidr,fidd2n,nqx,nl,nr,nz,ns,QHLA_opt, &
                       FI2D_opt,FI1Z_opt,FI1R_opt,varmasNZ,dnfactor,drhoi,ala1,inner,pjmassNZ,do_print

    Complex(pr) :: v2ig,vig_uig,dig,sumsum
    Complex(pr), Allocatable :: ank1(:,:),pfiun1(:,:),pfidn1(:,:)
    Complex(pr), Allocatable :: pakapj(:),propj(:), pdjpj(:), ptaupj(:),pdropj(:)
    Complex(pr), Allocatable :: pszfipj(:),psfizpj(:),psrfipj(:),psfirpj(:)
    Complex(pr), Pointer:: ppjk(:),pcpj(:,:),prpj(:,:),pypj(:,:)
    Complex(pr), Allocatable :: pupj(:,:)
    Real(pr), Pointer :: xxg_norm(:),xxg_phase(:)
    Real(pr) :: f,s,sd,su,sud,y,y2,sml2,cnzaa,cnraa,u,v2,tauin,xxx,yyy
    Real(pr) :: aav,anik,fi1r,fi1z,fi2d,qhla
    Real(pr) :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,xlampy
    Real(pr) :: tfiu,tfid,tfiur,tfidr,tfiuz,tfidz,tfiud2,tfidd2,TW_T
    Real(pr), Allocatable :: an2(:),ank2(:),pfiun2(:),pfidn2(:)
    Integer(ipr) :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus,n12
    Integer(ipr) :: imen,ib,m,im,ig,jg,it,j,jj,ja,jn,ILIHLI,k1,k2,kkk,kky
    Integer(ipr) :: k0(2),ky(2),kk(nqx),kyk(nqx)

    Allocate(ank1(nqx,ilpjmax),pfiun1(ndx,ilpjmax),pfidn1(ndx,ilpjmax))
    Allocate(pfiun2(ndx),pfidn2(ndx),an2(nqx),ank2(nqx))
    Allocate(pakapj(ilnghl),propj(ilnghl), pdjpj(ilnghl), ptaupj(ilnghl),pdropj(ilnghl),  &
             pSZFIpj(ilnghl),pSFIZpj(ilnghl),pSRFIpj(ilnghl),pSFIRpj(ilnghl))
    Allocate(pupj(nkx,ilpjmax))

    ! initialize parameters
    varmas = zero

    Do it=itmin,itmax
       If(it == 1) ilpj=ilpjn
       If(it == 2) ilpj=ilpjp
       ! Projection grid points:  keypj=max(1,keypj); ilpj=keypj;  ilpj2=ilpj**2 !all
       ! When a value two*pi is used the results are precisely the same but the accuracy for even L is slow with
       ! increasing L. When 'pi' is used it gives regular and better convergence with respect to both odd and even L.
       xxx = number_parity*pi/Real(ilpj,Kind=pr) ! equivalent to xxx = two*pi/Real(ilpj)
       Do i=1,ilpj
          yyy          = Real(i-1,Kind=pr)*xxx
          phypj(i)     = unit_i*yyy            ! phi_l = pi*l/L
          sinphy(i)    = unit_i*Sin(yyy)       ! i sin(phi)
          exp1iphy(i)  = Exp(unit_i*yyy)       ! exp(+i*phi)
          exp1iphym(i) = Exp(-unit_i*yyy)      ! exp(-i*phi)
          exp2iphy(i)  = Exp(two*unit_i*yyy)   ! exp(+2i*phi)
          exp2iphym(i) = Exp(-two*unit_i*yyy)  ! exp(-2i*phi)
       End Do

       ! zero for densities
       Do J=1,ilnghl
          pakapj(J)=zero; propj(J)=zero; pdjpj(J)=zero; ptaupj(J)=zero; pdropj(J)=zero
          pszfipj(J)=zero; psfizpj(J)=zero; psrfipj(J)=zero; psfirpj(J)=zero
       End Do

       ! it-pointers
       prpj => rpj(:,:,it); pcpj => cpj(:,:,it)
       pypj => ypj(:,:,it)
       ppjk => pjk(:,it)
       xxg_norm => xg_norm(:,it); xxg_phase => xg_phase(:,it)

       ! null for all pointers
       prpj=zero; pupj=zero; pcpj=zero; pypj=zero; ppjk=one

       ! particle-init (kkk-even: 2 x number of pairs)
       kkk=npr(it); If(kkk /= 2*(kkk/2)) kkk=npr(it)-1
       ppjk(1:ilpj) = exp1iphym(1:ilpj)**kkk !  exp(-i*N*phi_l)
       ! Initialization of the norm overlap x(phi) calculation
       Do ig=1,ilpj
          xxg_norm(ig)  = one
          xxg_phase(ig) = Atan2(Aimag(ppjk(ig)),Real(ppjk(ig)))
       End do

       ! start blocks
       k0(it)=0; ky(it)=0
       Do ib=1,nb
          nd=id(ib); im=ia(ib)
          If(Parity) Then
             LAPLUS=(ib+1)/2 !Yesp
          Else
             LAPLUS=ib       !Nop
          End If
          xlap=laplus; xlap2=xlap*xlap; xlam=xlap-one; xlam2=xlam*xlam

          ! charge block can quantities
          m=ib+(it-1)*nbx; k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=0
          ! Below the pwi cut-off
          If(k1 <= k2) Then
             imen = nd
             Do k = 1,nd
                k0(it) = k0(it) + 1; kk(k)  = k0(it); kkk = k0(it)
                ky(it) = ky(it) + 1; kyk(k) = ky(it); kky = ky(it)
                aav    = vk(kkk,it) ! v_k^2
                If(aav > 1.0e-12_pr) Then
                   Do ig=1,ilpj
                      v2ig = exp2iphy(ig)*aav ! = v_k^2 * exp(+2i*phi)
                      dig  = one - aav + v2ig ! = u_k^2 + v_k^2 * exp(+2i*phi)
                      If(kkk /= blocanon(it)) Then
                         ppjk(ig) = ppjk(ig)*dig  ! y(ig,it) <<<<<
                         xxg_norm(ig)  = xxg_norm(ig)*Abs(dig)
                         xxg_phase(ig) = xxg_phase(ig) + Atan2(Aimag(dig),Real(dig))
                      End If
                      prpj(kkk,ig) = v2ig/dig ! rho_k(phi), Eq.(74) in CPC 167 (2005) 43–63
                      pupj(kkk,ig) = exp1iphy(ig)*Sqrt(aav*(one-aav))/dig ! \tilde{rho}_k(phi), Eq.(74) in CPC 167 (2005) 43–63
                      pcpj(kky,ig) = exp2iphy(ig)/dig ! C(phi)
                      pypj(kky,ig) = exp1iphy(ig)/dig*unit_i*Sin(phypj(ig)/unit_i) ! sinphy(ig) !Y(mu,ig,it)
                   End Do
                End If
             End Do
             ! Cryptic left-over comment from Mario
             !   At this point density (and related) are strictly equivalent in qp- and can-representation
             !   (up to 10^-14). Pairing density is not so strict (up to 10^-5) due to uv from v^2 but
             !   pairing density is taken directly in qp representation so both representations
             !   qp and can are strictly exact (up to 10^-14)
             !
             ! Getting canonical wave function
             j=0
             Do jj = 1,nd
                Do k = 1,nd
                   j=j+1; n12 = jj+(k-1)*nd
                   an2(j) = ddc(jj,kk(k),it)
                End Do
             End Do
          Else
             ! Above the pwi cut-off (NB! Attention)
             ! here imem=0 and the contribution does not enter the densities but the Hamiltonian matrix
             ! used only in VAP regime
             ky(it)=ky(it)+1; kky = ky(it)
             Do ig=1,ilpj
                pcpj(kky,ig) = exp2iphy(ig)
                pypj(kky,ig) = exp1iphy(ig)*sinphy(ig)
             End Do
          End If

          ! calculate the densities only below the PWI cutoff
          If(imen > 0) Then
             ! gauss integration points
             Do il=1,ngl
                v2 = half/xl(il)
                Do ih=1,ngh
                   ihil = ih + (il-1)*ngh; ilihli=(ihil-1)*ilpj
                   u = xh(ih); y = y_opt(ihil); y2=y*y
                   xlamy=xlam*y; xlamy2=xlam2*y2;
                   xlapy=xlap*y; xlapy2=xlap2*y2;
                   xlampy=xlamy+xlapy

                   ! initialize spin up/down funct
                   Do k=1,nd
                      fiu(k)=zero; fiuz(k)=zero; fiur(k)=zero; fiud2n(k)=zero; pfiun2(k)=zero
                      fid(k)=zero; fidz(k)=zero; fidr(k)=zero; fidd2n(k)=zero; pfidn2(k)=zero
                   End Do

                   ! scan over basis states
                   jn=0
                   Do i=1,nd
                      ja = i+im; nla = nl(ja); nra = nr(ja); nza = nz(ja); nsa = ns(ja);
                      sml2  = nla*nla; cnzaa = nza+nza+1; cnraa = nra+nra+nla+1
                      QHLA=QHLA_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
                      FI1Z=FI1Z_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil)
                      !qha   = qh(nza,ih); qla = ql(nra,nla,il); qhla = qha*qla
                      !qhl1a = qha*ql1(nra,nla,il)*v2; qh1la = qh1(nza,ih)*qla
                      !fi1z  = fp1(ihil)*qhla+fp2(ihil)*qh1la+fp3(ihil)*qhl1a
                      !fi1r  = fp4(ihil)*qhla+fp5(ihil)*qh1la+fp6(ihil)*qhl1a
                      !fi2d  = (fs1(ihil)*qh1la**2 + four*fs4(ihil)*qh1la*qhl1a        &
                      !      +  fs2(ihil)*qhl1a**2 + two*(fs5(ihil)*qh1la              &
                      !      +  fs6(ihil)*qhl1a)*qhla + ((u*u - cnzaa)*fs1(ihil)       &
                      !      +  (p14-cnraa*v2+sml2*v2*v2)*fs2(ihil)+fs3(ihil))*qhla**2 &
                      !      -  two*(fi1r**2+fi1z**2))/(two*qhla)
                      !
                      ! wave function(spin:up,down; grad:r,z,d2)
                      If(nsa > 0) Then
                         Do k=1,nd
                            jn = jn+1; anik = an2(jn)
                            fiu(k)    = fiu(k)    + anik*qhla
                            fiur(k)   = fiur(k)   + anik*fi1r
                            fiuz(k)   = fiuz(k)   + anik*fi1z
                            fiud2n(k) = fiud2n(k) + anik*fi2d
                         End Do
                      Else
                         Do k=1,nd
                            jn = jn+1; anik = an2(jn)
                            fid(k)    = fid(k)    + anik*qhla
                            fidr(k)   = fidr(k)   + anik*fi1r
                            fidz(k)   = fidz(k)   + anik*fi1z
                            fidd2n(k) = fidd2n(k) + anik*fi2d
                         End Do
                      End If
                   End Do ! i

                   ! calculate densities
                   Do k=1,nd
                      kkk =kk(k)
                      tfiu=fiu(k); tfiuz=fiuz(k); tfiur=fiur(k); tfiud2=fiud2n(k)
                      tfid=fid(k); tfidz=fidz(k); tfidr=fidr(k); tfidd2=fidd2n(k)
                      Do ig=1,ilpj
                         v2ig=prpj(kkk,ig); vig_uig=pupj(kkk,ig)
                         I=ig+ilihli
                         pakapj(I)  = pakapj(I)  +  (tfiu**2+tfid**2)*vig_uig
                         propj(I)   = propj(I)   +  (tfiu**2+tfid**2)*v2ig
                         pdjpj(I)   = pdjpj(I)   +  (tfiur*tfidz-tfidr*tfiuz+xlamy*tfiu*(tfiur-tfidz) &
                                                 -   xlapy*tfid*(tfidr+tfiuz))*v2ig
                         TW_T=(tfiur**2+tfidr**2+tfiuz**2+tfidz**2)
                         tauin      = (xlamy2*tfiu**2+xlapy2*tfid**2+TW_T)
                         ptaupj(I)  = ptaupj(I)  +   tauin*v2ig
                         pdropj(I)  = pdropj(I)  +  (TW_T + tfiu*tfiud2 + tfid*tfidd2)*v2ig
                         psrfipj(I) = psrfipj(I) + (tfiur*tfid - tfidr*tfiu)*v2ig
                         psfirpj(I) = psfirpj(I) + (tfiu*tfid*xlampy)*v2ig
                         psfizpj(I) = psfizpj(I) + (xlamy*tfiu**2 - xlapy*tfid**2)*v2ig
                         pszfipj(I) = pszfipj(I) + (tfiuz*tfid - tfidz*tfiu)*v2ig
                      End Do !ig
                   End Do !k

                End Do !ih
             End Do !il
          End If
       End Do !ib

       ! normalized pjk
       sumsum = Sum(ppjk(1:ilpj))
       ppjk(1:ilpj) = ppjk(1:ilpj)/sumsum
       ! Y minus second term of Y
       Do k=1,ky(it)
          sumsum = Sum(ppjk(1:ilpj)*pypj(k,1:ilpj))
          pypj(k,1:ilpj) = pypj(k,1:ilpj) - sumsum
       End Do

       ! norm of the projected/unprojected density
       s = zero; sd = zero; su = zero; sud = zero;
       Do ihil=1,nghl
          ilihli=(ihil-1)*ilpj
          Do ig=1,ilpj
             I=ig+ilihli
             s=s+Real(two*propj(I)*ppjk(ig)); sd=sd+Real(four*pdropj(I)*ppjk(ig))
          End Do
          I=1+ilihli
          su=su+Real(two*propj(I)); sud=sud+Real(four*pdropj(I))
       End Do
       ! print unprojected normalization
       If(do_print == 1) Then
          Do iw=lout,lfile
             Write(iw,'(2(a,2(2x,D15.8)),(a,D15.8),a,i3)') &
                  '   pj/unpj  s= ',s,su,'   pj/unpj sd= ',sd,sud,' ala1= ',ala1(it),' inner= ',inner(it)
          End Do
       End If
       varmas = varmas + su
       varmasNZ(it) = su; pjmassNZ(it) = s

       s = Real(npr(it),Kind=pr)/s; dnfactor(it) = s; drhoi(it) = sd

       Do ihil = 1,nghl
          ilihli=(ihil-1)*ilpj
          ! wdcor moves out the int.weight and multiply by the jacobian
          f = two*wdcori(ihil)
          Do ig=1,ilpj
             I=ig+ilihli
             ropj (ihil,ig,it)  = f*propj (I)
             taupj(ihil,ig,it)  = f*ptaupj(I)
             dropj(ihil,ig,it)  = f*pdropj(I)*two
             djpj (ihil,ig,it)  = f*pdjpj (I)*two
             akapj(ihil,ig,it)  = f*pakapj(I)*half
             SRFIpj(ihil,ig,it) = f*psrfipj(I)
             SFIRpj(ihil,ig,it) = f*psfirpj(I)
             SFIZpj(ihil,ig,it) = f*psfizpj(I)
             SZFIpj(ihil,ig,it) = f*pszfipj(I)
          End Do !ig
       End Do !ihil

    End Do !it

    dnfactor(3)=dnfactor(1)+dnfactor(2)

    ! Calculate and print total norm overlap
    sumsum=c_zero
    Do ig=1,ilpjn
       Do jg=1,ilpjp
          sumsum = sumsum + xg_norm(ig,1)*xg_norm(jg,2)*Exp(unit_i*(xg_phase(ig,1)+xg_phase(jg,2)))
       End Do
    End Do
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'("Norm overlap: ",2f20.14)') sumsum/ilpjn/ilpjp
       End Do
    End If

    Deallocate(ank1,pfiun1,pfidn1)
    Deallocate(pakapj,propj,pdjpj,ptaupj,pdropj,pSZFIpj,pSFIZpj,pSRFIpj,pSFIRpj)

    Call coulompj() !complex coulomb fields

  End Subroutine densitpj
  !=======================================================================
  !> This routine computes the Coulomb field (direct part) with Gauge-dependent densities
  !=======================================================================
  Subroutine coulompj()
    Use bessik, Only : besei0
    Use HFBTHO_Gauss, Only : nghl,nleg,xleg,wleg,fl,fh,wdcor
    Use HFBTHO, Only : icacoupj,Parity,vc,chargee2,debug_solver

    Integer(ipr) :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,rrr,rz1,rz2,rrz1,rrz2,xx,alpha,beta,xxx

    If(debug_solver >= 1) Call get_CPU_time('coulompj',0)

    If(icacoupj == 0) Then

       icacoupj=1

       ! For parity-breaking shapes, the Coulomb potential was incorrectly
       ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
       ! were equivalent (see also below). This bug was corrected in version
       ! 139a
       If(Parity) Then
          fac1 = one;  fac2 = one
       Else
          fac1 = zero; fac2 = two
       End If
       ! Notes:
       !   - Missing factor 2 compared to Eq. (58) CPC paper because the density
       !     ro(:,it) already contains it (see routine DENSIT) due to T-invariance
       !   - Missing factor 1/2 when applying Gauss-Legendre quadrature (from [0,1]
       !     to the proper [-1,1] interval because it will be put back in subroutine
       !     expect() and is cancelled by a factor 2 in the HF field
       !   - For conserved parity, Gauss-Hermite points are all positive, the full
       !     integral over z' is split in z'<0 and z'>0, values of z and z1 below
       !     refer to the absolute values of z' (=-z' if z'<0)
       !
       bb=50.0_pr          ! Length scale L
       beta=2.00_pr
       alpha=one/beta
       f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)

       !$OMP PARALLEL DO        &
       !$OMP& DEFAULT(NONE)     &
       !$OMP& SCHEDULE(DYNAMIC) &
       !$OMP& SHARED(nghl,fl,fh,nleg,xleg,bb,fac1,fac2,wleg,wdcor,vc,f,alpha,beta) &
       !$OMP& PRIVATE(i,r,z,k,r1,z1,rrr,rr2,zd1,zd2,rz1,rz2,rrz1,rrz2, &
       !$OMP&         xx1,j,xx,y1,s1,t,y2,vik,xxx)
       Do i=1,nghl
          r = fl(i); z = fh(i)
          Do k=1,i

             r1 = fl(k); z1 = fh(k)
             rrr = two*r*r1; rr2 = (r - r1)**2
             ! z>0 part
             zd1 = (z - z1)**2
             rz1 = rr2 + zd1
             ! z<0 part
             zd2 = (z + z1)**2
             rz2 = rr2 + zd2
             ! Gauss-Legendre integration over u from 0 to D
             xx1=zero
             Do j=1,nleg
                xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
                xxx=(one-xleg(j)**beta)**(alpha+one)
                y1=(xleg(j)/(bb*xx))**2 ! u^2
                s1=y1*rrr               ! 2 u^2 r r'
                y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
                xx1=xx1+fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
             End Do
             vik=f*xx1/bb

             vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
             vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp

          End Do  !k
       End Do  !i
       !$OMP End Parallel Do

    End If

    ! Calculation of the Coulomb field
    coupj=zero
    Do i = 1,nghl
       Do k=1,ilpj
          coupj(:,k) = coupj(:,k) + vc(:,i)*ropj(i,k,2)
       End Do
    End Do

    If(debug_solver >= 1) Call get_CPU_time('coulompj',1)

  End Subroutine coulompj
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_PNP
