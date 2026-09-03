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
!              RESTORATION OF BROKEN SYMMETRIES PACKAGE                !
!                                                                      !
! ==================================================================== !

!-----------------------------------------------------------------------
!> This module restores broken symmetries (particle number, angular
!> momentum, parity) of HFB states by employing projection techniques.
!>
!>  @author
!>  Petar Marevic
!-----------------------------------------------------------------------
Module HFBTHO_projections

  Use HFBTHO_utilities

  Implicit None

  Public HFBTHO_restore

  Private initialize_projections,initialize_wavefunctions,initialize_angle,calculate_rotation_matrix, &
          calculate_overlaps,calculate_densities,calculate_energies,finalize_angle,project,print_project, &
          finalize_projections,simplex_basis,calculate_ry

  ! New types
  Type :: ptr_to_cmplx2darray
     Complex(pr), Dimension(:,:), Allocatable :: arr
  End Type ptr_to_cmplx2darray

  ! Stride size
  Integer(ipr), Parameter :: stride_size=23
  ! Symmetry restoration meshes
  Integer(ipr), Private, Save :: maxbet,maxphi,maxj
  Integer(ipr), Private, Save :: maxphi_eff,jjstep,maxN,maxP
  Real(pr), Private, Save :: phicyl_integration_step,integration_prefactor
  Real(pr), Allocatable, Private, Save :: betabs(:),betaweight(:),phiabs(:)
  Complex(pr), Allocatable, Private, Save :: rotated_overlap(:,:,:),rotation_matrix(:,:),inverse_rotation_matrix(:,:)
  Complex(pr), Allocatable, Private, Save :: ephi(:),ephic(:),ephicN(:,:),ephicP(:,:)
  Complex(pr), Allocatable, Private, Save :: detR(:)
  ! Cylindrical coordinate phi
  Integer(ipr), Parameter :: ngphi = 31
  Real(pr), Allocatable, Private, Save  :: phicyl(:)
  ! (rperp, z) x (iphicyl) mesh
  Integer(ipr), Allocatable, Private, Save :: ihil_convert(:),iphicyl_convert(:),ihil_iphicyl_convert(:,:)
  ! V and U Matrices
  Integer(ipr), Private, Save :: kdim(4)
  Complex(pr), Allocatable, Private, Save :: VmatrixN1(:,:),VmatrixN2(:,:)
  Complex(pr), Allocatable, Private, Save :: UmatrixN1(:,:),UmatrixN2(:,:)
  Complex(pr), Allocatable, Private, Save :: VmatrixP1(:,:),VmatrixP2(:,:)
  Complex(pr), Allocatable, Private, Save :: UmatrixP1(:,:),UmatrixP2(:,:)
  ! Rotated densities
  Real(pr), Private, Save :: piu(2)
  Complex(pr), Allocatable, Private, Save :: rotated_ro(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_tau(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_dj(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_dro(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_aka(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_akac(:,:,:)
  Complex(pr), Allocatable, Private, Save :: ro_projected(:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_sroz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_sror(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_srofi(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_jz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_jr(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_jfi(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_dsroz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_dsror(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_dsrofi(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_csroz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_csror(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_csrofi(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_sfiz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_sfir(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_srfi(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_szfi(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_szz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_srz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_srr(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_szr(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_sfifi(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_stauz(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_staur(:,:,:)
  Complex(pr), Allocatable, Private, Save :: rotated_staufi(:,:,:)
  !Real(pr), Allocatable, Public, Save  :: ro_normalization(:,:)
  Type(ptr_to_cmplx2darray), Allocatable, Private, Save ::  rotated_density(:,:,:),rotated_kappa(:,:,:),rotated_kappac(:,:,:)
  Integer(ipr), Private, Save :: nt_pair
  Integer(ipr), Allocatable, Private, Save :: iosc1_contributing(:,:),iosc2_contributing(:,:),nt_contributing(:)
  Integer(ipr), Allocatable, Private, Save :: iosc1_pair(:),iosc2_pair(:),itiphi_pair1(:),itiphi_pair2(:)
  ! Rotated and projected energies
  Complex(pr), Allocatable, Private, Save  :: all_overlaps(:),all_overlaps_gthr(:),all_energies(:),all_energies_gthr(:)
  Integer(ipr), Allocatable, Private, Save  :: beta_active(:),beta_active_gthr(:)
  Real(pr), Allocatable, Private, Save :: projected_overlap(:)
  Real(pr), Allocatable, Private, Save :: projected_ekinN(:),projected_ekinP(:)
  Real(pr), Allocatable, Private, Save :: projected_ecodi(:),projected_ecoex(:)
  Real(pr), Allocatable, Private, Save :: projected_EVOL_rho_tau(:),projected_EVOL_rho_rho(:)
  Real(pr), Allocatable, Private, Save :: projected_ESURF_rho_drho(:)
  Real(pr), Allocatable, Private, Save :: projected_ESO_rho_nablaj(:)
  Real(pr), Allocatable, Private, Save :: projected_eptN(:),projected_eptP(:)
  Real(pr), Allocatable, Private, Save :: projected_ess0(:),projected_essD(:)
  Real(pr), Allocatable, Private, Save :: projected_ecurr(:),projected_esDs(:)
  Real(pr), Allocatable, Private, Save :: projected_ejNxs(:)
  Real(pr), Allocatable, Private, Save :: projected_ETENS(:),projected_esT(:)
  Real(pr), Allocatable, Private, Save :: projected_xn1(:),projected_xn2(:)
  Real(pr), Allocatable, Private, Save :: projected_rms1(:),projected_rms2(:)
  Real(pr), Allocatable, Private, Save :: projected_delN(:),projected_delP(:)
  Real(pr), Allocatable, Private, Save :: projected_NP(:,:,:),projected_NP_norm(:)
  ! Simplex basis
  Integer(ipr), Allocatable, Private, Save :: nz_sim(:),nr_sim(:),nl_sim(:)
  ! Rperp nad z in fm, with (ih,il) indices
  Real(pr), Allocatable, Private, Save :: xl_ihil(:),xh_ihil(:)
  ! Timing variables
  Integer(ipr), Private, Save :: tstart,tfinish,clock_rate
  Real(pr), Private, Save :: time_ibet(6),time_fragments(4)
  ! Coulomb
  Complex(pr), Allocatable, Private, Save :: cou_rotated(:,:)
  ! MPI
  Integer(ipr), Private, Save :: beta_size,betaphi_size

  Integer(ipr), Allocatable, Private, Save :: block_vec(:),offset(:)
#if(USE_MPI>0)
  Integer(ipr), Allocatable, Private, Save :: beta_start(:)
#endif

Contains
  !==============================================================
  !> This routine performs the restoration of symmetries.
  !==============================================================
  Subroutine HFBTHO_restore()
#if(USE_MPI>0)
    Use mpi
    Use HFBTHO, Only : team_rank,team_size,ierr_mpi,COMM_team
#else
    Use HFBTHO, Only : team_rank
#endif

    Integer(ipr) :: ibet
#if(USE_MPI>0)
    Integer(ipr) :: i,ierr
    Integer(ipr), Allocatable :: rdispls(:),recvcnts(:)
#endif

    Call initialize_projections()
    Call initialize_wavefunctions()

#if(USE_MPI>0)
    If(maxbet > 1) Then
       Do ibet=beta_start(team_rank)+1,beta_start(team_rank)+block_vec(team_rank)
          Call initialize_angle()
          Call calculate_rotation_matrix(ibet)
          Call calculate_overlaps(ibet)
          Call calculate_densities(ibet)
          Call calculate_energies(ibet)
          Call finalize_angle()
       End Do
    Else
       ! In case of only 1 angle, only team_rank == 0 executes
       If(team_rank == 0) Then
          ibet = 1
          Call initialize_angle()
          Call calculate_rotation_matrix(ibet)
          Call calculate_overlaps(ibet)
          Call calculate_densities(ibet)
          Call calculate_energies(ibet)
          Call finalize_angle()
       End If
    End If

    If(team_size == 1) Then
       all_overlaps_gthr = all_overlaps
       all_energies_gthr = all_energies
       beta_active_gthr  = beta_active
    Else
       If(maxbet > 1) Then
          ! Team leader gathers rotated overlaps and energies. Attention, depending on the number of beta angles
          ! and/or MPI tasks requested, the number of angles processed by each task may be different: each task
          ! broadcasts vectors of different lengths, proportional to beta_active (different for each task), which
          ! are concatenated in arrays []_gthr. This is why we need to use mpi_gatherv and specify the offset, or
          ! displacement, occuring from task to task. Vector recvcnts() contains the size of data for each task,
          ! and vector rdispls() the offset for each task
          Allocate(recvcnts(0:team_size-1),rdispls(0:team_size-1))
          ! Sending out overlaps: size for each task is 2*beta_size*maxphi
          Do i=0,team_size-1
             recvcnts(i) = 2*block_vec(i)*maxphi
          End Do
          rdispls(0)=0
          Do i=1,team_size-1
             rdispls(i) = rdispls(i-1) + recvcnts(i-1)
          End Do
          Call MPI_Gatherv(all_overlaps,2*beta_size*maxphi,MPI_DOUBLE_COMPLEX,all_overlaps_gthr,recvcnts,rdispls, &
                                                                         MPI_DOUBLE_COMPLEX,0,COMM_team,ierr_mpi)
          ! Sending out energies: size for each task is stride_size*beta_size*maxphi**2
          Do i=0,team_size-1
             recvcnts(i) = stride_size*block_vec(i)*maxphi**2
          End Do
          rdispls(0)=0
          Do i=1,team_size-1
             rdispls(i) = rdispls(i-1) + recvcnts(i-1)
          End Do
          Call MPI_Gatherv(all_energies,stride_size*betaphi_size,MPI_DOUBLE_COMPLEX,all_energies_gthr,recvcnts,rdispls, &
                                                                      MPI_DOUBLE_COMPLEX,0,COMM_team,ierr_mpi)
          ! Sending out angles: size for each task is beta_size
          Do i=0,team_size-1
             recvcnts(i) = block_vec(i)
          End Do
          rdispls(0)=0
          Do i=1,team_size-1
             rdispls(i) = rdispls(i-1) + recvcnts(i-1)
          End Do
          Call MPI_Gatherv(beta_active,beta_size,MPI_INTEGER,beta_active_gthr,recvcnts,rdispls, &
                                                              MPI_INTEGER,0,COMM_team,ierr_mpi)
       Else
          ! If only one angle is requested, then only team member 0 had something to do
          If(team_rank == 0) Then
             all_overlaps_gthr = all_overlaps
             all_energies_gthr = all_energies
             beta_active_gthr  = beta_active
          End If
       End If
    End If
#else
    Do ibet=1,maxbet
       Call initialize_angle()
       Call calculate_rotation_matrix(ibet)
       Call calculate_overlaps(ibet)
       Call calculate_densities(ibet)
       Call calculate_energies(ibet)
       Call finalize_angle()
    End Do

    all_overlaps_gthr = all_overlaps
    all_energies_gthr = all_energies
    beta_active_gthr  = beta_active
#endif

    If(team_rank == 0) Then
       Call project()
       Call print_project()
    End If

    Call finalize_projections()

  End Subroutine HFBTHO_restore
  !============================================================
  !> This routine initializes the symmetry restoration procedure.
  !============================================================
  Subroutine initialize_projections()
    Use HFBTHO_Gauss, Only : gauleg,ngh,ngl,nghl,xl,xh
    Use HFBTHO_multipole_moments, Only : qmoment
#if(USE_MPI==0)
    Use HFBTHO, Only : number_of_gauge_points,number_of_rotational_angles,maximal_angular_momentum, &
                       delta_neutrons,delta_protons, force_parity,AMP_is_on,PNP_is_on,do_print,nb,id, &
                       ka,kd,nt,ro_normalization,team_rank,tz,n00,beta0,skyrme,blocking_mode
#else
    Use mpi
    Use HFBTHO, Only : number_of_gauge_points,number_of_rotational_angles,maximal_angular_momentum, &
                       delta_neutrons,delta_protons,force_parity,AMP_is_on,PNP_is_on,do_print,nb,id, &
                       ka,kd,nt,ro_normalization,team_rank,tz,ia,id,ka,kd,nbx,nqp,nuv,ierr_mpi,COMM_team, &
                       team_size,row_string,RVqpN,RVqpP,RUqpN,RUqpP,KpwiN,KpwiP,n00,beta0,skyrme,blocking_mode
#endif

    Integer(ipr) :: ibet,it,kcount,ib,nd,kaib,kdib,ki,kf,k1,iphicyl,ih,il,ihil,ihil_iphicyl,iosc1,iosc2,iosc_pair,iit,iphi,iit_iphi,nn_neut,nn_prot
    Integer(ipr) :: maxbet_old
#if(USE_MPI>0)
    Integer(ipr) :: block_size,block_mod,buffer,i
#endif
    Integer(ipr) :: tstart_ibet,tfinish_ibet
    Real(pr) :: xmin,xmax,s
    Real(pr),Allocatable :: betamesh(:)
    Character(50) :: projfname

    ! We consider only projections for even-even nuclei
    If(blocking_mode(1) /= 0 .Or. blocking_mode(2) /= 0) Stop 'STOP: Blocking not allowed in projections.'

    ! Initialize timing variables
    Call system_clock(tstart,clock_rate)
    tstart_ibet=tstart
    time_ibet=zero

    ! Setting file names and opening files. In parallel mode, only open the file
    ! for the team leader, since it is the only one that will access it
#if(USE_MPI>1)
    Write(projfname, '("projections",a7,a4)') row_string,'.out'
#else
    Write(projfname, '("projections",a4)') '.out'
#endif
#if(USE_MPI==0)
    Open(lproj,file=projfname,status='unknown')
#else
    If(do_print == 1) Open(lproj,file=projfname,status='unknown')
#endif

    maxphi = number_of_gauge_points
    maxbet = number_of_rotational_angles
    maxj   = maximal_angular_momentum
    maxN   = delta_neutrons
    maxP   = delta_protons

    ! Readjustment of the number of rotational angles and gauge angles based on symmetries
    ! For AMP, there are multiple options:
    !    (i)   no parity enforcement - the code will always consider and calculate the given
    !          number of rotational angles N, independently of whether N is even or odd
    !    (ii)  parity enforcement, even number of angles N - the code will consider N angles
    !          and use symmetries to calculate N/2 of them
    !    (iii) parity enforcement, odd number of angles N - the code will consider N+1 angles
    !          and use symmetries to calculate (N+1)/2 of them
    ! For PNP, there are multiple options:
    !    (i)   if the given number of gauge angles N is odd, the code will consider N angles
    !          and use symmetries to explicitly calculate (N+1)/2 of them
    !    (ii)  if the given number of gauge angles N is even, the code will consider N+1 angles
    !          and use symmetries to explicitly calculate (N+2)/2 of them
    maxbet_old = maxbet
    If(force_parity) Then
       If(Mod(maxbet,2) == 0) Then
          maxbet = maxbet/2
       Else
          maxbet = (maxbet+1)/2
       End If
    End If
    If(AMP_is_on /= 1) Then
       maxbet = 1
    End If
    If(PNP_is_on  >  0) Then
       If(Mod(maxphi,2) == 0) maxphi=maxphi+1
       maxphi_eff=(maxphi+1)/2
    Else
       maxphi=1; maxphi_eff=1
    End If

    ! Printing (conditional on rank within team as encoded in do_print)
    If(do_print == 1) Then
       Write(lproj,'(a)')     '              ===================================================='
       Write(lproj,'(a)')     '              Restoration of broken symmetries module, HFBTHO code'
       Write(lproj,'(a)')     '              ===================================================='
       Write(lproj,'(a44,2x,i2)') '  Number of oscillator shells:              ', n00
       Write(lproj,'(a41,2x,f5.3)') '  Basis deformation beta0:               ',  beta0
       Write(lproj,'(a)') '  ---'
       If(AMP_is_on == 1) Then
          Write(lproj, '(a)') '  Angular momentum projection is on.'
          Write(lproj,'(a44,1x,i3)')   '  Number of rotational angles given:         ', number_of_rotational_angles
          Write(lproj,'(a44,1x,i3)')   '  Number of rotational angles considered:    ', maxbet_old
          Write(lproj,'(a44,1x,i3)')   '  Number of rotational angles calculated:    ', maxbet
          Write(lproj,'(a44,1x,i3)')   '  Maximal angular momentum considered:       ', maxj
       Else
          Write(lproj, '(a)') '  Angular momentum projection is off.'
       End If
       Write(lproj,'(a)') '  ---'
       If(PNP_is_on  >  0) Then
          Write(lproj, '(a)') '  Particle number projection is on.'
          If(PNP_is_on == 1) Then
             Write(lproj, '(a)') '  Mixed prescription in PNP for non-integer powers of density.'
          Else If(PNP_is_on == 2) Then
             Write(lproj, '(a)') '  Projected prescription in PNP for non-integer powers of density.'
          End If
          Write(lproj,'(a44,1x,i3)')   '  Number of gauge angles given:             ', number_of_gauge_points
          Write(lproj,'(a44,1x,i3)')   '  Number of gauge angles considered:        ', maxphi
          Write(lproj,'(a44,1x,i3)')   '  Number of gauge angles calculated:        ', maxphi_eff
       Else
          Write(lproj, '(a)') '  Particle number projection is off.'
       End If
       Write(lproj,'(a)') '  ---'
       Write(lproj,'(a44,2x,l2)')   '  Parity enforced:                          ', force_parity
       Write(lproj,'(a42,2x,a)')      '  Energy functional:                      ', skyrme
    End If

    ! Counting the number of qp states for neutrons and protons
    Do it=1,2
       kcount=0
       Do ib=1,nb
          nd=id(ib); kaib=ka(ib,it); kdib=kd(ib,it); ki=kaib+1; kf=kaib+kdib
          Do k1=ki,kf
             kcount=kcount+1
          End Do
       End Do
       kdim(it)=kcount
    End Do
    kdim(3)=nt; kdim(4)=nt

    ! Setting up AMP mesh in beta angle
    If(maxbet == 1) Then
       Allocate(betabs(1),betaweight(1))
       betabs(1)=zero; betaweight(1)=one
       integration_prefactor = 1.0_pr
    Else
       If(force_parity) Then
          integration_prefactor = 2.0_pr ! Integration from 0 to pi/2
          xmin=0.0_pr; xmax=1.0_pr
          Allocate(betabs(maxbet),betamesh(maxbet),betaweight(maxbet))
          Call gauleg(xmin,xmax,betamesh,betaweight,maxbet)
          Do ibet=1,maxbet
             betabs(ibet)=Acos(betamesh(ibet))
          End Do
          Deallocate(betamesh)
       Else
          integration_prefactor = 1.0_pr ! Integration from 0 to pi
          xmin=-1.0_pr; xmax=1.0_pr
          Allocate(betabs(maxbet),betamesh(maxbet),betaweight(maxbet))
          Call gauleg(xmin,xmax,betamesh,betaweight,maxbet)
          Do ibet=1,maxbet
             betabs(ibet)=Acos(betamesh(ibet))
          End Do
          Deallocate(betamesh)
       End If
    End If

    ! Setting up PNP mesh in phi angle
    If(Mod(maxN,2) /= 0) maxN=maxN+1 ! Number of neutrons is even
    If(Mod(maxP,2) /= 0) maxP=maxP+1 ! Number of protons is even
    Allocate(phiabs(maxphi),ephi(maxphi),ephic(maxphi),ephicN(maxphi,-maxN/2:maxN/2),ephicP(maxphi,-maxP/2:maxP/2))
    Do iphi=1,maxphi
       phiabs(iphi)=(iphi-1)*pi/maxphi
       ephi(iphi)  = Cmplx(Cos(phiabs(iphi)), Sin(phiabs(iphi)),Kind=pr)
       ephic(iphi) = Cmplx(Cos(phiabs(iphi)),-Sin(phiabs(iphi)),Kind=pr)
       Do nn_neut=-maxN/2,maxN/2
          ephicN(iphi,nn_neut)=Cmplx(Cos(phiabs(iphi)*(tz(1)+two*nn_neut)),-Sin(phiabs(iphi)*(tz(1)+two*nn_neut)),Kind=pr)
       End Do
       Do nn_prot=-maxP/2,maxP/2
          ephicP(iphi,nn_prot)=Cmplx(Cos(phiabs(iphi)*(tz(2)+two*nn_prot)),-Sin(phiabs(iphi)*(tz(2)+two*nn_prot)),Kind=pr)
       End Do
    End Do

    ! Allocating memory for different global arrays, initializing them
    Allocate(rotated_overlap(maxbet,maxphi,2),detR(maxbet))
    Allocate(projected_overlap(0:maxj),projected_ekinN(0:maxj),projected_ekinP(0:maxj),projected_ecodi(0:maxj),projected_ecoex(0:maxj),&
             projected_EVOL_rho_tau(0:maxj),projected_EVOL_rho_rho(0:maxj),projected_ESURF_rho_drho(0:maxj),&
             projected_ESO_rho_nablaj(0:maxj),projected_eptN(0:maxj),projected_eptP(0:maxj),projected_ess0(0:maxj),projected_essD(0:maxj),&
             projected_ecurr(0:maxj),projected_esDs(0:maxj),projected_ejNxs(0:maxj),projected_ETENS(0:maxj),projected_esT(0:maxj),&
             projected_xn1(0:maxj),projected_xn2(0:maxj),projected_rms1(0:maxj),projected_rms2(0:maxj),&
             projected_delN(0:maxj),projected_delP(0:maxj),projected_NP(0:maxj,-maxN/2:maxN/2,-maxP/2:maxP/2),projected_NP_norm(0:maxj))
    Allocate(nz_sim(nt),nr_sim(nt),nl_sim(nt))
    Allocate(phicyl(ngphi))
    Allocate(ihil_convert(nghl*ngphi),iphicyl_convert(nghl*ngphi),ihil_iphicyl_convert(nghl,ngphi))
    Allocate(iosc1_pair(nt*(nt+1)/2),iosc2_pair(nt*(nt+1)/2))
    Allocate(itiphi_pair1(4*maxphi_eff),itiphi_pair2(4*maxphi_eff))

    rotated_overlap=c_zero;detR=c_zero
    projected_overlap=zero; projected_ekinN=zero; projected_ekinP=zero; projected_ecodi=zero; projected_ecoex=zero; projected_EVOL_rho_tau=zero
    projected_EVOL_rho_rho=zero; projected_ESURF_rho_drho=zero; projected_ESO_rho_nablaj=zero; projected_eptN=zero; projected_eptP=zero
    projected_ess0=zero; projected_essD=zero; projected_ecurr=zero; projected_esDs=zero; projected_ejNxs=zero; projected_ETENS=zero; projected_esT=zero
    projected_xn1=zero; projected_xn2=zero; projected_rms1=zero; projected_rms2=zero
    projected_delN=zero; projected_delP=zero; projected_NP=zero; projected_NP_norm=zero

    ! Setting up the simplex basis
    Call simplex_basis(.False.)

    ! Setting up mesh in cylindrical angle phi
    If(ngphi == 1) Then
       phicyl(1)=zero
    Else
       Do iphicyl=1,ngphi
          phicyl(iphicyl)=(iphicyl-1)*two*pi/(dble(ngphi)-one)
       End Do
    End If
    phicyl_integration_step=two*pi/(dble(ngphi)-one)

    ! Setting up conversion arrays for (z,perp) x (phicyl) coordinates
    ihil_iphicyl=0
    Do iphicyl=1,ngphi
       Do ihil=1,nghl
          ihil_iphicyl=ihil_iphicyl+1
          ihil_convert(ihil_iphicyl)=ihil
          iphicyl_convert(ihil_iphicyl)=iphicyl
          ihil_iphicyl_convert(ihil,iphicyl)=ihil_iphicyl
       End Do
    End Do

    ! Determining relevant angular momenta
    If(force_parity) Then
       jjstep=2
    Else
       jjstep=1
       If(Abs(qmoment(3,3)) < 1.d-6) jjstep=2 ! No octupole deformation
    End If

    ! Determining pairs of oscillator states
    iosc_pair=0
    Do iosc1=1,nt
       Do iosc2=iosc1,nt
          iosc_pair=iosc_pair+1
          iosc1_pair(iosc_pair)=iosc1;iosc2_pair(iosc_pair)=iosc2
       End Do
    End Do
    nt_pair=iosc_pair

    ! Determining iit-iphi pairs
    iit_iphi=0
    Do iit=1,4
       Do iphi=1,maxphi_eff
          iit_iphi=iit_iphi+1
          itiphi_pair1(iit_iphi)=iit
          itiphi_pair2(iit_iphi)=iphi
       End Do
    End Do

    ! Setting up xl_ihil, xh_ihil meshes
    If(Allocated(xl_ihil)) Deallocate(xl_ihil,xh_ihil)
    Allocate(xl_ihil(nghl),xh_ihil(nghl))
    Do ih=1,ngh
       Do il=1,ngl
          ihil=ih+(il-1)*ngh
          xl_ihil(ihil)=xl(il); xh_ihil(ihil)=xh(ih)
       End Do
    End Do

     ! MPI and broadcasting
#if(USE_MPI>0)
    ! This constructs gets the proper number of beta angles for each MPI task (not necessarily the same!)
    If(maxbet > 1) Then
       block_size = maxbet/team_size
       block_mod  = Mod(maxbet,team_size)
       Allocate(block_vec(0:team_size-1))
       If(block_mod > 0) Then
          buffer = block_mod
          Do i=0,team_size-1
             If(buffer > 0) Then
                block_vec(i) = block_size + 1
                buffer=buffer-1
             Else
                block_vec(i) = block_size
             End If
          End Do
       Else
          Do i=0,team_size-1
             block_vec(i) = block_size
          End Do
       End If
       Allocate(beta_start(0:team_size-1),offset(0:team_size-1))
       beta_start(0)=0; offset(0)=0
       Do i=1,team_size-1
          beta_start(i) = beta_start(i-1) + block_vec(i-1)
          offset(i) = offset(i-1) + block_vec(i-1)
       End Do
       beta_size    = block_vec(team_rank)
       betaphi_size = beta_size*maxphi**2 ! # of beta and phi angles calculated by each single process
    Else
       Allocate(block_vec(0:team_size-1),beta_start(0:team_size-1),offset(0:team_size-1))
       block_vec = 0; beta_start = 0; offset = 0
       block_vec(0) = 1 ! AMP with one angle, one team members processes it
       beta_size = block_vec(team_rank)
       betaphi_size = maxphi**2 ! # of beta and phi angles calculated by each single process
    End If
#else
    Allocate(offset(0:0),block_vec(0:0))
    offset       = 0
    beta_size    = maxbet
    block_vec    = beta_size
    betaphi_size = maxbet*maxphi**2
#endif

    If(team_rank == 0) Then
       Do it=1,2
          s=two*sum(ro_normalization(:,it))
          piu(it)=tz(it)/s
       End Do
    End If

#if(USE_MPI>0)
    Call mpi_bcast(piu,2,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(RVqpN,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(RVqpP,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(RUqpN,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(RUqpP,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(KpwiN,nqp,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(KpwiP,nqp,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(ia,nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(id,nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(ka,2*nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(kd,2*nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(kdim,4,mpi_integer,0,COMM_team,ierr_mpi)
#endif

    Allocate(all_overlaps(2*beta_size*maxphi),all_energies(stride_size*betaphi_size),beta_active(beta_size))
    all_overlaps=c_zero; all_energies=c_zero; beta_active=1
    If(team_rank == 0) Then
       Allocate(all_overlaps_gthr(2*maxbet*maxphi)); all_overlaps_gthr=c_zero
       Allocate(all_energies_gthr(stride_size*maxbet*maxphi*maxphi)); all_energies_gthr=c_zero
       Allocate(beta_active_gthr(maxbet))
    Else
       Allocate(all_overlaps_gthr(1:1)); all_overlaps_gthr=c_zero
       Allocate(all_energies_gthr(1:1)); all_energies_gthr=c_zero
       Allocate(beta_active_gthr(1:1))
    End If

    ! Timing
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(1)=time_ibet(1)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)

  End Subroutine initialize_projections
  !=================================================================================
  !> This routine defines the non-rotated U and V matrices in the simplex-y basis.
  !=================================================================================
  Subroutine initialize_wavefunctions()
    Use HFBTHO, Only : nt,nb,id,ia,ka,kd,RVqpN,RVqpP,RUqpN,RUqpP,KpwiN,KpwiP

    Integer(ipr) :: it,kcount,kcount_full,i_uv,ib,nd,i0,kaib,kdib,ki,kf,k1,nn,nn1,iosc
    Integer(ipr) :: tstart_ibet,tfinish_ibet

    ! Timing
    Call system_clock(tstart_ibet,clock_rate)

    Allocate(VmatrixN1(1:nt,1:kdim(1)),UmatrixN1(1:nt,1:kdim(1)),VmatrixN2(1:nt,1:kdim(3)),UmatrixN2(1:nt,1:kdim(3)),&
             VmatrixP1(1:nt,1:kdim(2)),UmatrixP1(1:nt,1:kdim(2)),VmatrixP2(1:nt,1:kdim(4)),UmatrixP2(1:nt,1:kdim(4)))
    VmatrixN1=c_zero; UmatrixN1=c_zero; VmatrixN2=c_zero; UmatrixN2=c_zero
    VmatrixP1=c_zero; UmatrixP1=c_zero; VmatrixP2=c_zero; UmatrixP2=c_zero

    Do it=1,2
       If(it == 1) Then ! Neutrons
          kcount=0; kcount_full=0; i_uv=0
          Do ib=1,nb ! Loop over blocks
             nd=id(ib); i0=ia(ib); kaib=ka(ib,it); kdib=kd(ib,it); ki=kaib+1; kf=kaib+kdib
             ! Bogoliubov matrices with a cut-off
             Do k1=ki,kf ! Loop over quasiparticle states
                kcount=kcount+1
                Do nn=1,nd ! Loop over basis states within the block ib
                   nn1=KpwiN(k1)+nn; iosc=i0+nn
                   If(nl_sim(iosc) >= 0) Then ! Simplex-y components of U and V Bogoliubov matrices
                      UmatrixN1(iosc,kcount)=Cmplx(RUqpN(nn1),zero,Kind=pr); VmatrixN1(iosc,kcount)=-Cmplx(RVqpN(nn1),zero,Kind=pr)
                   Else
                      UmatrixN1(iosc,kcount)=Cmplx(zero,RUqpN(nn1),Kind=pr); VmatrixN1(iosc,kcount)=-Cmplx(zero,RVqpN(nn1),Kind=pr)
                   End If
                End Do ! nn
             End Do ! k1
             ! Full Bogoliubov matrices (for overlap calculation)
             Do k1=1,nd
                kcount_full=kcount_full+1
                Do nn=1,nd
                   i_uv=i_uv+1; iosc=i0+nn
                   If(nl_sim(iosc) >= 0) Then ! Simplex-y components of U and V Bogoliubov matrices
                      UmatrixN2(iosc,kcount_full)=Cmplx(RUqpN(i_uv),zero,Kind=pr); VmatrixN2(iosc,kcount_full)=-Cmplx(RVqpN(i_uv),zero,Kind=pr)
                   Else
                      UmatrixN2(iosc,kcount_full)=Cmplx(zero,RUqpN(i_uv),Kind=pr); VmatrixN2(iosc,kcount_full)=-Cmplx(zero,RVqpN(i_uv),Kind=pr)
                   End If
                End Do
             End Do
          End Do ! ib
       Else ! Protons
          kcount=0; kcount_full=0; i_uv=0
          Do ib=1,nb ! Loop over blocks
             nd=id(ib); i0=ia(ib); kaib=ka(ib,it); kdib=kd(ib,it); ki=kaib+1; kf=kaib+kdib
             ! Bogoliubov matrices with a cut-off
             Do k1=ki,kf ! Loop over quasiparticle states
                kcount=kcount+1
                Do nn=1,nd ! Loop over basis states within the block ib
                   nn1=KpwiP(k1)+nn; iosc=i0+nn
                   If(nl_sim(iosc) >= 0) Then ! Simplex-y components of U and V Bogoliubov matrices
                      UmatrixP1(iosc,kcount)=Cmplx(RUqpP(nn1),zero,Kind=pr); VmatrixP1(iosc,kcount)=-Cmplx(RVqpP(nn1),zero,Kind=pr)
                   Else
                      UmatrixP1(iosc,kcount)=Cmplx(zero,RUqpP(nn1),Kind=pr); VmatrixP1(iosc,kcount)=-Cmplx(zero,RVqpP(nn1),Kind=pr)
                   End If
                End Do ! nn
             End Do ! k1
             ! Full Bogoliubov matrices (for overlap calculation)
             Do k1=1,nd
                kcount_full=kcount_full+1
                Do nn=1,nd
                   i_uv=i_uv+1; iosc=i0+nn
                   If(nl_sim(iosc) >= 0) Then ! Simplex-y components of U and V Bogoliubov matrices
                      UmatrixP2(iosc,kcount_full)=Cmplx(RUqpP(i_uv),zero,Kind=pr); VmatrixP2(iosc,kcount_full)=-Cmplx(RVqpP(i_uv),zero,Kind=pr)
                   Else
                      UmatrixP2(iosc,kcount_full)=Cmplx(zero,RUqpP(i_uv),Kind=pr); VmatrixP2(iosc,kcount_full)=-Cmplx(zero,RVqpP(i_uv),Kind=pr)
                   End If
                End Do
             End Do
          End Do ! ib
       End If
    End Do ! it

    ! Timing
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(1)=time_ibet(1)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)

  End Subroutine initialize_wavefunctions
  !===============================================================
  !> This routine initializes quantities specific to each angle.
  !===============================================================
  Subroutine initialize_angle()
    Use HFBTHO_Gauss, Only : nghl
    Use HFBTHO, Only : nt

    Integer(ipr) :: it,iphi,isimplex

    ! Allocating rotation matrices
    Allocate(rotation_matrix(nt,nt),inverse_rotation_matrix(nt,nt))
    Allocate(rotated_density(2,maxphi_eff,2),rotated_kappa(2,maxphi_eff,2),rotated_kappac(2,maxphi_eff,2))
    Do it=1,2
       Do iphi=1,maxphi_eff
          Do isimplex=1,2
             Allocate(rotated_density(it,iphi,isimplex)%arr(1:nt,1:nt)); rotated_density(it,iphi,isimplex)%arr(1:nt,1:nt)=c_zero
             Allocate(rotated_kappa(it,iphi,isimplex)%arr(1:nt,1:nt)); rotated_kappa(it,iphi,isimplex)%arr(1:nt,1:nt)=c_zero
             Allocate(rotated_kappac(it,iphi,isimplex)%arr(1:nt,1:nt)); rotated_kappac(it,iphi,isimplex)%arr(1:nt,1:nt)=c_zero
          End Do
       End Do
    End Do

    Allocate(rotated_ro(maxphi_eff,nghl*ngphi,2),rotated_tau(maxphi_eff,nghl*ngphi,2),rotated_dj(maxphi_eff,nghl*ngphi,2),rotated_dro(maxphi_eff,nghl*ngphi,2),&
             rotated_aka(maxphi_eff,nghl*ngphi,2),rotated_akac(maxphi_eff,nghl*ngphi,2),ro_projected(nghl*ngphi,2),&
             rotated_sroz(maxphi_eff,nghl*ngphi,2),rotated_sror(maxphi_eff,nghl*ngphi,2),rotated_srofi(maxphi_eff,nghl*ngphi,2),&
             rotated_jz(maxphi_eff,nghl*ngphi,2),rotated_jr(maxphi_eff,nghl*ngphi,2),rotated_jfi(maxphi_eff,nghl*ngphi,2),&
             rotated_dsroz(maxphi_eff,nghl*ngphi,2),rotated_dsror(maxphi_eff,nghl*ngphi,2),rotated_dsrofi(maxphi_eff,nghl*ngphi,2),&
             rotated_csroz(maxphi_eff,nghl*ngphi,2),rotated_csror(maxphi_eff,nghl*ngphi,2),rotated_csrofi(maxphi_eff,nghl*ngphi,2),&
             rotated_sfiz(maxphi_eff,nghl*ngphi,2),rotated_sfir(maxphi_eff,nghl*ngphi,2),rotated_srfi(maxphi_eff,nghl*ngphi,2),&
             rotated_szfi(maxphi_eff,nghl*ngphi,2),rotated_szz(maxphi_eff,nghl*ngphi,2),rotated_srz(maxphi_eff,nghl*ngphi,2),&
             rotated_srr(maxphi_eff,nghl*ngphi,2),rotated_szr(maxphi_eff,nghl*ngphi,2),rotated_sfifi(maxphi_eff,nghl*ngphi,2),&
             rotated_stauz(maxphi_eff,nghl*ngphi,2),rotated_staur(maxphi_eff,nghl*ngphi,2),rotated_staufi(maxphi_eff,nghl*ngphi,2))
    Allocate(iosc1_contributing(nt*nt,maxphi_eff),iosc2_contributing(nt*nt,maxphi_eff),nt_contributing(maxphi_eff))
    Allocate(cou_rotated(maxphi_eff,nghl))

    rotation_matrix=c_zero; inverse_rotation_matrix=c_zero
    rotated_ro=c_zero; rotated_tau=c_zero; rotated_dj=c_zero; rotated_dro=c_zero
    rotated_aka=c_zero; rotated_akac=c_zero; ro_projected=c_zero
    rotated_sroz=c_zero; rotated_sror=c_zero; rotated_srofi=c_zero
    rotated_jz=c_zero; rotated_jr=c_zero; rotated_jfi=c_zero
    rotated_dsroz=c_zero; rotated_dsror=c_zero; rotated_dsrofi=c_zero
    rotated_csroz=c_zero; rotated_csror=c_zero; rotated_csrofi=c_zero
    rotated_sfiz=c_zero; rotated_sfir=c_zero; rotated_srfi=c_zero
    rotated_szfi=c_zero; rotated_szz=c_zero; rotated_srz=c_zero
    rotated_srr=c_zero; rotated_szr=c_zero; rotated_sfifi=c_zero
    rotated_stauz=c_zero; rotated_staur=c_zero; rotated_staufi=c_zero
    iosc1_contributing=0; iosc2_contributing=0; nt_contributing=0
    cou_rotated=c_zero

  End Subroutine initialize_angle
  !=======================================================================================
  !> This routine calculates the rotation matrix \f$ \exp(-i\beta\hat{j}_y) \f$ in the simplex-y basis.
  !=======================================================================================
  Subroutine calculate_rotation_matrix(ibet)
    Use math, Only : calculate_inverse_cmplx,iv
    Use HFBTHO, Only : bp,bz,nt,team_rank

    Integer(ipr) :: ibet !< - Index of the current \f$ \beta \f$ angle

    Integer(ipr) :: iosc1,iosc2,iosc_pair,nza1,nra1,nla1,nza2,nra2,nla2
    Integer(ipr) :: tstart_ibet,tfinish_ibet
    Real(pr) :: beta,rotel1,rotel3
    Real(pr), Parameter :: tol_one=one-1.d-8, tol_zero=1.d-8
    Complex(pr) :: det
    Integer(ipr) :: info

    ! Timing
    Call system_clock(tstart_ibet,clock_rate)

    ! Defining rotational angle
    beta=betabs(ibet)

    ! Loop over oscillator states
    !$OMP  PARALLEL DO       &
    !$OMP& DEFAULT(NONE)     &
    !$OMP& SCHEDULE(DYNAMIC) &
    !$OMP& SHARED(nt_pair,beta,bp,bz,nz_sim,nr_sim,nl_sim,iv,rotation_matrix,iosc1_pair,iosc2_pair) &
    !$OMP& PRIVATE(iosc_pair,iosc1,nza1,nra1,nla1,iosc2,nza2,nra2,nla2,rotel1,rotel3)
    Do iosc_pair=1,nt_pair
       iosc1=iosc1_pair(iosc_pair); nza1=nz_sim(iosc1); nra1=nr_sim(iosc1); nla1=nl_sim(iosc1)
       iosc2=iosc2_pair(iosc_pair); nza2=nz_sim(iosc2); nra2=nr_sim(iosc2); nla2=nl_sim(iosc2)
       rotel1=zero; rotel3=zero
       call calculate_ry(beta,bp,bz,bp,bz,nza1,nra1,nla1,nza2,nra2, nla2,rotel1)
       call calculate_ry(beta,bp,bz,bp,bz,nza1,nra1,nla1,nza2,nra2,-nla2,rotel3)
       rotation_matrix(iosc1,iosc2)=Cmplx(Cos(half*beta)*rotel1,Sin(half*beta)*rotel3,Kind=pr)
       rotation_matrix(iosc2,iosc1)=iv(nla1+nla2)*Cmplx(Cos(half*beta)*rotel1,Sin(half*beta)*rotel3,Kind=pr) ! R(b,a)=iv(la+lb)*R(a,b)
    End Do ! iosc_pair
    !$OMP END PARALLEL DO

    ! Tests
    ! ------
    ! Some of the properties that the rotation matrix should satisfy
    ! when the basis is closed under rotation (spherical basis)
    !  1) det(R)=1
    !  2) R(alpha)*R(beta)=R(alpha+beta)
    !  3) (R^{-1})^T = R^* (Robledo, 1994)
    inverse_rotation_matrix=rotation_matrix

    Call calculate_inverse_cmplx(nt,inverse_rotation_matrix,det,info)
    If(info == 0 .And. Abs(det) > 1.d-300) Then
       detR(ibet)=det
       beta_active(ibet-offset(team_rank))=1
    Else
       detR(ibet)=c_zero
       beta_active(ibet-offset(team_rank))=0
    End If

    ! Test if R*R^-1 = identity
    !Call Zgemm('N','N',nt,nt,nt,c_unit,rotation_matrix,nt,inverse_rotation_matrix,nt,c_zero,identity,nt)
    !Do iosc1=1,nt
    !Do iosc2=1,nt
    !   If(iosc1 /= iosc2 .And. Zabs(identity(iosc1,iosc2)) > tol_zero) Then
    !      Write(lproj,'(2x,a29,i3,x,i3,a10,f12.8,a1,f13.9)') 'Error in R x R^{-1} element (',iosc1,iosc2,') at angle',betabs(ibet),':',Zabs(identity(iosc1,iosc2))
    !   End If
    !   If(iosc1 == iosc2 .And. Zabs(identity(iosc1,iosc2)) < tol_one) Then
    !      Write(lproj,'(2x,a29,i3,x,i3,a10,f12.8,a1,f13.9)') 'Error in R x R^{-1} element (',iosc1,iosc2,') at angle',betabs(ibet),':',Zabs(identity(iosc1,iosc2))
    !   End If
    !End Do
    !End Do

    ! Timing
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(2)=time_ibet(2)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)

  End Subroutine calculate_rotation_matrix
  !=======================================================================================
  !> This routine calculates overlaps between rotated states and rotated density matrix.
  !=======================================================================================
  Subroutine calculate_overlaps(ibet)
    Use HFBTHO, Only : AMP_is_on,PNP_is_on,nt,team_rank
    Use math, Only : calculate_inverse_cmplx

    Integer(ipr), Intent(In) :: ibet !< - Index of the current \f$ \beta \f$ angle

    Integer(ipr) :: it,ifl,iosc1,iosc2,iosc_contributing,iit,iphi,iit_iphi,isimplex
    Integer(ipr) :: tstart_ibet,tfinish_ibet
    Real(pr), Parameter :: eps=1.d-15
    Complex(pr), Allocatable :: OverlapMatrixTmp(:,:),AuxiliaryMatrix(:,:),AuxiliaryRotUMatrix(:,:),AuxiliaryDensity(:,:),aUmatrix(:,:),aVmatrix(:,:),AuxiliaryConjg(:,:)
    Type(ptr_to_cmplx2darray), Allocatable :: rotated_Vmatrix(:,:),rotated_invUmatrix(:,:),InverseOverlapMatrix(:,:,:),Vmatrix(:),Umatrix(:)
    Complex(pr) :: prefac,detA(2)

    If(beta_active(ibet-offset(team_rank)) == 0) Return

    ! Timing
    Call system_clock(tstart_ibet,clock_rate)

    !----------------------------------------- !
    !                AMP only                  !
    !----------------------------------------- !
    If(AMP_is_on == 1 .And. PNP_is_on == 0) Then

       ! Allocating local overlap arrays
       Allocate(rotated_Vmatrix(4,maxphi_eff),rotated_invUmatrix(4,maxphi_eff),InverseOverlapMatrix(4,maxphi_eff,1),Vmatrix(4),Umatrix(4))
       Do iit=1,4
          Allocate(Vmatrix(iit)%arr(1:nt,1:kdim(iit)),Umatrix(iit)%arr(1:nt,1:kdim(iit)))
          Do iphi=1,maxphi_eff
             Allocate(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit))); Allocate(rotated_invUmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)))
             isimplex=1
             Allocate(InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit)))
          End Do
       End Do
       Vmatrix(1)%arr(1:nt,1:kdim(1))=VmatrixN1(1:nt,1:kdim(1))
       Vmatrix(2)%arr(1:nt,1:kdim(2))=VmatrixP1(1:nt,1:kdim(2))
       Vmatrix(3)%arr(1:nt,1:kdim(3))=VmatrixN2(1:nt,1:kdim(3))
       Vmatrix(4)%arr(1:nt,1:kdim(4))=VmatrixP2(1:nt,1:kdim(4))
       Umatrix(1)%arr(1:nt,1:kdim(1))=UmatrixN1(1:nt,1:kdim(1))
       Umatrix(2)%arr(1:nt,1:kdim(2))=UmatrixP1(1:nt,1:kdim(2))
       Umatrix(3)%arr(1:nt,1:kdim(3))=UmatrixN2(1:nt,1:kdim(3))
       Umatrix(4)%arr(1:nt,1:kdim(4))=UmatrixP2(1:nt,1:kdim(4))
       ! Calculating overlaps
       !$OMP  PARALLEL DO       &
       !$OMP& DEFAULT(NONE)     &
       !$OMP& SCHEDULE(DYNAMIC) &
       !$OMP& SHARED(itiphi_pair1,itiphi_pair2,nt,kdim,ibet,rotation_matrix,inverse_rotation_matrix,Vmatrix,Umatrix,rotated_Vmatrix,rotated_invUmatrix,&
       !$OMP&        InverseOverlapMatrix,detR,maxphi_eff,maxphi,all_overlaps,team_rank,beta_size,offset) &
       !$OMP& PRIVATE(iit,iphi,ifl,detA,aUmatrix,aVmatrix,AuxiliaryConjg,OverlapMatrixTmp,isimplex)
       Do iit_iphi=1,4*maxphi_eff
          iit=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(aUmatrix(kdim(iit),kdim(iit)),aVmatrix(kdim(iit),kdim(iit))); aUmatrix=c_zero; aVmatrix=c_zero
          Allocate(AuxiliaryConjg(1:nt,1:kdim(iit)))
          ! Rotating V matrix
          Call Zgemm('N','N',nt,kdim(iit),nt,c_unit,rotation_matrix,nt,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,c_zero,rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt)
          AuxiliaryConjg(1:nt,1:kdim(iit))=conjg(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)))
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,c_unit,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,AuxiliaryConjg,nt,c_zero,aVmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          ! Rotating U matrix
          AuxiliaryConjg(1:nt,1:kdim(iit))=conjg(Umatrix(iit)%arr(1:nt,1:kdim(iit)))
          Call Zgemm('T','N',nt,kdim(iit),nt,c_unit,inverse_rotation_matrix,nt,AuxiliaryConjg,nt,c_zero,rotated_invUmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt)
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,c_unit,Umatrix(iit)%arr(1:nt,1:kdim(iit)),nt,rotated_invUmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt,c_unit,aUmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          !
          Allocate(OverlapMatrixTmp(kdim(iit),kdim(iit))); OverlapMatrixTmp=c_zero
          isimplex=1
          OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))=aUmatrix(1:kdim(iit),1:kdim(iit))+aVmatrix(1:kdim(iit),1:kdim(iit))
          ! Calculating inverse and determinant of the A matrix
          Call calculate_inverse_cmplx(kdim(iit),OverlapMatrixTmp,detA(isimplex),ifl)
          InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit))=OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))
          Deallocate(OverlapMatrixTmp)
          If(iit ==  3 .Or. iit == 4) Then
             all_overlaps((iit-3)*maxphi*beta_size+(iphi-1)*beta_size+ibet-offset(team_rank))=detA(1)*detR(ibet) ! Saving overlap
          End If
          Deallocate(aUmatrix,aVmatrix,AuxiliaryConjg)
       End Do
       !$OMP  End PARALLEL DO

       ! Calculating rotated densities in configuration space
       !$OMP  PARALLEL DO       &
       !$OMP& DEFAULT(NONE)     &
       !$OMP& SCHEDULE(DYNAMIC) &
       !$OMP& SHARED(nt,kdim,rotation_matrix,rotated_Vmatrix,InverseOverlapMatrix,Vmatrix,rotated_density,Umatrix,rotated_kappa,rotated_kappac,itiphi_pair1,itiphi_pair2,maxphi_eff) &
       !$OMP& PRIVATE(it,iphi,AuxiliaryMatrix,AuxiliaryRotUMatrix,AuxiliaryDensity,AuxiliaryConjg)
       Do iit_iphi=1,2*maxphi_eff
          it=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(AuxiliaryMatrix(nt,kdim(it)),AuxiliaryRotUMatrix(nt,kdim(it)),AuxiliaryDensity(nt,nt),AuxiliaryConjg(1:kdim(it),1:kdim(it)))
          ! rho density, simplex ++
          AuxiliaryConjg=conjg(InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it)))  ! a^{--} = conjg(a^{++})
          Call Zgemm('N','N',nt,kdim(it),kdim(it),c_unit,rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it)),nt,AuxiliaryConjg,kdim(it),c_zero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),c_unit,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_density(it,iphi,1)%arr(1:nt,1:nt)=AuxiliaryDensity(1:nt,1:nt)
          ! kappa density, simplex +-
          Call Zgemm('N','C',nt,nt,kdim(it),-c_unit,AuxiliaryMatrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_kappa(it,iphi,1)%arr(1:nt,1:nt)=AuxiliaryDensity(1:nt,1:nt)
          ! kappaC density, simplex -+
          Call Zgemm('N','N',nt,kdim(it),nt,c_unit,rotation_matrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryRotUMatrix,nt)
          Call Zgemm('N','N',nt,kdim(it),kdim(it),c_unit,AuxiliaryRotUMatrix,nt,AuxiliaryConjg,kdim(it),c_zero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),c_unit,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_kappac(it,iphi,2)%arr(1:nt,1:nt)=AuxiliaryDensity(1:nt,1:nt)
          ! rho density, simplex --
          rotated_density(it,iphi,2)%arr(1:nt,1:nt)=conjg(rotated_density(it,iphi,1)%arr(1:nt,1:nt))
          ! kappa density, simplex -+
          rotated_kappa(it,iphi,2)%arr(1:nt,1:nt)=-conjg(rotated_kappa(it,iphi,1)%arr(1:nt,1:nt))
          ! kappaC density, simplex +-
          rotated_kappac(it,iphi,1)%arr(1:nt,1:nt)=-conjg(rotated_kappac(it,iphi,2)%arr(1:nt,1:nt))
          Deallocate(AuxiliaryMatrix,AuxiliaryRotUMatrix,AuxiliaryDensity,AuxiliaryConjg)
       End Do
       !$OMP  End PARALLEL DO

       ! Determining non-zero density and kappa contributions
       Do iphi=1,maxphi_eff
          iosc_contributing=0
          Do iosc2=1,nt
             Do iosc1=1,nt
                If(Abs(rotated_density(1,iphi,1)%arr(iosc1,iosc2)) > eps .Or. Abs(rotated_density(2,iphi,1)%arr(iosc1,iosc2)) > eps .Or. &
                   Abs(rotated_kappa(1,iphi,1)%arr(iosc1,iosc2))   > eps .Or. Abs(rotated_kappa(2,iphi,1)%arr(iosc1,iosc2))   > eps .Or. &
                   Abs(rotated_kappac(1,iphi,1)%arr(iosc1,iosc2))  > eps .Or. Abs(rotated_kappac(2,iphi,1)%arr(iosc1,iosc2))  > eps) Then
                   iosc_contributing=iosc_contributing+1
                   iosc1_contributing(iosc_contributing,iphi)=iosc1
                   iosc2_contributing(iosc_contributing,iphi)=iosc2
                End If
             End Do
          End Do
          nt_contributing(iphi)=iosc_contributing
       End Do
       Deallocate(rotated_Vmatrix,rotated_invUmatrix,InverseOverlapMatrix,Umatrix,Vmatrix)
    ! ---------------------------------------- !
    !             PNP or PNP&AMP               !
    ! ---------------------------------------- !
    Else

       ! Allocating local overlap arrays
       Allocate(rotated_Vmatrix(4,maxphi_eff),rotated_invUmatrix(4,maxphi_eff),InverseOverlapMatrix(4,maxphi_eff,2),Vmatrix(4),Umatrix(4))
       Do iit=1,4
          Allocate(Vmatrix(iit)%arr(1:nt,1:kdim(iit)),Umatrix(iit)%arr(1:nt,1:kdim(iit)))
          Do iphi=1,maxphi_eff
             Allocate(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit))); Allocate(rotated_invUmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)))
             Do isimplex=1,2
                Allocate(InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit)))
             End Do
          End Do
       End Do
       Vmatrix(1)%arr(1:nt,1:kdim(1))=VmatrixN1(1:nt,1:kdim(1))
       Vmatrix(2)%arr(1:nt,1:kdim(2))=VmatrixP1(1:nt,1:kdim(2))
       Vmatrix(3)%arr(1:nt,1:kdim(3))=VmatrixN2(1:nt,1:kdim(3))
       Vmatrix(4)%arr(1:nt,1:kdim(4))=VmatrixP2(1:nt,1:kdim(4))
       Umatrix(1)%arr(1:nt,1:kdim(1))=UmatrixN1(1:nt,1:kdim(1))
       Umatrix(2)%arr(1:nt,1:kdim(2))=UmatrixP1(1:nt,1:kdim(2))
       Umatrix(3)%arr(1:nt,1:kdim(3))=UmatrixN2(1:nt,1:kdim(3))
       Umatrix(4)%arr(1:nt,1:kdim(4))=UmatrixP2(1:nt,1:kdim(4))

       ! Calculating overlaps
       !$OMP  PARALLEL DO       &
       !$OMP& DEFAULT(NONE)     &
       !$OMP& SCHEDULE(DYNAMIC) &
       !$OMP& SHARED(itiphi_pair1,itiphi_pair2,nt,kdim,ibet,rotation_matrix,inverse_rotation_matrix,Vmatrix,Umatrix,rotated_Vmatrix,rotated_invUmatrix,&
       !$OMP&        InverseOverlapMatrix,detR,phiabs,ephi,ephic,maxphi_eff,maxphi,all_overlaps,beta_size,team_rank,offset) &
       !$OMP& PRIVATE(iit,iphi,ifl,detA,OverlapMatrixTmp,aUmatrix,aVmatrix,AuxiliaryConjg,isimplex,prefac)
       Do iit_iphi=1,4*maxphi_eff
          iit=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(aUmatrix(kdim(iit),kdim(iit)),aVmatrix(kdim(iit),kdim(iit))); aUmatrix=c_zero; aVmatrix=c_zero
          Allocate(AuxiliaryConjg(1:nt,1:kdim(iit)))
          ! Rotating V matrix
          Call Zgemm('N','N',nt,kdim(iit),nt,c_unit,rotation_matrix,nt,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,c_zero,rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt)
          AuxiliaryConjg(1:nt,1:kdim(iit))=conjg(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)))
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,c_unit,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,AuxiliaryConjg,nt,c_zero,aVmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          ! Rotating U matrix
          AuxiliaryConjg(1:nt,1:kdim(iit))=conjg(Umatrix(iit)%arr(1:nt,1:kdim(iit)))
          Call Zgemm('T','N',nt,kdim(iit),nt,c_unit,inverse_rotation_matrix,nt,AuxiliaryConjg,nt,c_zero,rotated_invUmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt)
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,c_unit,Umatrix(iit)%arr(1:nt,1:kdim(iit)),nt,rotated_invUmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt,c_unit,aUmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          ! Loop over two simplex blocks
          Allocate(OverlapMatrixTmp(kdim(iit),kdim(iit)))
          Do isimplex=1,2
             OverlapMatrixTmp=c_zero
             If(isimplex == 1) Then
                OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))=ephic(iphi)*aUmatrix(1:kdim(iit),1:kdim(iit))+ephi(iphi)*aVmatrix(1:kdim(iit),1:kdim(iit))
             Else
                OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))=ephic(iphi)*conjg(aUmatrix(1:kdim(iit),1:kdim(iit)))+ephi(iphi)*conjg(aVmatrix(1:kdim(iit),1:kdim(iit)))
             End If
             ! Calculating inverse and determinant of the A matrix
             Call calculate_inverse_cmplx(kdim(iit),OverlapMatrixTmp,detA(isimplex),ifl)
             InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit))=OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))
          End Do
          Deallocate(OverlapMatrixTmp)
          If(iit ==  3 .Or. iit == 4) Then
             prefac=Cmplx(Cos(phiabs(iphi)*kdim(iit)),Sin(phiabs(iphi)*kdim(iit)),Kind=pr)
             all_overlaps((iit-3)*maxphi*beta_size+(iphi-1)*beta_size+ibet-offset(team_rank))=detA(1)*detR(ibet)*prefac ! saving overlap
             If(iphi > 1) Then
                all_overlaps((iit-3)*maxphi*beta_size+(maxphi-iphi+2-1)*beta_size+ibet-offset(team_rank))=conjg(detA(1)*detR(ibet)*prefac)
             End If
          End If
          Deallocate(aUmatrix,aVmatrix,AuxiliaryConjg)
       End Do
       !$OMP End PARALLEL DO

       ! Calculating rotated densities in configuration space
       !$OMP  PARALLEL DO       &
       !$OMP& DEFAULT(NONE)     &
       !$OMP& SCHEDULE(DYNAMIC) &
       !$OMP& SHARED(nt,kdim,rotation_matrix,rotated_Vmatrix,InverseOverlapMatrix,Vmatrix,rotated_density,Umatrix,rotated_kappa,rotated_kappac,itiphi_pair1,itiphi_pair2,ephi,ephic,maxphi_eff,maxphi) &
       !$OMP& PRIVATE(it,iit_iphi,iphi,AuxiliaryMatrix,AuxiliaryRotUMatrix,AuxiliaryDensity,AuxiliaryConjg)
       Do iit_iphi=1,2*maxphi_eff
          it=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(AuxiliaryMatrix(nt,kdim(it)),AuxiliaryRotUMatrix(nt,kdim(it)),AuxiliaryDensity(nt,nt),AuxiliaryConjg(1:nt,1:kdim(it)))
          ! rho density, simplex ++
          Call Zgemm('N','N',nt,kdim(it),kdim(it),c_unit,rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it)),nt,InverseOverlapMatrix(it,iphi,2)%arr(1:kdim(it),1:kdim(it)),kdim(it),c_zero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),c_unit,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_density(it,iphi,1)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! kappa density, simplex +-
          Call Zgemm('N','C',nt,nt,kdim(it),-c_unit,AuxiliaryMatrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_kappa(it,iphi,1)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! kappaC density, simplex -+
          Call Zgemm('N','N',nt,kdim(it),nt,c_unit,rotation_matrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryRotUMatrix,nt)
          Call Zgemm('N','N',nt,kdim(it),kdim(it),c_unit,AuxiliaryRotUMatrix,nt,InverseOverlapMatrix(it,iphi,2)%arr(1:kdim(it),1:kdim(it)),kdim(it),c_zero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),c_unit,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_kappac(it,iphi,2)%arr(1:nt,1:nt)=ephic(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! rho density, simplex --
          AuxiliaryConjg(1:nt,1:kdim(it))=conjg(rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it)))
          Call Zgemm('N','N',nt,kdim(it),kdim(it),c_unit,AuxiliaryConjg,nt,InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it)),kdim(it),c_zero,AuxiliaryMatrix,nt)
          Call Zgemm('N','T',nt,nt,kdim(it),c_unit,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_density(it,iphi,2)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! kappa density, simplex -+
          !Call Zgemm('N','N',nt,kdim(it),kdim(it),c_unit,AuxiliaryConjg,nt,InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it)),kdim(it),c_zero,AuxiliaryMatrix,nt)
          Call Zgemm('N','T',nt,nt,kdim(it),c_unit,AuxiliaryMatrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_kappa(it,iphi,2)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! kappaC density, simplex +-
          AuxiliaryRotUMatrix(1:nt,1:kdim(it))=conjg(AuxiliaryRotUMatrix(1:nt,1:kdim(it)))
          Call Zgemm('N','N',nt,kdim(it),kdim(it),c_unit,AuxiliaryRotUMatrix,nt,InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it)),kdim(it),c_zero,AuxiliaryMatrix,nt)
          Call Zgemm('N','T',nt,nt,kdim(it),c_unit,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,c_zero,AuxiliaryDensity,nt)
          rotated_kappac(it,iphi,1)%arr(1:nt,1:nt)=-ephic(iphi)*AuxiliaryDensity(1:nt,1:nt)
          Deallocate(AuxiliaryMatrix,AuxiliaryRotUMatrix,AuxiliaryDensity,AuxiliaryConjg)
       End Do
       !$OMP End PARALLEL DO

       ! Determining non-zero density and kappa contributions
       Do iphi=1,maxphi_eff
          iosc_contributing=0
          Do iosc2=1,nt
             Do iosc1=1,nt
                If(Abs(rotated_density(1,iphi,1)%arr(iosc1,iosc2)) > eps .Or. Abs(rotated_density(1,iphi,2)%arr(iosc1,iosc2)) > eps .Or. &
                   Abs(rotated_density(2,iphi,1)%arr(iosc1,iosc2)) > eps .Or. Abs(rotated_density(2,iphi,2)%arr(iosc1,iosc2)) > eps .Or. &
                   Abs(rotated_kappa(1,iphi,1)%arr(iosc1,iosc2))   > eps .Or. Abs(rotated_kappa(1,iphi,2)%arr(iosc1,iosc2))   > eps .Or. &
                   Abs(rotated_kappa(2,iphi,1)%arr(iosc1,iosc2))   > eps .Or. Abs(rotated_kappa(2,iphi,2)%arr(iosc1,iosc2))   > eps .Or. &
                   Abs(rotated_kappac(1,iphi,1)%arr(iosc1,iosc2))  > eps .Or. Abs(rotated_kappac(1,iphi,2)%arr(iosc1,iosc2))  > eps .Or. &
                   Abs(rotated_kappac(2,iphi,1)%arr(iosc1,iosc2))  > eps .Or. Abs(rotated_kappac(2,iphi,2)%arr(iosc1,iosc2))  > eps) Then
                   iosc_contributing=iosc_contributing+1
                   iosc1_contributing(iosc_contributing,iphi)=iosc1
                   iosc2_contributing(iosc_contributing,iphi)=iosc2
                End If
             End Do
          End Do
          nt_contributing(iphi)=iosc_contributing
       End Do
       Deallocate(rotated_Vmatrix,rotated_invUmatrix,InverseOverlapMatrix)
    End If

    ! Timing
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(3)=time_ibet(3)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)

  End Subroutine calculate_overlaps
  !==============================================
  !> This routine calculates rotated densities.
  !==============================================
  Subroutine calculate_densities(ibet)
    Use HFBTHO_Gauss, Only : nghl,wdcori
    Use HFBTHO, Only : team_rank,bp,bz,y_opt,qhla_opt,fi1z_opt,fi1r_opt,PNP_is_on,nz,nr,vc
    Use UNEDF, Only : use_j2terms

    Integer(ipr), Intent(In) :: ibet !< - Index of the current \f$ \beta \f$ angle

    Integer(ipr) :: ihil_iphicyl,ihil,iphicyl,iosc,iosc1,iosc2,it,nla1,nla2,iphi,ihil1,ihil2
    Integer(ipr) :: tstart_ibet,tfinish_ibet,istart,iend
    Complex(pr) :: ro_sumN,ro_sumP,tau_sumN,tau_sumP,dj_sumN,dj_sumP,dro_sumN,dro_sumP,aka_sumN,aka_sumP,akac_sumN,akac_sumP,&
                   sroz_sumN,sroz_sumP,sror_sumN,sror_sumP,srofi_sumN,srofi_sumP,jz_sumN,jz_sumP,jr_sumN,jr_sumP,jfi_sumN,jfi_sumP,&
                   dsroz_sumN,dsroz_sumP,dsror_sumN,dsror_sumP,dsrofi_sumN,dsrofi_sumP,csroz_sumN,csroz_sumP,csror_sumN,csror_sumP,&
                   csrofi_sumN,csrofi_sumP,sfiz_sumN,sfiz_sumP,sfir_sumN,sfir_sumP,srfi_sumN,srfi_sumP,szfi_sumN,szfi_sumP,&
                   szz_sumN,szz_sumP,srz_sumN,srz_sumP,srr_sumN,srr_sumP,szr_sumN,szr_sumP,sfifi_sumN,sfifi_sumP,&
                   stauz_sumN,stauz_sumP,staur_sumN,staur_sumP,staufi_sumN,staufi_sumP
    Complex(pr) :: fac_ro,fac_tau,fac_dj1,fac_dj2,fac_dro,fac_aka,fac_sroz,fac_sror,fac_srofi,fac_jz,fac_jr,fac_jfi,&
                   fac_dsroz,fac_dsror,fac_dsrofi,fac_csroz1,fac_csror1,fac_csror2,fac_csrofi1,fac_csrofi2,fac_sfiz,fac_sfir,&
                   fac_srfi,fac_szfi,fac_szz,fac_srz,fac_srr,fac_szr,fac_sfifi,fac_stauz,fac_staur,fac_staufi,facmul
    Complex(pr) :: ro_projN,ro_projP,ro_normN,ro_normP,overlapN,overlapP,rotated_roN,rotated_roP
    Real(pr) :: y,y2,phy
    Real(pr) :: q_h0l0a_h0l0b,q_h1l0a_h1l0b,q_h0l1a_h0l1b,q_h0l0a_h0l1b,q_h0l1a_h0l0b,q_h0l0a_h1l0b,q_h1l0a_h0l0b,q_h1l0a_h0l1b,q_h0l1a_h1l0b,&
                fac1,fac2,fac3,fac4,fac5
    Real(pr) :: bpi,bpi2,bzi,bzi2
    Complex(pr) :: densNpl,densPpl,densNmi,densPmi,kappNmi,kappPmi,kappcNmi,kappcPmi
    Complex(pr), Allocatable :: vc_pr(:,:),cou_rotated_tmp(:),rotated_ro_tmp(:)

    If(beta_active(ibet-offset(team_rank)) == 0) Return

    ! Timing
    Call system_clock(tstart_ibet,clock_rate)

    bpi=one/bp; bpi2=bpi**2; bzi=one/bz; bzi2=bzi**2

    ! Calculating various density components in the coordinate space
    !$OMP  PARALLEL DO       &
    !$OMP& DEFAULT(NONE)     &
    !$OMP& SCHEDULE(DYNAMIC) &
    !$OMP& SHARED(bpi2,bzi2,nghl,ihil_convert,iphicyl_convert,phicyl,y_opt,nl_sim,rotated_density,rotated_kappa,rotated_kappac,&
    !$OMP&       qhla_opt,fi1z_opt,fi1r_opt,rotated_ro,rotated_tau,rotated_dj,rotated_dro,rotated_aka,rotated_akac,&
    !$OMP&       rotated_sroz,rotated_sror,rotated_srofi,rotated_jz,rotated_jr,rotated_jfi,rotated_dsroz,rotated_dsror,rotated_dsrofi,&
    !$OMP&       rotated_csroz,rotated_csror,rotated_csrofi,rotated_sfiz,rotated_sfir,rotated_srfi,rotated_szfi,rotated_szz,rotated_srz,&
    !$OMP&       rotated_srr,rotated_szr,rotated_sfifi,rotated_stauz,rotated_staur,rotated_staufi,&
    !$OMP&       xl_ihil,xh_ihil,nr,nz,nt_contributing,iosc1_contributing,iosc2_contributing,maxphi,maxphi_eff,ibet,use_j2terms) &
    !$OMP& PRIVATE(ihil_iphicyl,ihil,iphicyl,phy,y,y2,iphi,ro_sumN,ro_sumP,tau_sumN,tau_sumP,dj_sumN,dj_sumP,dro_sumN,dro_sumP,&
    !$OMP&         aka_sumN,aka_sumP,akac_sumN,akac_sumP,sroz_sumN,sroz_sumP,sror_sumN,sror_sumP,srofi_sumN,srofi_sumP,&
    !$OMP&         jz_sumN,jz_sumP,jr_sumN,jr_sumP,jfi_sumN,jfi_sumP,dsroz_sumN,dsroz_sumP,dsror_sumN,dsror_sumP,dsrofi_sumN,dsrofi_sumP,&
    !$OMP&         csroz_sumN,csroz_sumP,csror_sumN,csror_sumP,csrofi_sumN,csrofi_sumP,sfiz_sumN,sfiz_sumP,sfir_sumN,sfir_sumP,srfi_sumN,srfi_sumP,&
    !$OMP&         szfi_sumN,szfi_sumP,szz_sumN,szz_sumP,srz_sumN,srz_sumP,srr_sumN,srr_sumP,szr_sumN,szr_sumP,sfifi_sumN,sfifi_sumP,&
    !$OMP&         stauz_sumN,stauz_sumP,staur_sumN,staur_sumP,staufi_sumN,staufi_sumP,iosc,iosc1,iosc2,nla1,nla2,&
    !$OMP&         q_h0l0a_h0l0b,q_h1l0a_h1l0b,q_h0l1a_h0l1b,q_h0l0a_h0l1b,q_h0l1a_h0l0b,q_h0l0a_h1l0b,q_h1l0a_h0l0b,q_h1l0a_h0l1b,&
    !$OMP&         q_h0l1a_h1l0b,fac1,fac2,fac3,fac4,fac5,fac_ro,fac_tau,fac_dj1,fac_dj2,fac_dro,fac_aka,fac_sroz,fac_sror,fac_srofi,&
    !$OMP&         fac_jz,fac_jr,fac_jfi,fac_dsroz,fac_dsror,fac_dsrofi,fac_csroz1,fac_csror1,fac_csror2,fac_csrofi1,fac_csrofi2,&
    !$OMP&         fac_sfiz,fac_sfir,fac_srfi,fac_szfi,fac_szz,fac_srz,fac_srr,fac_szr,fac_sfifi,fac_stauz,fac_staur,fac_staufi,&
    !$OMP&         densNpl,densPpl,densNmi,densPmi,kappNmi,kappPmi,kappcNmi,kappcPmi)
    Do ihil_iphicyl=1,nghl*ngphi
       ihil=ihil_convert(ihil_iphicyl); iphicyl=iphicyl_convert(ihil_iphicyl)

       phy=phicyl(iphicyl)
       y=y_opt(ihil);y2=y*y

       Do iphi=1,maxphi_eff
          ro_sumN     = c_zero; ro_sumP     = c_zero
          tau_sumN    = c_zero; tau_sumP    = c_zero
          dj_sumN     = c_zero; dj_sumP     = c_zero
          dro_sumN    = c_zero; dro_sumP    = c_zero
          aka_sumN    = c_zero; aka_sumP    = c_zero
          akac_sumN   = c_zero; akac_sumP   = c_zero

          sroz_sumN   = c_zero; sroz_sumP   = c_zero
          sror_sumN   = c_zero; sror_sumP   = c_zero
          srofi_sumN  = c_zero; srofi_sumP  = c_zero
          jz_sumN     = c_zero; jz_sumP     = c_zero
          jr_sumN     = c_zero; jr_sumP     = c_zero
          jfi_sumN    = c_zero; jfi_sumP    = c_zero
          dsroz_sumN  = c_zero; dsroz_sumP  = c_zero
          dsror_sumN  = c_zero; dsror_sumP  = c_zero
          dsrofi_sumN = c_zero; dsrofi_sumP = c_zero
          csroz_sumN  = c_zero; csroz_sumP  = c_zero
          csror_sumN  = c_zero; csror_sumP  = c_zero
          csrofi_sumN = c_zero; csrofi_sumP = c_zero

          sfiz_sumN   = c_zero; sfiz_sumP   = c_zero
          sfir_sumN   = c_zero; sfir_sumP   = c_zero
          srfi_sumN   = c_zero; srfi_sumP   = c_zero
          szfi_sumN   = c_zero; szfi_sumP   = c_zero
          szz_sumN    = c_zero; szz_sumP    = c_zero
          srz_sumN    = c_zero; srz_sumP    = c_zero
          srr_sumN    = c_zero; srr_sumP    = c_zero
          szr_sumN    = c_zero; szr_sumP    = c_zero
          sfifi_sumN  = c_zero; sfifi_sumP  = c_zero
          stauz_sumN  = c_zero; stauz_sumP  = c_zero
          staur_sumN  = c_zero; staur_sumP  = c_zero
          staufi_sumN = c_zero; staufi_sumP = c_zero
          !
          Do iosc=1,nt_contributing(iphi)
             iosc1=iosc1_contributing(iosc,iphi); iosc2=iosc2_contributing(iosc,iphi)
             nla1=nl_sim(iosc1); nla2=nl_sim(iosc2)

             densNpl  = rotated_density(1,iphi,1)%arr(iosc1,iosc2) + rotated_density(1,iphi,2)%arr(iosc1,iosc2)
             densNmi  = rotated_density(1,iphi,1)%arr(iosc1,iosc2) - rotated_density(1,iphi,2)%arr(iosc1,iosc2)
             densPpl  = rotated_density(2,iphi,1)%arr(iosc1,iosc2) + rotated_density(2,iphi,2)%arr(iosc1,iosc2)
             densPmi  = rotated_density(2,iphi,1)%arr(iosc1,iosc2) - rotated_density(2,iphi,2)%arr(iosc1,iosc2)
             kappNmi  = rotated_kappa(1,iphi,1)%arr(iosc1,iosc2)   - rotated_kappa(1,iphi,2)%arr(iosc1,iosc2)
             kappPmi  = rotated_kappa(2,iphi,1)%arr(iosc1,iosc2)   - rotated_kappa(2,iphi,2)%arr(iosc1,iosc2)
             kappcNmi = rotated_kappac(1,iphi,1)%arr(iosc1,iosc2)  - rotated_kappac(1,iphi,2)%arr(iosc1,iosc2)
             kappcPmi = rotated_kappac(2,iphi,1)%arr(iosc1,iosc2)  - rotated_kappac(2,iphi,2)%arr(iosc1,iosc2)

             q_h0l0a_h0l0b=qhla_opt(iosc1,ihil)*qhla_opt(iosc2,ihil)
             q_h1l0a_h1l0b=fi1z_opt(iosc1,ihil)*fi1z_opt(iosc2,ihil); q_h0l1a_h0l1b=fi1r_opt(iosc1,ihil)*fi1r_opt(iosc2,ihil)
             q_h0l0a_h0l1b=qhla_opt(iosc1,ihil)*fi1r_opt(iosc2,ihil); q_h0l1a_h0l0b=fi1r_opt(iosc1,ihil)*qhla_opt(iosc2,ihil)
             q_h0l0a_h1l0b=qhla_opt(iosc1,ihil)*fi1z_opt(iosc2,ihil); q_h1l0a_h0l0b=fi1z_opt(iosc1,ihil)*qhla_opt(iosc2,ihil)
             q_h1l0a_h0l1b=fi1z_opt(iosc1,ihil)*fi1r_opt(iosc2,ihil); q_h0l1a_h1l0b=fi1r_opt(iosc1,ihil)*fi1z_opt(iosc2,ihil)

             fac1=q_h1l0a_h1l0b+q_h0l1a_h0l1b+nla1*nla2*y2*q_h0l0a_h0l0b
             fac2=two*(xl_ihil(ihil)*bpi2+(xh_ihil(ihil))**2*bzi2-two*bpi2*(nr(iosc1)+nr(iosc2)+0.5_pr*(abs(nla1)+abs(nla2))+1.0_pr)&
                  -bzi2*(nz(iosc1)+nz(iosc2)+1.0_pr)+nla1*nla2*y2)*q_h0l0a_h0l0b+two*(q_h0l1a_h0l1b+q_h1l0a_h1l0b)
             fac3=two*(xl_ihil(ihil)*bpi2+(xh_ihil(ihil))**2*bzi2-two*bpi2*(nr(iosc1)+nr(iosc2)+0.5_pr*(abs(nla1)+abs(nla2))+1.0_pr)&
                  -bzi2*(nz(iosc1)+nz(iosc2)+1.0_pr)-nla1*nla2*y2)*q_h0l0a_h0l0b+two*(q_h0l1a_h0l1b+q_h1l0a_h1l0b)
             fac4=y*(nla1+nla2)*q_h0l0a_h0l0b
             fac5=y*(nla1-nla2)*q_h0l0a_h0l0b

             ! --------------------
             ! Time-even components
             ! --------------------
             ! ro density
             fac_ro=q_h0l0a_h0l0b*Cos((nla2-nla1)*phy)
             ro_sumN=ro_sumN+densNpl*fac_ro; ro_sumP=ro_sumP+densPpl*fac_ro
             ! tau density
             fac_tau=fac1*Cos((nla2-nla1)*phy)
             tau_sumN=tau_sumN+densNpl*fac_tau; tau_sumP=tau_sumP+densPpl*fac_tau
             ! dj density
             fac_dj1=y*(nla1*q_h0l0a_h0l1b+nla2*q_h0l1a_h0l0b)*Cos((nla2-nla1)*phy)
             fac_dj2=(y*(nla2*q_h1l0a_h0l0b-nla1*q_h0l0a_h1l0b)-q_h1l0a_h0l1b+q_h0l1a_h1l0b)*Cos((nla1+nla2+1)*phy)*unit_i
             dj_sumN=dj_sumN+densNpl*fac_dj1+densNmi*fac_dj2; dj_sumP=dj_sumP+densPpl*fac_dj1+densPmi*fac_dj2
             ! laplacian of ro
             fac_dro=fac2*half*Cos((nla2-nla1)*phy)
             dro_sumN=dro_sumN+densNpl*fac_dro; dro_sumP=dro_sumP+densPpl*fac_dro
             ! aka density
             fac_aka=half*q_h0l0a_h0l0b*Cos((nla2-nla1)*phy)
             aka_sumN=aka_sumN+kappNmi*fac_aka; aka_sumP=aka_sumP+kappPmi*fac_aka
             akac_sumN=akac_sumN+kappcNmi*fac_aka; akac_sumP=akac_sumP+kappcPmi*fac_aka
             ! -------------------
             ! Time-odd components
             ! -------------------
             ! sro (spin) density
             ! z component
             fac_sroz=unit_i*q_h0l0a_h0l0b*Sin((nla1-nla2)*phy)
             sroz_sumN=sroz_sumN+densNpl*fac_sroz; sroz_sumP=sroz_sumP+densPpl*fac_sroz
             ! rperp component
             fac_sror=-q_h0l0a_h0l0b*Sin((nla1+nla2+1)*phy)
             sror_sumN=sror_sumN+densNmi*fac_sror; sror_sumP=sror_sumP+densPmi*fac_sror
             ! fi component
             fac_srofi=-q_h0l0a_h0l0b*Cos((nla1+nla2+1)*phy)
             srofi_sumN=srofi_sumN+densNmi*fac_srofi; srofi_sumP=srofi_sumP+densPmi*fac_srofi
             ! current
             ! z component
             fac_jz=-half*unit_i*(q_h1l0a_h0l0b-q_h0l0a_h1l0b)*Cos((nla1-nla2)*phy)
             jz_sumN=jz_sumN+densNpl*fac_jz; jz_sumP=jz_sumP+densPpl*fac_jz
             ! rperp component
             fac_jr=-half*unit_i*(q_h0l1a_h0l0b-q_h0l0a_h0l1b)*Cos((nla1-nla2)*phy)
             jr_sumN=jr_sumN+densNpl*fac_jr; jr_sumP=jr_sumP+densPpl*fac_jr
             ! fi component
             fac_jfi=-half*unit_i*fac4*Sin((nla2-nla1)*phy)
             jfi_sumN=jfi_sumN+densNpl*fac_jfi; jfi_sumP=jfi_sumP+densPpl*fac_jfi
             ! laplacian of sro
             ! z component
             fac_dsroz=fac2*Sin((nla1-nla2)*phy)*unit_i
             dsroz_sumN=dsroz_sumN+densNpl*fac_dsroz; dsroz_sumP=dsroz_sumP+densPpl*fac_dsroz
             ! rperp component
             fac_dsror=fac3*Sin((nla1+nla2+1)*phy)*(-1.0_pr)
             dsror_sumN=dsror_sumN+densNmi*fac_dsror; dsror_sumP=dsror_sumP+densPmi*fac_dsror
             ! fi component
             fac_dsrofi=fac3*Cos((nla1+nla2+1)*phy)*(-1.0_pr)
             dsrofi_sumN=dsrofi_sumN+densNmi*fac_dsrofi; dsrofi_sumP=dsrofi_sumP+densPmi*fac_dsrofi
             ! curl sro
             ! z component
             fac_csroz1=(fac4-q_h0l1a_h0l0b-q_h0l0a_h0l1b)*Cos((nla1+nla2+1.0_pr)*phy)
             csroz_sumN=csroz_sumN+densNmi*fac_csroz1; csroz_sumP=csroz_sumP+densPmi*fac_csroz1
             ! rperp component
             fac_csror1=fac5*Cos((nla1-nla2)*phy)*unit_i
             fac_csror2=(q_h1l0a_h0l0b+q_h0l0a_h1l0b)*Cos((nla1+nla2+1.0_pr)*phy)
             csror_sumN=csror_sumN+densNpl*fac_csror1+densNmi*fac_csror2
             csror_sumP=csror_sumP+densPpl*fac_csror1+densPmi*fac_csror2
             ! fi component
             fac_csrofi1=(-q_h0l1a_h0l0b-q_h0l0a_h0l1b)*Sin((nla1-nla2)*phy)*unit_i
             fac_csrofi2=(-q_h1l0a_h0l0b-q_h0l0a_h1l0b)*Sin((nla1+nla2+1.0_pr)*phy)
             csrofi_sumN=csrofi_sumN+densNpl*fac_csrofi1+densNmi*fac_csrofi2
             csrofi_sumP=csrofi_sumP+densPpl*fac_csrofi1+densPmi*fac_csrofi2
             ! -----------------
             ! Tensor components
             ! -----------------
             If(use_j2terms .Eqv. .True.) Then
                ! tensor density components
                ! sfiz
                   fac_sfiz=half*fac4*Cos((nla2-nla1)*phy)
                sfiz_sumN=sfiz_sumN+densNpl*fac_sfiz; sfiz_sumP=sfiz_sumP+densPpl*fac_sfiz
                ! sfir
                   fac_sfir=half*unit_i*fac5*Cos((nla2+nla1+1)*phy)
                sfir_sumN=sfir_sumN+densNmi*fac_sfir; sfir_sumP=sfir_sumP+densPmi*fac_sfir
                ! srfi
                fac_srfi=half*unit_i*(q_h0l1a_h0l0b-q_h0l0a_h0l1b)*Cos((nla1+nla2+1)*phy)
                srfi_sumN=srfi_sumN+densNmi*fac_srfi; srfi_sumP=srfi_sumP+densPmi*fac_srfi
                ! szfi
                fac_szfi=half*unit_i*(q_h1l0a_h0l0b-q_h0l0a_h1l0b)*Cos((nla1+nla2+1)*phy)
                szfi_sumN=szfi_sumN+densNmi*fac_szfi; szfi_sumP=szfi_sumP+densPmi*fac_szfi
                ! szz
                fac_szz=half*(q_h0l0a_h1l0b-q_h1l0a_h0l0b)*Sin((nla2-nla1)*phy)
                szz_sumN=szz_sumN+densNpl*fac_szz; szz_sumP=szz_sumP+densPpl*fac_szz
                ! srz
                fac_srz=half*(q_h0l0a_h0l1b-q_h0l1a_h0l0b)*Sin((nla2-nla1)*phy)
                srz_sumN=srz_sumN+densNpl*fac_srz; srz_sumP=srz_sumP+densPpl*fac_srz
                ! srr
                fac_srr=half*unit_i*(q_h0l1a_h0l0b-q_h0l0a_h0l1b)*Sin((nla1+nla2+1)*phy)
                srr_sumN=srr_sumN+densNmi*fac_srr; srr_sumP=srr_sumP+densPmi*fac_srr
                ! szr
                fac_szr=half*unit_i*(q_h1l0a_h0l0b-q_h0l0a_h1l0b)*Sin((nla1+nla2+1)*phy)
                szr_sumN=szr_sumN+densNmi*fac_szr; szr_sumP=szr_sumP+densPmi*fac_szr
                ! sfifi
                fac_sfifi=-half*unit_i*fac5*Sin((nla1+nla2+1)*phy)
                sfifi_sumN=sfifi_sumN+densNmi*fac_sfifi; sfifi_sumP=sfifi_sumP+densPmi*fac_sfifi
                ! stau (spin) density components
                ! z
                fac_stauz=fac1*Sin((nla1-nla2)*phy)*unit_i
                stauz_sumN=stauz_sumN+densNpl*fac_stauz; stauz_sumP=stauz_sumP+densPpl*fac_stauz
                ! rperp
                fac_staur=(q_h0l0a_h0l0b*nla1*nla2*y2-q_h0l1a_h0l1b-q_h1l0a_h1l0b)*Sin((nla1+nla2+1.0_pr)*phy)
                staur_sumN=staur_sumN+densNmi*fac_staur; staur_sumP=staur_sumP+densPmi*fac_staur
                ! fi
                fac_staufi=(q_h0l0a_h0l0b*nla1*nla2*y2-q_h0l1a_h0l1b-q_h1l0a_h1l0b)*Cos((nla1+nla2+1.0_pr)*phy)
                staufi_sumN=staufi_sumN+densNmi*fac_staufi; staufi_sumP=staufi_sumP+densPmi*fac_staufi
             End If

          End Do ! iosc

          rotated_ro(iphi,ihil_iphicyl,1)   = ro_sumN;    rotated_ro(iphi,ihil_iphicyl,2)    = ro_sumP
          rotated_tau(iphi,ihil_iphicyl,1)  = tau_sumN;   rotated_tau(iphi,ihil_iphicyl,2)   = tau_sumP
          rotated_dj(iphi,ihil_iphicyl,1)   = dj_sumN;    rotated_dj(iphi,ihil_iphicyl,2)    = dj_sumP
          rotated_dro(iphi,ihil_iphicyl,1)  = dro_sumN;   rotated_dro(iphi,ihil_iphicyl,2)   = dro_sumP
          rotated_aka(iphi,ihil_iphicyl,1)  = aka_sumN;   rotated_aka(iphi,ihil_iphicyl,2)   = aka_sumP
          rotated_akac(iphi,ihil_iphicyl,1) = akac_sumN;  rotated_akac(iphi,ihil_iphicyl,2)  = akac_sumP

          rotated_sroz(iphi,ihil_iphicyl,1)  = sroz_sumN;   rotated_sroz(iphi,ihil_iphicyl,2)  = sroz_sumP
          rotated_sror(iphi,ihil_iphicyl,1)  = sror_sumN;   rotated_sror(iphi,ihil_iphicyl,2)  = sror_sumP
          rotated_srofi(iphi,ihil_iphicyl,1) = srofi_sumN;  rotated_srofi(iphi,ihil_iphicyl,2) = srofi_sumP
          rotated_jz(iphi,ihil_iphicyl,1)    = jz_sumN;     rotated_jz(iphi,ihil_iphicyl,2)    = jz_sumP
          rotated_jr(iphi,ihil_iphicyl,1)    = jr_sumN;     rotated_jr(iphi,ihil_iphicyl,2)    = jr_sumP
          rotated_jfi(iphi,ihil_iphicyl,1)   = jfi_sumN;    rotated_jfi(iphi,ihil_iphicyl,2)   = jfi_sumP
          rotated_dsroz(iphi,ihil_iphicyl,1) = dsroz_sumN;  rotated_dsroz(iphi,ihil_iphicyl,2)  = dsroz_sumP
          rotated_dsror(iphi,ihil_iphicyl,1) = dsror_sumN;  rotated_dsror(iphi,ihil_iphicyl,2)  = dsror_sumP
          rotated_dsrofi(iphi,ihil_iphicyl,1)= dsrofi_sumN; rotated_dsrofi(iphi,ihil_iphicyl,2) = dsrofi_sumP
          rotated_csroz(iphi,ihil_iphicyl,1) = csroz_sumN;  rotated_csroz(iphi,ihil_iphicyl,2)  = csroz_sumP
          rotated_csror(iphi,ihil_iphicyl,1) = csror_sumN;  rotated_csror(iphi,ihil_iphicyl,2)  = csror_sumP
          rotated_csrofi(iphi,ihil_iphicyl,1)= csrofi_sumN; rotated_csrofi(iphi,ihil_iphicyl,2) = csrofi_sumP

          If(use_j2terms .Eqv. .True.) Then
             rotated_sfiz(iphi,ihil_iphicyl,1) = sfiz_sumN;  rotated_sfiz(iphi,ihil_iphicyl,2)  = sfiz_sumP
             rotated_sfir(iphi,ihil_iphicyl,1) = sfir_sumN;  rotated_sfir(iphi,ihil_iphicyl,2)  = sfir_sumP
             rotated_srfi(iphi,ihil_iphicyl,1) = srfi_sumN;  rotated_srfi(iphi,ihil_iphicyl,2)  = srfi_sumP
             rotated_szfi(iphi,ihil_iphicyl,1) = szfi_sumN;  rotated_szfi(iphi,ihil_iphicyl,2)  = szfi_sumP
             rotated_szz(iphi,ihil_iphicyl,1)  = szz_sumN;   rotated_szz(iphi,ihil_iphicyl,2)   = szz_sumP
             rotated_srz(iphi,ihil_iphicyl,1)  = srz_sumN;   rotated_srz(iphi,ihil_iphicyl,2)   = srz_sumP
             rotated_srr(iphi,ihil_iphicyl,1)  = srr_sumN;   rotated_srr(iphi,ihil_iphicyl,2)   = srr_sumP
             rotated_szr(iphi,ihil_iphicyl,1)  = szr_sumN;   rotated_szr(iphi,ihil_iphicyl,2)   = szr_sumP
             rotated_sfifi(iphi,ihil_iphicyl,1)= sfifi_sumN; rotated_sfifi(iphi,ihil_iphicyl,2) = sfifi_sumP
             rotated_stauz(iphi,ihil_iphicyl,1) = stauz_sumN;  rotated_stauz(iphi,ihil_iphicyl,2)  = stauz_sumP
             rotated_staur(iphi,ihil_iphicyl,1) = staur_sumN;  rotated_staur(iphi,ihil_iphicyl,2)  = staur_sumP
             rotated_staufi(iphi,ihil_iphicyl,1)= staufi_sumN; rotated_staufi(iphi,ihil_iphicyl,2) = staufi_sumP
          End If

       End Do ! iphi

    End Do ! ihil_iphicyl
    !$OMP End PARALLEL DO

    ! Calculating projected density
    If(PNP_is_on == 2) Then
       !$OMP  PARALLEL DO       &
       !$OMP& DEFAULT(NONE)     &
       !$OMP& SCHEDULE(DYNAMIC) &
       !$OMP& SHARED(rotated_ro,nghl,maxphi,maxphi_eff,ibet,ephicN,ephicP,ro_projected,all_overlaps,beta_size,team_rank,offset) &
       !$OMP& PRIVATE(ihil_iphicyl,iphi,ro_projN,ro_projP,ro_normN,ro_normP,overlapN,overlapP,rotated_roN,rotated_roP)
       Do ihil_iphicyl=1,nghl*ngphi
          ro_projN=c_zero; ro_projP=c_zero; ro_normN=c_zero; ro_normP=c_zero
          Do iphi=1,maxphi
             If(iphi  <=  maxphi_eff) Then
                rotated_roN=rotated_ro(iphi,ihil_iphicyl,1)
                rotated_roP=rotated_ro(iphi,ihil_iphicyl,2)
             Else
                rotated_roN=conjg(rotated_ro(maxphi-iphi+2,ihil_iphicyl,1))
                rotated_roP=conjg(rotated_ro(maxphi-iphi+2,ihil_iphicyl,2))
             End If
             overlapN = all_overlaps((iphi-1)*beta_size+ibet-offset(team_rank))
             overlapP = all_overlaps(maxphi*beta_size+(iphi-1)*beta_size+ibet-offset(team_rank))
             ro_projN = ro_projN + ephicN(iphi,0)*overlapN*rotated_roN
             ro_projP = ro_projP + ephicP(iphi,0)*overlapP*rotated_roP
             ro_normN = ro_normN + ephicN(iphi,0)*overlapN
             ro_normP = ro_normP + ephicP(iphi,0)*overlapP
          End Do
          ro_projected(ihil_iphicyl,1) = ro_projN/ro_normN; ro_projected(ihil_iphicyl,2) = ro_projP/ro_normP
       End Do
       !$OMP End PARALLEL DO
    End If

    ! Normalization of densities
    Do it=1,2
      Do ihil_iphicyl=1,nghl*ngphi
         Do iphi=1,maxphi_eff
            rotated_ro(iphi,ihil_iphicyl,it)   = rotated_ro(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_tau(iphi,ihil_iphicyl,it ) = rotated_tau(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_dj(iphi,ihil_iphicyl,it)   = rotated_dj(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_dro(iphi,ihil_iphicyl,it)  = rotated_dro(iphi,ihil_iphicyl,it)*two*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_aka(iphi,ihil_iphicyl,it)  = rotated_aka(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_akac(iphi,ihil_iphicyl,it) = rotated_akac(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))

            rotated_sroz(iphi,ihil_iphicyl,it)  = rotated_sroz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_sror(iphi,ihil_iphicyl,it)  = rotated_sror(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_srofi(iphi,ihil_iphicyl,it) = rotated_srofi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_jz(iphi,ihil_iphicyl,it)    = rotated_jz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_jr(iphi,ihil_iphicyl,it)    = rotated_jr(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_jfi(iphi,ihil_iphicyl,it)   = rotated_jfi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_dsroz(iphi,ihil_iphicyl,it) = rotated_dsroz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_dsror(iphi,ihil_iphicyl,it) = rotated_dsror(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_dsrofi(iphi,ihil_iphicyl,it)= rotated_dsrofi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_csroz(iphi,ihil_iphicyl,it) = rotated_csroz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_csror(iphi,ihil_iphicyl,it) = rotated_csror(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            rotated_csrofi(iphi,ihil_iphicyl,it)= rotated_csrofi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))

            If(use_j2terms .Eqv. .True.) Then
               rotated_sfiz(iphi,ihil_iphicyl,it) = rotated_sfiz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_sfir(iphi,ihil_iphicyl,it) = rotated_sfir(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_srfi(iphi,ihil_iphicyl,it) = rotated_srfi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_szfi(iphi,ihil_iphicyl,it) = rotated_szfi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_szz(iphi,ihil_iphicyl,it)  = rotated_szz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_srz(iphi,ihil_iphicyl,it)  = rotated_srz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_srr(iphi,ihil_iphicyl,it)  = rotated_srr(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_szr(iphi,ihil_iphicyl,it)  = rotated_szr(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_sfifi(iphi,ihil_iphicyl,it)= rotated_sfifi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_stauz(iphi,ihil_iphicyl,it) = rotated_stauz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_staur(iphi,ihil_iphicyl,it) = rotated_staur(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
               rotated_staufi(iphi,ihil_iphicyl,it)= rotated_staufi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
            End If
         End Do
         If(PNP_is_on == 2) ro_projected(ihil_iphicyl,it)=ro_projected(ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
      End Do
    End Do

    ! Calculation of the Coulomb field
    Allocate(vc_pr(nghl,nghl)); vc_pr = c_zero
    !$OMP  PARALLEL DO       &
    !$OMP& DEFAULT(NONE)     &
    !$OMP& SCHEDULE(DYNAMIC) &
    !$OMP& SHARED(nghl,vc_pr,vc) &
    !$OMP& PRIVATE(ihil2,ihil1)
    Do ihil2=1,nghl
       Do ihil1=1,nghl
          vc_pr(ihil1,ihil2)=vc(ihil1,ihil2)*c_unit
       End Do
    End Do
    !$OMP End PARALLEL DO

    cou_rotated=c_zero; facmul = c_unit/Real(ngphi,Kind=pr)
    Allocate(cou_rotated_tmp(nghl),rotated_ro_tmp(nghl))
    !$OMP  PARALLEL DO       &
    !$OMP& DEFAULT(NONE)     &
    !$OMP& SCHEDULE(DYNAMIC) &
    !$OMP& SHARED(maxphi_eff,nghl,facmul,vc_pr,rotated_ro,cou_rotated) &
    !$OMP& PRIVATE(iphi,cou_rotated_tmp,istart,iend,iphicyl,rotated_ro_tmp)
    Do iphi=1,maxphi_eff
       cou_rotated_tmp=c_zero
       Do iphicyl=1,ngphi
          istart = (iphicyl-1)*nghl+1; iend = iphicyl*nghl
          rotated_ro_tmp = rotated_ro(iphi,istart:iend,2)
          Call zgemv('N',nghl,nghl,facmul,vc_pr,nghl,rotated_ro_tmp,1,c_unit,cou_rotated_tmp,1)
       End Do
       cou_rotated(iphi,1:nghl)=cou_rotated_tmp
    End Do
    !$OMP End PARALLEL DO
    Deallocate(vc_pr)
    Deallocate(cou_rotated_tmp,rotated_ro_tmp)

    ! Timing
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(4)=time_ibet(4)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)

  End Subroutine calculate_densities
  !=============================================
  !> This routine calculates rotated (gauge-dependent) energies.
  !=============================================
  Subroutine calculate_energies(ibet)
    Use HFBTHO_Gauss, Only : nghl,wdcor,fh,fl
    Use HFBTHO, Only : team_rank,PNP_is_on,facECM,hb0n,hb0p,icou,cex,rho_c
    Use UNEDF, Only : sigma,CpV0,CpV1,CExPar,use_j2terms,Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,&
                      Cs,Cds,Csj,CsDs,CsdJ,cT
    Use math, Only : integrate_complex_function

    Integer(ipr), Intent(In) ::  ibet !< - Index of the current \f$ \beta \f$ angle

    Integer(ipr) :: ihil,iphicyl,ihil_iphicyl,iphin,iphip,iphin_iphip,lloc
    Integer(ipr) :: tstart_ibet,tfinish_ibet
    Real(pr) :: z,zz,rrr,p2,twopii,whl
    Complex(pr) :: ekinN,ekinP,ecodi,ecoex,EVOL_rho_tau,EVOL_rho_rho,ESURF_rho_drho,ESO_rho_nablaj,eptN,eptP,&
                   ess0,essD,ecurr,esDs,ejNxs,ETENS,esT,xn1,rms1,xn2,rms2,delN,delP
    Complex(pr) :: rn,rp,tnt,tpt,dn,dp,djn,djp,akn,akp,aknc,akpc,akn2,akp2,rprojn,rprojp,adn,adp,srozn,srorn,srofin,srozp,srorp,srofip,&
                   jzn,jzp,jrn,jrp,jfin,jfip,dsrozn,dsrorn,dsrofin,dsrozp,dsrorp,dsrofip,csrozn,csrozp,csrorn,csrorp,csrofin,csrofip,&
                   sfizn,sfizp,sfirn,sfirp,szfin,szfip,srfin,srfip,szzn,szzp,srzn,srzp,srrn,srrp,szrn,szrp,sfifin,sfifip,&
                   stauzn,stauzp,staurn,staurp,staufin,staufip
    Complex(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,RHOPROJ_0,RHOPROJ_1,S2_0,S2_1,J2_0,J2_1,SDS_0,SDS_1,JDS_0,JDS_1,&
                   SFIZ_0,SFIZ_1,SFIR_0,SFIR_1,SZFI_0,SZFI_1,SRFI_0,SRFI_1,SZZ_0,SZZ_1,SRZ_0,SRZ_1,SRR_0,SRR_1,SZR_0,SZR_1,SFIFI_0,SFIFI_1,&
                   TENS_0,TENS_1,ST_0,ST_1
    Complex(pr) :: ekinN_phicyl(ngphi),ekinP_phicyl(ngphi),ecodi_phicyl(ngphi),ecoex_phicyl(ngphi),EVOL_rho_tau_phicyl(ngphi),&
                   EVOL_rho_rho_phicyl(ngphi),ESURF_rho_drho_phicyl(ngphi),ESO_rho_nablaj_phicyl(ngphi),eptN_phicyl(ngphi),eptP_phicyl(ngphi),&
                   ess0_phicyl(ngphi),essD_phicyl(ngphi),ecurr_phicyl(ngphi),esDs_phicyl(ngphi),ejNxs_phicyl(ngphi),ETENS_phicyl(ngphi),esT_phicyl(ngphi),&
                   row1_phicyl(ngphi),row2_phicyl(ngphi),delN_phicyl(ngphi),delP_phicyl(ngphi)
    Complex(pr) :: ekinN_integral,ekinP_integral,ecodi_integral,ecoex_integral,EVOL_rho_tau_integral,EVOL_rho_rho_integral,ESURF_rho_drho_integral,&
                   ESO_rho_nablaj_integral,eptN_integral,eptP_integral,ess0_integral,essD_integral,ecurr_integral,esDs_integral,ejNxs_integral,&
                   ETENS_integral,esT_integral,row1_integral,row2_integral,delN_integral,delP_integral
    Complex(pr) :: rho_coupling,rho_exchange,ccou

    If(beta_active(ibet-offset(team_rank)) == 0) Return

    ! Timing
    Call system_clock(tstart_ibet,clock_rate)

    twopii=one/(two*pi)
    iphin_iphip=0
    Do iphin=1,maxphi
       Do iphip=1,maxphi
          iphin_iphip=iphin_iphip+1
          ekinN=c_zero; ekinP=c_zero; ecodi=c_zero; ecoex=c_zero; EVOL_rho_tau=c_zero; EVOL_rho_rho=c_zero
          ESURF_rho_drho=c_zero; ESO_rho_nablaj=c_zero; eptN=c_zero; eptP=c_zero
          ess0=c_zero; essD=c_zero; ecurr=c_zero; esDs=c_zero; ejNxs=c_zero; ETENS=c_zero; esT=c_zero
          xn1=c_zero; rms1=c_zero; xn2=c_zero; rms2=c_zero; delN=c_zero; delP=c_zero
          Do ihil=1,nghl
             whl=wdcor(ihil)
             z=fh(ihil); zz=z*z; rrr=zz+fl(ihil)**2; p2=p32*zz-half*rrr

             Do iphicyl=1,ngphi

                ihil_iphicyl=ihil_iphicyl_convert(ihil,iphicyl)

                ! Projected density
                rprojn=ro_projected(ihil_iphicyl,1); rprojp=ro_projected(ihil_iphicyl,2)
                RHOPROJ_0=rprojn+rprojp; RHOPROJ_1=rprojn-rprojp

                ! rms and deformations
                z=fh(ihil); zz=z*z; rrr=zz+fl(ihil)**2
                p2=p32*zz-half*rrr

                ! np-representation
                If(iphin  <=  maxphi_eff) Then
                   rn      = rotated_ro(iphin,ihil_iphicyl,1)
                   tnt     = rotated_tau(iphin,ihil_iphicyl,1)
                   dn      = rotated_dro(iphin,ihil_iphicyl,1)
                   djn     = rotated_dj(iphin,ihil_iphicyl,1)
                   akn     = rotated_aka(iphin,ihil_iphicyl,1)
                   aknc    = rotated_akac(iphin,ihil_iphicyl,1)
                   akn2=akn*aknc; adn=akn*rn

                   srozn   = rotated_sroz(iphin,ihil_iphicyl,1)
                   srorn   = rotated_sror(iphin,ihil_iphicyl,1)
                   srofin  = rotated_srofi(iphin,ihil_iphicyl,1)
                   jzn     = rotated_jz(iphin,ihil_iphicyl,1)
                   jrn     = rotated_jr(iphin,ihil_iphicyl,1)
                   jfin    = rotated_jfi(iphin,ihil_iphicyl,1)
                   dsrozn  = rotated_dsroz(iphin,ihil_iphicyl,1)
                   dsrorn  = rotated_dsror(iphin,ihil_iphicyl,1)
                   dsrofin = rotated_dsrofi(iphin,ihil_iphicyl,1)
                   csrozn  = rotated_csroz(iphin,ihil_iphicyl,1)
                   csrorn  = rotated_csror(iphin,ihil_iphicyl,1)
                   csrofin = rotated_csrofi(iphin,ihil_iphicyl,1)

                   sfizn   = rotated_sfiz(iphin,ihil_iphicyl,1)
                   sfirn   = rotated_sfir(iphin,ihil_iphicyl,1)
                   szfin   = rotated_szfi(iphin,ihil_iphicyl,1)
                   srfin   = rotated_srfi(iphin,ihil_iphicyl,1)
                   szzn    = rotated_szz(iphin,ihil_iphicyl,1)
                   srzn    = rotated_srz(iphin,ihil_iphicyl,1)
                   srrn    = rotated_srr(iphin,ihil_iphicyl,1)
                   szrn    = rotated_szr(iphin,ihil_iphicyl,1)
                   sfifin  = rotated_sfifi(iphin,ihil_iphicyl,1)
                   stauzn  = rotated_stauz(iphin,ihil_iphicyl,1)
                   staurn  = rotated_staur(iphin,ihil_iphicyl,1)
                   staufin = rotated_staufi(iphin,ihil_iphicyl,1)
                Else
                   rn      = conjg(rotated_ro(maxphi-iphin+2,ihil_iphicyl,1))
                   tnt     = conjg(rotated_tau(maxphi-iphin+2,ihil_iphicyl,1))
                   dn      = conjg(rotated_dro(maxphi-iphin+2,ihil_iphicyl,1))
                   djn     = conjg(rotated_dj(maxphi-iphin+2,ihil_iphicyl,1))
                   akn     = conjg(rotated_aka(maxphi-iphin+2,ihil_iphicyl,1))
                   aknc    = conjg(rotated_akac(maxphi-iphin+2,ihil_iphicyl,1))
                   akn2=akn*aknc; adn=akn*rn

                   srozn   = conjg(rotated_sroz(maxphi-iphin+2,ihil_iphicyl,1))
                   srorn   = conjg(rotated_sror(maxphi-iphin+2,ihil_iphicyl,1))
                   srofin  = conjg(rotated_srofi(maxphi-iphin+2,ihil_iphicyl,1))
                   jzn     = conjg(rotated_jz(maxphi-iphin+2,ihil_iphicyl,1))
                   jrn     = conjg(rotated_jr(maxphi-iphin+2,ihil_iphicyl,1))
                   jfin    = conjg(rotated_jfi(maxphi-iphin+2,ihil_iphicyl,1))
                   dsrozn  = conjg(rotated_dsroz(maxphi-iphin+2,ihil_iphicyl,1))
                   dsrorn  = conjg(rotated_dsror(maxphi-iphin+2,ihil_iphicyl,1))
                   dsrofin = conjg(rotated_dsrofi(maxphi-iphin+2,ihil_iphicyl,1))
                   csrozn  = conjg(rotated_csroz(maxphi-iphin+2,ihil_iphicyl,1))
                   csrorn  = conjg(rotated_csror(maxphi-iphin+2,ihil_iphicyl,1))
                   csrofin = conjg(rotated_csrofi(maxphi-iphin+2,ihil_iphicyl,1))

                   sfizn   = conjg(rotated_sfiz(maxphi-iphin+2,ihil_iphicyl,1))
                   sfirn   = conjg(rotated_sfir(maxphi-iphin+2,ihil_iphicyl,1))
                   szfin   = conjg(rotated_szfi(maxphi-iphin+2,ihil_iphicyl,1))
                   srfin   = conjg(rotated_srfi(maxphi-iphin+2,ihil_iphicyl,1))
                   szzn    = conjg(rotated_szz(maxphi-iphin+2,ihil_iphicyl,1))
                   srzn    = conjg(rotated_srz(maxphi-iphin+2,ihil_iphicyl,1))
                   srrn    = conjg(rotated_srr(maxphi-iphin+2,ihil_iphicyl,1))
                   szrn    = conjg(rotated_szr(maxphi-iphin+2,ihil_iphicyl,1))
                   sfifin  = conjg(rotated_sfifi(maxphi-iphin+2,ihil_iphicyl,1))
                   stauzn  = conjg(rotated_stauz(maxphi-iphin+2,ihil_iphicyl,1))
                   staurn  = conjg(rotated_staur(maxphi-iphin+2,ihil_iphicyl,1))
                   staufin = conjg(rotated_staufi(maxphi-iphin+2,ihil_iphicyl,1))
                End If
                If(iphip  <=  maxphi_eff) Then
                   rp      = rotated_ro(iphip,ihil_iphicyl,2)
                   tpt     = rotated_tau(iphip,ihil_iphicyl,2)
                   dp      = rotated_dro(iphip,ihil_iphicyl,2)
                   djp     = rotated_dj(iphip,ihil_iphicyl,2)
                   akp     = rotated_aka(iphip,ihil_iphicyl,2)
                   akpc    = rotated_akac(iphip,ihil_iphicyl,2)
                   akp2=akp*akpc; adp=akp*rp
                   ccou    = cou_rotated(iphip,ihil)

                   srozp   = rotated_sroz(iphip,ihil_iphicyl,2)
                   srorp   = rotated_sror(iphip,ihil_iphicyl,2)
                   srofip  = rotated_srofi(iphip,ihil_iphicyl,2)
                   jzp     = rotated_jz(iphip,ihil_iphicyl,2)
                   jrp     = rotated_jr(iphip,ihil_iphicyl,2)
                   jfip    = rotated_jfi(iphip,ihil_iphicyl,2)
                   dsrozp  = rotated_dsroz(iphip,ihil_iphicyl,2)
                   dsrorp  = rotated_dsror(iphip,ihil_iphicyl,2)
                   dsrofip = rotated_dsrofi(iphip,ihil_iphicyl,2)
                   csrozp  = rotated_csroz(iphip,ihil_iphicyl,2)
                   csrorp  = rotated_csror(iphip,ihil_iphicyl,2)
                   csrofip = rotated_csrofi(iphip,ihil_iphicyl,2)

                   sfizp   = rotated_sfiz(iphip,ihil_iphicyl,2)
                   sfirp   = rotated_sfir(iphip,ihil_iphicyl,2)
                   szfip   = rotated_szfi(iphip,ihil_iphicyl,2)
                   srfip   = rotated_srfi(iphip,ihil_iphicyl,2)
                   szzp    = rotated_szz(iphip,ihil_iphicyl,2)
                   srzp    = rotated_srz(iphip,ihil_iphicyl,2)
                   srrp    = rotated_srr(iphip,ihil_iphicyl,2)
                   szrp    = rotated_szr(iphip,ihil_iphicyl,2)
                   sfifip  = rotated_sfifi(iphip,ihil_iphicyl,2)
                   stauzp  = rotated_stauz(iphip,ihil_iphicyl,2)
                   staurp  = rotated_staur(iphip,ihil_iphicyl,2)
                   staufip = rotated_staufi(iphip,ihil_iphicyl,2)
                Else
                   rp      = conjg(rotated_ro(maxphi-iphip+2,ihil_iphicyl,2))
                   tpt     = conjg(rotated_tau(maxphi-iphip+2,ihil_iphicyl,2))
                   dp      = conjg(rotated_dro(maxphi-iphip+2,ihil_iphicyl,2))
                   djp     = conjg(rotated_dj(maxphi-iphip+2,ihil_iphicyl,2))
                   akp     = conjg(rotated_aka(maxphi-iphip+2,ihil_iphicyl,2))
                   akpc    = conjg(rotated_akac(maxphi-iphip+2,ihil_iphicyl,2))
                   akp2=akp*akpc; adp=akp*rp
                   ccou    = conjg(cou_rotated(maxphi-iphip+2,ihil))

                   srozp   = conjg(rotated_sroz(maxphi-iphip+2,ihil_iphicyl,2))
                   srorp   = conjg(rotated_sror(maxphi-iphip+2,ihil_iphicyl,2))
                   srofip  = conjg(rotated_srofi(maxphi-iphip+2,ihil_iphicyl,2))
                   jzp     = conjg(rotated_jz(maxphi-iphip+2,ihil_iphicyl,2))
                   jrp     = conjg(rotated_jr(maxphi-iphip+2,ihil_iphicyl,2))
                   jfip    = conjg(rotated_jfi(maxphi-iphip+2,ihil_iphicyl,2))
                   dsrozp  = conjg(rotated_dsroz(maxphi-iphip+2,ihil_iphicyl,2))
                   dsrorp  = conjg(rotated_dsror(maxphi-iphip+2,ihil_iphicyl,2))
                   dsrofip = conjg(rotated_dsrofi(maxphi-iphip+2,ihil_iphicyl,2))
                   csrozp  = conjg(rotated_csroz(maxphi-iphip+2,ihil_iphicyl,2))
                   csrorp  = conjg(rotated_csror(maxphi-iphip+2,ihil_iphicyl,2))
                   csrofip = conjg(rotated_csrofi(maxphi-iphip+2,ihil_iphicyl,2))

                   sfizp   = conjg(rotated_sfiz(maxphi-iphip+2,ihil_iphicyl,2))
                   sfirp   = conjg(rotated_sfir(maxphi-iphip+2,ihil_iphicyl,2))
                   szfip   = conjg(rotated_szfi(maxphi-iphip+2,ihil_iphicyl,2))
                   srfip   = conjg(rotated_srfi(maxphi-iphip+2,ihil_iphicyl,2))
                   szzp    = conjg(rotated_szz(maxphi-iphip+2,ihil_iphicyl,2))
                   srzp    = conjg(rotated_srz(maxphi-iphip+2,ihil_iphicyl,2))
                   srrp    = conjg(rotated_srr(maxphi-iphip+2,ihil_iphicyl,2))
                   szrp    = conjg(rotated_szr(maxphi-iphip+2,ihil_iphicyl,2))
                   sfifip  = conjg(rotated_sfifi(maxphi-iphip+2,ihil_iphicyl,2))
                   stauzp  = conjg(rotated_stauz(maxphi-iphip+2,ihil_iphicyl,2))
                   staurp  = conjg(rotated_staur(maxphi-iphip+2,ihil_iphicyl,2))
                   staufip = conjg(rotated_staufi(maxphi-iphip+2,ihil_iphicyl,2))
                End If

                ! P-M representation densities
                RHO_0  = rn+rp;         RHO_1  = rn-rp
                TAU_0  = tnt+tpt;       TAU_1  = tnt-tpt
                DRHO_0 = dn+dp;         DRHO_1 = dn-dp
                DJ_0   = djn+djp;       DJ_1   = djn-djp

                If(PNP_is_on == 0 .Or. PNP_is_on == 1) Then
                   ! Mixed prescription
                   rho_coupling=RHO_0
                   rho_exchange=rp
                Else If(PNP_is_on == 2) Then
                   ! Projected prescription
                   If(int(sigma) == 0) rho_coupling=RHOPROJ_0
                   If(int(sigma) == 1) rho_coupling=RHO_0
                   rho_exchange=rprojp
                End If

                S2_0  = (srozn+srozp)**2+(srorn+srorp)**2+(srofin+srofip)**2
                S2_1  = (srozn-srozp)**2+(srorn-srorp)**2+(srofin-srofip)**2
                J2_0  = (jzn+jzp)**2+(jrn+jrp)**2+(jfin+jfip)**2
                J2_1  = (jzn-jzp)**2+(jrn-jrp)**2+(jfin-jfip)**2
                SDS_0 = (srozn+srozp)*(dsrozn+dsrozp)+(srorn+srorp)*(dsrorn+dsrorp)+(srofin+srofip)*(dsrofin+dsrofip)
                SDS_1 = (srozn-srozp)*(dsrozn-dsrozp)+(srorn-srorp)*(dsrorn-dsrorp)+(srofin-srofip)*(dsrofin-dsrofip)
                JDS_0 = (jzn+jzp)*(csrozn+csrozp)+(jrn+jrp)*(csrorn+csrorp)+(jfin+jfip)*(csrofin+csrofip)
                JDS_1 = (jzn-jzp)*(csrozn-csrozp)+(jrn-jrp)*(csrorn-csrorp)+(jfin-jfip)*(csrofin-csrofip)

                If(use_j2terms .Eqv. .True.) Then
                   SFIZ_0 = sfizn+sfizp;   SFIZ_1 = sfizn-sfizp
                   SFIR_0 = sfirn+sfirp;   SFIR_1 = sfirn-sfirp
                   SZFI_0 = szfin+szfip;   SZFI_1 = szfin-szfip;
                   SRFI_0 = srfin+srfip;   SRFI_1 = srfin-srfip;
                   SZZ_0  = szzn+szzp;     SZZ_1  = szzn-szzp
                   SRZ_0  = srzn+srzp;     SRZ_1  = srzn-srzp
                   SRR_0  = srrn+srrp;     SRR_1  = srrn-srrp
                   SZR_0  = szrn+szrp;     SZR_1  = szrn-szrp
                   SFIFI_0= sfifin+sfifip; SFIFI_1= sfifin-sfifip
                   TENS_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2+SZZ_0**2+SRZ_0**2+SRR_0**2+SZR_0**2+SFIFI_0**2
                   TENS_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2+SZZ_1**2+SRZ_1**2+SRR_1**2+SZR_1**2+SFIFI_1**2
                   ST_0  = (srozn+srozp)*(stauzn+stauzp)+(srorn+srorp)*(staurn+staurp)+(srofin+srofip)*(staufin+staufip)
                   ST_1  = (srozn-srozp)*(stauzn-stauzp)+(srorn-srorp)*(staurn-staurp)+(srofin-srofip)*(staufin-staufip)
                End If

                ! -----------------------------------------------------
                ! Energies for fixed (z,rperp) point and varying phicyl
                ! -----------------------------------------------------
                ! Time-even components
                ! Kinetic energy
                ekinN_phicyl(iphicyl)=hb0n*(TAU_0+TAU_1)*half*facECM                                   ! kinetic neutrons
                ekinP_phicyl(iphicyl)=hb0p*(TAU_0-TAU_1)*half*facECM                                   ! kinetic protons
                ! Coulomb energy
                ecodi_phicyl(iphicyl)=c_zero; ecoex_phicyl(iphicyl)=c_zero
                If (icou >= 1) ecodi_phicyl(iphicyl)=half*ccou*rp                                      ! Coulomb direct
                If (icou == 2.Or.icou == -4) ecoex_phicyl(iphicyl)=CExPar*cex*rho_exchange**p43        ! Coulomb exchange, Slater approximation
                ! Volume energy
                EVOL_rho_rho_phicyl(iphicyl)=(Crho(0)+Cdrho(0)*rho_coupling**sigma)*RHO_0**2+&         ! volume rho-rho
                                             (Crho(1)+Cdrho(1)*rho_coupling**sigma)*RHO_1**2
                EVOL_rho_tau_phicyl(iphicyl)=Ctau(0)*RHO_0*TAU_0+Ctau(1)*RHO_1*TAU_1                   ! volume rho tau
                ! Surface energy
                ESURF_rho_drho_phicyl(iphicyl)=CrDr(0)*RHO_0*DRHO_0+CrDr(1)*RHO_1*DRHO_1               ! surface rho delta rho
                ! Spin orbit
                ESO_rho_nablaj_phicyl(iphicyl)=CrdJ(0)*RHO_0*DJ_0+CrdJ(1)*RHO_1*DJ_1
                ! Pairing energy
                eptN_phicyl(iphicyl) = CpV0(0)*(one - (RHO_0/rho_c)*CpV1(0))*akn2
                eptP_phicyl(iphicyl) = CpV0(1)*(one - (RHO_0/rho_c)*CpV1(1))*akp2

                ! Time-odd components
                ess0_phicyl(iphicyl)=Cs(0)*S2_0+Cs(1)*S2_1
                essD_phicyl(iphicyl)=Cds(0)*rho_coupling**sigma*S2_0+Cds(1)*rho_coupling**sigma*S2_1
                ecurr_phicyl(iphicyl)=Csj(0)*J2_0+Csj(1)*J2_1
                esDs_phicyl(iphicyl)=Csds(0)*SDS_0+Csds(1)*SDS_1
                ejNxs_phicyl(iphicyl)=CsdJ(0)*JDS_0+CsdJ(1)*JDS_1

                ! Tensor components
                If(use_j2terms .Eqv. .True.) Then
                   ETENS_phicyl(iphicyl)=CJ(0)*TENS_0+CJ(1)*TENS_1
                   esT_phicyl(iphicyl)=CT(0)*ST_0+CT(1)*ST_1
                End If

                ! Rms and pairing gaps
                row1_phicyl(iphicyl)=rn; row2_phicyl(iphicyl)=rp
                delN_phicyl(iphicyl) = CpV0(0)*(one - (RHO_0/rho_c)*CpV1(0))*adn
                delP_phicyl(iphicyl) = CpV0(1)*(one - (RHO_0/rho_c)*CpV1(1))*adp
                !
             End Do ! iphicyl

             ! Integration over phicyl
             Call integrate_complex_function(ekinN_phicyl,ngphi,phicyl_integration_step,ekinN_integral)
             Call integrate_complex_function(ekinP_phicyl,ngphi,phicyl_integration_step,ekinP_integral)
             Call integrate_complex_function(ecodi_phicyl,ngphi,phicyl_integration_step,ecodi_integral)
             Call integrate_complex_function(ecoex_phicyl,ngphi,phicyl_integration_step,ecoex_integral)
             Call integrate_complex_function(EVOL_rho_tau_phicyl,ngphi,phicyl_integration_step,EVOL_rho_tau_integral)
             Call integrate_complex_function(EVOL_rho_rho_phicyl,ngphi,phicyl_integration_step,EVOL_rho_rho_integral)
             Call integrate_complex_function(ESURF_rho_drho_phicyl,ngphi,phicyl_integration_step,ESURF_rho_drho_integral)
             Call integrate_complex_function(ESO_rho_nablaj_phicyl,ngphi,phicyl_integration_step,ESO_rho_nablaj_integral)
             Call integrate_complex_function(eptN_phicyl,ngphi,phicyl_integration_step,eptN_integral)
             Call integrate_complex_function(eptP_phicyl,ngphi,phicyl_integration_step,eptP_integral)

             Call integrate_complex_function(ess0_phicyl,ngphi,phicyl_integration_step,ess0_integral)
             Call integrate_complex_function(essD_phicyl,ngphi,phicyl_integration_step,essD_integral)
             Call integrate_complex_function(ecurr_phicyl,ngphi,phicyl_integration_step,ecurr_integral)
             Call integrate_complex_function(esDs_phicyl,ngphi,phicyl_integration_step,esDs_integral)
             Call integrate_complex_function(ejNxs_phicyl,ngphi,phicyl_integration_step,ejNxs_integral)

             Call integrate_complex_function(row1_phicyl,ngphi,phicyl_integration_step,row1_integral)
             Call integrate_complex_function(row2_phicyl,ngphi,phicyl_integration_step,row2_integral)
             Call integrate_complex_function(delN_phicyl,ngphi,phicyl_integration_step,delN_integral)
             Call integrate_complex_function(delP_phicyl,ngphi,phicyl_integration_step,delP_integral)

             ! Integration over (z,rperp)
             ekinN          = ekinN + ekinN_integral*whl
             ekinP          = ekinP + ekinP_integral*whl
             ecodi          = ecodi + ecodi_integral*whl
             ecoex          = ecoex - ecoex_integral*whl
             EVOL_rho_tau   = EVOL_rho_tau + EVOL_rho_tau_integral*whl
             EVOL_rho_rho   = EVOL_rho_rho + EVOL_rho_rho_integral*whl
             ESURF_rho_drho = ESURF_rho_drho + ESURF_rho_drho_integral*whl
             ESO_rho_nablaj = ESO_rho_nablaj + ESO_rho_nablaj_integral*whl
             eptN           = eptN + eptN_integral*whl
             eptP           = eptP + eptP_integral*whl

             ess0           = ess0  + ess0_integral*whl
             essD           = essD  + essD_integral*whl
             ecurr          = ecurr + ecurr_integral*whl
             esDs           = esDs  + esDs_integral*whl
             ejNxs          = ejNxs + ejNxs_integral*whl

             xn1            = xn1 + row1_integral*whl
             xn2            = xn2 + row2_integral*whl
             rms1           = rms1+ row1_integral*rrr*whl
             rms2           = rms2+ row2_integral*rrr*whl
             delN           = delN + delN_integral*whl
             delP           = delP + delP_integral*whl

             If(use_j2terms .Eqv. .True.) Then
                Call integrate_complex_function(ETENS_phicyl,ngphi,phicyl_integration_step,ETENS_integral)
                Call integrate_complex_function(esT_phicyl,ngphi,phicyl_integration_step,esT_integral)

                ETENS = ETENS + ETENS_integral*whl
                esT   = esT + esT_integral*whl
             End If


          End Do ! ihil

          lloc = (ibet-offset(team_rank)-1)*maxphi*maxphi*stride_size+(iphin_iphip-1)*stride_size
          all_energies(lloc+ 1) = ekinN*twopii
          all_energies(lloc+ 2) = ekinP*twopii
          all_energies(lloc+ 3) = ecodi*twopii
          all_energies(lloc+ 4) = ecoex*twopii
          all_energies(lloc+ 5) = EVOL_rho_tau*twopii
          all_energies(lloc+ 6) = EVOL_rho_rho*twopii
          all_energies(lloc+ 7) = ESURF_rho_drho*twopii
          all_energies(lloc+ 8) = ESO_rho_nablaj*twopii
          all_energies(lloc+ 9) = eptN*twopii
          all_energies(lloc+10) = eptP*twopii
          all_energies(lloc+11) = ess0*twopii
          all_energies(lloc+12) = essD*twopii
          all_energies(lloc+13) = ecurr*twopii
          all_energies(lloc+14) = esDs*twopii
          all_energies(lloc+15) = ejNxs*twopii
          all_energies(lloc+16) = ETENS*twopii
          all_energies(lloc+17) = esT*twopii
          all_energies(lloc+18) = xn1*twopii
          all_energies(lloc+19) = xn2*twopii
          all_energies(lloc+20) = Sqrt(rms1/xn1)
          all_energies(lloc+21) = Sqrt(rms2/xn2)
          all_energies(lloc+22) = delN/xn1
          all_energies(lloc+23) = delP/xn2

       End Do ! iphip
    End Do ! iphin

    ! Timing
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(5)=time_ibet(5)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)

  End Subroutine calculate_energies
  !===============================================================
  !> This routine deallocates quantities specific to each angle.
  !===============================================================
  Subroutine finalize_angle()

    Deallocate(rotation_matrix,inverse_rotation_matrix,rotated_density,rotated_kappa,rotated_kappac)
    Deallocate(rotated_ro,rotated_tau,rotated_dj,rotated_dro,rotated_aka,rotated_akac,ro_projected,&
               rotated_sroz,rotated_sror,rotated_srofi,rotated_jz,rotated_jr,rotated_jfi,&
               rotated_dsroz,rotated_dsror,rotated_dsrofi,rotated_csroz,rotated_csror,rotated_csrofi,&
               rotated_sfiz,rotated_sfir,rotated_srfi,rotated_szfi,rotated_szz,rotated_srz,rotated_srr,&
               rotated_szr,rotated_sfifi,rotated_stauz,rotated_staur,rotated_staufi)
    Deallocate(iosc1_contributing,iosc2_contributing,nt_contributing)
    Deallocate(cou_rotated)

  End Subroutine finalize_angle
  !==============================================================
  !> This routine calculates the symmetry-projected quantities.
  !==============================================================
  Subroutine project()
    Use math, Only : wigner,iv
    Use HFBTHO, Only : AMP_is_on,PNP_is_on,force_parity,team_size

    Integer(ipr) :: ibet,iphi,jj,iphin,iphip,jjmax,nn_neut,nn_prot,icount,itask,i_beta,it
    Real(pr) :: beta,facj,betafac,parity_prefac,parity_prefac_pi
    Complex(pr) :: fac_overlap,fac_overlap_pi,cj,ekinN,ekinP,ecodi,ecoex,EVOL_rho_tau,EVOL_rho_rho,ESURF_rho_drho,ESO_rho_nablaj,&
                   eptN,eptP,ess0,essD,ecurr,esDs,ejNxs,ETENS,esT,xn1,xn2,rms1,rms2,delN,delP
    Complex(pr) :: cj_phi,ekinN_phi,ekinP_phi,ecodi_phi,ecoex_phi,EVOL_rho_tau_phi,EVOL_rho_rho_phi,ESURF_rho_drho_phi,ESO_rho_nablaj_phi,eptN_phi,eptP_phi,&
                   ess0_phi,essD_phi,ecurr_phi,esDs_phi,ejNxs_phi,ETENS_phi,esT_phi,xn1_phi,xn2_phi,rms1_phi,rms2_phi,delN_phi,delP_phi,cj_jj
    Complex(pr), Allocatable :: cj_NP(:,:)
    Integer(ipr) :: tstart_ibet,tfinish_ibet
    Complex(pr), Allocatable :: ekinN_rotated(:,:,:),ekinP_rotated(:,:,:),ecodi_rotated(:,:,:),ecoex_rotated(:,:,:),EVOL_rho_tau_rotated(:,:,:),&
                                EVOL_rho_rho_rotated(:,:,:),ESURF_rho_drho_rotated(:,:,:),ESO_rho_nablaj_rotated(:,:,:),eptN_rotated(:,:,:),eptP_rotated(:,:,:),&
                                ess0_rotated(:,:,:),essD_rotated(:,:,:),ecurr_rotated(:,:,:),esDs_rotated(:,:,:), ejNxs_rotated(:,:,:),ETENS_rotated(:,:,:),&
                                esT_rotated(:,:,:),xn1_rotated(:,:,:),xn2_rotated(:,:,:),rms1_rotated(:,:,:),rms2_rotated(:,:,:),delN_rotated(:,:,:),delP_rotated(:,:,:)
    ! Timing
    Call system_clock(tstart_ibet,clock_rate)

    ! Maximal angular momentum considered
    If(AMP_is_on == 1) Then
       jjmax=maxj
    Else
       jjmax=0
    End If

    ! Unwrapping rotated overlaps and energies
    Allocate(ekinN_rotated(maxbet,maxphi,maxphi),ekinP_rotated(maxbet,maxphi,maxphi),ecodi_rotated(maxbet,maxphi,maxphi),ecoex_rotated(maxbet,maxphi,maxphi),&
             EVOL_rho_tau_rotated(maxbet,maxphi,maxphi),EVOL_rho_rho_rotated(maxbet,maxphi,maxphi),ESURF_rho_drho_rotated(maxbet,maxphi,maxphi),&
             ESO_rho_nablaj_rotated(maxbet,maxphi,maxphi),eptN_rotated(maxbet,maxphi,maxphi),eptP_rotated(maxbet,maxphi,maxphi),ess0_rotated(maxbet,maxphi,maxphi),&
             essD_rotated(maxbet,maxphi,maxphi),ecurr_rotated(maxbet,maxphi,maxphi),esDs_rotated(maxbet,maxphi,maxphi),ejNxs_rotated(maxbet,maxphi,maxphi),&
             ETENS_rotated(maxbet,maxphi,maxphi),esT_rotated(maxbet,maxphi,maxphi),xn1_rotated(maxbet,maxphi,maxphi),xn2_rotated(maxbet,maxphi,maxphi),&
             rms1_rotated(maxbet,maxphi,maxphi),rms2_rotated(maxbet,maxphi,maxphi),delN_rotated(maxbet,maxphi,maxphi),delP_rotated(maxbet,maxphi,maxphi))

    icount=1
    Do itask=1,team_size
       Do it=1,2
          Do iphi=1,maxphi
             Do ibet=1,block_vec(itask-1)
                i_beta = ibet + offset(itask-1)
                rotated_overlap(i_beta,iphi,it)=all_overlaps_gthr(icount); icount=icount+1
             End Do
          End Do
       End Do
    End Do

    icount=1
    Do itask=1,team_size
       Do ibet=1,block_vec(itask-1)
          i_beta = ibet + offset(itask-1)
          Do iphin=1,maxphi
             Do iphip=1,maxphi
                ekinN_rotated(i_beta,iphin,iphip)         = all_energies_gthr(icount); icount=icount+1
                ekinP_rotated(i_beta,iphin,iphip)         = all_energies_gthr(icount); icount=icount+1
                ecodi_rotated(i_beta,iphin,iphip)         = all_energies_gthr(icount); icount=icount+1
                ecoex_rotated(i_beta,iphin,iphip)         = all_energies_gthr(icount); icount=icount+1
                EVOL_rho_tau_rotated(i_beta,iphin,iphip)  = all_energies_gthr(icount); icount=icount+1
                EVOL_rho_rho_rotated(i_beta,iphin,iphip)  = all_energies_gthr(icount); icount=icount+1
                ESURF_rho_drho_rotated(i_beta,iphin,iphip)= all_energies_gthr(icount); icount=icount+1
                ESO_rho_nablaj_rotated(i_beta,iphin,iphip)= all_energies_gthr(icount); icount=icount+1
                eptN_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
                eptP_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1

                ess0_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
                essD_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
                ecurr_rotated(i_beta,iphin,iphip)         = all_energies_gthr(icount); icount=icount+1
                esDs_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
                ejNxs_rotated(i_beta,iphin,iphip)         = all_energies_gthr(icount); icount=icount+1

                ETENS_rotated(i_beta,iphin,iphip)         = all_energies_gthr(icount); icount=icount+1
                esT_rotated(i_beta,iphin,iphip)           = all_energies_gthr(icount); icount=icount+1

                xn1_rotated(i_beta,iphin,iphip)           = all_energies_gthr(icount); icount=icount+1
                xn2_rotated(i_beta,iphin,iphip)           = all_energies_gthr(icount); icount=icount+1
                rms1_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
                rms2_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
                delN_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
                delP_rotated(i_beta,iphin,iphip)          = all_energies_gthr(icount); icount=icount+1
             End Do
          End Do
       End Do
    End Do

    ! Projecting on desired number of particles and on different angular momenta
    parity_prefac=0.5_pr; parity_prefac_pi=0.5_pr
    If(force_parity) Then
       parity_prefac=1.0_pr
       parity_prefac_pi=0.0_pr
    End If

    Do jj=0,jjmax,jjstep
       facj=half*(two*jj+one)*integration_prefactor
       cj=c_zero; ekinN=c_zero; ekinP=c_zero; ecodi=c_zero; ecoex=c_zero; EVOL_rho_tau=c_zero; EVOL_rho_rho=c_zero; ESURF_rho_drho=c_zero
       ESO_rho_nablaj=c_zero; eptN=c_zero; eptP=c_zero; ess0=c_zero; essD=c_zero; ecurr=c_zero; esDs=c_zero; ejNxs=c_zero
       ETENS=c_zero; esT=c_zero; xn1=c_zero; xn2=c_zero; rms1=c_zero; rms2=c_zero; delN=c_zero; delP=c_zero
       ! Angular momentum projection
       Do ibet=1,maxbet
          beta=betabs(ibet);betafac=betaweight(ibet)*wigner(jj,0,0,beta)
          ! Particle number projection
          cj_phi=c_zero;
          ekinN_phi = c_zero; ekinP_phi = c_zero; ecodi_phi = c_zero; ecoex_phi = c_zero; EVOL_rho_tau_phi = c_zero
          EVOL_rho_rho_phi = c_zero; ESURF_rho_drho_phi = c_zero; ESO_rho_nablaj_phi = c_zero; eptN_phi = c_zero; eptP_phi = c_zero;
          ess0_phi = c_zero; essD_phi = c_zero; ecurr_phi = c_zero; esDs_phi = c_zero; ejNxs_phi = c_zero; ETENS_phi = c_zero; esT_phi = c_zero
          xn1_phi = c_zero; xn2_phi = c_zero; rms1_phi = c_zero; rms2_phi = c_zero; delN_phi = c_zero; delP_phi = c_zero
          Do iphin=1,maxphi
             Do iphip=1,maxphi
                fac_overlap        = parity_prefac*rotated_overlap(ibet,iphin,1)*rotated_overlap(ibet,iphip,2)*ephicN(iphin,0)*ephicP(iphip,0)
                fac_overlap_pi     = parity_prefac_pi*iv(jj)*rotated_overlap(maxbet+1-ibet,iphin,1)*rotated_overlap(maxbet+1-ibet,iphip,2)*ephicN(iphin,0)*ephicP(iphip,0)
                cj_phi             = cj_phi + fac_overlap+fac_overlap_pi

                ekinN_phi          = ekinN_phi + fac_overlap*ekinN_rotated(ibet,iphin,iphip)+fac_overlap_pi*ekinN_rotated(maxbet+1-ibet,iphin,iphip)
                ekinP_phi          = ekinP_phi + fac_overlap*ekinP_rotated(ibet,iphin,iphip)+fac_overlap_pi*ekinP_rotated(maxbet+1-ibet,iphin,iphip)
                ecodi_phi          = ecodi_phi + fac_overlap*ecodi_rotated(ibet,iphin,iphip)+fac_overlap_pi*ecodi_rotated(maxbet+1-ibet,iphin,iphip)
                ecoex_phi          = ecoex_phi + fac_overlap*ecoex_rotated(ibet,iphin,iphip)+fac_overlap_pi*ecoex_rotated(maxbet+1-ibet,iphin,iphip)
                EVOL_rho_tau_phi   = EVOL_rho_tau_phi + fac_overlap*EVOL_rho_tau_rotated(ibet,iphin,iphip)+fac_overlap_pi*EVOL_rho_tau_rotated(maxbet+1-ibet,iphin,iphip)
                EVOL_rho_rho_phi   = EVOL_rho_rho_phi + fac_overlap*EVOL_rho_rho_rotated(ibet,iphin,iphip)+fac_overlap_pi*EVOL_rho_rho_rotated(maxbet+1-ibet,iphin,iphip)
                ESURF_rho_drho_phi = ESURF_rho_drho_phi + fac_overlap*ESURF_rho_drho_rotated(ibet,iphin,iphip) &
                                                         +fac_overlap_pi*ESURF_rho_drho_rotated(maxbet+1-ibet,iphin,iphip)
                ESO_rho_nablaj_phi = ESO_rho_nablaj_phi + fac_overlap*ESO_rho_nablaj_rotated(ibet,iphin,iphip) &
                                                        + fac_overlap_pi*ESO_rho_nablaj_rotated(maxbet+1-ibet,iphin,iphip)
                eptN_phi           = eptN_phi  + fac_overlap*eptN_rotated(ibet,iphin,iphip)  + fac_overlap_pi*eptN_rotated(maxbet+1-ibet,iphin,iphip)
                eptP_phi           = eptP_phi  + fac_overlap*eptP_rotated(ibet,iphin,iphip)  + fac_overlap_pi*eptP_rotated(maxbet+1-ibet,iphin,iphip)

                ess0_phi           = ess0_phi  + fac_overlap*ess0_rotated(ibet,iphin,iphip)  + fac_overlap_pi*ess0_rotated(maxbet+1-ibet,iphin,iphip)
                essD_phi           = essD_phi  + fac_overlap*essD_rotated(ibet,iphin,iphip)  + fac_overlap_pi*essD_rotated(maxbet+1-ibet,iphin,iphip)
                ecurr_phi          = ecurr_phi + fac_overlap*ecurr_rotated(ibet,iphin,iphip) + fac_overlap_pi*ecurr_rotated(maxbet+1-ibet,iphin,iphip)
                esDs_phi           = esDs_phi  + fac_overlap*esDs_rotated(ibet,iphin,iphip)  + fac_overlap_pi*esDs_rotated(maxbet+1-ibet,iphin,iphip)
                ejNxs_phi          = ejNxs_phi + fac_overlap*ejNxs_rotated(ibet,iphin,iphip) + fac_overlap_pi*ejNxs_rotated(maxbet+1-ibet,iphin,iphip)

                ETENS_phi          = ETENS_phi + fac_overlap*ETENS_rotated(ibet,iphin,iphip) + fac_overlap_pi*ETENS_rotated(maxbet+1-ibet,iphin,iphip)
                esT_phi            = esT_phi   + fac_overlap*esT_rotated(ibet,iphin,iphip)   + fac_overlap_pi*esT_rotated(maxbet+1-ibet,iphin,iphip)

                xn1_phi            = xn1_phi   + fac_overlap*xn1_rotated(ibet,iphin,iphip)   + fac_overlap_pi*xn1_rotated(maxbet+1-ibet,iphin,iphip)
                xn2_phi            = xn2_phi   + fac_overlap*xn2_rotated(ibet,iphin,iphip)   + fac_overlap_pi*xn2_rotated(maxbet+1-ibet,iphin,iphip)
                rms1_phi           = rms1_phi  + fac_overlap*rms1_rotated(ibet,iphin,iphip)  + fac_overlap_pi*rms1_rotated(maxbet+1-ibet,iphin,iphip)
                rms2_phi           = rms2_phi  + fac_overlap*rms2_rotated(ibet,iphin,iphip)  + fac_overlap_pi*rms2_rotated(maxbet+1-ibet,iphin,iphip)
                delN_phi           = delN_phi  + fac_overlap*delN_rotated(ibet,iphin,iphip)  + fac_overlap_pi*delN_rotated(maxbet+1-ibet,iphin,iphip)
                delP_phi           = delP_phi  + fac_overlap*delP_rotated(ibet,iphin,iphip)  + fac_overlap_pi*delP_rotated(maxbet+1-ibet,iphin,iphip)
             End Do ! iphip
          End Do ! iphin
          cj_phi=cj_phi/maxphi**2; ekinN_phi=ekinN_phi/maxphi**2; ekinP_phi=ekinP_phi/maxphi**2; ecodi_phi=ecodi_phi/maxphi**2; ecoex_phi=ecoex_phi/maxphi**2
          EVOL_rho_tau_phi=EVOL_rho_tau_phi/maxphi**2; EVOL_rho_rho_phi=EVOL_rho_rho_phi/maxphi**2; ESURF_rho_drho_phi=ESURF_rho_drho_phi/maxphi**2
          ESO_rho_nablaj_phi=ESO_rho_nablaj_phi/maxphi**2; eptN_phi=eptN_phi/maxphi**2; eptP_phi=eptP_phi/maxphi**2
          ess0_phi=ess0_phi/maxphi**2; essD_phi=essD_phi/maxphi**2; ecurr_phi=ecurr_phi/maxphi**2; esDs_phi=esDs_phi/maxphi**2; ejNxs_phi=ejNxs_phi/maxphi**2
          ETENS_phi=ETENS_phi/maxphi**2; esT_phi=esT_phi/maxphi**2; xn1_phi=xn1_phi/maxphi**2; xn2_phi=xn2_phi/maxphi**2;
          rms1_phi=rms1_phi/maxphi**2; rms2_phi=rms2_phi/maxphi**2; delN_phi=delN_phi/maxphi**2; delP_phi=delP_phi/maxphi**2

          If(AMP_is_on == 1) Then
             cj=cj+cj_phi*betafac
             ekinN=ekinN+betafac*ekinN_phi; ekinP=ekinP+betafac*ekinP_phi
             ecodi=ecodi+betafac*ecodi_phi; ecoex=ecoex+betafac*ecoex_phi
             EVOL_rho_tau=EVOL_rho_tau+betafac*EVOL_rho_tau_phi
             EVOL_rho_rho=EVOL_rho_rho+betafac*EVOL_rho_rho_phi
             ESURF_rho_drho=ESURF_rho_drho+betafac*ESURF_rho_drho_phi
             ESO_rho_nablaj=ESO_rho_nablaj+betafac*ESO_rho_nablaj_phi
             eptN=eptN+betafac*eptN_phi; eptP=eptP+betafac*eptP_phi

             ess0=ess0+betafac*ess0_phi; essD=essD+betafac*essD_phi
             ecurr=ecurr+betafac*ecurr_phi; esDs=esDs+betafac*esDs_phi
             ejNxs=ejNxs+betafac*ejNxs_phi

             ETENS=ETENS+betafac*ETENS_phi; esT=esT+betafac*esT_phi

             xn1=xn1+betafac*xn1_phi; xn2=xn2+betafac*xn2_phi
             rms1=rms1+betafac*rms1_phi; rms2=rms2+betafac*rms2_phi
             delN=delN+betafac*delN_phi; delP=delP+betafac*delP_phi
          Else
             projected_overlap(jj)=Dble(cj_phi)
             projected_ekinN(jj)=Dble(ekinN_phi)/projected_overlap(jj); projected_ekinP(jj)=Dble(ekinP_phi)/projected_overlap(jj)
             projected_ecodi(jj)=Dble(ecodi_phi)/projected_overlap(jj); projected_ecoex(jj)=Dble(ecoex_phi)/projected_overlap(jj)
             projected_EVOL_rho_tau(jj)=Dble(EVOL_rho_tau_phi)/projected_overlap(jj)
             projected_EVOL_rho_rho(jj)=Dble(EVOL_rho_rho_phi)/projected_overlap(jj)
             projected_ESURF_rho_drho(jj)=Dble(ESURF_rho_drho_phi)/projected_overlap(jj)
             projected_ESO_rho_nablaj(jj)=Dble(ESO_rho_nablaj_phi)/projected_overlap(jj)
             projected_eptN(jj)=Dble(eptN_phi)/projected_overlap(jj); projected_eptP(jj)=Dble(eptP_phi)/projected_overlap(jj)
             projected_ess0(jj)=Dble(ess0_phi)/projected_overlap(jj); projected_essD(jj)=Dble(essD_phi)/projected_overlap(jj)
             projected_ecurr(jj)=Dble(ecurr_phi)/projected_overlap(jj); projected_esDs(jj)=Dble(esDs_phi)/projected_overlap(jj)
             projected_ejNxs(jj)=Dble(ejNxs_phi)/projected_overlap(jj)
             projected_ETENS(jj)=Dble(ETENS_phi)/projected_overlap(jj); projected_esT(jj)=Dble(esT_phi)/projected_overlap(jj)
             projected_xn1(jj)=Dble(xn1_phi)/projected_overlap(jj); projected_xn2(jj)=Dble(xn2_phi)/projected_overlap(jj)
             projected_rms1(jj)=Dble(rms1_phi)/projected_overlap(jj); projected_rms2(jj)=Dble(rms2_phi)/projected_overlap(jj)
             projected_delN(jj)=Dble(delN_phi)/projected_overlap(jj); projected_delP(jj)=Dble(delP_phi)/projected_overlap(jj)
          End If

       End Do ! ibet

       If(AMP_is_on == 1) Then
          projected_overlap(jj)=Dble(cj)*facj
          projected_ekinN(jj)=Dble(ekinN)*facj/projected_overlap(jj); projected_ekinP(jj)=Dble(ekinP)*facj/projected_overlap(jj)
          projected_ecodi(jj)=Dble(ecodi)*facj/projected_overlap(jj); projected_ecoex(jj)=Dble(ecoex)*facj/projected_overlap(jj)
          projected_EVOL_rho_tau(jj)=Dble(EVOL_rho_tau)*facj/projected_overlap(jj)
          projected_EVOL_rho_rho(jj)=Dble(EVOL_rho_rho)*facj/projected_overlap(jj)
          projected_ESURF_rho_drho(jj)=Dble(ESURF_rho_drho)*facj/projected_overlap(jj)
          projected_ESO_rho_nablaj(jj)=Dble(ESO_rho_nablaj)*facj/projected_overlap(jj)
          projected_eptN(jj)=Dble(eptN)*facj/projected_overlap(jj); projected_eptP(jj)=Dble(eptP)*facj/projected_overlap(jj)
          projected_ess0(jj)=Dble(ess0)*facj/projected_overlap(jj); projected_essD(jj)=Dble(essD)*facj/projected_overlap(jj)
          projected_ecurr(jj)=Dble(ecurr)*facj/projected_overlap(jj); projected_esDs(jj)=Dble(esDs)*facj/projected_overlap(jj)
          projected_ejNxs(jj)=Dble(ejNxs)*facj/projected_overlap(jj)
          projected_ETENS(jj)=Dble(ETENS)*facj/projected_overlap(jj); projected_esT(jj)=Dble(esT)*facj/projected_overlap(jj)
          projected_xn1(jj)=Dble(xn1)*facj/projected_overlap(jj); projected_xn2(jj)=Dble(xn2)*facj/projected_overlap(jj)
          projected_rms1(jj)=Dble(rms1)*facj/projected_overlap(jj); projected_rms2(jj)=Dble(rms2)*facj/projected_overlap(jj)
          projected_delN(jj)=Dble(delN)*facj/projected_overlap(jj); projected_delP(jj)=Dble(delP)*facj/projected_overlap(jj)
       End If

    End Do ! jj

    ! Projecting on different numbers of particles
    If(PNP_is_on  >  0) Then
       Allocate(cj_NP(-maxN/2:maxN/2,-maxP/2:maxP/2))
       Do jj=0,jjmax,jjstep
          facj=half*(two*jj+one)*integration_prefactor
          cj_NP=c_zero
          Do ibet=1,maxbet
             beta=betabs(ibet);betafac=betaweight(ibet)*wigner(jj,0,0,beta)
             Do nn_neut=-maxN/2,maxN/2
                Do nn_prot=-maxP/2,maxP/2
                   cj_phi=c_zero
                   Do iphin=1,maxphi
                      Do iphip=1,maxphi
                         fac_overlap=parity_prefac*rotated_overlap(ibet,iphin,1)*rotated_overlap(ibet,iphip,2)*ephicN(iphin,nn_neut)*ephicP(iphip,nn_prot)
                         fac_overlap_pi=parity_prefac_pi*iv(jj)*rotated_overlap(maxbet+1-ibet,iphin,1)*rotated_overlap(maxbet+1-ibet,iphip,2)*ephicN(iphin,nn_neut)*ephicP(iphip,nn_prot)
                         cj_phi=cj_phi+fac_overlap+fac_overlap_pi
                      End Do
                   End Do
                   cj_phi=cj_phi/maxphi**2
                   If(AMP_is_on == 1) Then
                      cj_NP(nn_neut,nn_prot)=cj_NP(nn_neut,nn_prot)+cj_phi*betafac
                   Else
                      projected_NP(jj,nn_neut,nn_prot)=Dble(cj_phi)
                   End If
                End Do ! nn_neut
             End Do ! nn_prot
          End Do ! ibet
          If(AMP_is_on == 1) Then
             Do nn_neut=-maxN/2,maxN/2
                Do nn_prot=-maxP/2,maxP/2
                   projected_NP(jj,nn_neut,nn_prot)=Dble(cj_NP(nn_neut,nn_prot))*facj
                End Do
             End Do
          End If
       End Do ! jj
       Do jj=0,jjmax,jjstep
          cj_jj=c_zero
          Do nn_neut=-maxN/2,maxN/2
             Do nn_prot=-maxP/2,maxP/2
                cj_jj=cj_jj+projected_NP(jj,nn_neut,nn_prot)
             End Do
          End Do
          projected_NP_norm(jj)=Dble(cj_jj)
       End Do
       Deallocate(cj_NP)
    End If

    Deallocate(ekinN_rotated,ekinP_rotated,ecodi_rotated,ecoex_rotated,EVOL_rho_tau_rotated,EVOL_rho_rho_rotated,ESURF_rho_drho_rotated,&
               ESO_rho_nablaj_rotated,eptN_rotated,eptP_rotated,ess0_rotated,essD_rotated,ecurr_rotated,esDs_rotated,ejNxs_rotated,&
               ETENS_rotated,esT_rotated,xn1_rotated,xn2_rotated,rms1_rotated,rms2_rotated,delN_rotated,delP_rotated)

    ! Timing
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(6)=time_ibet(6)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)

  End Subroutine project
  !==================================================================
  !> This routine prints to file the symmetry-projected quantities.
  !==================================================================
  Subroutine print_project()
    Use HFBTHO, Only : AMP_is_on,PNP_is_on,amas,r00,tz
    Use UNEDF, Only : use_j2terms

    Integer(ipr) :: ibet,jj,beta_active_sum,jjmax,iphi,nn_neut,nn_prot
    Real(pr) :: cjtot,cnztot,ekt(3),ecoul,evol,esurf,eso,ept(3),ess0,essD,ecurr,esDs,ejNxs,etens,esT,eteven,etodd,eproj,&
                cjnorm,cjnormtot,cnznorm,cnznormtot
    Real(pr) :: xn(3),rms(3),del(2),r212,r222,RpRMSsq,RnRMSsq,DarwinFoldy,rc

    ! Angular momentum projection only
    If(AMP_is_on == 1 .And. PNP_is_on == 0) Then
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =================       Mesh in rotational angle beta      ================='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(3x,a,4x,a,2x,a,6x,a,12x,a)') 'ibet','angle','active','overlap neutrons','overlap protons'
       Write(lproj,'(3x,a,4x,a,2x,a,5x,a,5x,a)') '---','-------','---','-----------------------','-----------------------'
       beta_active_sum=0
       Do ibet=1,maxbet
          Write(lproj,'(2x,i3,3x,f9.6,2x,i2,3x,f12.9,1x,f12.9,"i",2x,f12.9,1x,f12.9,"i")') ibet,betabs(ibet),beta_active_gthr(ibet), &
                        rotated_overlap(ibet,1,1),rotated_overlap(ibet,1,2)
          beta_active_sum=beta_active_sum+beta_active_gthr(ibet)
       End Do
       Write(lproj,'(19x,a)') '---'
       Write(lproj,'(18x,i3)') beta_active_sum
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =============     Angular momentum content of an HFB state     ============='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(22x,a,10x,a,7x,a)') 'J','|cJ|^2','partial sum'
       Write(lproj,'(21x,a,8x,a,7x,a)') '---','--------','----------'
       cjtot=zero
       Do jj=0,maxj,jjstep
          cjtot=cjtot+abs(projected_overlap(jj))
          Write(lproj,'(20x,i3,7x,f10.7,7x,f10.7)')  jj,projected_overlap(jj),cjtot
       End Do
       Write(lproj,'(47x,a)') ' ---------'
       Write(lproj,'(47x,f10.7)') cjtot
    End If

    ! Particle number projection only
    If(AMP_is_on == 0 .And. PNP_is_on > 0) Then
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  ===================       Mesh in gauge angle phi        ==================='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(3x,a,4x,a,10x,a,12x,a)') 'iphi','angle','overlap neutrons','overlap protons'
       Write(lproj,'(3x,a,3x,a,6x,a,5x,a)') '---','--------','-----------------------','-----------------------'
       Do iphi=1,maxphi
          Write(lproj,'(2x,i3,3x,f9.6,3x,f12.9,1x,f12.9,"i",2x,f12.9,1x,f12.9,"i")') iphi,phiabs(iphi),rotated_overlap(1,iphi,1),&
                                                                                                       rotated_overlap(1,iphi,2)
       End Do
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  ===========     Particle-number decomposition of an HFB state    ==========='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(20x,a,7x,a,7x,a,5x,a)') 'N', 'Z', '|cNZ|^2', 'partial sum'
       Write(lproj,'(19x,a,5x,a,6x,a,5x,a)') '---', '---', '---------', '---------'
       cnztot=zero
       Do nn_neut=-maxN/2,maxN/2
          Do nn_prot=-maxP/2,maxP/2
             cnztot=cnztot+projected_NP(0,nn_neut,nn_prot)
             Write(lproj,'(16x,f6.2,2x,f6.2,3x,f12.9,2x,f12.9)') tz(1)+nn_neut*2, tz(2)+nn_prot*2, projected_NP(0,nn_neut,nn_prot),cnztot
          End Do
       End Do
       Write(lproj,'(50x,a)') '---------'
       Write(lproj,'(47x,f12.9)') cnztot
    End If

    ! Angular momentum and Particle number projection
    If(AMP_is_on == 1 .And. PNP_is_on > 0) Then
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  ==========   Mesh in rotational angle beta and gauge angle phi    =========='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(3x,a,4x,a,3x,a,2x,a,4x,a,10x,a,13x,a)') 'ibet','beta','active','iphi','phi','overlap neutrons','overlap protons'
       Write(lproj,'(3x,a,4x,a,2x,a,4x,a,3x,a,5x,a,5x,a)') '---','-------','---','---','-------','-----------------------','-----------------------'
       beta_active_sum=0
       Do ibet=1,maxbet
          beta_active_sum=beta_active_sum+beta_active_gthr(ibet)
          Do iphi=1,maxphi
             Write(lproj,'(2x,i3,3x,f9.6,2x,i2,4x,i3,3x,f9.6,2x,f12.9,1x,f12.9,"i",2x,f12.9,1x,f12.9,"i")') &
                           ibet,betabs(ibet),beta_active_gthr(ibet),&
                           iphi,phiabs(iphi),rotated_overlap(ibet,iphi,1),rotated_overlap(ibet,iphi,2)
          End Do
       End Do
       Write(lproj,'(19x,a)') '---'
       Write(lproj,'(18x,i3)') beta_active_sum
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =============    Angular momentum content of PNP-HFB state     ============='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(7x,a,10x,a,7x,a,6x,a,4x,a)') 'J','|cJ|^2','partial sum','|cJ|^2 (norm)','partial sum'
       Write(lproj,'(6x,a,8x,a,7x,a,7x,a,7x,a)') '---','--------','----------','----------','----------'
       cjnorm=zero
       Do jj=0,maxj,jjstep
          cjnorm=cjnorm+abs(projected_overlap(jj))
       End Do
       cjtot=zero;cjnormtot=zero
       Do jj=0,maxj,jjstep
          cjtot=cjtot+abs(projected_overlap(jj))
          cjnormtot=cjnormtot+abs(projected_overlap(jj))/cjnorm
          Write(lproj,'(5x,i3,7x,f10.7,7x,f10.7,7x,f10.7,7x,f10.7)')  jj,projected_overlap(jj),cjtot,projected_overlap(jj)/cjnorm,cjnormtot
       End Do
       Write(lproj,'(32x,a,24x,a)') ' ---------',' ---------'
       Write(lproj,'(32x,f10.7,24x,f10.7)') cjtot,cjnormtot
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =======     Particle-number decomposition of AMP-projected states    ======='
       Write(lproj,'(a)')     '  ============================================================================'
       jjmax=maxj
       If(maxj > 40) jjmax=40
       Do jj=0,jjmax,jjstep
          Write(lproj, '(32x,a10)') '----------'
          Write(lproj,'(30x,a7,1x,i3,a2)') '| J = ',jj,' |'
          Write(lproj, '(32x,a10)') '----------'
          Write(lproj,'(6x,a)') '--------------------------------------------------------------------'
          Write(lproj,'(7x,a,7x,a,7x,a,5x,a,2x,a,1x,a)') 'N', 'Z', '|cNZ|^2', 'partial sum','|cNZ|^2 (norm)','partial sum'
          Write(lproj,'(6x,a,5x,a,6x,a,5x,a,5x,a,5x,a)') '---', '---', '---------', '---------','---------', '---------'
          cnznorm=zero
          Do nn_neut=-maxN/2,maxN/2
             Do nn_prot=-maxP/2,maxP/2
                cnznorm=cnznorm+projected_NP(jj,nn_neut,nn_prot)
             End Do
          End Do
          cnztot=zero;cnznormtot=zero
          Do nn_neut=-maxN/2,maxN/2
             Do nn_prot=-maxP/2,maxP/2
                cnztot=cnztot+projected_NP(jj,nn_neut,nn_prot)
                cnznormtot=cnznormtot+projected_NP(jj,nn_neut,nn_prot)/cnznorm
                Write(lproj,'(3x,f6.2,2x,f6.2,3x,f12.9,2x,f12.9,2x,f12.9,2x,f12.9)') tz(1)+nn_neut*2, tz(2)+nn_prot*2,&
                      projected_NP(jj,nn_neut,nn_prot),cnztot,projected_NP(jj,nn_neut,nn_prot)/cnznorm,cnznormtot
             End Do
          End Do
          Write(lproj,'(37x,a,19x,a)') '---------','---------'
          Write(lproj,'(34x,f12.9,16x,f12.9)') cnztot,cnznormtot
       End Do
    End If

    Write(lproj,'(a)')     '  ============================================================================'
    Write(lproj,'(a)')     '  ====================     Symmetry-projected energies    ===================='
    Write(lproj,'(a)')     '  ============================================================================'
    If(AMP_is_on == 1) Then
       jjmax=maxj
       If(maxj > 60) jjmax=60
    Else
       jjmax=0
    End If

    Do jj=0,jjmax,jjstep
       ! Energies
       ekt(1)=projected_ekinN(jj); ekt(2)=projected_ekinP(jj); ekt(3)=ekt(1)+ekt(2)
       ecoul=projected_ecodi(jj)+projected_ecoex(jj)
       evol=projected_EVOL_rho_tau(jj)+projected_EVOL_rho_rho(jj); esurf=projected_ESURF_rho_drho(jj)
       eso=projected_ESO_rho_nablaj(jj)
       ept(1)=projected_eptN(jj); ept(2)=projected_eptP(jj); ept(3)=ept(1)+ept(2)
       ess0=projected_ess0(jj); essD=projected_essD(jj); ecurr=projected_ecurr(jj); esDs=projected_esDs(jj)
       ejNxs=projected_ejNxs(jj); etens=projected_ETENS(jj); esT=projected_esT(jj)

       eteven=ekt(3)+ecoul+evol+esurf+eso+ept(3)
       etodd=ess0+essD+ecurr+esDs+ejNxs
       If(use_j2terms .Eqv. .True.) Then
          eteven=eteven+etens
          etodd=etodd+esT
       End If

       eproj=eteven+etodd
       ! Radii
       xn(1)=projected_xn1(jj); xn(2)=projected_xn2(jj); xn(3)=xn(1)+xn(2)
       rms(1)=projected_rms1(jj); rms(2)=projected_rms2(jj)
       r212=rms(1)**2; r222=rms(2)**2
       rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
       del(1)=projected_delN(jj); del(2)=projected_delP(jj)
       ! Charge radius, from Adv. Nucl. Phys. 8, 219 (1975)
       RpRMSsq=0.769_pr
       RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
       DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
       rc=Sqrt(r222+RpRMSsq+(xn(1)/xn(2))*RnRMSsq+DarwinFoldy)
       If(AMP_is_on == 1) Then
          Write(lproj, '(12x,a10)') '----------'
          Write(lproj,'(10x,a7,1x,i3,a2)') '| J = ',jj,' |'
          Write(lproj, '(12x,a10)') '----------'
       End If
       Write(lproj,'(a,3f15.6)') '  rms-radius ..........',rms
       Write(lproj,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
       Write(lproj,'(a,2f15.6)')    '  delta(n,p) ..........',del(1),del(2)
       Write(lproj,'(a,3f15.6)')    '  kinetic energy ......',ekt(1),ekt(2),ekt(3)
       Write(lproj,'(a,3f15.6)')    '  pairing energy ......',ept
       Write(lproj,'(a,30x,f15.6)') '  volume energy .......',evol
       Write(lproj,'(a,30x,f15.6)') '        rho_tau .......',projected_EVOL_rho_tau(jj)
       Write(lproj,'(a,30x,f15.6)') '        rho_rho .......',projected_EVOL_rho_rho(jj)
       Write(lproj,'(a,30x,f15.6)') '  surface energy ......',esurf
       Write(lproj,'(a,30x,f15.6)') '   rho_DELTA_rho ......',projected_ESURF_rho_drho(jj)
       Write(lproj,'(a,30x,f15.6)') '  spin-orbit energy ...',eso
       Write(lproj,'(a,30x,f15.6)') '        rho_NABLA_J ...',projected_ESO_rho_nablaj(jj)
       Write(lproj,'(a,30x,f15.6)') '  coulomb energy ......',ecoul
       Write(lproj,'(a,30x,f15.6)') '          direct ......',projected_ecodi(jj)
       Write(lproj,'(a,30x,f15.6)') '          exchange ....',projected_ecoex(jj)
       If(use_j2terms .Eqv. .True.) Then
          Write(lproj,'(a,30x,f15.6)') '  tensor energy .......',etens
       End If
       Write(lproj,'(a,30x,f15.6)') '  time-odd terms ......',etodd
       Write(lproj,'(a,30x,f15.6)') '             s_s ------',ess0
       Write(lproj,'(a,30x,f15.6)') '             s_s (DD) -',essD
       Write(lproj,'(a,30x,f15.6)') '             j_j .....-',ecurr
       Write(lproj,'(a,30x,f15.6)') '             s_Ds -----',esDs
       Write(lproj,'(a,30x,f15.6)') '             j_Vxs ----',ejNxs
       If(use_j2terms .Eqv. .True.) Then
          Write(lproj,'(a,30x,f15.6)') '             s_T ------',esT
       End If
       Write(lproj,'(a,30x,f15.6)') '  Projected Energy ....',eproj
    End Do

    ! Timing
    Call system_clock(tfinish,clock_rate)
    Write(lproj,'(a)')     '  ============================================================================'
    Write(lproj,'(a)')     '  =======================      System clock times      ======================='
    Write(lproj,'(a)')     '  ============================================================================'
    Write(lproj,'(43x,a,7x,a,8x,a)') 'seconds','minutes','hours'
    Write(lproj,'(43x,a,7x,a,7x,a)') '-------','-------','-------'
    Write(lproj,'(a37,1x,f12.2,2x,f12.2,2x,f12.2)') '  Total clock time:                  ', Real(tfinish-tstart)/ Real(clock_rate),Real(tfinish-tstart)/(Real(clock_rate)*60_pr),&
                                                                                             Real(tfinish-tstart)/(Real(clock_rate)*60_pr*60_pr)
    Write(lproj,'(43x,a)') '-----------------------------------'
    Write(lproj,'(a37,1x,f12.2,2x,f12.2,2x,f12.2)') '  Time in initialization routines:   ',time_ibet(1),time_ibet(1)/60_pr,time_ibet(1)/(60_pr*60_pr)
    Write(lproj,'(a37,1x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_rotation_matrix: ',time_ibet(2),time_ibet(2)/60_pr,time_ibet(2)/(60_pr*60_pr)
    Write(lproj,'(a37,1x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_overlaps:        ',time_ibet(3),time_ibet(3)/60_pr,time_ibet(3)/(60_pr*60_pr)
    Write(lproj,'(a37,1x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_densities:       ',time_ibet(4),time_ibet(4)/60_pr,time_ibet(4)/(60_pr*60_pr)
    Write(lproj,'(a37,1x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_energies:        ',time_ibet(5),time_ibet(5)/60_pr,time_ibet(5)/(60_pr*60_pr)
    Write(lproj,'(a37,1x,f12.2,2x,f12.2,2x,f12.2)') '  Time in project:                   ',time_ibet(6),time_ibet(6)/60_pr,time_ibet(6)/(60_pr*60_pr)

  End Subroutine print_project
  !==========================================================
  !> This routine finalizes symmetry restoration procedure.
  !==========================================================
  Subroutine finalize_projections()
    Use HFBTHO, Only : ro_normalization

    Close(lproj)
    Deallocate(betabs,betaweight)
    Deallocate(phiabs,ephi,ephic,ephicN,ephicP)
    Deallocate(rotated_overlap,detR)
    Deallocate(projected_overlap,projected_ekinN,projected_ekinP,projected_ecodi,projected_ecoex,projected_EVOL_rho_tau,&
               projected_EVOL_rho_rho,projected_ESURF_rho_drho,projected_ESO_rho_nablaj,projected_eptN,projected_eptP,&
               projected_ess0,projected_essD,projected_ecurr,projected_esDs,projected_ejNxs,&
               projected_ETENS,projected_esT,projected_xn1,projected_xn2,projected_rms1,projected_rms2,&
               projected_delN,projected_delP,projected_NP,projected_NP_norm)
    Deallocate(nz_sim,nr_sim,nl_sim)
    Deallocate(phicyl)
    Deallocate(ihil_convert,iphicyl_convert,ihil_iphicyl_convert)
    Deallocate(iosc1_pair,iosc2_pair)
    Deallocate(itiphi_pair1,itiphi_pair2)
    Deallocate(ro_normalization,xl_ihil,xh_ihil)
    Deallocate(VmatrixN1,UmatrixN1,VmatrixN2,UmatrixN2,VmatrixP1,UmatrixP1,VmatrixP2,UmatrixP2)
    Deallocate(all_overlaps,all_energies,beta_active)
    Deallocate(offset,block_vec)
#if(USE_MPI>0)
    Deallocate(beta_start)
#endif
    Deallocate(all_overlaps_gthr,all_energies_gthr,beta_active_gthr)

  End Subroutine finalize_projections
  !===========================================================================
  !> This routine calculates the states \f$ |n_z n_{\perp} \Lambda \rangle \f$ used to represent
  !> different operators in the basis of eigenstates of the simplex-y operator.
  !===========================================================================
  Subroutine simplex_basis(lpr)
    Use HFBTHO, Only : nt,nz,nr,ns,nl

    Logical, Intent(In)  :: lpr

    Integer(ipr) :: iosc

    Do iosc=1,nt
       nz_sim(iosc) = nz(iosc)
       nr_sim(iosc) = nr(iosc)
       If(ns(iosc) == 1) Then
          nl_sim(iosc) =  nl(iosc)
       Else
          nl_sim(iosc) = -nl(iosc)
       End If
    End Do

    If(lpr) Then
       Do iosc=1,nt
          If(iosc == 1) Then
             Write(lproj,'(a)')     '  ================================================================'
             Write(lproj,'(a)')     '  ======= Correspondence between HO and simplex-y bases  ========='
             Write(lproj,'(a)')     '  ================================================================'
             Write(lproj,'(2x,a,9x,a8,9x,a,9x,a15,5x,a)')   '|', 'HO basis', '|','Simplex-y basis','|'
             Write(lproj,'(2x,a)')   '----------------------------------------------------------'
             Write(lproj,'(4x,a,4x,a2,3x,a2,3x,a2,3x,a2,12x,a,3x,a2,3x,a2,3x,a2)') 'N','nz', &
                          'nr','nl','ns','N','nz','nr','nl'
          End If
          Write(lproj,'(i5,3x,i3,2x,i3,2x,i3,2x,i3,8x,i5,2x,i3,2x,i3,2x,i3)') iosc,nz(iosc),nr(iosc),nl(iosc), &
                       ns(iosc),iosc,nz_sim(iosc),nr_sim(iosc),nl_sim(iosc)
       End Do
    End If

  End Subroutine simplex_basis
  !========================================================================================
  !> This routine  calculates the rotation matrix element \f$ \langle n_z n_r \Lambda | \hat{R}(\beta) |
  !> n'_z n'_r \Lambda' \rangle \f$ between two spin-less HO states, based on Eq.(47) in Ref.
  !> \cite nazmitdinov1996representation.
  !>  @param[in] bb - Value of the angle \f$\beta \f$
  !>  @param[in] bp1 - Radial oscillator length \f$ b_{\perp} \f$ for the bra
  !>  @param[in] bz1 - Longitudinal oscillator length \f$ b_{z} \f$ for the bra
  !>  @param[in] bp2 - Radial oscillator length \f$ b'_{\perp} \f$ for the ket
  !>  @param[in] bz2 - Longitudinal oscillator length \f$ b'_{z} \f$ for the ket
  !>  @param[in] nz1 - Quantum number \f$ n_z \f$ for the bra
  !>  @param[in] nr1 - Quantum number \f$ n_r \f$ for the bra
  !>  @param[in] nl1 - Quantum number \f$ \Lambda \f$ for the bra
  !>  @param[in] nz2 - Quantum number \f$ n'_z \f$ for the ket
  !>  @param[in] nr2 - Quantum number \f$ n'_r \f$ for the ket
  !>  @param[in] nl2 - Quantum number \f$ \Lambda' \f$ for the ket
  !>  @param[inout] rotel - Actual value fo the matrix element
  !========================================================================================
  Subroutine calculate_ry(bb,bp1,bz1,bp2,bz2,nz1,nr1,nl1,nz2,nr2,nl2,rotel)
    Use math, Only : calculate_inverse_real,iv,fak,fi

    Real(pr), Intent(In) :: bb,bp1,bz1,bp2,bz2
    Real(pr), Intent(Inout) :: rotel
    Integer(ipr), Intent(In) :: nz1,nr1,nl1,nz2,nr2,nl2

    Real(pr):: eta,mub,muc,omega,deltap,deltam,detm
    Real(pr), Allocatable:: M(:,:)
    Real(pr):: fmat11,fmat13,fmat22,fmat33,ffmat11,ffmat22,ffmat33,ffmat12,ffmat13,ffmat23,gmat11,gmat13,gmat22,gmat33, &
               ggmat11,ggmat22,ggmat33,ggmat12,ggmat13,ggmat23,kmat11,kmat13,kmat22,kmat31,kmat33,kkmat11,kkmat12,kkmat13, &
               kkmat21,kkmat22,kkmat23,kkmat31,kkmat32,kkmat33
    Integer(ipr) :: ifl,np1,nm1,nn1,np2,nm2,nn2,ip1,mm1,mm2,mm3,mm4,mm5,mm6,lim1,lim2,lim3,ii1,jj1,kk1,ll1
    Real(pr) :: p1,prefac,bfac,bsum1,bsum2,bsum3,bsum4,cfac,csum1,csum2,csum3,csum4,afac,asum1,asum2,asum3,asum4,asum5,asum6

    rotel=zero
    nn1=2*nr1+Abs(nl1)+nz1; nn2=2*nr2+Abs(nl2)+nz2
    If(nn1 /= nn2) Return

    ! Defining F, G, and K matrices
    eta=(bp2/bp1)**2; mub=(bp1/bz1)**2; muc=(bp2/bz2)**2; omega=mub*eta+muc/eta+(one+mub*muc)*(sin(bb))**2+(mub+muc)*(cos(bb))**2
    deltap=(one-mub*muc)*(sin(bb))**2+(mub-muc)*(cos(bb))**2; deltam=(one-mub*muc)*(sin(bb))**2-(mub-muc)*(cos(bb))**2

    fmat11=one/omega*(mub*eta+deltam-muc/eta); fmat13=-one/omega*sqrt(mub)*(muc-one)*sin(two*bb)
    fmat22=(eta-one)/(eta+one); fmat33=one/omega*(mub*eta-deltam-muc/eta)

    ffmat11=half*(fmat11-fmat22); ffmat22=half*(fmat11-fmat22); ffmat33=fmat33
    ffmat12=-half*(fmat11+fmat22); ffmat13=-sqrt(half)*fmat13; ffmat23=sqrt(half)*fmat13

    gmat11=-one/omega*(mub*eta-deltap-muc/eta); gmat13=one/omega*sqrt(muc)*(mub-one)*sin(two*bb)
    gmat22=-fmat22; gmat33=-one/omega*(mub*eta+deltap-muc/eta)

    ggmat11=half*(gmat11-gmat22); ggmat22=half*(gmat11-gmat22); ggmat33=gmat33
    ggmat12=-half*(gmat11+gmat22); ggmat13=-sqrt(half)*gmat13; ggmat23=sqrt(half)*gmat13

    kmat11=two/omega*(mub*sqrt(eta)+muc/sqrt(eta))*cos(bb); kmat13=-two/omega*sqrt(mub)*(muc/sqrt(eta)+sqrt(eta))*sin(bb)
    kmat22=two*sqrt(eta)/(eta+one); kmat31=two/omega*sqrt(muc)*(one/sqrt(eta)+mub*sqrt(eta))*sin(bb)
    kmat33=two/omega*sqrt(mub*muc)*(sqrt(eta)+one/sqrt(eta))*cos(bb)

    kkmat11=half*(kmat11+kmat22); kkmat12=half*(-kmat11+kmat22); kkmat13=-sqrt(half)*kmat13
    kkmat21=half*(-kmat11+kmat22); kkmat22=half*(kmat11+kmat22); kkmat23=sqrt(half)*kmat13
    kkmat31=-sqrt(half)*kmat31; kkmat32=sqrt(half)*kmat31; kkmat33=kmat33

    ! Calculating determinant of the M matrix
    Allocate(M(1:3,1:3))
    M=zero
    M(1,1)=half*(((bp1/bp2)*cos(bb))**2+((bp1/bz2)*sin(bb))**2+one)
    M(1,3)=half*(-(bz1/bp2)*(bp1/bp2)+(bz1/bz2)*(bp1/bz2))*cos(bb)*sin(bb)
    M(2,2)=half*((bp1/bp2)**2+one)
    M(3,1)=half*(-(bz1/bp2)*(bp1/bp2)+(bz1/bz2)*(bp1/bz2))*cos(bb)*sin(bb)
    M(3,3)=half*(((bz1/bz2)*cos(bb))**2+((bz1/bp2)*sin(bb))**2+one)
    Call calculate_inverse_real(3,M,detm,ifl)
    Deallocate(M)

    ! Establishing quantum numbers np1,nm1,np2,nm2,nn1,nn2
    If(nl1 >= 0) Then
       np1=nr1+nl1; nm1=nr1
    Else
       np1=nr1; nm1=nr1+Abs(nl1)
    End If

    If(nl2 >= 0) Then
       np2=nr2+nl2; nm2=nr2
    Else
       np2=nr2; nm2=nr2+Abs(nl2)
    End If

    ! Calculating prefactor
    ip1=iv(np1+np2+nr1+nr2)
    p1=(bp1/bp2)*Sqrt(bz1/bz2)/Sqrt(detm*two**(nn1+nn2))*Sqrt(fak(np1)*fak(nm1)*fak(nz1)*fak(np2)*fak(nm2)*fak(nz2))
    prefac=ip1*p1

    Do mm1=0,np1  ! sum over n1
       Do mm2=0,nm1  ! sum over n2
          Do mm3=0,nz1  ! sum over n3

             bfac=zero
             If(Mod(np1+nm1+nz1-mm1-mm2-mm3,2) == 0) Then ! N-N_tilde is even
                lim1=(np1+nm1+nz1-mm1-mm2-mm3)/2-np1+mm1
                lim2=(np1+nm1+nz1-mm1-mm2-mm3)/2-nm1+mm2
                lim3=(np1+nm1+nz1-mm1-mm2-mm3)/2-nz1+mm3
                Do ii1=0,(lim2+lim3)/2
                   Do jj1=0,(lim1+lim3)/2
                      Do kk1=0,(lim1+lim2)/2
                         If(lim3-ii1-jj1+kk1 >= 0 .And. lim2-ii1-kk1+jj1 >= 0 .And. lim1-jj1-kk1+ii1 >= 0) Then
                            bsum1=(ffmat11**ii1)*(ffmat22**jj1)*(ffmat33**kk1)*fi(ii1)*fi(jj1)*fi(kk1)*two**(-ii1-jj1-kk1)
                            bsum2=ffmat12**(lim3-ii1-jj1+kk1)*fi(lim3-ii1-jj1+kk1)
                            bsum3=ffmat13**(lim2-ii1-kk1+jj1)*fi(lim2-ii1-kk1+jj1)
                            bsum4=ffmat23**(lim1-jj1-kk1+ii1)*fi(lim1-jj1-kk1+ii1)
                            bfac=bfac+bsum1*bsum2*bsum3*bsum4
                         End If
                      End Do ! kk1
                   End Do ! jj1
                End Do ! ii1
                bfac=bfac*two**((np1+nm1+nz1-mm1-mm2-mm3)/2)
             End If ! Mod(np1+nm1+nz1-mm1-mm2-mm3,2) == 0

             Do mm4=0,np2  ! sum over n1'
                Do mm5=0,nm2  ! sum over n2'
                   Do mm6=0,nz2 ! sum over n3'
                      cfac=zero;afac=zero
                      If(Mod(np2+nm2+nz2-mm4-mm5-mm6,2) == 0) Then
                         If(mm1+mm2+mm3 == mm4+mm5+mm6) Then ! Condition on A, that N=N'. Otherwise, A component will be zero, and the entire term as well.
                            lim1=(np2+nm2+nz2-mm4-mm5-mm6)/2-np2+mm4
                            lim2=(np2+nm2+nz2-mm4-mm5-mm6)/2-nm2+mm5
                            lim3=(np2+nm2+nz2-mm4-mm5-mm6)/2-nz2+mm6
                            Do ii1=0,(lim2+lim3)/2
                               Do jj1=0,(lim1+lim3)/2
                                  Do kk1=0,(lim1+lim2)/2
                                     If(lim3-ii1-jj1+kk1 >= 0 .And. lim2-ii1-kk1+jj1 >= 0 .And. lim1-jj1-kk1+ii1 >= 0) Then
                                        csum1=(ggmat11**ii1)*(ggmat22**jj1)*(ggmat33**kk1)*fi(ii1)*fi(jj1)*fi(kk1)*two**(-ii1-jj1-kk1)
                                        csum2=ggmat12**(lim3-ii1-jj1+kk1)*fi(lim3-ii1-jj1+kk1)
                                        csum3=ggmat13**(lim2-ii1-kk1+jj1)*fi(lim2-ii1-kk1+jj1)
                                        csum4=ggmat23**(lim1-jj1-kk1+ii1)*fi(lim1-jj1-kk1+ii1)
                                        cfac=cfac+csum1*csum2*csum3*csum4
                                     End If
                                  End Do
                               End Do
                            End Do
                            cfac=cfac*two**((np2+nm2+nz2-mm4-mm5-mm6)/2)
                            Do ii1=0,mm1
                               Do jj1=0,mm2
                                  Do kk1=0,mm1
                                     Do ll1=0,mm2
                                        If(mm4-ii1-jj1 >= 0 .And. mm5-kk1-ll1 >= 0 .And. mm1-ii1-kk1 >= 0 .And. mm2-jj1-ll1 >= 0 &
                                                                                   .And. mm6-mm1-mm2+ii1+jj1+kk1+ll1 >= 0) Then
                                           asum1=(kkmat11**ii1)*(kkmat12**jj1)*(kkmat21**kk1)*(kkmat22**ll1)*fi(ii1)*fi(jj1)*fi(kk1)*fi(ll1)
                                           asum2=kkmat13**(mm4-ii1-jj1)*fi(mm4-ii1-jj1)
                                           asum3=kkmat23**(mm5-kk1-ll1)*fi(mm5-kk1-ll1)
                                           asum4=kkmat31**(mm1-ii1-kk1)*fi(mm1-ii1-kk1)
                                           asum5=kkmat32**(mm2-jj1-ll1)*fi(mm2-jj1-ll1)
                                           asum6=kkmat33**(mm6-mm1-mm2+ii1+jj1+kk1+ll1)*fi(mm6-mm1-mm2+ii1+jj1+kk1+ll1)
                                           afac=afac+asum1*asum2*asum3*asum4*asum5*asum6
                                        End If
                                     End Do ! ll1
                                  End Do ! kk1
                               End Do ! jj1
                            End Do ! ii1
                            afac=afac*two**(mm1+mm2+mm3)
                         End If ! N=N'
                      End If ! Mod(np2+nm2+nz2-mm4-mm5-mm6,2) == 0
                      !
                      rotel=rotel+bfac*afac*cfac
                   End Do ! mm6
                End Do ! mm5
             End Do ! mm4

          End Do ! mm3
       End Do ! mm2
    End Do ! mm1

    rotel=prefac*rotel

  End Subroutine calculate_ry
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_projections
