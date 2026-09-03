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
!                        QRPA-pnFAM INTERFACE PACKAGE                  !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides an interface between HFBTHO and the QRPA-pnFAM
!> code developed by the Chapel Hill group.
!>
!>  @author
!>    Mikka Mustonen, Thomas Shafer
!----------------------------------------------------------------------
!  Subroutines: - save_HFBTHO_solution(plaintext)
!----------------------------------------------------------------------
Module HFBTHO_storage

  Use HFBTHO_utilities

  Implicit None

  Public save_HFBTHO_solution

  Private check_and_open

  ! Version number to check against
  Integer(ipr), Parameter :: FAM_VERSION = 10

Contains
  !=======================================================================
  ! Saves the relevant HFBTHO solution information to file for re-use by FAM.
  !=======================================================================
  Subroutine save_HFBTHO_solution(plaintext)
    Use UNEDF, Only: cr0 => crho, crr => cdrho, ctau, cdrho => crdr, ctj => cj, crdj, use_j2terms, CpV0, CpV1
    Use HFBTHO_Gauss, Only : nghl,wdcori,fh
    Use HFBTHO, only: nb, id, REqpP, RUqpP, RVqpP, REqpN, RUqpN, RVqpN, nr, nz, nl, ns, npar,     &
                      y_opt, qhla_opt, fi1r_opt, fi1z_opt, fi2d_opt, npr, ala, alast, del, bet,   &
                      ehfb, pwi, fn_T, fp_T, temper, entropy, switch_on_temperature, nbx, ka, kd, &
                      KqpP, KpwiP, KqpN, KpwiN, keyblo, blok1k2d, ntx

    Logical, Intent(In), optional :: plaintext
    Logical :: save_as_text
    Integer(ipr) :: ifh, ierr, ib, is
    Integer(ipr) :: i, iexit
    Real(pr), Dimension(nghl, ntx) :: qhla_opt_t, fi1r_opt_t, fi1z_opt_t, fi2d_opt_t

    ! Option to save in plain text
    save_as_text = .False.
    If(present(plaintext)) Then
       If(plaintext .eqv. .True.) Then
          save_as_text = .True.
       End If
    End If

    ! Check and open file
    ifh = 77
    iexit = check_and_open(ifh,save_as_text)
    If(iexit /= 0) Then
       Write(*,'(a,i0)') ' Error saving HFB solution: could not open file to write. Error', ierr
       stop
    End If

    ! ------------------------------------------------------------------------
    ! Store the HFB solution details to 'solution.hfb(.txt if plain text)'
    ! ------------------------------------------------------------------------
    If(.not. save_as_text) Then
       ! The solution version number
       Write(ifh) FAM_VERSION

       ! Basic HFB quantities
       ! N, Z, A where applicable
       Write(ifh) npr(:)             ! Particle number
       Write(ifh) ala(:)             ! Lambdas ala (quasiparticles measured w.r.t. these, I think)
       Write(ifh) alast(:)           ! Lambdas alast (last-bound s.p. energy)
       Write(ifh) del(:)             ! Pairing gaps
       Write(ifh) pwi                ! Pairing window
       Write(ifh) CpV0(:), CpV1(:)   ! Pairing strengths
       Write(ifh) bet                ! Total deformation
       Write(ifh) ehfb               ! Binding energy

       ! The HFB quasiparticle energies and amplitudes
       Write(ifh) nbx
       Write(ifh) nb
       Write(ifh) id(:)
       Write(ifh) REqpP(:)
       Write(ifh) RVqpP(:)
       Write(ifh) RUqpP(:)
       Write(ifh) REqpN(:)
       Write(ifh) RVqpN(:)
       Write(ifh) RUqpN(:)

       ! Pairing window active q.p. levels
       Write(ifh) ka(:,:)
       Write(ifh) kd(:,:)
       Write(ifh) KqpP(:)
       Write(ifh) KpwiP(:)
       Write(ifh) KqpN(:)
       Write(ifh) KpwiN(:)

       ! Basis state quantum numbers
       Write(ifh) nr(:)
       Write(ifh) nz(:)
       Write(ifh) nl(:)
       Write(ifh) ns(:)
       Write(ifh) npar(:)

       ! Wave functions and integration data
       Write(ifh) nghl                      ! number of integration points
       Write(ifh) wdcori(:)                 ! inverse of integration weights
       Write(ifh) y_opt(:)                  ! 1/rho in 'fm^(-1)'
       Write(ifh) fh(:)                     ! z in 'fm'
       Write(ifh) transpose(qhla_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
       Write(ifh) transpose(fi1r_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
       Write(ifh) transpose(fi1z_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
       Write(ifh) transpose(fi2d_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)

       ! Time-even isovector coupling constants from HFB mean field
       Write(ifh) cr0(1), crr(1), cdrho(1), ctau(1), ctj(1), crdj(1)
       ! Store .TRUE. => 1, .FALSE. => 0
       Write(ifh) merge(1, 0, use_j2terms)

       ! Temperature-dependence of the HFB solution
       ! Store .TRUE. => 1, .FALSE. => 0
       Write(ifh) merge(1, 0, switch_on_temperature)
       Write(ifh) temper
       Write(ifh) entropy
       Write(ifh) fp_T(:)
       Write(ifh) fn_T(:)

       ! Neutron Blocking
       If(keyblo(1) == 0) Then
          Write(ifh) 0, 0, 0, 0, 0.0d0
       Else
          ! Write down the block and state
          Do ib=0, nb
             If(sum(id(1:ib)) >= KqpN(blok1k2d(1))) exit
          End Do
          is = KqpN(blok1k2d(1)) - sum(id(1:ib-1))
          Write(ifh) keyblo(1), KqpN(blok1k2d(1)), ib, is, REqpN(KqpN(blok1k2d(1)))
       End If
       ! Proton Blocking
       If(keyblo(2) == 0) Then
          Write(ifh) 0, 0, 0, 0, 0.0d0
       Else
          Do ib=0, nb
             If(sum(id(1:ib)) >= KqpP(blok1k2d(2))) exit
          End Do
          is = KqpP(blok1k2d(2)) - sum(id(1:ib-1))
          Write(ifh) keyblo(2), KqpP(blok1k2d(2)), ib, is, REqpP(KqpP(blok1k2d(2)))
       End If

    Else

       ! The solution version number
       Write(ifh,*) 'FAM_VERSION'
       Write(ifh,*) FAM_VERSION
       Write(ifh,*)

       ! Basic HFB quantities
       ! N, Z, A where applicable
       Write(ifh,*) 'npr(:)'
       Write(ifh,*) npr(:)     ! Particle number
       Write(ifh,*) 'ala(:)'
       Write(ifh,*) ala(:)     ! Lambdas ala (quasiparticles measured w.r.t. these, I think)
       Write(ifh,*) 'alast(:)'
       Write(ifh,*) alast(:)   ! Lambdas alast (last-bound s.p. energy)
       Write(ifh,*) 'del(:)'
       Write(ifh,*) del(:)     ! Pairing gaps
       Write(ifh,*) 'pwi'
       Write(ifh,*) pwi        ! Pairing window
       Write(ifh,*) 'CpV0(:)'
       Write(ifh,*) CpV0(:)    ! Pairing strength
       Write(ifh,*) 'CpV1(:)'
       Write(ifh,*) CpV1(:)    ! Pairing strength
       Write(ifh,*) 'bet'
       Write(ifh,*) bet        ! Total deformation
       Write(ifh,*) 'ehfb'
       Write(ifh,*) ehfb       ! Binding energy
       Write(ifh,*)

       ! The HFB quasiparticle energies and amplitudes
       Write(ifh,*) 'nbx'
       Write(ifh,*) nbx
       Write(ifh,*) 'nb'
       Write(ifh,*) nb
       Write(ifh,*) 'id(:)'
       Write(ifh,*) id(:)
       Write(ifh,'(2x,6(A,2x))') 'REqpP(:)', 'RVqpP(:)', 'RUqpP(:)', 'REqpN(:)', 'RVqpN(:)', 'RUqpN(:)'
       Do i=0,size(REqpP)
          Write(ifh,*) REqpP(i), RVqpP(i), RUqpP(i), REqpN(i), RVqpN(i), RUqpN(i)
       End Do
       Write(ifh,*)

       ! Pairing window active q.p. levels
       Write(ifh,*) 'ka(:,:)'
       Do i=1,2
          Write(ifh,*) ka(:,i)
       End Do
       Write(ifh,*) 'kd(:,:)'
       Do i=1,2
          Write(ifh,*) kd(:,i)
       End Do
       Write(ifh,'(2x,4(A,2x))')  'KqpP(:)', 'KpwiP(:)', 'KqpN(:)', 'KpwiN(:)'
       Do i=1,size(KqpP)
           Write(ifh,*) KqpP(i), KpwiP(i), KqpN(i), KpwiN(i)
       End Do
       Write(ifh,*)

       ! Basis state quantum numbers
       Write(ifh,'(2x,5(A,2x))') 'nr(:)', 'nz(:)', 'nl(:)', 'ns(:)', 'npar(:)'
       Do i=1,size(nr)
           Write(ifh,*) nr(i), nz(i), nl(i), ns(i), npar(i)
       End Do
       Write(ifh,*)

       ! Wave functions and integration data
       Write(ifh,*) 'nghl'
       Write(ifh,*) nghl  ! number of integration points
       Write(ifh,'(2x,3(A,2x))') 'wdcori(:)', 'y_opt(:)', 'fh(:)'
       Do i=1,nghl
          Write(ifh,*) wdcori(i), y_opt(i), fh(i)  ! (inverse of integration weights), (1/rho in 'fm^(-1)'), (z in 'fm')
       End Do
       Write(ifh,*) 'trans(qhla_opt(:,:))'
       qhla_opt_t = transpose(qhla_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
       Do i=1,nghl
          Write(ifh,*) qhla_opt_t(i,:)
       End Do
       Write(ifh,*) 'trans(fi1r_opt(:,:))'
       fi1r_opt_t = transpose(fi1r_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
       Do i=1,nghl
          Write(ifh,*) fi1r_opt_t(i,:)
       End Do
       Write(ifh,*) 'trans(fi1z_opt(:,:))'
       fi1z_opt_t = transpose(fi1z_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
       Do i=1,nghl
          Write(ifh,*) fi1z_opt_t(i,:)
       End Do
       Write(ifh,*) 'trans(fi2d_opt(:,:))'
       fi2d_opt_t = transpose(fi2d_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
       Do i=1,nghl
          Write(ifh,*) fi2d_opt_t(i,:)
       End Do
       Write(ifh,*)

       ! Time-even isovector coupling constants from HFB mean field
       Write(ifh,'(2x,6(A,2x))') 'cr0(1)', 'crr(1)', 'cdrho(1)', 'ctau(1)', 'ctj(1)', 'crdj(1)'
       Write(ifh,*) cr0(1), crr(1), cdrho(1), ctau(1), ctj(1), crdj(1)
       ! Store .TRUE. => 1, .FALSE. => 0
       Write(ifh,*) 'Use_J^2_terms?'
       Write(ifh,*) merge(1, 0, use_j2terms)
       Write(ifh,*)

       ! Temperature-dependence of the HFB solution
       ! Store .TRUE. => 1, .FALSE. => 0
       Write(ifh,*) 'temp_on?'
       Write(ifh,*) merge(1, 0, switch_on_temperature)
       Write(ifh,*) 'temper'
       Write(ifh,*) temper
       Write(ifh,*) 'entropy'
       Write(ifh,*) entropy
       Write(ifh,'(2x,2(A,2x))') 'fp_T(:)', 'fn_T(:)'
       Do i=1,size(fp_T)
          Write(ifh,*) fp_T(i), fn_T(i)
       End Do
       Write(ifh,*)

       ! Neutron Blocking
       Write(ifh,*) 'Neutron_Blocking'
       Write(ifh,'(2x,5(A,2x))') 'keyblo(1)', 'KqpN', 'ib', 'is', 'REqpN'
       If(keyblo(1) == 0) Then
          Write(ifh,*) 0, 0, 0, 0, 0.0d0
       Else
          ! Write down the block and state
          Do ib=0, nb
             If(sum(id(1:ib)) >= KqpN(blok1k2d(1))) exit
          End Do
          is = KqpN(blok1k2d(1)) - sum(id(1:ib-1))
          Write(ifh, *) keyblo(1), KqpN(blok1k2d(1)), ib, is, REqpN(KqpN(blok1k2d(1)))
       End If
       ! Proton Blocking
       Write(ifh,*) 'Proton_Blocking'
       Write(ifh,'(2x,5(A,2x))') 'keyblo(2)', 'KqpP', 'ib', 'is', 'REqpP'
       If(keyblo(2) == 0) Then
          Write(ifh, *) 0, 0, 0, 0, 0.0d0
       Else
          Do ib=0, nb
             If(sum(id(1:ib)) >= KqpP(blok1k2d(2))) exit
          End Do
          is = KqpP(blok1k2d(2)) - sum(id(1:ib-1))
          Write(ifh, *) keyblo(2), KqpP(blok1k2d(2)), ib, is, REqpP(KqpP(blok1k2d(2)))
       End If
    End If

    Close(ifh)
    Write(*,'(a)') ' Storage completed.'

  End Subroutine save_HFBTHO_solution
  !=======================================================================
  ! This function checks status of solution.hfb file (mirroring 'check_file'
  ! in hfbtho_io). If it exists and can be opened, it is opened and the
  ! function returns exit status 0. If it does not exit, a new file is opened
  ! the function returns exit status 0. Otherwise exit status is 1.
  !=======================================================================
  Integer(ipr) Function check_and_open(fileunit,save_as_text)
    Logical, Intent(In) :: save_as_text
    Integer(ipr), Intent(In) :: fileunit

    Logical :: file_exists,file_opened
    Integer(ipr) :: ierr, iexit

    iexit = 0
    If(save_as_text) Then
       file_exists=.False.; inquire(file='solution.hfb.txt', exist=file_exists); ierr=0
       Write(*,'(/a)') ' ### STORING HFB SOLUTION (PLAINTEXT format)'
       Write(*,'(a)')  ' Filename: "solution.hfb.txt"'
       If(file_exists) Then
          file_opened=.False.; inquire(unit=fileunit, opened=file_opened)
          If(file_opened) Then
             Close(fileunit)
          End If
          open(unit=fileunit, file='solution.hfb.txt', status='old', iostat=ierr)
       Else
          open(unit=fileunit, file='solution.hfb.txt', status='new', iostat=ierr)
       End If
    Else
       file_exists=.False.; inquire(file='solution.hfb', exist=file_exists); ierr=0
       Write(*,'(/a)') ' ### STORING HFB SOLUTION (BINARY format)'
       Write(*,'(a)')  ' Filename: "solution.hfb"'
       If(file_exists) Then
          file_opened=.False.; inquire(unit=fileunit, opened=file_opened)
          If(file_opened) Then
             Close(fileunit)
          End If
          open(unit=fileunit, file='solution.hfb', status='old', form='unformatted', iostat=ierr)
       Else
          open(unit=fileunit, file='solution.hfb', status='new', form='unformatted', iostat=ierr)
       End If
    End If
    If(ierr /= 0) Then
       iexit = 1
    End If
    check_and_open = iexit

  End Function check_and_open
  !=======================================================================
  !
  !=======================================================================
end module HFBTHO_storage
