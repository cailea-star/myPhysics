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
!                      INPUT/OUTPUT PACKAGE                            !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module contains all routines dealing with input (reading
!> data from disk) and output (writing data on disk).
!>
!> @authors
!> Nicolas Schunck, Evan Ney
!----------------------------------------------------------------------
!  Subroutines: - set_ASCII_filenames(toggle_output,filename_output,filename_dat)
!               - set_BINARY_filenames(filename_binary)
!               - inout(is,iexit,filename_binary)
!               - write_version()
!               - read_data_old(iexit)
!               - read_data(iexit)
!               - write_data(is)
!               - write_data_old()
!               - blosort(it,n)
!               - set_grids()
!               - cartesian_grid()
!               - test_pairing(Nqp)
!               - test_normality_bogo(iwave,jwave,summ)
!               - test_normality(iwave,jwave)
!               - basis_rspace()
!               - basis_functions()
!               - constraining_field_mesh()
!  Functions: - version_number()
!             - check_file()
!----------------------------------------------------------------------
Module HFBTHO_io

  Use HFBTHO_utilities

  Implicit None

  Public set_ASCII_filenames,set_BINARY_filenames,inout

  Private write_version,version_number,check_file,read_data_old,read_data,write_data,write_data_old,blosort

  Integer(ipr), Public, Save :: VERSION_DATA = 5 !< Version number of the output binary file.
                                                 ! VERSION_DATA=4 -> version without nqp, nuv,Z_NECK, qmoment
                                                 ! VERSION_DATA=3 -> version 3.00 of HFBTHO, without rk and ak
                                                 ! VERSION_DATA=2 -> version 2.00d of HFBTHO
                                                 ! VERSION_DATA=1 -> 1.66 (not implemented)
  Integer(ipr), Public, Save :: VERSION_READ = 5 !< Version number of the input binary file.
  Character(Len=256), Public, Save :: welfile !< Name of the binary file

  ! Characteristics of the basis read from file
  Integer(ipr), Public, Save :: Z_r,N_r,n00_r,nb_r,nt_r,nqp_r,nuv_r
  Integer(ipr), Public, Save :: ngh_r,ngl_r,nleg_r
  Integer(ipr), Allocatable, Public, Save :: nr_r(:),nz_r(:),nl_r(:),ns_r(:),ID_r(:)
  Real(pr), Public, Save :: b0_r,bz_r,bp_r,neckLag_r,neckValue_r,neckRequested_r,pwi_r,varmas_r,bet_r
  Real(pr), Public, Save :: del_r(2),ala_r(2),ept_r(3),ala2_r(2),varmasNZ_r(2)
  Real(pr), Dimension(0:1), Public, Save :: CpV0_r,CpV1_r
  Real(pr), Allocatable, Public, Save :: xh_r(:),xl_r(:),wh_r(:),wl_r(:)

  Integer(ipr), Private, Save :: debug_io = 0

Contains
  !=======================================================================
  !> This subroutine defines the names of the ASCII file containining the 'screen' output. This routine is called
  !> at the beginning of \ref hfbtho_dft_solver(). The user has the option to either use the preset names hfbtho.out
  !> and thoout.dat (possibly modified to contain the job number), or to pass specific names as arguments of the
  !> routine.
  !>   - USE_MPI=0: no MPI, default file names are hfbtho.out and thoout.dat
  !>   - USE_MPI=1: multi-core HFBTHO execution and/or additional parallelism from calling application, e.g. Python
  !>                script. Execution is supposed to be in a separate directory, e.g., run_000001/, hence filenames
  !>                don't need to carry the task number. Only task leader writes.
  !>   - USE_MPI=2: team-based parallelism in Fortran, all filenames are written in the same directory, hence should carry
  !>                their task number. Onl task leader writes. Default names are hfbtho_000001.out and thoout_000001.dat.
  !=======================================================================
  Subroutine set_ASCII_filenames(toggle_output,filename_output,filename_dat)
#if(USE_MPI>0)
    Use HFBTHO, Only : do_print,team_rank,row_string
#else
    Use HFBTHO, Only : do_print
#endif

    Logical, Intent(In) :: toggle_output !< - Activates opening of output files. If 0, no data is written out
    Character(Len=256), Intent(In), Optional :: filename_output !< - User-defined name for the 'condensed' output file
    Character(Len=256), Intent(In), Optional :: filename_dat !< - User-defined name for the 'extended' output file

    Character(Len=256) :: output_name, dat_name

    ! Possibly overwrite the name of the file containing the basic output
    If(Present(filename_output)) Then
       output_name = Trim(filename_output)
    Else
#if(USE_MPI==2)
       output_name = 'hfbtho'//row_string//'.out'
#else
       output_name = 'hfbtho.out'
#endif
    End If

    ! Possibly overwrite the name of the file containing the detailed output (with qp information)
    If(Present(filename_dat)) Then
       dat_name = Trim(filename_dat)
    Else
#if(USE_MPI==2)
       dat_name = 'thoout'//row_string//'.out'
#else
       dat_name = 'thoout.dat'
#endif
    End If

#if(USE_MPI==0)
    do_print = 1
    If(toggle_output) Then
       Open(lout,file=output_name,status='unknown')
       If(lout < lfile) Open(lfile,file=dat_name,status='unknown')
    End If
#else
    If(team_rank == 0) Then
       do_print = 1
       If(toggle_output) Then
          Open(lout,file=output_name,status='unknown')
          If(lout < lfile) Open(lfile,file=dat_name,status='unknown')
       End If
    Else
       do_print = 0
    End If
#endif

  End Subroutine set_ASCII_filenames
  !=======================================================================
  !> This subroutine defines the filenames for the binary file. This routine is called at the
  !> beginning of \ref inout().
  !=======================================================================
  Subroutine set_BINARY_filenames(filename_binary)
#if(USE_MPI==2)
    Use HFBTHO, Only : iLST1,row_string
#else
    Use HFBTHO, Only : iLST1
#endif

    Character(Len=256), Intent(In), Optional :: filename_binary !< - Name of the user-defined binary file

    ! Possibly overwrite the name of the file containing the basic output
    If(Present(filename_binary)) Then
       welfile = Trim(filename_binary)
    Else
#if(USE_MPI==2)
       If(iLST1 <= 0) Write(welfile,'("hfbtho_output",a7,a4)')  row_string,'.hel'
       If(iLST1 >  0) Write(welfile,'("hfbtho_output",a7,a4)')  row_string,'.tel'
#else
       If(iLST1 <= 0) Write(welfile,'("hfbtho_output.hel")')
       If(iLST1 >  0) Write(welfile,'("hfbtho_output.tel")')
#endif
    End If

  End Subroutine set_BINARY_filenames
  !=======================================================================
  !> This subroutine is the central interface to all input/output operations in HFBTHO. Depending on the
  !> value of its input parameter it will either read a binary file to get input data, or it will write
  !> the current data in a binary file. In the reading phase, the code tries to read a version number. If
  !> this read is successful, it is assumed the binary file is in version VERSION_DATA.
  !=======================================================================
  Subroutine inout(is,iexit,filename_binary)
#if(USE_MPI>0)
    Use mpi
    Use HFBTHO_mpi_communication, Only : pack, broadcast_binary_to_team,unpack
    Use HFBTHO, Only : inin,do_print,ierror_flag,iasswrong,HFB_cores,COMM_team,ierr_mpi
#else
    Use HFBTHO, Only : inin,do_print,ierror_flag,iasswrong
#endif

    Integer(ipr), Intent(In) :: is !< - Integer specifying if the file is read (=1) or written (=2,3)
    Integer(ipr), Intent(Inout) :: iexit !< - Integer giving the exit status of the routine (0: OK, >0: not OK)
    Character(Len=256), Intent(In), Optional :: filename_binary !< - Name of the binary file

    ! Local variables
    Integer(ipr) :: checked,checked_again,iw
#if(USE_MPI>0)
    Integer(ipr) :: n_sizes
#endif
    Character(Len=50) :: action

    ! Defining filename for the binary file
    Call set_BINARY_filenames(filename_binary)
    If(ierror_flag /= 0) Return
    !---------------------------------------------------------------------
    ! Read data to start the calculation
    !---------------------------------------------------------------------
    If(is == 1) Then
       action = 'Read'
       ! Start from scratch
       If(inin > 0) Then
          iexit=1
       Else
          If(do_print == 1) Then
             ! Check the file is valid
             checked = check_file(welfile, action, is)
             If(checked==0) Then
                ! Read the version number of the file
                VERSION_READ = version_number()
                ! If version is current, try to read the data
                If(VERSION_READ >= 3) Then
                   iexit=0
                   Call read_data(iexit)
                ! Try to read the file assuming it is in the old format
                Else
                   checked_again = check_file(welfile, action, is)
                   iexit=0
                   Call read_data_old(iexit)
                End If
             Else
                iexit=1
             End If
          End If
       End If
       ! If more than 1 MPI process is used per HFB calculation, only 1 core reads the data, and we need to
       ! broadcast it to all afterwards
#if(USE_MPI>0)
       If(HFB_cores > 1) Then
          Call mpi_bcast(iexit, 1, mpi_integer, 0, COMM_team, ierr_mpi)
          If(iexit == 0) Then
             n_sizes = 3
             If(do_print == 1) Call pack()
             Call broadcast_binary_to_team(n_sizes)
             If(do_print == 0) Call unpack()
          End If
       End If
#endif
    End If
    !---------------------------------------------------------------------
    ! Write data on disk
    !---------------------------------------------------------------------
    If(is >= 2 .And. iasswrong(3) == 0) Then
       If(do_print == 1) Then
          action = 'Write'
          ! Check the file is valid
          checked = check_file(welfile, action, is)
          If(checked==0) Then
             If(is == 2) Call write_version()
             If(is == 3 .And. VERSION_DATA == 2) Call write_data_old()
             If(VERSION_DATA == 5) Call write_data(is)
             iexit=0
          Else
             iexit=1
          End If
       End If
    End If

  End Subroutine inout
  !=======================================================================
  !> This subroutine just writes the version number on disk.
  !=======================================================================
  Subroutine write_version()
    Use HFBTHO, Only : lwou

    Write(lwou) VERSION_DATA

  End Subroutine write_version
  !=======================================================================
  !> This function reads and returns the version number from the input file. Returns -1 in case the
  !> version cannot be read.
  !=======================================================================
  Integer(ipr) Function version_number()
    Use HFBTHO, Only : lwin

    Integer(ipr) :: ierr,version

    version = 1
    Read(lwin,IOSTAT=ierr) version
    If(ierr /= 0) version = -1
    version_number = version

  End Function version_number
  !=======================================================================
  !> This function checks the status of the binary file. If the file exists and can be opened, it is
  !> opened and the functions returns an exit status of 0. If the file should be read but does not exist,
  !> the function returns 1; if the file should be written but does not exist, the function opens a new
  !> file and returns 0.
  !=======================================================================
  Integer(ipr) Function check_file(filename, action, is)
    Use HFBTHO, Only : lwin,lwou

    Character(Len=256), Intent(In) :: filename !< - Name of the file to check
    Character(Len=50), Intent(In) :: action !< - Equal to 'Read' or 'Write', defines what to do if the file does not exist
    Integer(ipr), Intent(In) :: is !< - 1=read, 2=write pre-hfbdiag quantities, 3=write post-hfbdiag quantities

    Logical :: file_exists,file_opened
    Integer(ipr) :: ierr,iexit

    iexit = 0
    file_exists=.False.; Inquire(file=filename, exist=file_exists); ierr=0
    If(Trim(action) == 'Read') Then
       If(file_exists) Then
          file_opened=.False.; Inquire(unit=lwin, opened=file_opened)
          If(file_opened) Close(lwin)
          Open(lwin,file=filename,status='old',form='unformatted',IOSTAT=ierr)
          If(ierr /= 0) iexit = 1
       Else
          iexit = 1
       End If
    End If
    If(Trim(action) == 'Write') Then
       If(file_exists) Then
          file_opened=.False.; Inquire(unit=lwou, opened=file_opened)
          If(file_opened) Close(lwou)
          If(is == 3) Then ! Append to existing file
             Open(lwou,file=filename,status='old',form='unformatted',position='append',IOSTAT=ierr)
          Else             ! Overwrite existing file
             Open(lwou,file=filename,status='old',form='unformatted',IOSTAT=ierr)
          End If
          If(ierr /= 0) iexit = 1
       Else
          If(is == 3) Then ! File should exist
             iexit = 1
          Else
             Open(lwou,file=filename,status='new',form='unformatted',IOSTAT=ierr)
             iexit = 0
          End If
       End If
    End If
    check_file = iexit

  End Function check_file
  !=======================================================================
  !> This subroutine reads the data for 'old' binary files corresponding to HFBTHO version 200d
  !=======================================================================
  Subroutine read_data_old(iexit)
    Use HFBTHO, Only : lwin,n00,nb,siold,etot,rms,xmix,ept,alast,tz,varmas,varmasNZ,pjmassNZ,ass,skass,brin,bbroyden,si, &
                       kindhfb,Add_Pairing,nhhdim2,nhhdim3,bloblo,blo123,blok1k2,blomax,bloqpdif,iLST,itass,iqqmax,do_print, &
                       bloall
    Use HFBTHO_THO, Only : decay,rmm3,amm3,bmm3,cmm3,fdsx,fdsy,fdsy1,fdsy2,fdsy3, &
                           fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
    Use HFBTHO_multipole_moments, Only : multLag

    Integer(ipr), Intent(Inout) :: iexit !< - Integer giving the exit status of the routine (0: OK, >0: not OK)

    Integer(ipr) :: iw,n1,n2,nd,ib,bloall1,lambdaMax1,counterLine
    Integer(ipr) :: npr1,npr11,n001,ibro,i,ibasis
    Integer(ipr) :: nhhdim1,NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
    Real(pr) :: beta1,v0r(2),v1r(2),pwir

    ! Read data
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,*) ' Reading from wel_file: ',Trim(welfile)
          Write(iw,*)
       End Do
    End If
    counterLine = 0; iexit = 0
    Read(lwin,Err=100,End=100) npr11,npr1,ngh_r,ngl_r,n001,nb_r,nt_r
    counterLine = counterLine+1
    If(Abs(n001) /= Abs(n00) .And. nb_r /= nb) go to 100
    Read(lwin,Err=100,End=100) b0_r,bz_r,bp_r,beta1,siold,etot,rms,bet_r,xmix,v0r,v1r,pwir, &
                               del_r,ept,ala_r,ala2_r,alast,tz,varmas,varmasNZ,pjmassNZ, &
                               ass,skass
    brin=zero; bbroyden='L'; si=siold
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) nt_r,nb_r,nhhdim1
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) lambdaMax1
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) multLag
    counterLine = counterLine+1
    If(Allocated(ID_r)) Deallocate(ID_r)
    Allocate(ID_r(1:nb_r))
    Read(lwin,Err=100,End=100) ID_r
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) brin
    counterLine = counterLine+1
    !
    ! Add small pairing de=de+0.1 in the no-LN case to prevent pairing collapse
    If(kindhfb == 1 .And. Add_Pairing) Then
       ibro=0
       Do ib=1,NB
          ND=ID_r(ib)
          I=ibro
          Do N1=1,ND
             Do N2=1,N1
                I=I+1
                brin(i+nhhdim2)=brin(i+nhhdim2)+0.10_pr
                brin(i+nhhdim3)=brin(i+nhhdim3)+0.10_pr
             End Do !N2
          End Do !N1
          ibro=i
       End Do !IB
    End If
    If(Allocated(nr_r)) Deallocate(nr_r,nz_r,nl_r,ns_r)
    Allocate(nr_r(1:nt_r),nz_r(1:nt_r),nl_r(1:nt_r),ns_r(1:nt_r))
    ibasis=0
    Do ib=1,NB
       ND=ID_r(ib)
       Do N1=1,ND
          ibasis=ibasis+1;
          Read(lwin,ERR=100,End=100) NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
          nr_r(ibasis)=(NZA2NRA-NLA1-NZA1)/2
          nz_r(ibasis)=NZA1
          nl_r(ibasis)=NLA1
          ns_r(ibasis)=NLANSA0-2*NLA1
       End Do
    End Do
    counterLine = counterLine+1
    ! blocking
    Read(lwin,ERR=100,End=100) bloall1
    counterLine = counterLine+1
    Read(lwin,ERR=100,End=100) bloblo,blo123,blok1k2,blomax,bloqpdif
    counterLine = counterLine+1
    If(bloall1 /= bloall) go to 100
    !tel
    If(iLST > 0) Then
       Read(lwin,ERR=100,End=100) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
       If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0, &
                                      fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3)
       Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),fdsy2(iqqmax),fdsy3(iqqmax), &
                fspb0(iqqmax),fspc0(iqqmax),fspd0(iqqmax),fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),  &
                fspb2(iqqmax),fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),fspc3(iqqmax),fspd3(iqqmax))
       Read(lwin,ERR=100,End=100) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0, &
                                  fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
    End If
    Close(lwin)
    Return

100 Continue
    iexit=1
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'(1x,a,a,a)')
          Write(iw,'(1x,a,a,a)')   ' The file ',Trim(welfile),' is corrupted!'
          Write(iw,'(1x,a,i2,a)')  ' Problem occurs at line ',counterLine,'        '
          Write(iw,'(1x,a,a,a)')   ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
          Write(iw,'(1x,a,a,a)')
       End Do
    End If

  End Subroutine read_data_old
  !=======================================================================
  !> This subroutine reads a binary file containing the results of a HFBTHO calculation. In the new
  !> format (VERSION_DATA>=3), the binary file is structured by keywords and contains the HF and
  !> pairing field on the Gauss quadrature mesh. This new format allows restarts even when (i) the
  !> basis is different (as long as the quadrature mesh is the same) (ii) the number of constraints
  !> is different (iii) metadata such as number of protons, neutrons, characteristics of the EDF, etc.
  !> are different. The price to pay for this flexibility is an increased size of the binary file.
  !=======================================================================
  Subroutine read_data(iexit)
    Use UNEDF, Only : use_INM,use_cm_cor,use_j2terms,force_is_dme,finite_range,hb0_charge_dependent, &
                      E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM, &
                      Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CpV0,CpV1,sigma,hbzero
    Use HFBTHO_Gauss, Only : ngh,ngl
    Use HFBTHO_fission_fragments, Only : Z_NECK
    Use HFBTHO_Gogny, Only : NumVz,NumVr,VrGogny,VzGogny
    Use HFBTHO_Lipkin, Only : ala2
    Use HFBTHO_THO, Only : decay,rmm3,amm3,bmm3,cmm3,fdsx,fdsy,fdsy1,fdsy2,fdsy3, &
                           fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
    Use HFBTHO_multipole_moments, Only : numberCons,multLambda,qmoment,multLag
    Use HFBTHO, Only : lwin,do_print,lambdaMax,ncons_max,pwi,hb0,hb0n,hb0p,b0,bp,bz,n00,nt,nb,si,etot,rms,ept,alast,pjmassNZ, &
                       siold,ala,del,tz,varmas,varmasNZ,qfield,neckLag,neckValue,ro,aka,rk,ak,nqp,fn_T,fp_T,ass,skass,bet,xmix, &
                       vn,vhbn,vrn,vzn,vdn,vsn,vSFIZn,vSZFIn,vSFIRn,vSRFIn,dvn,MEFFn, &
                       vp,vhbp,vrp,vzp,vdp,vsp,vSFIZp,vSZFIp,vSFIRp,vSRFIp,dvp,MEFFp, &
                       bloblo,blo123,blok1k2,blomax,bloqpdif,blocking_never_done,temper,entropy,geff_inv,pairing_regularization, &
                       brin,nhhdim,nhhdim2,nhhdim3,nhhdim4,functional,skyrme,bbroyden,itass,iqqmax,set_neck_constrain,record_HFB_matrix

    Integer(ipr), Intent(Inout) :: iexit !< - Integer giving the exit status of the routine (0: OK, >0: not OK)

    Character(Len=8) :: key
    Character(Len=30) :: skyrme_r
    Logical :: different_basis
    Integer(ipr) :: counterLine,iw,ib,nd,n1,icons,lambda,jcons,lambda_r,bloall_r
    Integer(ipr) :: switch_to_THO_r, projection_is_on_r
    Logical :: collective_inertia_r, fission_fragments_r, pairing_regularization_r, localization_functions_r, &
               set_temperature_r, set_neck_constrain_r
    Logical :: use_INM_r, use_cm_cor_r, use_j2terms_r, force_is_dme_r, finite_range_r, hb0_charge_dependent_r
    Integer(ipr) :: numberCons_r,lambdaMax_r,nhhdim_r
    Integer(ipr) :: ibasis
    Integer(ipr), Dimension(1:ncons_max) :: multLambda_r
    Real(pr) :: E_NM_r,K_NM_r,SMASS_NM_r,RHO_NM_r,ASS_NM_r,LASS_NM_r,VMASS_NM_r,P_NM_r,KA_NM_r,sigma_r
    Real(pr) :: hbzero_r,hb0_r,hb0n_r,hb0p_r,xmix_r
    Real(pr), Dimension(0:1) :: Crho_r,Cdrho_r,Ctau_r,CrDr_r,CrdJ_r,CJ_r
    Real(pr), Dimension(2) :: tz_r
    Real(pr), Dimension(0:lambdaMax) :: multRequested_r,multLag_r
    Real(pr), Dimension(0:lambdaMax,1:3) :: qmoment_r
    ! Arrays
    Real(pr), Allocatable, Target :: fn_T_r(:),fp_T_r(:)

    ! Read data
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,*) ' Reading from wel_file: ',Trim(welfile)
          Write(iw,*)
       End Do
    End If
    counterLine = 0; iexit = 0
    ! Loop over all keywords in the file
    Do
       Read(lwin,Err=100,End=99) key
       ! Metadata: N, Z, force, optional flags
       If(Trim(key)=='Metadata') Then
          Read(lwin,Err=100,End=100) Z_r,N_r
          Read(lwin,Err=100,End=100) collective_inertia_r, fission_fragments_r, pairing_regularization_r, localization_functions_r
          Read(lwin,Err=100,End=100) switch_to_THO_r, projection_is_on_r, set_temperature_r, set_neck_constrain_r
       End If
       ! EDF: all coupling constants, pairing cutoff, name
       If(Trim(key)=='SkyFunct') Then
          Read(lwin,Err=100,End=100) use_INM_r, use_cm_cor_r, use_j2terms_r, force_is_dme_r, finite_range_r, skyrme_r, &
                                     hb0_charge_dependent_r, pwi_r
          Read(lwin,Err=100,End=100) E_NM_r,K_NM_r,SMASS_NM_r,RHO_NM_r,ASS_NM_r,LASS_NM_r,VMASS_NM_r,P_NM_r,KA_NM_r
          Read(lwin,Err=100,End=100) Crho_r,Cdrho_r,Ctau_r,CrDr_r,CrdJ_r,CJ_r,CpV0_r,CpV1_r,sigma_r
          Read(lwin,Err=100,End=100) hbzero_r,hb0_r,hb0n_r,hb0p_r
          ! Overwrite (Skyrme) EDF characteristics if so requested
          If(Trim(functional) == 'READ') Then
             skyrme = skyrme_r
             ! Nuclear matter properties
             E_NM=E_NM_r; K_NM=K_NM_r; SMASS_NM=SMASS_NM_r; RHO_NM=RHO_NM_r; ASS_NM=ASS_NM_r
             LASS_NM=LASS_NM_r; VMASS_NM=VMASS_NM_r; P_NM=P_NM_r; KA_NM=KA_NM_r
             ! All Skyrme coupling constants and
             Crho=Crho_r; Cdrho=Cdrho_r; Ctau=Ctau_r; CrDr=CrDr_r; CrdJ=CrdJ_r; CJ=CJ_r; CpV0=CpV0_r; CpV1=CpV1_r
             sigma=sigma_r
             ! Pairing cut off
             pwi=pwi_r
             ! Logical flags that activate various terms
             use_INM=use_INM_r; use_cm_cor=use_cm_cor_r; use_j2terms=use_j2terms_r; force_is_dme=force_is_dme_r; finite_range=finite_range_r
             pairing_regularization = pairing_regularization_r
             ! Kinetic energy h^2/2m factors
             hb0_charge_dependent = hb0_charge_dependent_r; hbzero=hbzero_r;  hb0=hb0_r; hb0n=hb0n_r; hb0p=hb0p_r
          End If
       End If
       ! Basis: deformations, number of quanta, integration mesh
       If(Trim(key)=='HO-Basis') Then
          Read(lwin,Err=100,End=100) b0_r,bz_r,bp_r
          Read(lwin,Err=100,End=100) n00_r,nb_r,nt_r,ngh_r,ngl_r,nleg_r
          If(ngh_r /= ngh .Or. ngl_r /= ngl .And. do_print == 1) Then
             Do iw=lout,lfile
                Write(iw,'("Error in read_data - Inconsistent quadrature mesh!")')
                Write(iw,'("The code will start from scratch")')
             End Do
             iexit=1
          End If
          If(Allocated(xh_r)) Deallocate(xh_r,xl_r,wh_r,wl_r)
          Allocate(xh_r(1:ngh_r),xl_r(1:ngl_r),wh_r(1:ngh_r),wl_r(1:ngl_r))
          Read(lwin,Err=100,End=100) xh_r,xl_r,wh_r,wl_r
          If(Trim(functional) == 'READ') Then
             b0=b0_r; bp=bp_r; bz=bz_r
             n00 = n00_r; nb = nb_r; nt = nt_r
          End If
       End If
       ! Quantum numbers and \Omega-blocks information
       If(Trim(key)=='QuantNum') Then
          If(Allocated(ID_r)) Deallocate(ID_r)
          Allocate(ID_r(1:nb_r))
          Read(lwin,Err=100,End=100) ID_r
          If(Allocated(nr_r)) Deallocate(nr_r,nz_r,nl_r,ns_r)
          Allocate(nr_r(1:nt_r),nz_r(1:nt_r),nl_r(1:nt_r),ns_r(1:nt_r))
          ibasis=0
          Do ib=1,nb_r
             nd=ID_r(ib)
             Do n1=1,nd
                ibasis=ibasis+1;
                Read(lwin,Err=100,End=100) nr_r(ibasis),nz_r(ibasis),nl_r(ibasis),ns_r(ibasis)
             End Do
          End Do
          ! Check if the file has a different basis
          different_basis = nb /= nb_r .Or. Abs(b0_r-b0) > 1.e-6_pr .Or. Abs(bz_r-bz) > 1.e-6_pr &
                                       .Or. Abs(bp_r-bp) > 1.e-6_pr .Or. n00_r /= n00 .Or. nt_r /= nt
          If(different_basis .And. Trim(functional) /= 'READ' .And. do_print == 1) Then
             Do iw=lout,lfile
                Write(iw,'("nb  =",i4," nt  =",i4," n00  =",i4)') nb,nt,n00
                Write(iw,'("nb_r=",i4," nt_r=",i4," n00_r=",i4)') nb_r,nt_r,n00_r
                Write(iw,'("Warning in read_data - Inconsistent basis!")')
             End Do
          End If
       End If
       ! Various
       If(Trim(key)=='Various.') Then
          Read(lwin,Err=100,End=100) si,etot,rms,bet_r,xmix_r
          Read(lwin,Err=100,End=100) pwi_r,del_r,ept_r,ala_r,ala2_r,alast
          Read(lwin,Err=100,End=100) tz_r,varmas_r,varmasNZ_r,pjmassNZ,ass,skass
          siold=si; bet=bet_r; xmix=xmix_r
          pwi=pwi_r; del=del_r; ept=ept_r; ala=ala_r; ala2=ala2_r
          tz=tz_r; varmas=varmas_r; varmasNZ=varmasNZ_r
       End If
       ! Constraints: requested values, Lagrange parameters
       If(Trim(key)=='Constrai') Then
          Read(lwin,Err=100,End=100) numberCons_r,lambdaMax_r
          multLambda_r=0; multRequested_r=0.0_pr; qmoment_r=0.0_pr; multLag_r=0.0_pr
          If(VERSION_READ > 3) Then
             Read(lwin,Err=100,End=100) multLambda_r
             Read(lwin,Err=100,End=100) multRequested_r
             If(VERSION_READ > 4) Read(lwin,Err=100,End=100) qmoment_r
             Read(lwin,Err=100,End=100) multLag_r
          Else
             Read(lwin,Err=100,End=100) (multLambda_r(icons),icons=1,numberCons_r)
             Read(lwin,Err=100,End=100) (multRequested_r(lambda),lambda=0,lambdaMax_r)
             Read(lwin,Err=100,End=100) (multLag_r(lambda),lambda=1,lambdaMax_r)
          End If
          If(set_neck_constrain_r) Then
             Read(lwin,Err=100,End=100) neckRequested_r
             If(VERSION_READ > 4) Then
                Read(lwin,Err=100,End=100) neckValue_r,Z_NECK
             Else
                neckValue_r = 5.0_pr; Z_NECK = 0.0_pr
             End If
             Read(lwin,Err=100,End=100) neckLag_r
          End If
          Read(lwin,Err=100,End=100) qfield
          ! Reset values of Lagrange parameters based on the values read on disk
          qmoment = 0.0_pr
          Do icons=1,numberCons
             lambda=multLambda(icons)
             Do jcons=1,numberCons_r
                lambda_r=multLambda_r(jcons)
                If(lambda==lambda_r) Then
                   multLag(lambda) = multLag_r(lambda_r)
                   If(VERSION_READ > 4) qmoment(lambda,1:3) = qmoment_r(lambda,1:3)
                End If
             End Do
          End Do
          If(set_neck_constrain .And. set_neck_constrain_r) Then
             neckLag   = neckLag_r
             neckValue = neckValue_r
          End If
       End If
       ! Density matrix and pairing tensor in coordinate space
       If(Trim(key)=='Densits.') Then
          Read(lwin,Err=100,End=100) ro     ! 2*rho
          Read(lwin,Err=100,End=100) aka    ! Kappa
       End If
       ! HF and pairing field in coordinate space
       If(Trim(key)=='FieldsN.') Then
          Read(lwin,Err=100,End=100) vn     ! RHO_ij
          Read(lwin,Err=100,End=100) vhbn   ! TAU_ij
          Read(lwin,Err=100,End=100) vrn    ! NABLAr RHO__ij
          Read(lwin,Err=100,End=100) vzn    ! NABLAz RHO__ij
          Read(lwin,Err=100,End=100) vdn    ! DELTA RHO_ij
          Read(lwin,Err=100,End=100) vsn    ! NABLA . J__ij
          Read(lwin,Err=100,End=100) vSFIZn ! JFIZ_ij
          Read(lwin,Err=100,End=100) vSZFIn ! JZFI_ij
          Read(lwin,Err=100,End=100) vSFIRn ! JFIR_ij
          Read(lwin,Err=100,End=100) vSRFIn ! JRFI_ij
          Read(lwin,Err=100,End=100) dvn    ! \Delta_ij
       End If
       If(Trim(key)=='FieldsP.') Then
          Read(lwin,Err=100,End=100) vp
          Read(lwin,Err=100,End=100) vhbp
          Read(lwin,Err=100,End=100) vrp
          Read(lwin,Err=100,End=100) vzp
          Read(lwin,Err=100,End=100) vdp
          Read(lwin,Err=100,End=100) vsp
          Read(lwin,Err=100,End=100) vSFIZp
          Read(lwin,Err=100,End=100) vSZFIp
          Read(lwin,Err=100,End=100) vSFIRp
          Read(lwin,Err=100,End=100) vSRFIp
          Read(lwin,Err=100,End=100) dvp
       End If
       ! Blocking
       If(Trim(key)=='Blocking') Then
          Read(lwin,Err=100,End=100) bloall_r
          Read(lwin,Err=100,End=100) bloblo,blo123,blok1k2,blomax,bloqpdif
       End If
       ! Blocking
       If(Trim(key)=='Blk-rest') Then
          Read(lwin,Err=100,End=100) blocking_never_done
       End If
       ! Temperature
       If(Trim(key)=='Temperat') Then
          Read(lwin,Err=100,End=100) temper,entropy
          If(VERSION_READ > 4) Then
             Read(lwin,Err=100,End=100) nuv_r,nqp_r
          Else
             nqp_r=nqp
          End If
          If(Allocated(fn_T_r)) Deallocate(fn_T_r,fp_T_r)
          Allocate(fn_T_r(nqp_r),fp_T_r(nqp_r))
          Read(lwin,Err=100,End=100) fp_T_r
          Read(lwin,Err=100,End=100) fn_T_r
          If(nqp_r == nqp) Then
             fn_T = fn_T_r; fp_T = fp_T_r
          End If
       End If
       ! THO
       If(Trim(key)=='THObasis') Then
          If(Allocated(fdsx)) Then
             Read(lwin,Err=100,End=100) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
             Read(lwin,Err=100,End=100) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0,  &
                                        fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
          End If
       End If
       ! Pairing regularization
       If(Trim(key)=='Regular.') Then
          Read(lwin,Err=100,End=100) MEFFn
          Read(lwin,Err=100,End=100) MEFFp
          Read(lwin,Err=100,End=100) geff_inv
       End If
       ! Gogny force
       If(Trim(key)=='GognyVNN') Then
          ! Finite-range basis is implemented in configuration space and the
          ! basis on file must be the same as the current one for smooth restart
          If(different_basis) iexit=1
          Read(lwin,Err=100,End=100) NumVz,NumVr
          If(VERSION_READ > 3) Then
             Read(lwin,Err=100,End=100) rk
             Read(lwin,Err=100,End=100) ak
          End If
          Read(lwin,Err=100,End=100) vrGogny
          Read(lwin,Err=100,End=100) vzGogny
       End If
       ! HFB matrix
       If(Trim(key)=='HFBmatrX') Then
          Read(lwin,Err=100,End=100) nhhdim_r
          If(Allocated(brin)) Deallocate(brin)
          If(set_neck_constrain_r) Then
             Allocate(brin(4*nhhdim_r+lambdaMax+2))
          Else
             Allocate(brin(4*nhhdim_r+lambdaMax))
          End If
          brin = 0.0_pr
          Read(lwin,Err=100,End=100) brin
          If(nhhdim_r == nhhdim .And. .Not. different_basis) Then
             record_HFB_matrix = .True.; bbroyden='L'
             nhhdim=nhhdim_r; nhhdim2=2*nhhdim; nhhdim3=3*nhhdim; nhhdim4=4*nhhdim
          End If
       End If
    End Do

 99 Continue
    Close(lwin)
    Return

100 Continue
    iexit=1
    Close(lwin)
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'(1x,a,a,a)')
          Write(iw,'(1x,a,a,a)')   ' The file ',Trim(welfile),' is corrupted!'
          Write(iw,'(1x,a,a8,a)')  ' Problem occurs for key ',Trim(key),'        '
          Write(iw,'(1x,a,a,a)')   ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
          Write(iw,'(1x,a,a,a)')
       End Do
    End If

  End Subroutine read_data
  !=======================================================================
  !> This subroutine writes a binary file containing the results of a HFBTHO calculation. In the new
  !> format (VERSION_DATA>=3), the binary file is structured by keywords and contains the HF and pairing
  !> field on the Gauss quadrature mesh.
  !=======================================================================
  Subroutine write_data(is)
    Use UNEDF, Only : use_INM,use_cm_cor,use_j2terms,force_is_dme,finite_range,coulomb_gaussian,hb0_charge_dependent, &
                      E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM, &
                      Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CpV0,CpV1,sigma,hbzero
    Use HFBTHO_Gauss, Only : ngh,ngl,nleg,xh,xl,wh,wl
    Use HFBTHO_fission_fragments, Only : Z_NECK
    Use HFBTHO_Gogny, Only : NumVz,NumVr,VrGogny,VzGogny
    Use HFBTHO_Lipkin, Only : ala2
    Use HFBTHO_THO, Only : decay,rmm3,amm3,bmm3,cmm3,fdsx,fdsy,fdsy1,fdsy2,fdsy3, &
                           fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
    Use HFBTHO_multipole_moments, Only : numberCons,multLambda,multRequested,qmoment,multLag
    Use HFBTHO, Only : lwou,npr,collective_inertia,fission_fragments,pairing_regularization,localization_functions, &
                       switch_to_THO,projection_is_on,set_temperature,set_neck_constrain,full_HFB_matrix,skyrme,ro,aka,rk,ak, &
                       hb0,hb0n,hb0p,b0,bp,bz,n00,nb,nt,ID,nb,NL,NR,NZ,NS,si,etot,rms,bet,xmix,pwi,del,ept,ala,alast, &
                       tz,varmas,varmasNZ,pjmassNZ,ass,skass,lambdaMax,neckRequested,neckValue,neckLag,qfield,iLST,iqqmax,itass, &
                       vn,vhbn,vrn,vzn,vdn,vsn,vSFIZn,vSZFIn,vSFIRn,vSRFIn,dvn,MEFFn, &
                       vp,vhbp,vrp,vzp,vdp,vsp,vSFIZp,vSZFIp,vSFIRp,vSRFIp,dvp,MEFFp, &
                       bloall,bloblo,blo123,blok1k2,blomax,bloqpdif,blocking_never_done,temper,entropy,geff_inv, &
                       brin,nhhdim,nuv,nqp,fn_T,fp_T

    Integer(ipr), Intent(In) :: is !< - 1=read, 2=write pre-hfbdiag quantities, 3=write post-hfbdiag quantities

    Integer(ipr) :: it,ib,nd,n1,ibasis,nla,nra,nza,nsa

    ! is=2: Data from before diagonalization
    If(is == 2) Then
       ! Metadata: N, Z, force, optional flags
       Write(lwou) 'Metadata'
       Write(lwou) npr(2),npr(1)
       Write(lwou) collective_inertia, fission_fragments, pairing_regularization, localization_functions
       Write(lwou) switch_to_THO, projection_is_on, set_temperature, set_neck_constrain
       ! EDF: all coupling constants, pairing cutoff, name
       Write(lwou) 'SkyFunct'
       Write(lwou) use_INM, use_cm_cor, use_j2terms, force_is_dme, finite_range, skyrme, hb0_charge_dependent, pwi
       Write(lwou) E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM
       Write(lwou) Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CpV0,CpV1,sigma
       Write(lwou) hbzero,hb0,hb0n,hb0p
       ! Basis: deformations, number of quanta, integration mesh
       Write(lwou) 'HO-Basis'
       Write(lwou) b0,bz,bp
       Write(lwou) n00,nb,nt,ngh,ngl,nleg
       Write(lwou) xh,xl,wh,wl
       ! Quantum numbers and \Omega-blocks information
       Write(lwou) 'QuantNum'
       Write(lwou) ID
       ibasis=0
       Do ib=1,nb
          nd=ID(ib)
          Do n1=1,nd
             ibasis=ibasis+1
             nla=NL(ibasis); nra=NR(ibasis); nza=NZ(ibasis); nsa=NS(ibasis)
             Write(lwou) nra,nza,nla,nsa
          End Do
       End Do
       ! Pairing
       Write(lwou) 'Various.'
       Write(lwou) si,etot,rms,bet,xmix
       Write(lwou) pwi,del,ept,ala,ala2,alast
       Write(lwou) tz,varmas,varmasNZ,pjmassNZ,ass,skass
       ! Constraints: requested values, Lagrange parameters
       Write(lwou) 'Constrai'
       Write(lwou) numberCons,lambdaMax
       Write(lwou) multLambda
       Write(lwou) multRequested
       Write(lwou) qmoment
       Write(lwou) multLag
       If(set_neck_constrain) Then
          Write(lwou) neckRequested
          Write(lwou) neckValue,Z_NECK
          Write(lwou) neckLag
       End If
       Write(lwou) qfield
       ! Density matrix and pairing tensor in coordinate space
       Write(lwou) 'Densits.'
       Write(lwou) ro     ! Rho
       Write(lwou) aka    ! Kappa
       ! HF and pairing field in coordinate space
       Write(lwou) 'FieldsN.'
       Write(lwou) vn     ! RHO_ij
       Write(lwou) vhbn   ! TAU_ij
       Write(lwou) vrn    ! NABLAr RHO__ij
       Write(lwou) vzn    ! NABLAz RHO__ij
       Write(lwou) vdn    ! DELTA RHO_ij
       Write(lwou) vsn    ! NABLA . J__ij
       Write(lwou) vSFIZn ! JFIZ_ij
       Write(lwou) vSZFIn ! JZFI_ij
       Write(lwou) vSFIRn ! JFIR_ij
       Write(lwou) vSRFIn ! JRFI_ij
       Write(lwou) dvn    ! \Delta_ij
       Write(lwou) 'FieldsP.'
       Write(lwou) vp
       Write(lwou) vhbp
       Write(lwou) vrp
       Write(lwou) vzp
       Write(lwou) vdp
       Write(lwou) vsp
       Write(lwou) vSFIZp
       Write(lwou) vSZFIp
       Write(lwou) vSFIRp
       Write(lwou) vSRFIp
       Write(lwou) dvp
       ! THO
       If(switch_to_THO /= 0) Then
          Write(lwou) 'THObasis'
          If(iLST > 0) Then
             If(Allocated(fdsx)) Then
                Write(lwou) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
                Write(lwou) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0, &
                            fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
             End If
          End If
       End If
       ! Pairing regularization
       If(pairing_regularization) Then
          Write(lwou) 'Regular.'
          Write(lwou) MEFFn
          Write(lwou) MEFFp
          Write(lwou) geff_inv
       End If
       ! Gogny force
       If(finite_range .Or. coulomb_gaussian) Then
          Write(lwou) 'GognyVNN'
          Write(lwou) NumVz,NumVr
          Write(lwou) rk
          Write(lwou) ak
          Write(lwou) vrGogny
          Write(lwou) vzGogny
       End If
       ! Full HFB matrix (optional)
       If(full_HFB_matrix) Then
          Write(lwou) 'HFBmatrX'
          Write(lwou) nhhdim
          Write(lwou) brin
       End If
    End If

    ! is=3: Data from after diagonalization
    If(is == 3) Then
       ! Blocking
       Write(lwou) 'Blocking'
       Do it=1,2
          Call blosort(it,blomax(it))
       End Do
       Write(lwou) bloall
       Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
       Write(lwou) 'Blk-rest'
       Write(lwou) blocking_never_done
       ! Temperature
       If(set_temperature) Then
          Write(lwou) 'Temperat'
          Write(lwou) temper,entropy
          Write(lwou) nuv,nqp
          Write(lwou) fp_T
          Write(lwou) fn_T
       End If
    End If

    Close(lwou)

  End Subroutine write_data
  !=======================================================================
  !> This subroutine writes the binary file using the old convention of HFBTHO version 200d. It is only used
  !> for debugging purposes.
  !=======================================================================
  Subroutine write_data_old()
    Use UNEDF, Only : CpV0,CpV1
    Use HFBTHO, Only : n00,nb,nt,b0,bz,bp,beta0,si,etot,rms,bet,xmix,pwi,del,ept,ala,alast,tz,varmas,varmasNZ, &
                       pjmassNZ,ass,skass,ntx,NB,nhhdim,lambdaMax,ID,brin,NL,NR,NZ,NS,iLST,lwou,npr, &
                       bloall,bloblo,blo123,blok1k2,blomax,bloqpdif,itass,iqqmax,do_print
    Use HFBTHO_Gauss, Only : ngh,ngl
    Use HFBTHO_Lipkin, Only : ala2
    Use HFBTHO_THO, Only : decay,rmm3,amm3,bmm3,cmm3,fdsx,fdsy,fdsy1,fdsy2,fdsy3, &
                           fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
    Use HFBTHO_multipole_moments, Only : multLag

    Integer(ipr) :: iw,N1,ND,it,ib,ibasis
    Integer(ipr) :: npr1,npr11,NLANSA1,NLA,NRA,NZA,NSA

    npr11=npr(1); npr1=npr(2)
    Write(lwou) npr11,npr1,ngh,ngl,n00,nb,nt
    Write(lwou) b0,bz,bp,beta0,si,etot,rms,bet,xmix,CpV0,CpV1,pwi,  &
                del,ept,ala,ala2,alast,tz,varmas,varmasNZ,pjmassNZ, &
                ass,skass
    Write(lwou) ntx,NB,nhhdim
    Write(lwou) lambdaMax
    Write(lwou) multLag
    Write(lwou) ID
    Write(lwou) brin
    ibasis=0
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
          ibasis=ibasis+1
          NLA=NL(ibasis); NRA=NR(ibasis); NZA=NZ(ibasis); NSA=NS(ibasis); NLANSA1=(-1)**(NZA+NLA)
          Write(lwou) 2*NLA+NSA,NLANSA1,NZA+2*NRA+NLA,NZA,NLA
       End Do
    End Do
    ! blocking: sort blocking candidates first
    Do it=1,2
       Call blosort(it,blomax(it))
    End Do
    Write(lwou) bloall
    Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
    ! THO
    If(iLST > 0) Then
       If(Allocated(fdsx)) Then
          Write(lwou) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
          Write(lwou) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0, &
                      fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
       End If
    End If
    Close(lwou)
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'(a,a,a)')
          Write(iw,'(a,a,a)') '  Writing to wel_file: ',welfile
          Write(iw,'(a,a,a)') ' __________________________________  '
          Write(iw,'(a,a,a)') '  The tape ',welfile,' recorded:     '
          Write(iw,'(a,a,a)') '  nucname,npr,ngh,ngl,n00,nb,nt      '
          Write(iw,'(a,a,a)') '  b0,beta0,si,etot,rms,bet,xmix      '
          Write(iw,'(a,a,a)') '  pairing:     CpV0,CpV1,pwi         '
          Write(iw,'(a,a,a)') '  delta:       del,ept               '
          Write(iw,'(a,a,a)') '  lambda:      ala,ala2,alast,tz     '
          Write(iw,'(a,a,a)') '  asymptotic:  varmas,ass,skass      '
          Write(iw,'(a,a,a)') '  ntx,nb,nhhdim,id,N_rz,n_r,n_z      '
          Write(iw,'(a,a,a)') '  Omega2,Sigma2,Parity,Lambda        '
          Write(iw,'(a,a,a)') '  matrices(inbro):    hh,de          '
          Write(iw,'(a,a,a)') '  *all blocking candidates           '
          If(Allocated(fdsx)) Write(iw,'(a,a,a)') '  *all THO arrays                    '
          Write(iw,'(a,a,a)') ' __________________________________  '
          Write(iw,'(a,a,a)')
       End Do
    End If

  End Subroutine write_data_old
  !=======================================================================
  !> This routine sorts blocking candidates before writing them to disk.
  !=======================================================================
  Subroutine blosort(it,n)
    Use HFBTHO, Only : bloqpdif,bloblo,blo123,blok1k2

    Integer(ipr), Intent(In) :: it !< - Isospin (=particle type): it=1 neutrons, it=2, protons
    Integer(ipr), Intent(In) :: n  !< - Number of blocking candidates

    Integer(ipr) :: ip,i,k,j
    Real(pr) :: p

    Do i=1,n
       k=i; p=bloqpdif(i,it)
       If(i < n) Then
          Do j=i+1,n
             If(bloqpdif(j,it) < p) Then
                k=j; p=bloqpdif(j,it)
             End If
          End Do
          If(k /= i) Then
             bloqpdif(k,it)=bloqpdif(i,it); bloqpdif(i,it)=p
             ip = bloblo(k,it);  bloblo(k,it)  = bloblo(i,it);  bloblo(i,it)  = ip
             ip = blo123(k,it);  blo123(k,it)  = blo123(i,it);  blo123(i,it)  = ip
             ip = blok1k2(k,it); blok1k2(k,it) = blok1k2(i,it); blok1k2(i,it) = ip
          End If
       End If
    End Do

  End Subroutine blosort
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_io
