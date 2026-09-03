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

!=======================================================================
!> This routine launches an instance of the HFBTHO program. This instance can correspond to a single
!> HFBTHO calculation, possibly distributed across different MPI ranks, or a full mass table, dripline
!> or PES calculation.
!
!> @param[in] filename_hfbtho - Name of the HFBTHO input file; default: hfbtho_NAMELIST.dat
!> @param[in] filename_unedf - Name of the UNEDF input file (optional); default: hfbtho_UNEDF.dat
!> @param[in] my_comm_world - MPI communicator that contains all MPI processes
!> @param[in] my_comm_team - MPI communicator of the process' team
!> @param[in] my_n_teams - Total number of teams in existence
!> @param[in] my_team_color - Unique ID of the process' team
!> @param[in] toggle_output - Integer that decides if results will be printed out to the files
!> @param[in] filename_output - Name of the 'condensed' output file; default: hfbtho.out
!> @param[in] filename_dat - Name of the 'extended' output file; default: thoout.dat
!> @param[in] filename_binary - Name of the binary file; default: hfbtho_output.hel
!=======================================================================
Subroutine Main_Program(filename_hfbtho,filename_unedf, &
                        my_comm_world,my_comm_team,my_n_teams,my_team_color,  &
                        toggle_output,filename_output,filename_dat,filename_binary)
  Use HFBTHO_utilities
#if(USE_QRPA==1)
  Use HFBTHO_storage
#endif
#if(USE_MPI>0)
  Use mpi
  Use HFBTHO_mpi_communication
  Use HFBTHO, Only : ID_string,mpi_taskid,mpi_size,lambda_active,automatic_basis,proton_number,neutron_number, &
                     expectation_values,do_print,task_error,ierror_flag,ierror_info,task_error_gthr,COMM_world, &
                     ierr_mpi,team_rank,team_size
#else
  Use HFBTHO, Only : ID_string,mpi_taskid,lambda_active,automatic_basis,proton_number,neutron_number, &
                     expectation_values,do_print,task_error,ierror_flag,ierror_info
#endif
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
  Use HFBTHO_large_scale
#endif
#if(READ_FUNCTIONAL==1)
  Use HFBTHO_read_functional
#endif
  Use HFBTHO_io, Only : set_ASCII_filenames
  Use HFBTHO_solver, Only : initialize_HFBTHO_NAMELIST,read_HFBTHO_NAMELIST,adjust_basis

  Implicit None

  ! Inputs
  Character(Len=256), Intent(In) :: filename_hfbtho,filename_unedf
  Logical, Intent(In) :: toggle_output
  Character(Len=256), Intent(In) :: filename_output,filename_dat,filename_binary
  Integer, Intent(In) :: my_comm_world
  Integer, Intent(In) :: my_comm_team
  Integer, Intent(In) :: my_n_teams
  Integer, Intent(In) :: my_team_color

  Integer :: ierr,i

  ierror_flag = 0

  !-------------------------------------------------------------------
  ! Setting up MPI parallel environment. All variables related to MPI, particularly ranks and size of current
  ! task in various subcommunicators, are set there
  !-------------------------------------------------------------------
#if(USE_MPI>0)
  Call configure_MPI_environment(my_comm_world, my_comm_team, my_n_teams, my_team_color, ierr)
#endif

  task_error = 0
  Write(ID_string,'(i6.6)') mpi_taskid
  !-------------------------------------------------------------------
  ! Preset namelist data for the requested nucleus
  !-------------------------------------------------------------------
  Call initialize_HFBTHO_NAMELIST()
  !-------------------------------------------------------------------
  ! USE_MPI=1: spreads a single HFBTHO calculations across multiple MPI ranks, possibly within
  !            a HFBTHO library. All communications are performed within a HFB team.
  ! Team leader reads general input data. When there is no team, team_rank = 0
  ! coincides with overall rank = 0
  !-------------------------------------------------------------------
#if(USE_MPI==1)
  If(team_rank == 0) Call read_HFBTHO_NAMELIST(filename_hfbtho)
  Call mpi_bcast(ierror_flag,1,mpi_integer,0,my_comm_team,ierr)
  If(ierror_flag /= 0) Then
     If(team_rank == 0) Then
        Do i=1,ierror_flag
           Write(*,'(a11,i2,2x,a)') ' error_flag=',i,ierror_info(i)
        End Do
     End If
     Return
  End If
  !-------------------------------------------------------------------
  ! Broadcast of process-independent input data in parallel mode
  !-------------------------------------------------------------------
  If(team_rank == 0) then
     Call Construct_Vectors()
     If(Allocated(task_error_gthr)) Deallocate(task_error_gthr)
     Allocate(task_error_gthr(0:team_size-1))
  Else
     If(Allocated(task_error_gthr)) Deallocate(task_error_gthr)
     Allocate(task_error_gthr(0:0))
  End If
  Call broadcast_vectors()
  If(team_rank /= 0) Call Deconstruct_Vectors()
  Deallocate(task_error_gthr)
#else
  !-------------------------------------------------------------------
  ! USE_MPI=0: no parallelism at all
  ! USE_MPI=2: hard-coded large-scale calculations, incompatible with HFBTHO library. The code
  !            can only be run as an external program. All communications are performed within
  !            MPI_COMM_WORLD.
  ! Process 0 reads general input data. In case of specific calculation modes (mass tables,
  ! drip lines, potential energy surfaces), additional data is also read by process 0
  !-------------------------------------------------------------------
  If(mpi_taskid == 0) Call read_HFBTHO_NAMELIST(filename_hfbtho)
#if(USE_MPI==2)
  Call mpi_bcast(ierror_flag,1,mpi_integer,0,my_comm_world,ierr)
#endif
  If(ierror_flag /= 0) Then
     If(mpi_taskid == 0) Then
        Do i=1,ierror_flag
           Write(*,'(a11,i2,2x,a)') ' error_flag=',i,ierror_info(i)
        End Do
     End If
     Return
  End If
  If(mpi_taskid == 0) Then
#if(DO_MASSTABLE==1)
     Call read_HFBTHO_MassTable()
#endif
#if(DO_PES==1)
     Call read_HFBTHO_PES()
#endif
#if(DRIP_LINES==1)
     Call read_HFBTHO_StableLine()
#endif
#if(READ_FUNCTIONAL==1)
     Call read_HFBTHO_Functional()
#endif
  End If
  ! Broadcast of process-independent input data in parallel mode
#if(USE_MPI==2)
  If(mpi_taskid == 0) then
     Call Construct_Vectors()
     If(Allocated(task_error_gthr)) Deallocate(task_error_gthr)
     Allocate(task_error_gthr(0:mpi_size-1))
  Else
     If(Allocated(task_error_gthr)) Deallocate(task_error_gthr)
     Allocate(task_error_gthr(0:0))
  End If
  Call broadcast_vectors()
  If(mpi_taskid > 0) Call Deconstruct_Vectors()
#if(READ_FUNCTIONAL==1)
  Call broadcast_functional()
#endif
#endif
#endif
  !-------------------------------------------------------------------
  ! Overwrite basis characteristics if so requested
  !-------------------------------------------------------------------
  If(lambda_active(2) > 0 .And. automatic_basis) Call adjust_basis(expectation_values(2),.False.,proton_number,neutron_number)
  !-------------------------------------------------------------------
  ! Calculation modes: mass table, driplines, PES or stand-alone
  !-------------------------------------------------------------------
#if(DO_MASSTABLE==1)
  Call compute_mass_table(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
#elif(DRIP_LINES==1)
  Call compute_driplines(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
#elif(DO_PES==1)
  Call compute_PES(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
#else
  !-------------------------------------------------------------------
  ! Defining filenames and opening the ASCII files for printing
  !-------------------------------------------------------------------
  Call set_ASCII_filenames(toggle_output,filename_output,filename_dat)
  !-------------------------------------------------------------------
  ! Perform single HFBTHO calculation
  !-------------------------------------------------------------------
  Call execute_HFBTHO(filename_unedf,filename_binary)
#endif
  !-------------------------------------------------------------------
  ! Collecting data and closing files
  !-------------------------------------------------------------------
#if(USE_MPI==2)
  Call mpi_gather(task_error,1,mpi_integer,task_error_gthr,1,mpi_integer,0,COMM_world,ierr_mpi)
  If(mpi_taskid == 0) Then
     If(sum(task_error_gthr) == 0) Then
        Write(*,*)
        Write(*,'(a)') ' Parallel execution ended without errors'
        Write(*,*)
     End If
  End If
#endif
#if(DO_PES==0 && DO_MASSTABLE==0 && DRIP_LINES==0)
  If(lout < lfile) Close(lfile) ! close the output
  If(do_print == 1) Close(lout) ! close the output
#endif

End Subroutine Main_Program
!=======================================================================
!> This routines runs a single HFB calculation (possibly follwoed by symmetry restoration).
!=======================================================================
Subroutine execute_HFBTHO(filename_unedf,filename_binary)
  Use HFBTHO_utilities
  Use HFBTHO, Only : n00_INI,b0_INI,q_INI,npr_INI,kindhfb_INI,b2_0,b3_0,b4_0,MAX_ITER_INI,epsi_INI,inin_INI,skyrme_INI, &
                     Add_Pairing_INI,icou_INI,set_pairing,V0n_INI,V0p_INI,pwi_INI,cpv1_INI,iLST_INI,keypjn_INI,keypjp_INI, &
                     iproj_INI,npr1pj_INI,npr2pj_INI,switch_on_temperature,temper,neck_constraints,neckRequested,ngh_INI, &
                     ngl_INI,nleg_INI,basis_HFODD_INI,nstate_INI,Parity_INI,IDEBUG_INI,DO_FITT_INI,Print_HFBTHO_Namelist_INI, &
                     number_of_shells,oscillator_length,basis_deformation,neutron_number,proton_number,type_of_calculation, &
                     beta2_deformation,beta3_deformation,beta4_deformation,number_iterations,accuracy,restart_file,functional, &
                     add_initial_pairing,type_of_coulomb,include_3N_force,user_pairing,vpair_n,vpair_p,pairing_cutoff, &
                     pairing_feature,neutron_blocking,proton_blocking,switch_to_THO,gauge_points,projection_is_on,delta_N, &
                     delta_Z,set_temperature,temperature,collective_inertia,fission_fragments,pairing_regularization, &
                     automatic_basis,localization_functions,set_neck_constrain,neck_value,number_Gauss,number_Laguerre, &
                     number_Legendre,compatibility_HFODD,number_states,force_parity,print_time,iiter,ierror_flag,mpi_taskid, &
                     lambdaMax,lambda_active,fragment_properties,tz_fragments,bloqpdif,nkblo_INI,blomax,PNP_is_on,AMP_is_on, &
                     do_print,ierror_info,real_N,real_Z,task_error,spinors_on_grid,Nx_lattice,Ny_lattice,Nz_lattice
  Use UNEDF, Only : read_UNEDF_NAMELIST,use_3N_couplings
  Use HFBTHO_localization, Only : localization
#if(READ_FUNCTIONAL==1)
  Use HFBTHO_read_functional, Only : replace_functional
#endif
  Use HFBTHO_io, Only : welfile
  Use HFBTHO_solver, Only : check_consistency,HFBTHO_DFT_SOLVER
  Use HFBTHO_projections, Only : HFBTHO_restore

  Implicit None

  Character(Len=256), Intent(In) :: filename_unedf !< - Name of the UNEDF file (optional)
  Character(Len=256), Intent(In) :: filename_binary !< - Name of the binary file

  ! Local variables
  Integer(ipr) :: iblocase(2),nkblocase(2,5)
  Integer(ipr) :: it,i,icount,jcount,l,noForce,iw

  n00_INI                   = number_of_shells       ! number of shells
  b0_INI                    = oscillator_length      ! oscillator length
  q_INI                     = basis_deformation      ! deformation beta_2 of the basis
  npr_INI(1)                = neutron_number         ! N
  npr_INI(2)                = proton_number          ! Z
  kindhfb_INI               = type_of_calculation    ! 1: HFB, -1: HFB+LN

  b2_0                      = beta2_deformation      ! beta2 parameter of the initial WS solution
  b3_0                      = beta3_deformation      ! beta3 parameter of the initial WS solution
  b4_0                      = beta4_deformation      ! beta4 parameter of the initial WS solution

  MAX_ITER_INI              = number_iterations      ! max number of iterations
  epsi_INI                  = accuracy               ! convergence of iterations
  inin_INI                  = restart_file           ! restart from file

  skyrme_INI                = TRIM(functional)       ! functional
  Add_Pairing_INI           = add_initial_pairing    ! add pairing starting from file
  icou_INI                  = type_of_coulomb        ! coul: no-(0), dir.only-(1), plus exchange-(2)
  use_3N_couplings          = include_3N_force       ! Include 3N force on certain DME functionals

  set_pairing               = user_pairing           ! pairing is defined by user if .True.
  V0n_INI                   = vpair_n                ! pairing strength for neutrons
  V0p_INI                   = vpair_p                ! pairing strength for protons
  pwi_INI                   = pairing_cutoff         ! pairing q.p. cutoff
  cpv1_INI                  = pairing_feature        ! Type of pairing: volume, surface, mixed

  nkblocase(1,:)            = neutron_blocking       ! config. of neutron blocked state
  nkblocase(2,:)            = proton_blocking        ! config. of proton blocked state

  iLST_INI                  = switch_to_THO          ! 0:HO, -1:HO->THO, 1:THO
  keypjn_INI                = gauge_points           ! PNP: number of gauge points
  keypjp_INI                = gauge_points           ! PNP: number of gauge points
  iproj_INI                 = projection_is_on       ! projecting on different nucleus
  npr1pj_INI                = delta_N                ! its neutron number
  npr2pj_INI                = delta_Z                ! its proton number

  switch_on_temperature     = set_temperature        ! switches on temperature mode
  temper                    = temperature            ! value of the temperature

  collective_inertia        = collective_inertia     ! calculate collective mass and zero-point energy
  fission_fragments         = fission_fragments      ! calculate fission fragment characteristics
  pairing_regularization    = pairing_regularization ! activates the regularization of the pairing force
  automatic_basis           = automatic_basis        ! computes localization functions
  localization_functions    = localization_functions ! computes localization functions

  neck_constraints          = set_neck_constrain     ! activate the constraint on the neck
  neckRequested             = neck_value             ! set the requested value for the neck

  ngh_INI                   = number_Gauss           ! number of Gauss-Hermite points for z-direction
  ngl_INI                   = number_Laguerre        ! number of Gauss-Laguerre points for rho-direction
  nleg_INI                  = number_Legendre        ! number of Gauss-Legendre points for Coulomb
  basis_HFODD_INI           = compatibility_HFODD    ! flag to enforce same basis as HFODD
  nstate_INI                = number_states          ! total number of states in basis
  Parity_INI                = force_parity           ! reflection symmetry
  IDEBUG_INI                = print_time             ! debug

  DO_FITT_INI               = .False.                ! calculates quantities for reg.optimization
  Print_HFBTHO_Namelist_INI = .False.                ! Print Namelist

  ! Reset iteration number to 0
  iiter = 0
  ! Checking consistency of *_INI variables
  Call check_consistency()
  If(ierror_flag /= 0) Then
     If(do_print == 1) Then
        Do iw=lout,lfile
           Write(iw,'(a33,i4)') ' ERRORS IN Main_Program, process',mpi_taskid
           Write(iw,'(2(a3,i4),a18,f6.1)')' Z=',npr_INI(2),' N=',npr_INI(1), ' basis deformation', q_INI
           Do i=1,ierror_flag
              Write(iw,'(a8,i4,a12,i2,2x,a)') ' process',mpi_taskid,', error_flag=',i,ierror_info(i)
           End Do
           Write(iw,'("Terminating very early...")')
        End Do
     End If
     Return
  End If
  ! Check if there is at least one constraint
  icount=0
  Do l=1,lambdaMax
     If(lambda_active(l) > 0) icount=icount+1
  End Do
  ! If there is at least one constraint, check if any breaks parity
  If(icount > 0) Then
     jcount=0
     Do l=1,lambdaMax,2
        If(lambda_active(l) > 0) jcount=jcount+1
     End Do
     If(jcount > 0) Parity_INI=.False.
  Else
     collective_inertia = .False.
  End If
  If(fission_fragments) Parity_INI=.False.
  ! Enforces no-temperature mode if T < 1.e-10
  If(set_temperature .And. Abs(temper) <= 1.e-10_pr) switch_on_temperature=.False.

  Call read_UNEDF_NAMELIST(skyrme_INI,noForce,filename_unedf)
  If(skyrme_INI == 'SeaLL1' .Or. skyrme_INI == 'SLY4mod') Then
     pairing_regularization = .True.; pwi_INI=100.0
  End If
  ! If functional is used, projection automatically switched off
  If(noForce == 0) iproj_INI=0
  ! Read parameters of the energy functionals from a file
#if(READ_FUNCTIONAL==1)
  Call replace_functional()
#endif
  !---------------------------------------------------------------------------
  !                                BLOCKING
  !  Blocking candidates are selected from HFB results in the even-even parent nucleus
  !  (with Z+/- 1 or/and N +/- 1 particles).
  !   - If nkblocase(it,2) = 0, we block all blocking candidates within an
  !     energy window of pwiblo = 1 MeV (set in preparer()) around the Fermi level of
  !     the even-even parent nucleus
  !   - If nkblocase(it,2)  /=  0, we block the single candidate identified by the quantum
  !     numbers contained in the array nkblocase
  !---------------------------------------------------------------------------
  iblocase=0; bloqpdif=zero ! blomax will be charged from the previous solution
  Do it=1,2
     ! Blocking all qp within an energy window
     If(nkblocase(it,1) /= 0 .And. nkblocase(it,2) == 0) Then
        If(it == 1) Then
           iblocase(1)=iblocase(1)+1
           If(iblocase(1) > blomax(1)) iblocase(1)=1
        Else
           If(iblocase(1) <= 1) iblocase(2)=iblocase(2)+1
        End If
        nkblo_INI(it,1)=Sign(iblocase(it),nkblocase(it,1))
        nkblo_INI(it,2)=0
     Else
        ! Blocking a single qp
        nkblo_INI(it,:)=nkblocase(it,:)
     End If
  End Do
  !--------------------------------------------------------------------
  ! Run the solver in all cases EVEN/ODDS, FITS/NO-FITS
  !--------------------------------------------------------------------
  Call HFBTHO_DFT_SOLVER(filename_binary)
  !--------------------------------------------------------------------
  ! Restore broken symmetries of HFB states
  !--------------------------------------------------------------------
  If((PNP_is_on >= 1 .Or. AMP_is_on >= 1) .And. .Not. fission_fragments) Call HFBTHO_restore()
  !--------------------------------------------------------------------
  ! Calculate localization functions
  !--------------------------------------------------------------------
  If(localization_functions) Call localization()
  !--------------------------------------------------------------------
  ! Display error messages in case of problems
  !--------------------------------------------------------------------
#if(USE_MPI>0)
  If(do_print == 1) Then
     Do iw=lout,lfile
        If(ierror_flag == 0) Then
           Write(iw,*)
           Write(iw,'(a)') ' HFBTHO_SOLVER ended without errors'
           Write(iw,*)
        Else
           task_error = 1
           Write(iw,*)
           Write(iw,'(a33,i4)') ' ERRORS IN HFBTHO_SOLVER, process',mpi_taskid
           Do i=1,ierror_flag
              Write(iw,'(a8,i4,a12,i2,2x,a)') ' process',mpi_taskid,', error_flag=',i,ierror_info(i)
           End Do
           Write(iw,*)
        End If
     End Do
     ierror_flag = 0
  End If
#else
  If(do_print == 1) Then
     Do iw=lout,lfile
        If(ierror_flag == 0) Then
           Write(iw,*)
           Write(iw,'(a)') ' HFBTHO_SOLVER ended without errors'
           Write(iw,*)
        Else
           task_error = 1
           Write(iw,*)
           Write(iw,'(a33,i4)') ' ERRORS IN HFBTHO_SOLVER, process',mpi_taskid
           Write(iw,'(2(a3,i4),a18,f6.1)')' Z=',npr_INI(2),' N=',npr_INI(1), ' basis deformation', q_INI
           Do i=1,ierror_flag
              Write(iw,'(a8,i4,a12,i2,2x,a)') ' process',mpi_taskid,', error_flag=',i,ierror_info(i)
           End Do
           Write(iw,*)
        End If
     End Do
     ierror_flag = 0
  End If
#endif

End Subroutine execute_HFBTHO
!=======================================================================
!> This routine computes a full mass table. It is only compiled and called if DO_MASSTABLE=1.
!>
!> @param[in] filename_unedf - Name of the UNEDF input file (optional); default: hfbtho_UNEDF.dat
!> @param[in] toggle_output - Integer that decides if results will be printed out to the files
!> @param[in] filename_output - Name of the 'condensed' output file; default: hfbtho.out
!> @param[in] filename_dat - Name of the 'extended' output file; default: thoout.dat
!> @param[in] filename_binary - Name of the binary file; default: hfbtho_output.hel
!=======================================================================
Subroutine compute_mass_table(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
#if(DO_MASSTABLE==1)
  Use HFBTHO_utilities
  Use HFBTHO, Only : mpi_taskid,iRow,Z_chain,N_chain,A_chain,beta_deformation,Q20,Z_masstable,N_masstable,beta_masstable, &
                     row_string,nRows,mpi_size,proton_number,neutron_number,expectation_values,basis_deformation,ID_string, &
                     team_rank,lfile,Q20_masstable
  Use HFBTHO_io, Only : set_ASCII_filenames
#if(USE_MPI==2)
  Use HFBTHO_large_scale
  Use HFBTHO_mpi_communication
#endif

  Character(Len=256), Intent(In) :: filename_unedf
  Logical, Intent(In) :: toggle_output
  Character(Len=256), Intent(In) :: filename_output,filename_dat,filename_binary

  Integer(ipr) :: icalc

  ! Allocation of vectors
  If(nRows > 0) Then
#if(USE_MPI==2)
     Call allocate_out_vectors()
#else
     Call allocate_mass_table()
#endif
  End If

  icalc = 0
  ! Team leader creates a file for bookkeeping
  If(mpi_taskid == 0) Open(127,file='TableLog.dat')
  ! Loop over elements of the mass table
  Do iRow = 0,nRows
     Z_chain = Z_masstable(iRow)
     N_chain = N_masstable(iRow)
     A_chain = Z_chain+N_chain
     beta_deformation = beta_masstable(iRow)
     Q20 = Q20_masstable(iRow)!beta_deformation*sqrt(5/pi)*(A_chain)**(5/3._pr)/100._pr
     Write(row_string,'("_",i6.6)') iRow
     If(iRow == 0 .And. nRows > 0) Cycle
     If(mpi_taskid == 0) Write(127,'(a7,2i5,2f15.8)') row_string,Z_chain,N_chain,Q20,beta_deformation
     !only do the calculations that correspond to your task id
     If(Mod(iRow,mpi_size) /= mpi_taskid) Cycle
     proton_number  = Z_chain
     neutron_number = N_chain
     expectation_values(2) = Q20
     basis_deformation = beta_deformation
     !-------------------------------------------------------------
     ! Defining filenames and opening the ASCII files for printing
     !-------------------------------------------------------------
     Call set_ASCII_filenames(toggle_output, filename_output, filename_dat)
     !-------------------------------------------------------------------
     ! Perform single HFBTHO calculation
     !-------------------------------------------------------------------
     Call execute_HFBTHO(filename_unedf,filename_binary)
     If(nRows > 0) Then
        Write(*,'("task ",a6," finished row ",a6)') ID_string, row_string(2:7)
#if(USE_MPI==2)
        Call fill_out_vectors(icalc)
#else
        Call fill_mass_table(icalc)
#endif
       icalc = icalc + 1
     End If
     Close(lfile) ! close the output
  End Do
  If(nRows > 0) Then
#if(USE_MPI==2)
     Call gather_results()
#endif
     Call print_mass_table()
  End If
  If(team_rank == 0) Close(127)
#else
  Character(Len=256), Intent(In) :: filename_unedf
  Logical, Intent(In) :: toggle_output
  Character(Len=256), Intent(In) :: filename_output, filename_dat, filename_binary
#endif

End Subroutine compute_mass_table
!=======================================================================
!> This routine computes a full mass table from dripline to dripline. It is only compiled and called
!> if DRIP_LINES=1.
!>
!> @param[in] filename_unedf - Name of the UNEDF input file (optional); default: hfbtho_UNEDF.dat
!> @param[in] toggle_output - Integer that decides if results will be printed out to the files
!> @param[in] filename_output - Name of the 'condensed' output file; default: hfbtho.out
!> @param[in] filename_dat - Name of the 'extended' output file; default: thoout.dat
!> @param[in] filename_binary - Name of the binary file; default: hfbtho_output.hel
!=======================================================================
Subroutine compute_driplines(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
#if(DRIP_LINES==1)
  Use HFBTHO_utilities
  Use HFBTHO, Only : mpi_taskid,iRow,Z_chain,N_chain,A_chain,beta_deformation,Q20,Z_stable_line,N_stable_line,i_deformation, &
                     row_string,nRows,proton_number,neutron_number,expectation_values,basis_deformation,kindhfb_INI,Energy_chain, &
                     team_rank,number_teams,team_color,lfile,pi,beta_step,team_string,separation_2N, Minimum_Energy_Prev, &
                     Minimum_Energy,direction_sl,Calc_counter,Energy_chain_gthr,direction_str,team_size
  Use HFBTHO_io, Only : set_ASCII_filenames
#if(USE_MPI==2)
  Use HFBTHO_large_scale
  Use HFBTHO_mpi_communication
#endif

  Character(Len=256), Intent(In) :: filename_unedf
  Logical, Intent(In) :: toggle_output
  Character(Len=256), Intent(In) :: filename_output,filename_dat,filename_binary

  beta_step = 1.0_pr/Real(number_deformations-1,kind=pr)
  Write(team_string,'(1i3.3)') team_color
  ! Team leader allocates array to recieve energies and creates a file for bookkeeping
  If(team_rank == 0) Then
     Allocate(Energy_chain_gthr(0:team_size-1))
     open(127,file='TeamTable'//team_string//'.dat')
  End If
  calc_counter = 0
  ! Loop over the nuclei in the "stable line"
  Do iRow = 1,nRows
     ! Only do the chains that correspond to your team
     If(Mod(iRow,number_teams) /= team_color) Cycle
     Minimum_Energy_Prev = 100
     Z_chain = Z_stable_line(iRow)
     N_chain = N_stable_line(iRow)
     ! Move along the isotopic (or isotonic) chain until the drip line is reached
     Do
        A_chain = Z_chain + N_chain
        beta_deformation = -0.5_pr - beta_step
        ! Loop over the different basis deformations of each nucleus
        Do i_deformation = 1,number_deformations
           beta_deformation = beta_deformation + beta_step
           calc_counter = calc_counter + 1
           Q20 = beta_deformation*Sqrt(5.0_pr/pi)*(A_chain)**(5.0_pr/3.0_pr)/100._pr
           Write(row_string,'("_",a3,"_",i6.6)') team_string,calc_counter
           ! Team leader writes type of calculation for bookkeeping
           If(team_rank == 0) Then
              If(direction_sl(iRow) == 1) Then
                 direction_str = ' isotopic'
              Else
                 direction_str = ' isotonic'
              End If
              Write(127,'(a11,2i5,2f15.8,a9)') row_string,z_chain,n_chain,Q20,beta_deformation,direction_str
           End If
           If(number_deformations == 1) beta_deformation = 0._pr
           ! Only calculate what corresponds to each process
           If(Mod(i_deformation,team_size) /= team_rank) cycle
           proton_number  = Z_chain
           neutron_number = N_chain
           expectation_values(2) =  Q20
           basis_deformation = beta_deformation
           !-------------------------------------------------------------
           ! Defining filenames and opening the ASCII files for printing
           !-------------------------------------------------------------
           Call set_ASCII_filenames(toggle_output, filename_output, filename_dat)
           !-------------------------------------------------------------------
           ! Perform single HFBTHO calculation
           !-------------------------------------------------------------------
           Call execute_HFBTHO(filename_unedf,filename_binary)
           ! Calculations without Lipkin-Nogami
           If(kindhfb_INI > 0) Then
              Energy_chain = ehfb
           ! Calculations with Lipkin-Nogami
           Else
              Energy_chain = etot
           End If
           Close(lfile) ! close the output
        ! end loop over deformations
        End Do
        Call find_minimum_energy()
        separation_2N = Minimum_Energy_Prev - Minimum_Energy
        Minimum_Energy_Prev = Minimum_Energy
        If(separation_2N < 0._pr) Then
           ! Drip line has been reached
           Exit
        Else
           ! Drip line has not been reached
           If(direction_sl(iRow) == 1) Then
              If(N_chain >= 310) Exit
              N_chain = N_chain + 2
           Else
              If(Z_chain >= 120) Exit
              Z_chain = Z_chain + 2
           End If
        End If
     ! end isotopic (or isotonic) chain
     End Do
     ! team leader announces isotopic (or isotonic) chain finished
     If(team_rank == 0) Then
        If(direction_sl(iRow) == 1) Then
           Write(*,'(a4,i3,a27,i4,a6,i4)') 'team',team_color,' finished isotopic chain Z=',Z_chain,' at N=',N_chain
        Else
           Write(*,'(a4,i3,a27,i4,a6,i4)') 'team',team_color,' finished isotonic chain N=',N_chain,' at Z=',Z_chain
        End If
     End If
     Close(lfile) ! close the output
  ! end loop over nuclei inside the "stable line
  End Do
  ! team leader closes bookkeeping file
  If(mpi_taskid == 0) Close(127)
#else
  Character(Len=256), Intent(In) :: filename_unedf
  Logical, Intent(In) :: toggle_output
  Character(Len=256), Intent(In) :: filename_output, filename_dat, filename_binary
#endif

End Subroutine compute_driplines
!=======================================================================
!> This routine computes a potential energy surface. It is only compiled and called  if DO_PES=1.
!>
!> @param[in] filename_unedf - Name of the UNEDF input file (optional); default: hfbtho_UNEDF.dat
!> @param[in] toggle_output - Integer that decides if results will be printed out to the files
!> @param[in] filename_output - Name of the 'condensed' output file; default: hfbtho.out
!> @param[in] filename_dat - Name of the 'extended' output file; default: thoout.dat
!> @param[in] filename_binary - Name of the binary file; default: hfbtho_output.hel
!=======================================================================
Subroutine compute_PES(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
#if(DO_PES==1)
  Use HFBTHO_utilities
  Use HFBTHO, Only : mpi_taskid,number_teams,team_color,iRow,Z_chain,N_chain,A_chain,Z_PES,N_PES,lambda_PES, &
                     beta2_deformation,beta3_deformation,beta4_deformation,bet2_PES,bet3_PES,bet4_PES,Q_PES, &
                     automatic_basis,row_string,npoints,proton_number,neutron_number,expectation_values,ndefs, &
                     basis_deformation,ID_string,lfile,lambda_active
  Use HFBTHO_io, Only : set_ASCII_filenames
#if(USE_MPI==2)
  Use HFBTHO_large_scale
  Use HFBTHO_mpi_communication
#endif
  Use HFBTHO_solver, Only : adjust_basis

  Character(Len=256), Intent(In) :: filename_unedf
  Logical, Intent(In) :: toggle_output
  Character(Len=256), Intent(In) :: filename_output,filename_dat,filename_binary

  Integer(ipr) :: icalc,j,lambda

  icalc = 0
  ! Team leader creates a file for bookkeeping
  If(mpi_taskid == 0) Open(127,file='TableLog.dat')
  ! Loop over points in the potential energy surface
  Do iRow = 1,npoints
     Z_chain = Z_PES(iRow)
     N_chain = N_PES(iRow)
     A_chain = Z_chain+N_chain
     Write(row_string,'("_",i6.6)') iRow
     If(mpi_taskid == 0) Write(127,'(a7,2i5,9f10.3)') row_string,Z_chain,N_chain,(Q_PES(iRow,j),j=1,ndefs)
     ! Only do the calculation that correspond to your team
     If(Mod(iRow,number_teams) /= team_color) Cycle
     proton_number  = Z_chain
     neutron_number = N_chain
     ! Default basis deformation and WS deformation based on input file
     If(bet2_PES(iRow) > -8.0) Then
        !basis_deformation = bet2_PES(iRow)
        beta2_deformation = bet2_PES(iRow)
     End If
     If(bet3_PES(iRow) > -8.0) beta3_deformation = bet3_PES(iRow)
     If(bet4_PES(iRow) > -8.0) beta4_deformation = bet4_PES(iRow)
     Do j=1,ndefs
        lambda = lambda_PES(j)
        ! More advanced fit based on value of Q2 only
        If(lambda == 2 .And. automatic_basis) Call adjust_basis(Q_PES(iRow,j),.False.,proton_number,neutron_number)
        expectation_values(lambda) = Q_PES(iRow,j)
        lambda_active(lambda) = 1
     End Do
     !-------------------------------------------------------------
     ! Defining filenames and opening the ASCII files for printing
     !-------------------------------------------------------------
     Call set_ASCII_filenames(toggle_output, filename_output, filename_dat)
     !-------------------------------------------------------------------
     ! Perform single HFBTHO calculation
     !-------------------------------------------------------------------
     Call execute_HFBTHO(filename_unedf,filename_binary)
     If(npoints > 0) Then
        Write(*,'("task ",a6," finished row ",a6)') ID_string, row_string(2:7)
        icalc = icalc + 1
     End If
     If(do_print == 1) Close(lfile) ! close the output
  End Do ! End of loop over points in the PES
  If(do_print == 1) Close(127)
#else
  Character(Len=256), Intent(In) :: filename_unedf
  Logical, Intent(In) :: toggle_output
  Character(Len=256), Intent(In) :: filename_output, filename_dat, filename_binary
#endif

End Subroutine compute_PES
