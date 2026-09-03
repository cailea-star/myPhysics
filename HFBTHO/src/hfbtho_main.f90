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

Program hfbthoprog

  Use HFBTHO_utilities
#if(USE_MPI>0)
  Use mpi
  Use HFBTHO_mpi_communication, Only : Create_MPI_Teams
  Use HFBTHO, Only : ierr_mpi,COMM_world,COMM_team,HFB_cores,mpi_size,number_teams,team_color, &
#if(DRIP_LINES==1)
                     number_deformations, &
#endif
                     mpi_taskid
#else
  Use HFBTHO, Only : COMM_world,COMM_team,HFB_cores,mpi_size, &
#if(DRIP_LINES==1)
                     number_deformations, &
#endif
                     mpi_taskid
#endif

  Implicit None

  Logical :: toggle_output
  Integer :: n_args
  Integer :: my_comm_world
  Integer :: my_comm_team
  Integer :: my_n_teams
  Integer :: my_team_color
  Character(Len=256) :: filename_hfbtho !< Name of the input file containing the namelists
  Character(Len=256) :: filename_unedf  !< Name of the input file containing the EDF parametrization
  Character(Len=256) :: filename_output !< Name of the ASCII output file (short version)
  Character(Len=256) :: filename_dat    !< Name of the ASCII output file (long version)
  Character(Len=256) :: filename_binary !< Name of the binary output file

  ! Initialize MPI environment.
  !   * USE_MPI=0: no MPI requested, presetting default values
  !   * USE_MPI=1: single, multi-core HFBTHO calculation spread across HFB_cores
  !   * USE_MPI=2: large-scale, possibly multi-core HFBTHO calculations
#if(USE_MPI>0)
  Call MPI_INIT(ierr_mpi)
  COMM_world = MPI_COMM_WORLD
  Call MPI_COMM_SIZE(COMM_world, mpi_size, ierr_mpi)
  Call MPI_COMM_RANK(COMM_world, mpi_taskid, ierr_mpi)
#if(USE_MPI==1)
  HFB_cores = mpi_size
#endif
#else
  COMM_world = 0
  COMM_team = 0
  HFB_cores = 1
  mpi_size = 1
  mpi_taskid = 0
#endif

  ! Set values for all filenames
  filename_hfbtho = 'hfbtho_NAMELIST.dat'
  filename_unedf  = 'UNEDF_NAMELIST.dat'
  filename_output = 'hfbtho.out'
  filename_dat    = 'thoout.dat'
  filename_binary = 'hfbtho_output.hel'

  ! Reading the input arguments to the code
  n_args = command_argument_count()
  If      (n_args == 1) Then
    Call get_command_argument(1, filename_hfbtho)
  Else If (n_args == 2) Then
    Call get_command_argument(1, filename_hfbtho)
    Call get_command_argument(2, filename_unedf)
  Else If (n_args /= 0) Then
    Write(*,*) ""
    Write(*,*) "Invalid number of command line arguments"
    Write(*,*) ""
    Write(*,*) "The program accepts either no, one, or two command line"
    Write(*,*) "arguments.  "
    Write(*,*) ""
    Write(*,*) "For no arguments, the program will use the configuration"
    Write(*,*) "files hfbtho_NAMELIST.dat and UNEDF_NAMELIST.dat."
    Write(*,*) ""
    Write(*,*) "For one arguments, the program will use the given argument"
    Write(*,*) "as the filename for the HFBTHO namelist and use "
    Write(*,*) "UNEDF_NAMELIST.dat for the UNEDF Namelist."
    Write(*,*) ""
    Write(*,*) "For two arguments, the first argument must be the filename"
    Write(*,*) "for the HFBTHO namelist and the second to the UNEDF "
    Write(*,*) "Namelist file."
    Write(*,*) ""
    Call Exit(1)
  End If

  ! In parallel mode (MPI activated), creating teams and setting values to pass to main program.
  ! If no MPI, setting dummy values to all these quantities (they won't be used)
#if(USE_MPI>0)
#if(DRIP_LINES==1)
  Call Create_MPI_Teams(number_deformations)
#else
  Call Create_MPI_Teams(HFB_cores)
#endif
  my_comm_world = COMM_world
  my_comm_team  = COMM_team
  my_n_teams    = number_teams
  my_team_color = team_color
#else
  my_comm_world = COMM_world
  my_comm_team  = COMM_team
  my_n_teams    = 1
  my_team_color = 0
#endif

  ! Calling main HFBTHO program: this will read all the data needed for the run and execute HFBTHO.
  ! If specific calculation modes (mass tables, drip lines, potential energy surfaces) are requested,
  ! this routine will call HFBTHO several times (serial mode) or in parallel to loop over all the
  ! configurations. Otherwise, it will simply execute HFBTHO once.
  toggle_output = .True.
  Call Main_Program(filename_hfbtho, filename_unedf, &
                    my_comm_world, my_comm_team, my_n_teams, my_team_color, &
                    toggle_output, filename_output, filename_dat, filename_binary)

#if(USE_MPI>0)
  ! Wait here until all processes are done
  Call mpi_barrier(COMM_world, ierr_mpi)
  ! Close MPI environment
  Call mpi_finalize(ierr_mpi)
#endif

End Program hfbthoprog
