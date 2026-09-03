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

#if(READ_FUNCTIONAL==1)
!----------------------------------------------------------------------
!>  This module contains subroutines to read the parameters of an
!>  arbitrary Skyrme-like functional from a file instead of using the
!>  ones that are hardcoded in the UNEDF module (see hfbtho_unedf.f90).
!>  This is particularly useful for making calculations with statistical
!>  ensembles of parameter sets.
!>
!> @author
!>    Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - read_HFBTHO_Functional
!               - replace_functional
!               - broadcast_functional
!  Functions: - n_items(line)
!----------------------------------------------------------------------
Module HFBTHO_read_functional

  Use HFBTHO_utilities

  Implicit None

#if(USE_MPI==2)
  Public read_HFBTHO_Functional,replace_functional,broadcast_functional
#else
  Public read_HFBTHO_Functional,replace_functional
#endif

  Private n_items

Contains
  !=======================================================================
  !>
  !=======================================================================
  Subroutine read_HFBTHO_Functional()
    Use HFBTHO, Only : n_func_param,functional_vector

    Integer(ipr) :: lfunctional=16,i
    Character(len=500) :: line

    Open(lfunctional,file='hfbtho_FUNCTIONAL.dat')
    Read(lfunctional,*)
    Read(lfunctional,'(a)') line
    Close(lfunctional)
    n_func_param =  n_items(line)
    If(n_func_param < 12 .Or. n_func_param > 15) Then
       stop 'The number of parameters in hfbtho_FUNCTIONAL.dat has to be between 12 and 15'
    End If
    If(Allocated(functional_vector)) Deallocate(functional_vector)
    Allocate(functional_vector(1:n_func_param))
    Read(line,*) functional_vector(1:n_func_param)

  End Subroutine read_HFBTHO_Functional
  !=======================================================================
  !>
  !=======================================================================
  Function n_items(line) result(nitems)
    Character, Intent(In) :: line*(*)
    Integer(ipr) :: nitems
    Integer(ipr) :: i,n,toks

    i = 1
    n = len_trim(line)
    toks = 0
    nitems = 0
    Do While(i <= n)
       Do While(line(i:i) == ' ')
          i = i + 1
          If(n < i) Return
       End Do
       toks = toks + 1
       nitems = toks
       Do
          i = i + 1
          If(n < i) Return
          If(line(i:i) == ' ') Exit
       End Do
    End Do

  End Function n_items
  !=======================================================================
  !>
  !=======================================================================
  Subroutine replace_functional()
    Use HFBTHO, Only : n_func_param,functional_vector
    Use UNEDF, Only : E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,CrDr,CpV0,CrdJ,CJ

    RHO_NM   = functional_vector( 1)
    E_NM     = functional_vector( 2)
    K_NM     = functional_vector( 3)
    ASS_NM   = functional_vector( 4)
    LASS_NM  = functional_vector( 5)
    SMASS_NM = functional_vector( 6)
    CrDr(0)  = functional_vector( 7)
    CrDr(1)  = functional_vector( 8)
    CpV0(0)  = functional_vector( 9)
    CpV0(1)  = functional_vector(10)
    CrdJ(0)  = functional_vector(11)
    CrdJ(1)  = functional_vector(12)
    If(n_func_param == 13) Then
       VMASS_NM = functional_vector(13)
    Else If(n_func_param == 14) Then
       CJ(0)    = functional_vector(13)
       CJ(1)    = functional_vector(14)
    Else If(n_func_param == 15) Then
       CJ(0)    = functional_vector(13)
       CJ(1)    = functional_vector(14)
       VMASS_NM = functional_vector(15)
    End If

  End Subroutine replace_functional
#if(USE_MPI==2)
  !=======================================================================
  !>
  !=======================================================================
  Subroutine broadcast_functional()
    Use mpi
    Use HFBTHO, Only : n_func_param,functional_vector,COMM_world,ierr_mpi

    Call mpi_bcast(n_func_param,1,MPI_INTEGER,0,COMM_world,ierr_mpi)
    If(mpi_taskid > 0) Then
       If(Allocated(functional_vector)) Deallocate(functional_vector)
       Allocate(functional_vector(1:n_func_param))
    End If
    Call mpi_bcast(functional_vector,n_func_param,MPI_DOUBLE_PRECISION,0,COMM_world,ierr_mpi)

  End Subroutine broadcast_functional
#endif
  !=======================================================================
  !>
  !=======================================================================
End Module HFBTHO_read_functional
#endif
