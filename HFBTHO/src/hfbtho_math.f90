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
!                    MATHEMATICAL ROUTINES PACKAGE                     !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module contains several routines implementing basic linear algebra
!> and mathematical operations
!>
!>  @author
!>  Nicolas Schunck, Petar Marevic, Mario Stoitsov
!----------------------------------------------------------------------
!  Subroutines: - lingd(ma,mx,n,m,a,x,d,Ifl)
!               - calculate_inverse_cmplx(ndim,mat,det,info)
!               - calculate_inverse_real(ndim,mat,det,info)
!               - csplin(n, x, y, b, c, d)
!               - cseval(n,u,x,y,b,c,d,splf0)
!               - deri(h,n,f1,dunl)
!               - gfv()
!               - wigner(j,mm1,mm2,beta)
!               - integrate_complex_function(complex_function,n,h,complex_integral)
!               - integrate_simpson(f,n,h,r)
!               - SPLINE(X,Y,N,YP1,YPN,Y2,IERROR)
!  Functions: - PROINT(ZPOINT)
!             - DEFLEG(LAMACT,MIUACT,ZVALUE)
!             - ZBRENT(FONCTI,XLOWER,XUPPER,TOLERA,IFOUND)
!             - SPLINT(XENTRY,YENTRY,Y2A,NPOINT,XARGUM)
!----------------------------------------------------------------------
Module math

  Use HFBTHO_utilities

  Implicit None

  Public lingd,calculate_inverse_cmplx,calculate_inverse_real,csplin,cseval,deri,gfv,wigner, &
         integrate_complex_function,integrate_simpson,PROINT,DEFLEG,ZBRENT,SPLINT,SPLINE

  Private bubble_sort

  Integer(ipr), Allocatable, Public, Save :: iv(:) !< \f$ \mathtt{iv(n) } = (-1)^n    \f$
  Real(pr), Allocatable, Public, Save :: fak(:)    !< \f$ \mathtt{sq(n) } = \sqrt{n}, \f$
  Real(pr), Allocatable, Public, Save :: fi(:)     !< \f$ \mathtt{sqi(n)} = 1/\sqrt{n}\f$
  Real(pr), Allocatable, Public, Save :: sq(:)     !< \f$ \mathtt{fak(n)} = n!        \f$
  Real(pr), Allocatable, Public, Save :: sqi(:)    !< \f$ \mathtt{fi(n)} = 1/n!        \f$
  Real(pr), Allocatable, Public, Save :: wf(:)     !< \f$ \mathtt{wf(n) } = \sqrt{n!} \f$
  Real(pr), Allocatable, Public, Save :: wfi(:)    !< \f$ \mathtt{wfi(n)} = 1/\sqrt{n!}\f$

Contains
  !=======================================================================
  !> This routine olves the system of linear equations A*X = B. At the beginning the matrix B is
  !> stored in X. During the calculation it will be overwritten. D is the determinant of A.
  !=======================================================================
  Subroutine lingd(ma,mx,n,m,a,x,d,Ifl)
    Integer(ipr), Intent(In) :: ma,mx,n,m
    Integer(ipr), Intent(Inout) :: Ifl
    Real(pr), Dimension(ma,m), Intent(Inout) :: a
    Real(pr), Dimension(mx,m), Intent(Inout) :: x
    Real(pr), Intent(Inout) :: d

    Integer(ipr), Save :: i,j,k,l,k1,n1
    Real(pr), Save :: p,q,tol,cp,cq
    Real(pr), Parameter :: tollim = 1.0e-10_pr

    Ifl = 1; p = zero
    Do i=1,n
       q = zero
       Do j=1,n
          q = q + Abs(a(i,j))
       End Do
       If(q > p) p = q
    End Do
    tol = tollim*p; d   = one
    Do k=1,n
       p = zero
       Do j=k,n
          q = Abs(a(j,k))
          If(q < p) Cycle
          p = q; i = j
       End Do
       If(p < tol) Then
          Write(*,'(80("-"))')
          Write(*,'(" *****  ERROR IN LINGD, TOLERANCE =",e11.4," VALUE OF A(",i3,",",i3,") IS ",e11.4)') &
                      tol,i,k,a(i,k)
          Write(*,'(80("-"))')
          Ifl = -1
          Return
       End If
       cp = one/a(i,k)
       If(i /= k) Then
          d = -d
          Do l=1,m
             cq = x(i,l); x(i,l) = x(k,l); x(k,l) = cq
          End Do
          Do l=k,n
             cq = a(i,l); a(i,l) = a(k,l); a(k,l) = cq
          End Do
       End If
       d = d*a(k,k)
       If(k == n) Exit
       k1 = k + 1
       Do i=k1,n
          cq=a(i,k)*cp
          Do l=1,m
             x(i,l)=x(i,l)-cq*x(k,l)
          End Do
          Do l=k1,n
             a(i,l)=a(i,l)-cq*a(k,l)
          End Do
       End Do
    End Do
    Do l=1,m
       x(n,l)=x(n,l)*cp
    End Do
    If(n == 1) Return
    n1=n-1
    Do k=1,n1
       cp = one/a(n-k,n-k)
       Do l=1,m
          cq = x(n-k,l)
          Do i=1,k
             cq = cq-a(n-k,n+1-i)*x(n+1-i,l)
          End Do
          x(n-k,l) = cq*cp
       End Do
    End Do

  End Subroutine lingd
  !============================================================================
  !> This routine calculates the inverse and determinant of a complex square matrix.
  !============================================================================
  Subroutine calculate_inverse_cmplx(ndim,mat,det,info)
    Integer(ipr), Intent(In) :: ndim !< - Actual dimension of the matrix
    Complex(pr), Allocatable, Intent(Inout) :: mat(:,:) !< - Complex-valued square matrix
    Complex(pr), Intent(Inout) :: det !< - Value of the determinant
    Integer(ipr), Intent(Inout) :: info !< - Error message

    Integer(ipr) :: info1,info2,n
    Real(pr) :: sgn
    Integer(ipr), Allocatable :: ipiv(:)
    Complex(pr), Allocatable :: work(:),mat_diag(:)

    Allocate(ipiv(ndim),work(ndim))
    ipiv=0; work=c_zero

    ! ZGETRF computes an LU factorization of a complex M-by-N matrix A
    ! using partial pivoting with row interchanges.
    Call Zgetrf(ndim,ndim,mat,ndim,ipiv,info1)

    ! Calculating determinant
    Allocate(mat_diag(ndim))
    Do n=1,ndim
       mat_diag(n)=mat(n,n)
    End Do

    Call bubble_sort(mat_diag,ndim)

    det=Cmplx(1.0_pr,0.0_pr,Kind=pr); sgn=1.0_pr
    If(Mod(ndim,2) == 0) Then
       Do n=1,ndim/2
          det=det*mat_diag(n)*mat_diag(ndim+1-n)
          If(ipiv(n) /= n) sgn=-sgn
          If(ipiv(ndim+1-n) /= ndim+1-n) sgn=-sgn
       End Do
       det=sgn*det
    Else
       Do n=1,(ndim-1)/2
          det=det*mat_diag(n)*mat_diag(ndim-n)
          If(ipiv(n) /= n) sgn=-sgn
          If(ipiv(ndim-n) /= ndim-n) sgn=-sgn
       End Do
       det=det*mat_diag(ndim)
       If(ipiv(ndim) /= ndim) sgn=-sgn
       det=sgn*det
    End If

    Deallocate(mat_diag)

    ! ZGETRI computes the inverse of a matrix using the LU factorization computed by ZGETRF.
    Call Zgetri(ndim,mat,ndim,ipiv,work,ndim,info2)

    info=info1+info2

    Deallocate(ipiv,work)

  End Subroutine calculate_inverse_cmplx
  !=================================================================
  !> This subroutine sorts matrix elements based on their module and the bubble sorting algorithm.
  !=================================================================
  Subroutine bubble_sort(array_to_sort,ndim)
    Complex(pr), Allocatable, Intent(Inout) :: array_to_sort(:) !< - Array to sort
    Integer(ipr), Intent(In) :: ndim !< - Actual number of elements in the array

    Integer(ipr) :: i,j
    Complex(pr) :: temp
    Logical :: swapped

    Do j=ndim-1, 1, -1
      swapped = .False.
      Do i=1,j
        If(Abs(array_to_sort(i)) > Abs(array_to_sort(i+1))) Then
          temp = array_to_sort(i)
          array_to_sort(i) = array_to_sort(i+1)
          array_to_sort(i+1) = temp
          swapped = .True.
        End If
      End Do
      If(.not. swapped) Exit
    End Do

  End Subroutine bubble_sort
  !==========================================================================
  !> This routine calculates the inverse and determinant of a real square matrix.
  !==========================================================================
  Subroutine calculate_inverse_real(ndim,mat,det,info)
    Integer(ipr), Intent(In) :: ndim !< - Actual dimension of the matrix
    Real(pr), Allocatable, Intent(In) :: mat(:,:) !< - Real-valued square matrix
    Real(pr), Intent(Inout) :: det !< - Value of the determinant
    Integer(ipr), Intent(Inout) :: info !< - Error message

    Integer(ipr) :: info1,info2,n
    Real(pr) :: sgn
      Integer(ipr), Allocatable :: ipiv(:)
    Real(pr), Allocatable :: work(:)

      Allocate(ipiv(ndim),work(ndim))
      ipiv=0; work=0.0_pr

      ! DGETRF computes an LU factorization of a real M-by-N matrix A
      ! using partial pivoting with row interchanges.
      Call Dgetrf(ndim,ndim,mat,ndim,ipiv,info1)

      ! Calculating determinant
      det=1.0_pr; sgn=1.0_pr
      Do n=1,ndim
         det=det*mat(n,n)
         If(ipiv(n) /= n) sgn=-sgn
      End Do
      det=sgn*det

      ! DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
      Call Dgetri(ndim,mat,ndim,ipiv,work,ndim,info2)

      info=info1+info2

      Deallocate(ipiv,work)

  End Subroutine calculate_inverse_real
  !=======================================================================
  !> This routine calculates the coefficients for cubic interpolation. The coefficients are  noted
  !> \f$ b_i, c_i, d_i, i=1,2,\dots,n \f$ and are computed for a cubic interpolating spline
  !>   \f[
  !>      s(x) = y_i + b_i (x-x_i) + c_i (x-x_i)^2 + d_i (x-x_i)^3, \quad x_i \leq x \leq x_{i+1}
  !>   \f]
  !> The coefficients are given by
  !>   \f[
  !>       y_i = s(x_i), \quad\quad
  !>       b_i = \frac{ds}{dx}(x_i),  \quad\quad
  !>       c_i = \frac{1}{2}\frac{d^{2}s}{dx^{2}}(x_i),  \quad\quad
  !>       d_i = \frac{1}{6}\frac{d^{3}s}{dx^{3}}(x_i)
  !>   \f]
  !> The accompanying function \ref cseval can be used to evaluate the spline, its derivative or even
  !> its 2nd derivative.
  !>
  !>  @param[in] n - the number of data points or knots (n >= 2)
  !>  @param[in] x - the abscissas of the knots in strictly increasing order
  !>  @param[in] y - the ordinates of the knots
  !>  @param[inout] b - arrays of spline coefficients as defined above
  !>  @param[inout] c - arrays of spline coefficients as defined above
  !>  @param[inout] d - arrays of spline coefficients as defined above
  !=======================================================================
  Subroutine csplin(n,x,y,b,c,d)
    Integer(ipr), Intent(In) :: n
    Real(pr), Dimension(1:n), Intent(In) :: x,y
    Real(pr), Dimension(1:n), Intent(Inout) :: b,c,d

    Integer(ipr) :: nm1,i,ib
    Real(pr) :: t

    ! check input for consistency
    If(n < 2) Stop '-n < 2 in csplin call--'
    nm1 = n-1
    Do i = 1, nm1
       If(x(i) >= x(i+1)) Stop 'x not strictly ascending in csplin call'
    End Do
    If(n /= 2) Then
       ! set up tridiagonal system
       ! b = diagonal, d = offdiagonal, c = right hand side.
       d(1) = x(2) - x(1); c(2) = (y(2) - y(1))/d(1)
       Do i = 2, nm1
          d(i) = x(i+1) - x(i); b(i) = two*(d(i-1) + d(i))
          c(i+1) = (y(i+1) - y(i))/d(i); c(i) = c(i+1) - c(i)
       End Do
       ! end conditions.  third derivatives at  x(1)  and  x(n)
       ! obtained from divided dIfferences
       b(1) = -d(1); b(n) = -d(n-1); c(1) = zero; c(n) = zero
       If(n /= 3) Then
          c(1) =  c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
          c(n) =  c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
          c(1) =  c(1)*d(1)**2/(x(4)-x(1))
          c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
          ! forward elimination
       Else
          Do i = 2, n
             t = d(i-1)/b(i-1); b(i) = b(i) - t*d(i-1); c(i) = c(i) - t*c(i-1)
          End Do
       End If
       ! back substitution
       c(n) = c(n)/b(n)
       Do ib = 1, nm1
          i = n-ib
          c(i) = (c(i) - d(i)*c(i+1))/b(i)
       End Do
       ! compute polynomial coefficients
       b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + two*c(n))
       Do i = 1, nm1
          b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + two*c(i))
          d(i) = (c(i+1) - c(i))/d(i); c(i) = three*c(i)
       End Do
       c(n) = three*c(n); d(n) = d(n-1)
       Return
    Else
       b(1) = (y(2)-y(1))/(x(2)-x(1)); c(1) = zero; d(1) = zero
       Return
    End If

  End Subroutine csplin
  !=======================================================================
  !> This routine evaluates the interpolating spline function at point x.
  !>
  !>  @param[in] n - the number of data points or knots (n >= 2)
  !>  @param[in] u - the point where to evaluate the interpolating function
  !>  @param[in] x - the abscissas of the knots in strictly increasing order
  !>  @param[in] y - the ordinates of the knots
  !>  @param[in] b - arrays of spline coefficients as defined above
  !>  @param[in] c - arrays of spline coefficients as defined above
  !>  @param[in] d - arrays of spline coefficients as defined above
  !>  @param[inout] splf0 - value of the interpolated function
  !>  @param[inout] i - index tracking the nearest point x(i)
  !=======================================================================
  Subroutine cseval(n,u,x,y,b,c,d,splf0,i)
    Integer(ipr), Intent(In) :: n
    Integer(ipr), Intent(Inout), Optional :: i
    Real(pr), Intent(In) :: x(n),y(n),b(n),c(n),d(n),u
    Real(pr), Intent(Inout) :: splf0

    Integer(ipr) :: j,k
    Real(pr) :: dx,xip1,xi

    If(.Not. Present(i)) i = 1

    j = 0
    If(i >= n) i = 1
    xi=x(i); xip1=x(i+1)

    If(u < xi .Or. (u >= xi .And. u > xip1)) Then
       i = 1
       j = n + 1
       Do While (j > i+1)
          k = (i+j)/2
          If(u <  x(k)) j = k
          If(u >= x(k)) i = k
       End Do
       ! evaluate splf0
       dx = u - x(i)
       splf0 = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    Else
       dx = u - x(i)
       splf0 = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    End If

  End Subroutine cseval
  !=======================================================================
  !> This routine computes the first derivative of 'f1' if the step is 'h'
  !=======================================================================
  Subroutine deri(h,n,f1,dunl)
    Integer(ipr), Intent(In) :: n !< - Number of points in the function
    Real(pr), Intent(In) :: h !< - Step size of the grid
    Real(pr), Dimension(1:n), Intent(In) :: f1 !< - Vector containing all values of the function
    Real(pr), Dimension(1:n), Intent(Inout) :: dunl !< - Vector contanining all values of the derivative

    Integer(ipr) :: k
    Real(pr) :: t60,t12

    t60 =1.0_pr/(h*60.0_pr); t12 =1.0_pr/(h*12.0_pr)

    dunl(1)  =(eight*f1(2)-f1(3)+f1(1))*t12
    dunl(2)  =(45.0_pr*(f1(3)-f1(1))-nine*f1(4)+f1(5)-f1(1))*t60
    dunl(3)  =(45.0_pr*(f1(4)-f1(2))-nine*(f1(5)-f1(1))+f1(6))*t60
    dunl(n)  =(-eight*f1(n-1)+f1(n)+f1(n-2))*t12
    dunl(n-1)=(45.0_pr*(f1(n)-f1(n-2))+nine*f1(n-3)-f1(n)-f1(n-4))*t60
    dunl(n-2)=(45.0_pr*(f1(n-1)-f1(n-3))-nine*(f1(n)-f1(n-4))-f1(n-5))*t60
    Do k=4,n-3
       dunl(k) =(45.0_pr*(f1(k+1)-f1(k-1))-nine*(f1(k+2)-f1(k-2))+f1(k+3)-f1(k-3))*t60
    End Do

  End Subroutine deri
  !=======================================================================
  !> This routine calculates sign, Sqrt, factorials, etc. of integers and half integers
  !>   - \f$ \mathtt{iv(n) } = (-1)^n    \f$
  !>   - \f$ \mathtt{sq(n) } = \sqrt{n}, \f$
  !>   - \f$ \mathtt{sqi(n)} = 1/\sqrt{n}\f$
  !>   - \f$ \mathtt{fak(n)} = n!        \f$
  !>   - \f$ \mathtt{fi(n)} = 1/n!        \f$
  !>   - \f$ \mathtt{wf(n) } = \sqrt{n!} \f$
  !>   - \f$ \mathtt{wfi(n)} = 1/\sqrt{n!}\f$
  !=======================================================================
  Subroutine gfv()
    Integer(ipr), Parameter :: igfv=170 !maximal number for GFV
    Integer(ipr) :: i

    If(Allocated(iv)) Deallocate(iv,fak,fi,sq,sqi,wf,wfi)
    Allocate(iv(-igfv:igfv),fak(0:igfv),fi(0:igfv),sq(0:igfv),sqi(0:igfv))
    Allocate(wf(0:igfv),wfi(0:igfv))

    iv(0)=1; sq(0)=zero; sqi(0)=1.0e30_pr
    fak(0)=one; fi(0)=one; wf(0)=one; wfi(0)=one
    Do i=1,igfv
       iv(i)=-iv(i-1)
       iv(-i) = iv(i)
       sq(i)=Sqrt(Real(i,Kind=pr)); sqi(i)=one/sq(i)
       fak(i)=Real(i,Kind=pr)*fak(i-1); fi(i)=one/fak(i)
       wf(i)=sq(i)*wf(i-1); wfi(i)=one/wf(i)
    End Do

  End Subroutine gfv
  !======================================================================================
  !> This function calculates the small Wigner function \f$ d_{mm'}^{j}(\beta) \f$ for an arbitrary
  !> value of the angle \f$ \beta \f$
  !>  @param[in] j - angular momentum j
  !>  @param[in] mm1 - Quantum number m
  !>  @param[in] mm2 - Quantum number m'
  !>  @param[in] beta - Euler angle \f$ \beta \f$
  !======================================================================================
  Real(pr) Function wigner(j,mm1,mm2,beta)
    Integer(ipr), Intent(In) :: j,mm1,mm2
    Real(pr), Intent(In) :: beta

    Integer(ipr) :: m,imax,i,i1,i2,i3
    Real(pr) :: fac1,fac2,sum1

    m=Max(mm1,mm2)
    imax=j+mm2

    sum1=zero
    Do i=0,imax
       i1=j+mm2-i
       i2=mm1-mm2+i
       i3=j-mm1-i
       If(i1 >= 0 .And. i2 >= 0 .And. i3 >= 0) Then
          fac1=fi(i)*wf(j+mm1)*fi(i1)*wf(j-mm1); fac2=fi(i2)*wf(j+mm2)*fi(i3)*wf(j-mm2)
          sum1=sum1+iv(i1)*fac1*fac2*(Cos(half*beta))**(2*j+mm2-mm1-2*i)*(Sin(half*beta))**(mm1-mm2+2*i)
       End If
    End Do
    wigner=sum1*iv(j)

  End Function wigner
  !=============================================================================
  !> This subroutine calculates the integral of a complex-valued function using Simspon rule.
  !=============================================================================
  Subroutine integrate_complex_function(complex_function,n,h,complex_integral)
    Integer(ipr), Intent(In) :: n
    Real(pr), Intent(In) :: h
    Complex(pr), Dimension(1:n) :: complex_function
    Complex(pr), Intent(Inout) :: complex_integral

    Real(pr) :: real_integral,imag_integral
    Real(pr), Dimension(1:n) :: real_function,imag_function

    real_function=Dble(complex_function); imag_function=Aimag(complex_function)
    Call integrate_simpson(real_function,n,h,real_integral)
    Call integrate_simpson(imag_function,n,h,imag_integral)

    complex_integral=Cmplx(real_integral,imag_integral,kind=pr)

  End Subroutine integrate_complex_function
  !=============================================================================
  !> This subroutine performs Simpson's rule integration of a real-valued function
  !> defined by an array of equidistant values.
  !=============================================================================
  Subroutine integrate_simpson(f,n,h,r)
    Integer(ipr), Intent(In) :: n !< - Number of points
    Real(pr), Intent(In) :: h !< - Uniform spacing between function values
    Real(pr), Dimension(1:n), Intent(In) :: f !< - Array of values of the function
    Real(pr), Intent(Inout) :: r !< - Estimate of the integral that is returned to caller

    Integer(ipr) :: npanels,nbegin,nend,nn
    Real(pr) :: x
    Real(pr), Parameter :: c3d8=0.3750_pr

    ! No integration for n=1
    If(n == 1) Then
       r=f(1)
       Return
    End If

    npanels=n-1
    nbegin=1
    r=zero
    If(Mod(npanels,2) /= 0) Then
    ! Number of panels is odd. Use 3/8 rule on the first three panels.
       r = h*c3d8*(f(1) + 3.0_pr*(f(2)+f(3)) + f(4))
       If(n == 3) Return
       nbegin = 4
    Endif
    ! Apply 1/3 rule for the remaining panels.
    r = r + h*p13*(f(nbegin)+4_pr*f(nbegin+1)+f(n))
    nbegin = nbegin+2
    If(nbegin == n) Return
    x = zero
    nend = n-2
    Do nn=nbegin,nend,2
       x = x + f(nn) + two*f(nn+1)
    End Do
    r = r + h*two*p13*x

  End Subroutine integrate_simpson
  !====================================================================
  !> This function computes the probability integral at point z > 0
  !>   \f[
  !>      \Phi(z) = \frac{2}{\sqrt{\pi}} \int_{0}^{z} e^{-t^2}dt
  !>   \f]
  !>  by using the Simpson 3/8 integration rule.
  !====================================================================
  Real(pr) Function PROINT(ZPOINT)
    Real(pr), Intent(In) :: ZPOINT !< - Value of the z argument

    Integer(ipr) :: i, Nx
    Real(pr) :: dx,x,VALRES
    Real(pr), Allocatable :: FUNCTI(:)

    Nx = 100000
    Allocate(FUNCTI(1:Nx))

    dx = Abs(ZPOINT/Real(Nx-1,Kind=pr))
    Do i=1,Nx
       x=dx*Real(i-1,Kind=pr); FUNCTI(i)=Exp(-x**2)
    End Do

    VALRES=zero; Call integrate_simpson(FUNCTI,Nx,dx,VALRES)

    Deallocate(FUNCTI)

    PROINT=two*VALRES/Sqrt(pi)

  End Function PROINT
  !====================================================================
  !> This subroutine computes the value of the associated Legendre
  !> polynomial \f$ P_{\ell m}(x) \f$
  !====================================================================
  Real(pr) Function DEFLEG(LAMACT,MIUACT,ZVALUE)
    Integer(ipr), Intent(In) :: LAMACT !< - \f$ \ell \f$ value
    Integer(ipr), Intent(In) :: MIUACT !< - m value
    Real(pr), Intent(In) :: ZVALUE !< - x value

    Integer(ipr) :: L,M
    Real(pr) :: ZVALPP,ZVAL_P,ZPOLYN,TMPVAL,epsilon,one

    ! Argument z of polynomial P_{l,m}(z) must be lower than 1
    epsilon=1.0e-14_pr; one=1.0_pr
    If(Abs(Abs(ZVALUE)-one) <= epsilon) Then
       Write(*,'("ZVALUE=",f20.16," one=",f20.16," epsilon=",f20.16)') ZVALUE,one,epsilon
       Stop 'Error in DEFLEG - ARGUMENT |z| = 1'
    End If
    ! Initialization
    ZVALPP = one; ZVAL_P = ZVALUE
    ! P{0,0}(z) = 1
    If(LAMACT == 0) Then
       DEFLEG=one; Return
    End If
    ! P{1,0}(z) = 1
    If(LAMACT == 1 .And. MIUACT == 0) Then
       DEFLEG=ZVALUE; Return
    End If
    ! P{1,1}(z) = 1
    If(LAMACT == 1 .And. MIUACT == 1) Then
       DEFLEG=Sqrt(one-ZVALUE**2); Return
    End If
    !
    ! Obtaining P_{l,0}(z) by recurrence
    If(LAMACT > 1) Then

       Do L = 1,LAMACT-1
          ZPOLYN = (Real(2*L+1,Kind=pr)*ZVALUE*ZVAL_P - Real(L,Kind=pr)*ZVALPP) /Real(L+1,Kind=pr)
          TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
       End Do

       If(MIUACT == 0) Then
          DEFLEG=ZVAL_P
       Else
          ! Obtaining P_{l,1} from P_{l,0}(z) and P_{l-1,0}(z)
          !   - ZVAL_P contains P_{l,0}(z)
          !   - ZVALPP contains P_{l-1,0}(z)
          ZPOLYN =( Real(LAMACT,Kind=pr)*ZVALUE*ZVAL_P - Real(LAMACT,Kind=pr)*ZVALPP)/Sqrt(one - ZVALUE**2)
          TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
          ! Obtaining P_{l,m+1} from P_{l,m}(z) and P_{l,m-1}(z)
          !   - ZVAL_P contains P_{l,1}(z)
          !   - ZVALPP contains P_{l,0}(z)
          Do M = 1,MIUACT-1
             ZPOLYN =-2.0_pr*Real(M,Kind=pr)*ZVALUE*ZVAL_P/Sqrt(one - ZVALUE**2) &
                           - Real((LAMACT+M)*(LAMACT-M+1),Kind=pr)*ZVALPP
             TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
          End Do
          DEFLEG=(-1)**MIUACT * ZVAL_P
       End If

    End If

  End Function DEFLEG
  !====================================================================
  !> This function implements the Brent method to find the root of a function and returns the
  !> estimate of the root
  !====================================================================
  Real(pr) Function ZBRENT(FONCTI,XLOWER,XUPPER,TOLERA,IFOUND)
    Real(pr), Intent(In) :: XLOWER !< - Lower bound for the root
    Real(pr), Intent(In) :: XUPPER !< - Upper bound for the root
    Real(pr), Intent(In) :: TOLERA !< - Numerical precision for the Brent method
    Integer(ipr), Intent(Inout) :: IFOUND !< - Integer equal to 1 if minimum found

    Integer(ipr), Parameter :: ITRMAX=100
    Integer(ipr) :: ITERAT
    Real(pr), Parameter :: EPSCPU=3.D-8
    Real(pr) :: A,B,C,D,E,FA,FB,FC,P,Q,R,S,EPSILO,XM

    Interface
       Real(pr) Function FONCTI(X)
         Use HFBTHO_utilities
         Use HFBTHO, Only : bz
         Implicit None
         Real(pr), Intent(In) :: X
       End Function FONCTI
    End Interface

    IFOUND=1; A=XLOWER; B=XUPPER; FA=FONCTI(A); FB=FONCTI(B)
    If((FA > 0.0_pr .And. FB > 0.0_pr) .Or. (FA < 0.0_pr .And. FB < 0.0_pr)) Then
       IFOUND=0
       Write(*,'("A = ",F20.14," FA = ",f20.14)') A,FA
       Write(*,'("B = ",F20.14," FB = ",f20.14)') B,FB
       Write(*,'("ROOT MUST BE BRACKETED FOR ZBRENT")')
       Return
    End If

    C=B; FC=FB

    Do ITERAT=1,ITRMAX

       If((FB > 0.0_pr .And. FC > 0.0_pr) .Or. (FB < 0.0_pr .And. FC < 0.0_pr)) Then
         C=A; FC=FA; D=B-A; E=D
       End If

       If(Abs(FC) < Abs(FB)) Then
         A=B; B=C; C=A; FA=FB; FB=FC; FC=FA
       End If

       EPSILO=2.0_pr*EPSCPU*Abs(B) + 0.5_pr*TOLERA; XM=0.5_pr*(C - B)

       If(Abs(XM) <= EPSILO .Or. Abs(FB) <= EPSCPU) Then
          ZBRENT=B; Return
       End If

       If(Abs(E) >= EPSILO .And. Abs(FA) >= Abs(FB)) Then

           S = FB/FA

           If(Abs(A-C) <= EPSCPU) Then
              P=2.0_pr*XM*S
              Q=1.0_pr - S
           Else
              ! Attempt inverse quadratic interpolation
              Q=FA/FC; R=FB/FC
              P=S*(2.0_pr*XM*Q*(Q - R) - (B - A)*(R - 1.0_pr))
              Q=(Q - 1.0_pr)*(R - 1.0_pr)*(S - 1.0_pr)
           End If

           ! Check whether in bounds
           If (P > 0.0_pr) Q = -Q

           P=Abs(P)

           ! Test quality of interpolation. If too bad, switch to bisection method
           If(2.0_pr*P  <  Min(3.0_pr*XM*Q-Abs(EPSILO*Q),Abs(E*Q))) Then
              E=D; D=P/Q
           Else
              D=XM; E=D
           End If

       Else ! Bounds decreasing too slowly, use bisection

          D=XM; E=D

       End If

       ! Move last best guess to A
       A=B; FA=FB

       ! Evaluate new trial root
       If(Abs(D)  >  EPSILO) Then
          B=B+D
       Else
          B=B+Sign(EPSILO,XM)
       End If

       FB=FONCTI(B)

    End Do

    IFOUND=2
    ZBRENT=B

  End Function ZBRENT
  !====================================================================
  !> This function returns the cubic spline interpolation at point x. It duplicates routine
  !> \ref cseval().
  !====================================================================
  Real(pr) Function SPLINT(XENTRY,YENTRY,Y2A,NPOINT,XARGUM)
    Integer(ipr), Intent(In) :: NPOINT !< - Number of data points N
    Real(pr), Intent(In) :: XARGUM !< - Point x at which the function f is interpolated
    Real(pr), Dimension(1:NPOINT), Intent(In) :: XENTRY !< - Vector of size 1:N of data alues \f$ x_{k} \f$
    Real(pr), Dimension(1:NPOINT), Intent(In) :: YENTRY !< - Vector of size 1:N of data values \f$ y_{k} = f(x_{k}) \f$
    Real(pr), Dimension(1:NPOINT), Intent(In) :: Y2A !< - Vector of size 1:N containing the second derivatives of f at point x

    Integer(ipr) :: KLOWER, KHIGHR, K
    Real(pr) :: HLNGTH,A,B

    KLOWER=1; KHIGHR=NPOINT
    Do While((KHIGHR-KLOWER) > 1)
       K=(KHIGHR+KLOWER)/2
       If(XENTRY(K) > XARGUM) Then
          KHIGHR=K
       Else
          KLOWER=K
       End If
    End Do

    HLNGTH=XENTRY(KHIGHR)-XENTRY(KLOWER)

    If(HLNGTH == 0) Stop 'SPLINT01'

    A=(XENTRY(KHIGHR)-XARGUM)/HLNGTH
    B=(XARGUM-XENTRY(KLOWER))/HLNGTH

    SPLINT = A*YENTRY(KLOWER) + B*YENTRY(KHIGHR) + ( (A**3-A)*Y2A(KLOWER) + (B**3-B)*Y2A(KHIGHR) )*(HLNGTH**2)/6.0_pr

  End Function SPLINT
  !====================================================================
  !> This routine calculates the vector of second derivatives needed for the cubic spline interpolation
  !> of the function f. It duplicates routine \ref csplin()
  !====================================================================
  Subroutine SPLINE(X,Y,N,YP1,YPN,Y2,IERROR)
    Integer(ipr), Intent(In) :: N !< - Number of data points N
    Real(pr), Intent(In) :: YP1 !< - Value of the first derivative of f at the first point \f$ x_{1} \f$
    Real(pr), Intent(In) :: YPN !< - Value of the first derivative of f at the last point \f$ x_{N} \f$
    Real(pr), Dimension(1:N), Intent(In) :: X !< - Vector of size 1:N of data alues \f$ x_{k} \f$
    Real(pr), Dimension(1:N), Intent(In) :: Y !< - Vector of size 1:N of data values \f$ y_{k} = f(x_{k}) \f$
    Integer(ipr), Intent(Inout) :: IERROR !< - Error flag is non-zero if something wrong happened
    Real(pr), Dimension(1:N), Intent(Out) :: Y2 !< - Vector of size 1:N of second derivatives of the function f

    Integer(ipr), Parameter :: NDSPLN=200
    Integer(ipr) :: I,K
    Real(pr) :: SIG,P,QN,UN
    Real(pr), Dimension(1:NDSPLN) :: U

                            IERROR=0
    If(N < 4.Or.N > NDSPLN) IERROR=1

    If(YP1 > 0.99D+30) Then
       Y2(1)=0.0_pr
       U(1)=0.0_pr
    Else
       Y2(1)=-0.5_pr
       U(1)=(3.0_pr/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    End If

    Do I=2,N-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2.0_pr
       Y2(I)=(SIG-1.0_pr)/P
       U(I)=(6.0_pr*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    End Do

    If(YPN > 0.99D+30) Then
       QN=0.0_pr
       UN=0.0_pr
    Else
       QN=0.5_pr
       UN=(3.0_pr/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
    End If

    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1)

    Do K=N-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
    End Do

  End Subroutine SPLINE
  !=======================================================================
  !
  !=======================================================================
End Module math
