!         ___
!       //   ) )
!      //         ___      ___
!     //  ____  //___) ) //   ) )
!    //    / / //       //   / /
!   ((____/ / ((____   ((___/ /  MATERIALS
!
!    Copyright (c) :  Vahid Galavi
!
!    Main author(s):  Vahid Galavi
!
module ModMathLibrary
use ModGlobalConstants
use ModDebug

public

interface swap
 module procedure swap_double
 module procedure swap_integer
end interface swap

contains

!-----------------------------------------------------------------------
subroutine MatrixToVectorProd(xMat, Vec, n, VecR)
implicit none
integer, intent(IN):: n
double precision, dimension(n,n),  intent(IN):: xMat
double precision, dimension(n),    intent(IN):: Vec
double precision, dimension(n), intent(INOUT):: VecR

integer i,j
double precision x

do i=1,n
  X=0
  do j=1,n
    X=X+xMat(i,j)*Vec(j)
  end do
  VecR(i)=X
end do

!VecR = reshape(matmul(xMat, reshape(Vec, (/n,1/))), (/n/))

end subroutine MatrixToVectorProd

!-----------------------------------------------------------------------
subroutine sumVectors(Vec1, Vec2, R1, R2, n, VecR)
implicit none
integer, intent(IN):: n
double precision, intent(IN):: R1, R2
double precision, dimension(n), intent(IN):: Vec1, Vec2
double precision, dimension(n), intent(OUT):: VecR

VecR(1:n) = R1*Vec1(1:n) + R2*Vec2(1:n)

end subroutine sumVectors

!----------------------------------------------------------------------------------
subroutine calInverseMatrix(AInput, B, n)
implicit none
integer, intent(IN):: n
double precision, dimension(n,n), intent(IN) :: AInput
double precision, dimension(n,n), intent(OUT):: B

integer:: i, j, k, iPiv
double precision:: T, X
double precision, dimension(n,n):: A

A = AInput

B = 0.0d0
do i=1,n
  B(i,i) = 1d0
end do

do i=1,n
  T=A(i,i)
  iPiv=i
  do j=i+1,n
    if ( Abs(A(j,i)) > Abs(A(iPiv,i))  ) iPiv=j
  end do
  if (iPiv /= i) then
    do j=1,n
      x         = A( i  ,j)
      A( i  ,j) = A(iPiv,j)
      A(iPiv,j) = x
      x         = B( i  ,j)
      B( i  ,j) = B(iPiv,j)
      B(iPiv,j) = x
    end do
    T=A(i,i)
  end if
  do j=1,n
    A(i,j)=A(i,j)/(T + TINY)
    B(i,j)=B(i,j)/(T + TINY)
  end do
  do k=1,n
    if (k /= i) then
      T=A(k,i)
      do j=1,n
        A(k,j)=A(k,j)-T*A(i,j)
        B(k,j)=B(k,j)-T*B(i,j)
      end do
    end if
  end do
end do

end subroutine calInverseMatrix

!-----------------------------------------------------------------------
subroutine calPrincipalStressesAndVectors(stress,xN1,xN2,xN3,prSig)
implicit none
double precision, dimension(N_STRESS_VECTOR), intent(IN) :: stress
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(OUT):: xN1,xN2,xN3
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(OUT):: prSig

call calEigenValuesAndVectors(stress, xN1, xN2, xN3, prSig)

end subroutine calPrincipalStressesAndVectors

!-----------------------------------------------------------------------
subroutine StressTensorToVector(tensor, vector)
implicit none

double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(IN) :: tensor
double precision, dimension(N_STRESS_VECTOR), intent(OUT) :: vector

integer :: i
! Stress vector XX, YY, ZZ, XY, YZ, ZX

do i=1,3
  vector(i) = tensor(i,i)
end do

vector(4) = tensor(2,1)
vector(5) = tensor(3,2)
vector(6) = tensor(3,1)

end subroutine StressTensorToVector

!-----------------------------------------------------------------------
subroutine StressVectorToTensor(vector, tensor)
implicit none

double precision, dimension(N_STRESS_VECTOR), intent(IN):: vector
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(OUT) :: tensor

! Stress vector XX, YY, ZZ, XY, YZ, ZX

tensor(1,1) = vector(1) ! xx
tensor(1,2) = vector(4) ! xy = yx
tensor(1,3) = vector(6) ! zx = xz

tensor(2,1) = vector(4) ! xy = yx
tensor(2,2) = vector(2) ! yy
tensor(2,3) = vector(5) ! zy = yz

tensor(3,1) = vector(6) ! zx = xz
tensor(3,2) = vector(5) ! zy = yz
tensor(3,3) = vector(3) ! zz

end subroutine StressVectorToTensor

!-----------------------------------------------------------------------
subroutine StrainVectorToTensor(vector, tensor)
implicit none

double precision, dimension(N_STRESS_VECTOR), intent(IN) :: vector
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(OUT) :: tensor

! Strain vector XX, YY, ZZ, XY, YZ, ZX

tensor(1,1) = vector(1)     ! xx
tensor(1,2) = vector(4)/2d0 ! xy = yx
tensor(1,3) = vector(6)/2d0 ! zx = xz

tensor(2,1) = vector(4)/2d0 ! xy = yx
tensor(2,2) = vector(2)     ! yy
tensor(2,3) = vector(5)/2d0 ! zy = yz

tensor(3,1) = vector(6)/2d0 ! zx = xz
tensor(3,2) = vector(5)/2d0 ! zy = yz
tensor(3,3) = vector(3)     ! zz

end subroutine StrainVectorToTensor

!-----------------------------------------------------------------------
! Get Eigenvalues/Eigenvectors 3D
subroutine calEigenValuesAndVectors(stress,xN1,xN2,xN3,prSig)
implicit none

double precision, dimension(N_STRESS_VECTOR), intent(IN) :: stress
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(OUT):: xN1,xN2,xN3
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(OUT):: prSig
!
! Applied on principal stresses, directions
!
integer :: i, j, k, it, itmax, ip, iq, is1, is2, is3
double precision :: c, s, a1p, a2p, a3p, v1p, v2p, v3p, ap1, ap2, ap3
double precision :: absMaxS, tol, tau, sign_tau, t
double precision :: A(3,3),V(3,3)

call StressVectorToTensor(stress, A)

! Set V to unity matrix
V(1,1) = 1
V(2,1) = 0
V(3,1) = 0

V(1,2) = 0
V(2,2) = 1
V(3,2) = 0

V(1,3) = 0
V(2,3) = 0
V(3,3) = 1

absMaxS=0.0
do i=1,3
  do j=1,3
    if (abs(a(i,j)) > absMaxS) absMaxS=abs(a(i,j))
  end do
end do
tol = 1d-16 * absMaxS
it = 0
itmax = 50
do while ( it < itMax .and. abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) > tol )
  it=it+1
  do k=1,3
    if (k == 1) then
      ip=1
      iq=2
    Else if (k == 2) then
      ip=2
      iq=3
    Else
      ip=1
      iq=3
    end if
    if (abs(a(ip,iq)) > TINY) then
      tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
      if (tau .Ge.0.0) then
        sign_tau=1.0
      Else
        sign_tau=-1.0
      end if
      t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
      c=1.0/sqrt(1.0+t*t)
      s=t*c
      a1p=c*a(1,ip)-s*a(1,iq)
      a2p=c*a(2,ip)-s*a(2,iq)
      a3p=c*a(3,ip)-s*a(3,iq)
      a(1,iq)=s*a(1,ip)+c*a(1,iq)
      a(2,iq)=s*a(2,ip)+c*a(2,iq)
      a(3,iq)=s*a(3,ip)+c*a(3,iq)
      a(1,ip)=a1p
      a(2,ip)=a2p
      a(3,ip)=a3p

      v1p=c*v(1,ip)-s*v(1,iq)
      v2p=c*v(2,ip)-s*v(2,iq)
      v3p=c*v(3,ip)-s*v(3,iq)
      v(1,iq)=s*v(1,ip)+c*v(1,iq)
      v(2,iq)=s*v(2,ip)+c*v(2,iq)
      v(3,iq)=s*v(3,ip)+c*v(3,iq)
      v(1,ip)=v1p
      v(2,ip)=v2p
      v(3,ip)=v3p

      ap1=c*a(ip,1)-s*a(iq,1)
      ap2=c*a(ip,2)-s*a(iq,2)
      ap3=c*a(ip,3)-s*a(iq,3)
      a(iq,1)=s*a(ip,1)+c*a(iq,1)
      a(iq,2)=s*a(ip,2)+c*a(iq,2)
      a(iq,3)=s*a(ip,3)+c*a(iq,3)
      a(ip,1)=ap1
      a(ip,2)=ap2
      a(ip,3)=ap3
    end if ! a(ip,iq)<>0
  end do ! k
end do ! while
! principal values on diagonal of a

do i=1,3
  prSig(i) = a(i,i)
enddo

! Sort eigenvalues S1 <= S2 <= S3
is1 = 1
is2 = 2
is3 = 3
if (prSig(1) > prSig(2)) then
  t   = prSig(2)
  prSig(2)  = prSig(1)
  prSig(1)  = t
  it  = is2
  is2 = is1
  is1 = it
end if
if (prSig(2) > prSig(3)) then
  t   = prSig(3)
  prSig(3)  = prSig(2)
  prSig(2)  = t
  it  = is3
  is3 = is2
  is2 = it
end if
if (prSig(1) > prSig(2)) then
  t   = prSig(2)
  prSig(2)  = prSig(1)
  prSig(1)  = t
  it  = is2
  is2 = is1
  is1 = it
end if
do i=1,3
  xN1(i) = v(i,is1)
  xN2(i) = v(i,is2)
  xN3(i) = v(i,is3)
end do

end subroutine calEigenValuesAndVectors

!-----------------------------------------------------------------------
subroutine calPrincipalStresses(stress,prSig)
implicit none

double precision, dimension(N_STRESS_VECTOR), intent(IN) :: stress
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(OUT):: prSig

integer :: i, j, k, it, itmax, ip, iq
double precision :: c, s, a1p, a2p, a3p, ap1, ap2, ap3
double precision :: absMaxS, tol, tau, sign_tau, t
double precision :: A(3,3)
!
!
! Applied on principal stresses, directions
!
call StressVectorToTensor(stress, A)

absMaxS=0.0
do i=1,3
  do j=1,3
    if (abs(a(i,j)) > absMaxS) absMaxS=abs(a(i,j))
  end do
end do
tol = 1d-16 * absMaxS

it=0
itmax = 50
do while ( it < itmax .and. abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) > tol )

  it=it+1
  do k=1,3
    if (k == 1) then
      ip=1
      iq=2
    Else if (k  == 2) then
      ip=2
      iq=3
    Else
      ip=1
      iq=3
    end if

    if (abs(a(ip,iq)) > TINY) then
      tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
      if (tau .Ge.0.0) then
        sign_tau=1.0
      Else
        sign_tau=-1.0
      end if
      t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
      c=1.0/sqrt(1.0+t*t)
      s=t*c
      a1p=c*a(1,ip)-s*a(1,iq)
      a2p=c*a(2,ip)-s*a(2,iq)
      a3p=c*a(3,ip)-s*a(3,iq)
      a(1,iq)=s*a(1,ip)+c*a(1,iq)
      a(2,iq)=s*a(2,ip)+c*a(2,iq)
      a(3,iq)=s*a(3,ip)+c*a(3,iq)
      a(1,ip)=a1p
      a(2,ip)=a2p
      a(3,ip)=a3p

      ap1=c*a(ip,1)-s*a(iq,1)
      ap2=c*a(ip,2)-s*a(iq,2)
      ap3=c*a(ip,3)-s*a(iq,3)
      a(iq,1)=s*a(ip,1)+c*a(iq,1)
      a(iq,2)=s*a(ip,2)+c*a(iq,2)
      a(iq,3)=s*a(ip,3)+c*a(iq,3)
      a(ip,1)=ap1
      a(ip,2)=ap2
      a(ip,3)=ap3
    end if ! a(ip,iq)<>0
  end do ! k
end do ! while
! principal values on diagonal of a
do i=1,3
  prSig(i) = a(i,i)
enddo

if (prSig(1) > prSig(2)) then
  t   = prSig(2)
  prSig(2)  = prSig(1)
  prSig(1)  = t
end if
if (prSig(2) > prSig(3)) then
  t   = prSig(3)
  prSig(3)  = prSig(2)
  prSig(2)  = t
end if
if (prSig(1) > prSig(2)) then
  t   = prSig(2)
  prSig(2)  = prSig(1)
  prSig(1)  = t
end if

end subroutine calPrincipalStresses

!-----------------------------------------------------------------------
! Returns the Cartesian stresses using the principal stresses
! and the principal directions
subroutine PrincipalStressToCartesian(prSig,xN1,xN2,xN3,stress)
implicit none
integer, parameter :: V_SIZE = N_PRINCIPAL_STRESS_VECTOR
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(IN) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(IN) :: xN1,xN2,xN3
double precision, dimension(N_STRESS_VECTOR), intent(OUT):: stress

integer :: i
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR):: StressTensor,T,TT,STT

! Fill transformation (rotation) matrix
do i=1,N_PRINCIPAL_STRESS_VECTOR
  T(i,1) = xN1(i)
  T(i,2) = xN2(i)
  T(i,3) = xN3(i)
  TT(1,i) = T(i,1)
  TT(2,i) = T(i,2)
  TT(3,i) = T(i,3)
end do

StressTensor = 0.0

do i=1,N_PRINCIPAL_STRESS_VECTOR
  StressTensor(i,i) = prSig(i)
enddo

! stress = T*StressTensor*TT
STT = matmul(StressTensor, TT)
StressTensor = matmul(T, STT)

call StressTensorToVector(StressTensor, stress)

end subroutine PrincipalStressToCartesian

!-----------------------------------------------------------------------
subroutine CalPrincipalStrainsAndVectors(strain,xN1,xN2,xN3,prStrain)
implicit none
double precision, dimension(N_STRESS_VECTOR), intent(IN):: strain
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(OUT) :: xN1,xN2,xN3
double precision, intent(OUT):: prStrain(3)

double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR) :: A,V

integer :: i, j, k, it, itmax, ip, iq
double precision :: c, s, a1p, a2p, a3p, v1p, v2p, v3p, ap1, ap2, ap3
double precision :: absMaxS, tol, tau, sign_tau, t

call StrainVectorToTensor(strain, A)

! Set V to unity matrix
V(1,1) = 1
V(2,1) = 0
V(3,1) = 0

V(1,2) = 0
V(2,2) = 1
V(3,2) = 0

V(1,3) = 0
V(2,3) = 0
V(3,3) = 1

absMaxS=0.0
do i=1,3
  do j=1,3
    if (abs(a(i,j)) > absMaxS) absMaxS=abs(a(i,j))
  end do
end do
tol = 1d-20 * absMaxS
it = 0
itmax = 50
do while ( it < itMax .and. abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) > tol )
  it=it+1
  do k=1,3
    if (k == 1) then
      ip=1
      iq=2
    Else if (k ==2) then
      ip=2
      iq=3
    Else
      ip=1
      iq=3
    end if
    if ( abs(a(ip,iq)) < 1d-50 ) a(ip,iq) = 0
    if (a(ip,iq)  /=  0.0) then         ! ongelijk nul ?
      tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
      if (tau .Ge.0.0) then
        sign_tau=1.0
      Else
        sign_tau=-1.0
      end if
      t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
      c=1.0/sqrt(1.0+t*t)
      s=t*c
      a1p=c*a(1,ip)-s*a(1,iq)
      a2p=c*a(2,ip)-s*a(2,iq)
      a3p=c*a(3,ip)-s*a(3,iq)
      a(1,iq)=s*a(1,ip)+c*a(1,iq)
      a(2,iq)=s*a(2,ip)+c*a(2,iq)
      a(3,iq)=s*a(3,ip)+c*a(3,iq)
      a(1,ip)=a1p
      a(2,ip)=a2p
      a(3,ip)=a3p

      v1p=c*v(1,ip)-s*v(1,iq)
      v2p=c*v(2,ip)-s*v(2,iq)
      v3p=c*v(3,ip)-s*v(3,iq)
      v(1,iq)=s*v(1,ip)+c*v(1,iq)
      v(2,iq)=s*v(2,ip)+c*v(2,iq)
      v(3,iq)=s*v(3,ip)+c*v(3,iq)
      v(1,ip)=v1p
      v(2,ip)=v2p
      v(3,ip)=v3p

      ap1=c*a(ip,1)-s*a(iq,1)
      ap2=c*a(ip,2)-s*a(iq,2)
      ap3=c*a(ip,3)-s*a(iq,3)
      a(iq,1)=s*a(ip,1)+c*a(iq,1)
      a(iq,2)=s*a(ip,2)+c*a(iq,2)
      a(iq,3)=s*a(ip,3)+c*a(iq,3)
      a(ip,1)=ap1
      a(ip,2)=ap2
      a(ip,3)=ap3
    end if ! a(ip,iq)<>0
  end do ! k
end do ! while

! principal values on diagonal of a
do i=1,3
  prStrain(i) = a(i,i)
enddo

do i=1,3
  xN1(i) = v(i,1) ! first  column
  xN2(i) = v(i,2) ! second column
  xN3(i) = v(i,3) ! third  column
end do

end subroutine CalPrincipalStrainsAndVectors

!-----------------------------------------------------------------------
double precision function lengthVector(A,n) result(length)
implicit none
integer, intent(in) :: n
double precision, dimension(n), intent(in)::A
integer i
double precision sum

sum = 0.0
do i=1,n
  sum = sum + A(i)*A(i)
enddo
length = sqrt(sum)

end function lengthVector

!-----------------------------------------------------------------------
subroutine swap_double(a,b)
implicit none
double precision, intent(inout) :: a,b
double precision :: temp

temp = a
a = b
b = temp

end subroutine swap_double

!-----------------------------------------------------------------------
subroutine swap_integer(a,b)
implicit none
integer, intent(inout) :: a,b
integer :: temp

temp = a
a = b
b = temp

end subroutine swap_integer

end module ModMathLibrary
