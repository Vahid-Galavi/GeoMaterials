!         ___
!       //   ) )
!      //         ___      ___
!     //  ____  //___) ) //   ) )
!    //    / / //       //   / /
!   ((____/ / ((____   ((___/ /  MATERIALS
!
!
!
!    Main authors:    Vahid Galavi
!
module ModStressAndStrainLibrary
use ModGlobalConstants
use ModMathLibrary
use ModDebug
use ModString

implicit none

interface getdQdSigma
  module procedure getdQdSigmaFromPrStress
  module procedure getdQdSigmaFromPQPrStress
end interface getdQdSigma

interface getLodeAngle
  module procedure getLodeAngleJ
  module procedure getLodeAngleS
end interface getLodeAngle

!> Stress indices
integer, parameter :: INDEX_XX = 1
integer, parameter :: INDEX_YY = 2
integer, parameter :: INDEX_ZZ = 3
integer, parameter :: INDEX_XY = 4
integer, parameter :: INDEX_YZ = 5
integer, parameter :: INDEX_ZX = 6

!> principal stress indices
integer, parameter :: INDEX_S1 = 1
integer, parameter :: INDEX_S2 = 2
integer, parameter :: INDEX_S3 = 3

public

contains
!--------------------------------------------------------------------------------------------------
subroutine getdEps3dEps(dEps3dEps)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dEps3dEps

dEps3dEps    = 0.0d0
dEps3dEps(3) = 1.0d0

end subroutine getdEps3dEps

!--------------------------------------------------------------------------------------------------
subroutine getdEpsVdEps(dEpsVdEps)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dEpsVdEps

dEpsVdEps = 1.0d0

end subroutine getdEpsVdEps

!--------------------------------------------------------------------------------------------------
! get derivative Von Mises strain with respect to principal strains:
! eps_q = (sqrt(2)/2) * eps_oct
! eps_q = (sqrt(2)/3) * sqrt[(eps_1 - eps_2)^2 + (eps_1 - eps_3)^2 + (eps_3 - eps_2)^2]
subroutine getdEpsQdEps(eps, dEpsQdEps)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)::  eps
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dEpsQdEps

double precision, parameter :: factor = 2.0d0 / 9.0d0
double precision :: epsQ, epsV

epsV = getEpsV(eps)
epsQ = getEpsQPr(eps)

if (epsQ > 0.0d0) then
  dEpsQdEps(:) = (factor / epsQ) * (3.0d0*eps(:) - epsV)
else
  dEpsQdEps    = 0.0d0
endif

end subroutine getdEpsQdEps

!--------------------------------------------------------------------------------------------------
! get derivative deviatoric strain, defined in HS model, with respect to principal strains:
! gamma = eps_1 - eps_2 - eps_3
subroutine getdGammadEps(dEpsQdEps)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dEpsQdEps

dEpsQdEps(1) =  1
dEpsQdEps(2) = -1
dEpsQdEps(3) = -1

end subroutine getdGammadEps

!--------------------------------------------------------------------------------------------------
! get derivative deviatoric strain, defined in HS model, with respect to principal strains:
! eps_q = (3/2) * (sqrt(2)/2) * eps_oct
! eps_q = sqrt{0.5 * [(eps_1 - eps_2)^2 + (eps_1 - eps_3)^2 + (eps_3 - eps_2)^2] }
subroutine getdEpsGammaSdEps(eps, dEpsQdEps)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)::  eps
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dEpsQdEps

double precision, parameter :: factor = 3.0d0 / 2.0d0

call getdEpsQdEps(eps, dEpsQdEps)

dEpsQdEps = factor * dEpsQdEps

end subroutine getdEpsGammaSdEps

!--------------------------------------------------------------------------------------------------
subroutine getdPdSigma(dPdSigma)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dPdSigma

dPdSigma = 1.0d0/3.0d0

end subroutine getdPdSigma

!--------------------------------------------------------------------------------------------------
subroutine getdQdSigmaFromPQPrStress(p, q, prSig, dqdSigma)
implicit none
double precision, intent(in):: p, q
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)::  prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dqdSigma

double precision, parameter :: factor = 3.0d0 / 2.0d0

dqdSigma(:) = (factor/(q+TINY)) * (prSig(:) - p)

end subroutine getdQdSigmaFromPQPrStress

!--------------------------------------------------------------------------------------------------
subroutine getdQdSigmaFromPrStress(prSig, dqdSigma)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)::  prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dqdSigma

double precision:: p, q

p = calP(prSig)
q = calQ(prSig)

call getdQdSigmaFromPQPrStress(p, q, prSig, dqdSigma)

end subroutine getdQdSigmaFromPrStress

!--------------------------------------------------------------------------------------------------
subroutine getdSigmaBardSigma(J2, S, dSigmaBardSigma)
implicit none
double precision, intent(in):: J2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dSigmaBardSigma

double precision :: factor
double precision:: dJ2dSigma(N_PRINCIPAL_STRESS_VECTOR)

call getdJ2dSigma(S, dJ2dSigma)

factor = 0.5d0 / (sqrt(J2) + TINY)
dSigmaBardSigma = factor * dJ2dSigma

end subroutine getdSigmaBardSigma

!--------------------------------------------------------------------------------------------------
subroutine getdJ2dSigma(S, dJ2dSigma)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dJ2dSigma

integer :: i, j
double precision :: sum, factor

do i=1,N_PRINCIPAL_STRESS_VECTOR
  sum = 0.0
  do j=1,N_PRINCIPAL_STRESS_VECTOR
    if (i==j) then
      factor = 2.0d0/3.0d0
    else
      factor = -1.0d0/3.0d0
    endif
    sum = sum + factor* S(j)
  enddo
  dJ2dSigma(i) = sum
enddo

end subroutine getdJ2dSigma

!--------------------------------------------------------------------------------------------------
subroutine getdJ3dSigma(S, dJ3dSigma)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dJ3dSigma

dJ3dSigma(1) = S(2)*S(3)
dJ3dSigma(2) = S(1)*S(3)
dJ3dSigma(3) = S(2)*S(1)

end subroutine getdJ3dSigma

!--------------------------------------------------------------------------------------------------
! lode angle derivative
subroutine getdSin3ThetadSigma(S, J2, J3, dSin3ThetadSigma)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S
double precision, intent(in):: J2, J3
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dSin3ThetadSigma

double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dJ3dSigma, dJ2dSigma
double precision:: factor, J23

factor = - 3.0d0 * sqrt(3.0d0) * 0.5d0
call getdJ2dSigma(S, dJ2dSigma)
call getdJ3dSigma(S, dJ3dSigma)

J23 = J2*J2*J2

if (J23 > 0.0d0) then
  dSin3ThetadSigma = factor * (1.0d0/ J23) *(sqrt(J23) * dJ3dSigma - sqrt(J2)*dJ2dSigma*J3)
else
  dSin3ThetadSigma = 0.0d0
endif


end subroutine getdSin3ThetadSigma

!--------------------------------------------------------------------------------------------------
subroutine sortPrincipalStress(prSig)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(inout) :: prSig

if (prSig(1) < prSig(2)) call swap(prSig(1), prSig(2))
if (prSig(1) < prSig(3)) call swap(prSig(1), prSig(3))
if (prSig(2) < prSig(3)) call swap(prSig(2), prSig(3))

end subroutine sortPrincipalStress

!--------------------------------------------------------------------------------------------------
double precision function getEpsV(eps) result (epsV)
implicit none
double precision, intent(in):: eps(*)

epsV = eps(1) + eps(2) + eps(3)

end function getEpsV

!--------------------------------------------------------------------------------------------------
double precision function getEpsq(eps) result (epsq)
implicit none
double precision, dimension(N_STRESS_VECTOR), intent(in)::  eps

epsq = (sqrt(2.0d0)/3.0d0)*sqrt(    (eps(1)-eps(2))**2 &
                                +   (eps(1)-eps(3))**2 &
                                +   (eps(2)-eps(3))**2 &
                                + 6*(eps(4)       )**2 &
                                + 6*(eps(5)       )**2 &
                                + 6*(eps(6)       )**2 )

end function getEpsq

!--------------------------------------------------------------------------------------------------
subroutine calEpsDev(eps, epsD)
implicit none
double precision, dimension(N_STRESS_VECTOR), intent(in)::  eps
double precision, dimension(N_STRESS_VECTOR), intent(out):: epsD

double precision:: epsVol

epsVol  = getEpsV(eps)
epsD(INDEX_XX) = eps(INDEX_XX) - epsVol/3d0
epsD(INDEX_YY) = eps(INDEX_YY) - epsVol/3d0
epsD(INDEX_ZZ) = eps(INDEX_ZZ) - epsVol/3d0
epsD(INDEX_XY) = eps(INDEX_XY)
epsD(INDEX_YZ) = eps(INDEX_YZ)
epsD(INDEX_ZX) = eps(INDEX_ZX)

end subroutine calEpsDev

!--------------------------------------------------------------------------------------------------
double precision function getEpsQPr(eps) result (epsq)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: eps

epsq = (sqrt(2.0d0)/3.0d0)*sqrt(  (eps(1)-eps(2))**2 &
                                + (eps(1)-eps(3))**2 &
                                + (eps(2)-eps(3))**2 )

end function getEpsQPr

!--------------------------------------------------------------------------------------------------
double precision function getPrStressNorm(prSig) result (sigNorm)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

sigNorm = sqrt( prSig(1)**2 + prSig(2)**2 + prSig(3)**2 )

end function getPrStressNorm

!--------------------------------------------------------------------------------------------------
double precision function getI1(prSig) result (aI)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

aI = (prSig(1) + prSig(2) + prSig(3))

end function getI1

!--------------------------------------------------------------------------------------------------
double precision function getI2(prSig) result (aI)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

aI =  prSig(1) * prSig(2)               &
    +            prSig(2) * prSig(3) &
    + prSig(1)            * prSig(3)

end function getI2

!--------------------------------------------------------------------------------------------------
double precision function getI3(prSig) result (aI)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

aI = (prSig(1) * prSig(2) * prSig(3))

end function getI3

!--------------------------------------------------------------------------------------------------
double precision function getP2D(prSig) result(p)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

p = (getSigma1(prSig) + getSigma3(prSig))*0.5d0

end function getP2D

!--------------------------------------------------------------------------------------------------
double precision function calP(prSig) result(p)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

p = (prSig(1) + prSig(2) + prSig(3))/3.0d0

end function calP

!--------------------------------------------------------------------------------------------------
double precision function getStressRatio(prSig) result(eta)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

double precision p, q

p = calP(prSig)
q = calQ(prSig)
eta = q / (p + TINY)

end function getStressRatio

!--------------------------------------------------------------------------------------------------
double precision function getSigma3(prSig) result(s3)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

s3 = prSig(3) !minVal(prSig)

end function getSigma3
!--------------------------------------------------------------------------------------------------
double precision function getSigma1(prSig) result(s1)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

s1 = prSig(1) !maxVal(prSig)

end function getSigma1

!--------------------------------------------------------------------------------------------------
double precision function calQ(prSig) result(q)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

q = sqrt( 0.5*(  (prSig(1)-prSig(2))**2  &
               + (prSig(2)-prSig(3))**2  &
               + (prSig(1)-prSig(3))**2) )

end function calQ

!--------------------------------------------------------------------------------------------------
double precision function getJ2(S) result (J2)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S

!J2 = -(s(1)*s(2) + s(2)*s(3) +s(1)*s(3))
! or
J2 = 0.5*(s(1)**2 + s(2)**2+ s(3)**2)

end function getJ2

!--------------------------------------------------------------------------------------------------
subroutine getPrincipalDeviatoricS(prSig, S)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)::  prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: S

double precision :: p

p = calP(prSig)
s(:) = prSig(:) - p

end subroutine getPrincipalDeviatoricS

!--------------------------------------------------------------------------------------------------
subroutine getPrincipalDeviatoricSWithP(p, prSig, S)
implicit none
double precision, intent(in):: p
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)::  prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: S

s(:) = prSig(:) - p

end subroutine getPrincipalDeviatoricSWithP

!--------------------------------------------------------------------------------------------------
double precision function getJ3(S) result (J3)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S

J3 = s(1) * s(2) * s(3)

end function getJ3

!--------------------------------------------------------------------------------------------
subroutine getdQSigmaKinematic(J2, S, prBackStress, dqdSigma)
implicit none
double precision, intent(in):: J2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S, prBackStress
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dqdSigma

double precision, parameter :: factor = sqrt(3.0d0)
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: dSigmaBardSigma

call getdSigmaBardSigmaKinematic(J2, S, prBackStress, dSigmaBardSigma)
dqdSigma = factor * dSigmaBardSigma

end subroutine getdQSigmaKinematic

!--------------------------------------------------------------------------------------------
subroutine getdSigmaBardSigmaKinematic(J2, S, prBackStress, dSigmaBardSigma)
implicit none
double precision, intent(in):: J2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S, prBackStress
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dSigmaBardSigma

double precision:: factor
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dJ2dSigma

call getdJ2dSigmaKinematic(S, prBackStress, dJ2dSigma)

factor = 0.5d0 / (sqrt(J2) + TINY)
dSigmaBardSigma = factor * dJ2dSigma

end subroutine getdSigmaBardSigmaKinematic

!--------------------------------------------------------------------------------------------
subroutine getdJ2dSigmaKinematic(S, prBackStress, dJ2dSigma)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: S, prBackStress
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dJ2dSigma

integer :: i, j
double precision :: sum, factor

dJ2dSigma(:) = factor * S(:) * (2.0d0 - prBackStress(:))

do i=1,N_PRINCIPAL_STRESS_VECTOR
  sum = 0.0
  do j=1,N_PRINCIPAL_STRESS_VECTOR
    if (i==j) then
      factor =  (2.0d0- prBackStress(j))/3.0d0
    else
      factor = -(1.0d0+prBackStress(j))/3.0d0
    endif
    sum = sum + factor* S(j)
  enddo
  dJ2dSigma(i) = sum
enddo

end subroutine getdJ2dSigmaKinematic

!--------------------------------------------------------------------------------------------------
double precision function getLodeAngleJ(J2, J3) result (theta)
implicit none
double precision, intent(in):: J2, J3

double precision:: factor, sin3Theta

factor = - 3.0d0 * sqrt(3.0d0) * 0.5d0
if (J2 > 0.0d0) then
  sin3Theta = factor * (J3 / sqrt(J2*J2*J2))
else
  sin3Theta = 0.0d0
endif

sin3Theta = min(sin3Theta,  1.0d0)
sin3Theta = max(sin3Theta, -1.0d0)

theta = (1.0d0/3.0d0) * asin(sin3Theta)

end function getLodeAngleJ

!----------------------------------------------------------------------------------------------------
double precision function getLodeAngleS(prSig) result (theta)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S
double precision ::J2, J3

call getPrincipalDeviatoricS(prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)
theta = getLodeAngleJ(J2,J3)

end function getLodeAngleS

!-----------------------------------------------------------------
double precision function getqTX(prSig) result(q)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

q = prSig(1) - prSig(3)

end function getqTX

!-----------------------------------------------------------------
subroutine getPrincipalDMatrix(G, xNu, D)
implicit none
double precision, intent(in)::  G, xNu
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(out):: D

integer i,j
double precision factor

factor = 2.0d0*G/(1.0-2.0*xNu)

do i=1,N_PRINCIPAL_STRESS_VECTOR
  do j=1,N_PRINCIPAL_STRESS_VECTOR
    if (i==j) then
      D(i,j) = factor*(1-xNu)
    else
      D(i,j) = factor*(  xNu)
    endif
  enddo
enddo

end subroutine getPrincipalDMatrix

!--------------------------------------------------------------------------------------------------
subroutine getPrincipalCMatrix(G, xNu, C)
implicit none
double precision, intent(in)::  G, xNu
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(out):: C

integer i,j
double precision factor

factor = 1.0 / (2.0d0*G*(1.0 + xNu))

do i=1,N_PRINCIPAL_STRESS_VECTOR
  do j=1,N_PRINCIPAL_STRESS_VECTOR
    if (i==j) then
      C(i,j) = factor*(1)
    else
      C(i,j) = factor*(-xNu)
    endif
  enddo
enddo

end subroutine getPrincipalCMatrix

!--------------------------------------------------------------------------------------------------
subroutine DMatrixDrained(G, xNu, D)
implicit none
double precision, intent(in)::  G, xNu
double precision, dimension(N_STRESS_VECTOR, N_STRESS_VECTOR), intent(out):: D

integer i,j
double precision fac1, fac2

fac1 = 2*G*(1-xNu)/(1-2*xNu)
fac2 = 2*G*( xNu )/(1-2*xNu)

D = 0.0

do i=1,3
  do j=1,3
    D(i,j) = fac2
  enddo
  D(i,i) = fac1
  D(i+3,i+3) = G
enddo

end subroutine DMatrixDrained
!--------------------------------------------------------------------------------------------------
subroutine AddBulkWToDMatrixDrained(bulkW, D)
implicit none
double precision, intent(in)::    bulkW
double precision, dimension(N_STRESS_VECTOR,N_STRESS_VECTOR), intent(inout):: D

integer i,j

do i=1,3
  do j=1,3
    D(i,j) = D(i,j) + bulkW
  enddo
enddo

end subroutine AddBulkWToDMatrixDrained

!--------------------------------------------------------------------------------------------------
double precision function prStressLength(vector) result(vectorLength)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: vector

vectorLength = lengthVector(vector, N_PRINCIPAL_STRESS_VECTOR)

end function prStressLength

end module ModStressAndStrainLibrary
