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
module ModSolversLibrary
use ModMathLibrary
use ModString
use ModDebug
use ModStressAndStrainLibrary
use ModGlobalConstants
use ModYieldSurfaceLibrary

implicit none

!-- solvers constants
! constant to convert input eta to mu parameter used in the formulation
double precision, parameter:: ETA_TO_MU_DP = 3.0d0

! loading regime
integer, parameter:: UNDEFINED_LOADING     = -1
integer, parameter:: ELASTIC_LOADING       = 0
integer, parameter:: FAILURE_LOADING       = 1
integer, parameter:: TENSILE_LOADING       = 2
integer, parameter:: CAP_LOADING           = 3
integer, parameter:: PRIMARY_CAP_LOADING   = 4
integer, parameter:: FRICTION_HARDENING_LOADING = 5
integer, parameter:: PRIMARY_TENSILE_LOADING  = 6
integer, parameter:: SECONDARY_LOADING        = 7
integer, parameter:: SECONDARY_CAP_LOADING    = 8

! loading redime (characters for debugging)
character(len =  8), parameter:: PRIMARY         = 'PRIMARY'
character(len =  9), parameter:: SECONDARY       = 'SECONDARY'
character(len =  7), parameter:: TENSILE         = 'TENSILE'
character(len =  3), parameter:: CAP             = 'CAP'
character(len = 11), parameter:: PRIMARY_CAP     = 'PRIMARY_CAP'
character(len = 15), parameter:: PRIMARY_TENSILE = 'PRIMARY_TENSILE'
character(len = 13), parameter:: SECONDARY_CAP   = 'SECONDARY_CAP'
character(len = 12), parameter:: APEX_TENSILE    = 'APEX_TENSILE'
character(len = 12), parameter:: MOHR_COULOMB    = 'MOHR_COULOMB'

! limits and constants
double precision, parameter:: LODE_ANGLE_TRIAXIAL_LIMIT   = 25.0d0 * RADIANS
double precision, parameter:: CRITICAL_SHEAR_STRESS_RATIO = 1.0d0
double precision, parameter:: MAX_SHEAR_STRESS_RATIO      = 1.0d0
double precision, parameter:: UNDRAINED_NU                = 0.495d0
double precision, parameter:: SIN_PHI_MAX =  0.999d0
double precision, parameter:: SIN_PHI_MIN = -0.999d0

! tolerated errors
double precision, parameter:: TOLERANCE_TOLER         = 0.0d0
double precision, parameter:: TOLERANCE_MN            = 1d-3
double precision, parameter:: TOLERANCE_MC            = 1d-11
double precision, parameter:: TOLERANCE_HS_MC         = 1d-11
double precision, parameter:: TOLERANCE_DP            = 1d-11
double precision, parameter:: TOLERANCE_TEN           = 1d-11
double precision, parameter:: TOLERANCE_CAP           = 1d-8
double precision, parameter:: TOLERANCE_COMPLEMENTARY = 1d-12
double precision, parameter:: TOLERANCE_CAM_CLAY      = 1d-5
double precision, parameter:: TOLERANCE_DEFAULT              = 1d-12
double precision, parameter:: THRESHOLD_HIGH_TOLERATED_ERROR = 5.0d0
double precision, parameter:: HIGH_TOLERANCE_FACTOR          = 1.0d0

! substepping limints and convergence
integer, parameter:: N_SUB_MAX_INITIAL = 10 !100
integer, parameter:: N_SUB_MAX = N_SUB_MAX_INITIAL * 10
integer, parameter:: N_SUB_FACTOR_CRITICAL = 5
integer, parameter:: UP_SCALE_FACTOR = 2
double precision, parameter:: MIN_EPS_VOLUMETRIC = 1d-2
double precision, parameter:: MIN_EPS_DEVIATORIC = 1d-4
double precision, parameter:: DIVERGENCE_LIMIT = 1d30
double precision, parameter:: MAX_SUB_SIG_DEVIATORIC_FACTOR = 5.0d0
double precision, parameter:: MAX_SUB_SIG_VOLUMETRIC_FACTOR = 1d-1

! yield/plastic potential functions
integer, parameter:: MOHR_COULOMB_YIELD              = 1
integer, parameter:: DRUCKER_PRAGER_YIELD            = 2
integer, parameter:: MATSUOKA_NAKAI_CONVEX_YIELD     = 3
integer, parameter:: MATSUOKA_NAKAI_YIELD            = 4
integer, parameter:: TENSION_YIELD                   = 5
integer, parameter:: COMPLEMENTARY_YIELD             = 6
integer, parameter:: DRUCKER_PRAGER_YIELD_EXTENSION  = 7
integer, parameter:: HARDENING_SOIL_MC_YIELD         = 8
integer, parameter:: CAP_HS_YIELD                    = 9
integer, parameter:: CAM_CLAY_YIELD                  = 10
integer, parameter:: DRUCKER_PRAGER_LODE_ANGLE_YIELD = 11

! overstress functions
integer, parameter:: BINGHAM_MODEL = 0
integer, parameter:: PERZYNA_MODEL = 1

! flow rules
integer, parameter:: FLOW_RULE_HARDENING_SOIL = 0
integer, parameter:: FLOW_RULE_HS_SMALL       = 1
integer, parameter:: FLOW_RULE_ROWE           = 2
integer, parameter:: FLOW_RULE_UBC            = 3
integer, parameter:: FLOW_RULE_WAN            = 4
integer, parameter:: FLOW_RULE_BAUER          = 5
integer, parameter:: FLOW_RULE_VARIABLE_BAUER = 6

integer, parameter:: N_CROSS_YIELD_SURFACES = 2

! IDTask for Plaxis User Defined Soil Models.
integer, parameter :: IDTASK_INITIALISATION = 1
integer, parameter :: IDTASK_STRESS         = 2
integer, parameter :: IDTASK_DEP            = 3
integer, parameter :: IDTASK_NSTAT          = 4
integer, parameter :: IDTASK_ATTRIBUTES     = 5
integer, parameter :: IDTASK_DE             = 6

contains

!--------------------------------------------------------------------------------------------------
double precision function CalculateBulkModulusWater(G, nu, nu_undrained) result(BulkW)
implicit none

double precision, intent(in) :: G, nu, nu_undrained

BulkW = (2.0/3.0) * G * ( (1.0d0 + nu_undrained) / (1.0 - 2.0d0*nu_undrained) &
                         -(1.0 +   nu          ) / (1.0 - 2.0 * nu         ) )

end function CalculateBulkModulusWater

!--------------------------------------------------------------------------------------------------
subroutine CalculateExcessPorePressure(G, nu, propUndrainedNu, &
                                       excessPP0, dEpsV, &
                                       BulkW, excessPP )
implicit none
double precision, intent(in):: G, nu, propUndrainedNu
double precision, intent(in):: excessPP0, dEpsV
double precision, intent(out):: BulkW, excessPP

! local variables:
double precision:: dExcessPP, nu_undrained

if (propUndrainedNu < TINY) then
  nu_undrained = UNDRAINED_NU
else
  nu_undrained = propUndrainedNu
endif

! excess pore pressure
BulkW = CalculateBulkModulusWater(G, nu, nu_undrained)
dExcessPP = BulkW * dEpsV

excessPP  = excessPP0 + dExcessPP

end subroutine CalculateExcessPorePressure

!--------------------------------------------------------------------------------------------------
integer function getnSubStepsBasedOnEpsq(dEps, MinEpsDev) result(nSub)
implicit none
double precision, dimension(N_STRESS_VECTOR),  intent(in) :: dEps
double precision, optional,  intent(in) :: MinEpsDev

! local variables:
integer :: nSubDev, nSubVol
double precision:: dEpsq, dEpsV

dEpsq = getEpsq(dEps)
if (present(MinEpsDev)) then
  nSubDev  = ceiling(abs(dEpsq)/(MinEpsDev+TINY))
else
  nSubDev  = ceiling(abs(dEpsq)/MIN_EPS_DEVIATORIC)
endif

dEpsV = getEpsV(dEps)
nSubVol = ceiling(abs(dEpsV)/MIN_EPS_VOLUMETRIC)

nSub = max(nSubDev, nSubVol)
nSub  = min(nSub, N_SUB_MAX)
nSub  = max(nSub, 1)

end function getnSubStepsBasedOnEpsq

!--------------------------------------------------------------------------------------------------
subroutine getdFdLamdaCrossPoint(dFdh,                       &
                                 Coupling,                   &
                                 dGdS, DdGdS, dFdS, dhdEpsP, &
                                 dFdLamda)
implicit none
double precision, dimension(N_CROSS_YIELD_SURFACES),                            intent(in) :: dFdh
double precision, dimension(N_CROSS_YIELD_SURFACES, N_CROSS_YIELD_SURFACES),    intent(in) :: Coupling
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_CROSS_YIELD_SURFACES), intent(in) :: dGdS, DdGdS, dFdS, dhdEpsP
double precision, dimension(N_CROSS_YIELD_SURFACES, N_CROSS_YIELD_SURFACES),    intent(out):: dFdLamda

! - dfdSigma*D*dGdSigma + dfdh*dhdEpsP*dGdSigma
dFdLamda(1,1) = dot_product(DdGdS(:,1), dFdS(:,1)) - dFdh(1) * dot_product(dGdS(:,1), dhdEpsP(:,1)) * Coupling(1,1)
dFdLamda(1,2) = dot_product(DdGdS(:,2), dFdS(:,1)) - dFdh(1) * dot_product(dGdS(:,2), dhdEpsP(:,1)) * Coupling(1,2)
dFdLamda(2,1) = dot_product(DdGdS(:,1), dFdS(:,2)) - dFdh(2) * dot_product(dGdS(:,1), dhdEpsP(:,2)) * Coupling(2,1)
dFdLamda(2,2) = dot_product(DdGdS(:,2), dFdS(:,2)) - dFdh(2) * dot_product(dGdS(:,2), dhdEpsP(:,2)) * Coupling(2,2)

end subroutine getdFdLamdaCrossPoint

!--------------------------------------------------------------------------------------------------
integer function getnSubStepsBasedOnSigCone(sinPhiPeak, pApex, prSig0, prSigElas, stressRatio) result(nSub)
implicit none
double precision, intent(in):: sinPhiPeak, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig0, prSigElas
double precision, optional, intent(in):: stressRatio

! local variables:
double precision:: q, qf, qStep

q  = calQTX_MC(prSigElas, FindActivePart_MC(prSigElas))
qf = calQf_MC(sinPhiPeak, pApex, prSig0, FindActivePart_MC(prSig0))

qStep = MAX_SUB_SIG_DEVIATORIC_FACTOR * qf

nSub  = ceiling(q/qStep)

if (present(stressRatio)) then
  if (stressRatio > CRITICAL_SHEAR_STRESS_RATIO .and. stressRatio < MAX_SHEAR_STRESS_RATIO) then
    nSub = nSub * N_SUB_FACTOR_CRITICAL
  endif
  if (stressRatio >= MAX_SHEAR_STRESS_RATIO) nSub = 1
endif

nSub  = min(nSub, N_SUB_MAX_INITIAL)
nSub  = max(nSub, 1)

end function getnSubStepsBasedOnSigCone

!--------------------------------------------------------------------------------------------------
double precision function calResidual(f, dLamda, dTime, viscosity) result(residual)
implicit none
double precision, intent(in):: f, dLamda, dTime, viscosity

residual = f - dLamda * viscosity / (dTime + TINY)

end function calResidual

!--------------------------------------------------------------------------------------------------
double precision function calPerzynaDFactor(f, f0, n) result(factor)
implicit none
double precision, intent(in):: f, f0, n

factor = n * (f**(n - 1.0d0)) / (f0**n)

end function calPerzynaDFactor

!--------------------------------------------------------------------------------------------------
double precision function calOverStressPhi(iOverSigFunc, f, F0, rateVisco) result(phiF)
implicit none
integer, intent(in):: iOverSigFunc
double precision, intent(in) :: f, F0, rateVisco

select case(iOverSigFunc)
case(BINGHAM_MODEL)
  phiF = f
case(PERZYNA_MODEL)
  phiF = (f/F0)**rateVisco
case default
  phiF = f
end select

end function calOverStressPhi

!--------------------------------------------------------------------------------------------------
double precision function calOverStressDFactor(iOverSigFunc, f, F0, rateVisco) result(OverStressDFactor)
implicit none
integer, intent(in):: iOverSigFunc
double precision, intent(in) :: f, F0, rateVisco

select case(iOverSigFunc)
case(PERZYNA_MODEL)
  OverStressDFactor = calPerzynaDFactor(f, F0, rateVisco)
case default
  ! Bigham model as a default model
  OverStressDFactor = 1.0d0
end select

end function calOverStressDFactor

!--------------------------------------------------------------------------------------------------
subroutine calOverStressdFdS(f, f0, rateVisco, dFdS)
implicit none
double precision, intent(in)::  f, f0, rateVisco
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(inout):: dFdS

dFdS = calPerzynaDFactor(f, f0, rateVisco) * dFdS

end subroutine calOverStressdFdS

!--------------------------------------------------------------------------------------------------
subroutine setprSigApex(pApex, prSigApex)
implicit none
double precision, intent(in)::  pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: prSigApex

prSigApex(1:N_PRINCIPAL_STRESS_VECTOR) = - pApex

end subroutine setprSigApex

!--------------------------------------------------------------------------------------------------
double precision function getTolerance(IDYieldFunc, prSig) result(res)
implicit none
integer, intent(in):: IDYieldFunc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), optional, intent(in):: prSig

! local variables:
double precision :: p, factor

factor = 1.0d0

if (present(prSig)) then
  p = calP(prSig)
  if (p < THRESHOLD_HIGH_TOLERATED_ERROR) then
    factor = HIGH_TOLERANCE_FACTOR
  endif
else
  p = 0.0d0
endif

select case(IDYieldFunc)
case(MOHR_COULOMB_YIELD)
  res = TOLERANCE_MC * factor
case(DRUCKER_PRAGER_YIELD, DRUCKER_PRAGER_YIELD_EXTENSION)
  res = TOLERANCE_DP * factor
case(MATSUOKA_NAKAI_YIELD, MATSUOKA_NAKAI_CONVEX_YIELD)
  res = TOLERANCE_MN * factor
case(TENSION_YIELD)
  res = TOLERANCE_TEN
case(COMPLEMENTARY_YIELD)
  res = TOLERANCE_COMPLEMENTARY
case(HARDENING_SOIL_MC_YIELD)
  res = TOLERANCE_HS_MC * factor
case(CAP_HS_YIELD)
  if (p > 0.0d0) then
    factor = max(1.0d0, log(p)*log10(p))
  endif
  res = TOLERANCE_CAP * factor
case(CAM_CLAY_YIELD)
  if (p > 0.0d0) then
    factor = max(1.0d0, log(p)*log10(p))
  endif
  res = TOLERANCE_CAM_CLAY * factor
case default
  res = TOLERANCE_DEFAULT
  call SetError('unknown yield function in getTolerance, IDYieldFunc: ' // trim(string(IDYieldFunc)))
end select

res = res + TOLERANCE_TOLER

end function getTolerance

!--------------------------------------------------------------------------------------------------
double precision function calCapStarValues(value, e0) result(valueStar)
implicit none
double precision, intent(in) :: value, e0

valueStar = value / (1.0d0 + e0)

end function calCapStarValues

!--------------------------------------------------------------------------------------------------
double precision function CalKCapCamClay(lambda, kappa, e0) result(K)
implicit none
double precision, intent(in) :: lambda, kappa, e0

! local variables:
double precision :: lambdaStar, kappaStar

lambdaStar = calCapStarValues(lambda, e0)
kappaStar  = calCapStarValues(kappa, e0)
K = 1.0d0 / (lambdaStar - kappaStar)

end function CalKCapCamClay

!--------------------------------------------------------------------------------------------------
double precision function calVoidRatio(e0,dEpsVi) result(eCurrent)
implicit none
double precision, intent(in) :: e0, dEpsVi

! local variables:
double precision :: dEpsV

! Soil mechanics sign convention
dEpsV = max(dEpsVi, -100.0d0)

if (.true.) then
  eCurrent = (1.0d0 + e0) * exp(-dEpsV) - 1.0d0
else
  eCurrent = -dEpsV*(1.0 + e0) + e0
endif

eCurrent = max(eCurrent, 0.0d0)

end function calVoidRatio

!--------------------------------------------------------------------------------------------------
subroutine setSubroutineCalled(IsDebug, subroutineCalled, text)
implicit none
logical, intent(in):: IsDebug
character(len=1023), intent(inout):: subroutineCalled
character(*), intent(in):: text

if (IsDebug .or. IS_MINOR_DEBUG) then
  subroutineCalled = trim(subroutineCalled) // '.' // trim(text)
endif

end subroutine setSubroutineCalled

!--------------------------------------------------------------------------------------------------
logical function isDiverged(prSig) result(res)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variable
integer:: i

res = .false.

do i=1,N_PRINCIPAL_STRESS_VECTOR
  if (abs(prSig(i)) > DIVERGENCE_LIMIT .or. isNAN(prSig(i))) then
    res = .true.
    RETURN
  endif
enddo

end function isDiverged

!--------------------------------------------------------------------------------------------------
subroutine CheckAndModifyPrincipalStressAtCornersMC(prSigElas, prSig, isStressModified, IsDebug)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)    :: prSigElas
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(inout) :: prSig
logical, intent(out):: isStressModified
logical, intent(in) :: IsDebug

! local variables:
integer :: indexActivePartMC
isStressModified = .false.

indexActivePartMC = FindActivePart_MC(prSig)
if (IsDebug) then
  call WriteInDebugFile('indexActivePartMC    : ' // trim(string(indexActivePartMC)))
  call WriteInDebugFile('indexActivePartMCElas: ' // trim(string(FindActivePart_MC(prSigElas))))
endif

select case(indexActivePartMC)
case(INDEX_COULOMB_TRIAXIAL_COMP)
  prSig = prSigElas
  prSig(2) = 0.5d0* (prSig(2) + prSig(3))
  prSig(3) = prSig(2)
  isStressModified = .true.
  if (IsDebug) call WriteInDebugFile('Triaxial compression: Stress is modified for sig2 < sig3')
case(INDEX_COULOMB_TRIAXIAL_EXTN)
  prSig = prSigElas
  prSig(2) = 0.5d0* (prSig(2) + prSig(1))
  prSig(1) = prSig(2)
  isStressModified = .true.
  if (IsDebug) call WriteInDebugFile('Triaxial extension: Stress is modified for sig1 < sig2')
end select

end subroutine CheckAndModifyPrincipalStressAtCornersMC

!--------------------------------------------------------------------------------------------------
subroutine normalisedGdS(dGdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(inout):: dGdS

dGdS = dGdS / (prStressLength(dGdS) + TINY)

end subroutine normalisedGdS

!---------------------------------------------------------------------------------------------------------
subroutine calElasticPrincipalStress(Sig0,           &
                                     dEps,           &
                                     Sig, prSigElas, &
                                     DE,             &
                                     isDebug)
 
implicit none
double precision, dimension(N_STRESS_VECTOR),                 intent(in)   :: dEps, sig0
double precision, dimension(N_STRESS_VECTOR),                 intent(inout):: sig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR),       intent(inout):: prSigElas
double precision, dimension(N_STRESS_VECTOR,N_STRESS_VECTOR), intent(in)   :: DE
logical,                                                      intent(in)   :: isDebug

! Local variables
double precision, dimension(N_STRESS_VECTOR):: dSig

! elastic stress sub_increment
call MatrixToVectorProd(DE, dEps, N_STRESS_VECTOR, dSig)
! elastic stress
call sumVectors(Sig0, dSig, 1d0, 1d0, N_STRESS_VECTOR, Sig)

! calculate principal stresses and directions
call calPrincipalStresses(Sig, prSigElas)

end subroutine calElasticPrincipalStress

end module ModSolversLibrary