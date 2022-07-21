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
module ModYieldSurfaceLibrary
use ModMathLibrary
use ModString
use ModDebug
use ModStressAndStrainLibrary
use ModGlobalConstants

implicit none

! factors
double precision, parameter:: SIGN_FACTOR_GMC = 1.0d0

public

! switches
logical, parameter:: IS_GRIFFITH_MATSUOKA_NAKAI        = .false.
logical, parameter:: USE_ALTERNATIVE_CAM_CLAY          = .false.
! Note: USE_P_TENSION_FUNCTION should be true for MN and MNC models:
logical, parameter:: USE_P_TENSION_FUNCTION            = .false.
logical, parameter:: USE_ALWAYS_SIG3_TENSION_FUNCTION  = .false.
logical, parameter:: USE_ALWAYS_MIN_SIG_TENSION_FUNCTION  = .false.

logical, parameter:: USE_SEPARATED_MOHR_COULOMB_ZONES  = .true.
logical, parameter:: USE_MOHR_COULOMB_SHAPE_CAP        = .true.
logical, parameter:: USE_NEW_F_HS_MC_FORMULA           = .false.

! indices constants for Mohr-Coulomb model
integer, parameter:: INDEX_COULOMB_YIELD_GENERAL = 1
integer, parameter:: INDEX_COULOMB_TRIAXIAL_COMP = 2
integer, parameter:: INDEX_COULOMB_TRIAXIAL_EXTN = 3

contains

!--------------------------------------------------------------------------------------------------
!--- General functions: ---------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
double precision function calM_TriaxComp(sinPhi) result(xM)
implicit none
double precision, intent(in):: sinPhi

xM = 6.0 * sinPhi / (3.0 - sinPhi)

end function calM_TriaxComp

!--------------------------------------------------------------------------------------------------
double precision function calM_TriaxExt(sinPhi) result(xM)
implicit none
double precision, intent(in):: sinPhi

xM = 6.0 * sinPhi / (3.0 + sinPhi)

end function calM_TriaxExt

!--------------------------------------------------------------------------------------------------
double precision function getSinPhiFromMTriaxialCompression(xM) result(sinPhi)
implicit none
double precision, intent(in):: xM

sinPhi = 3.0d0 * xM / (6.0 + xM)

end function getSinPhiFromMTriaxialCompression

!--------------------------------------------------------------------------------------------------
double precision function getSinPhiFromMTriaxialExtension(xM) result(sinPhi)
implicit none
double precision, intent(in):: xM

sinPhi = 3.0d0 * xM / (6.0 - xM)

end function getSinPhiFromMTriaxialExtension

!--------------------------------------------------------------------------------------------------
!--- Cap model functions: -------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! return f of cap yield functions
double precision function calF_CAP_HS(xM, pc, prSig, sinPhi, indexActivePartMC) result(f)
implicit none
double precision, intent(in):: xM, pc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, optional, intent(in):: sinPhi
integer, optional, intent(in) :: indexActivePartMC

if (USE_MOHR_COULOMB_SHAPE_CAP .and. present(sinPhi)) then
  f = calF_CAP_HS_MC(xM, pc, prSig, sinPhi, indexActivePartMC)
else
  f = calF_CAP_HS_DP(xM, pc, prSig)
endif

end function calF_CAP_HS

!--------------------------------------------------------------------------------------------------
double precision function calF_CAP_HS_DP(xM, pc, prSig) result(f)
implicit none
double precision, intent(in):: xM, pc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: p, q

p = calP(prSig)
q = calQ(prSig)
f = (q/xM)**2 + p**2 - pc**2

end function calF_CAP_HS_DP

!--------------------------------------------------------------------------------------------------
double precision function calF_CAP_HS_MC(xM, pc, prSig, sinPhi, indexActivePartMC_I) result(f)
implicit none
double precision, intent(in):: xM, pc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, intent(in):: sinPhi
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    ! Triaxial compression
    f = calF_CAP_HS_MC_TriaxComp(xM, pc, prSig)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    ! Triaxial extension
    f = calF_CAP_HS_MC_TriaxExt(xM, pc, prSig, sinPhi)
  case(INDEX_COULOMB_YIELD_GENERAL)
    ! between triaxial compression and extension
    f = calF_CAP_HS_MC_General(xM, pc, prSig, sinPhi)
  case default
    call SetError('unknown indexActivePartMC in calF_CAP_HS_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  ! between triaxial compression and extension
  f = calF_CAP_HS_MC_General(xM, pc, prSig, sinPhi)
endif

end function calF_CAP_HS_MC

!--------------------------------------------------------------------------------------------------
double precision function calF_CAP_HS_MC_General(xM, pc, prSig, sinPhi) result(f)
implicit none
double precision, intent(in):: xM, pc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, intent(in):: sinPhi

! local variables:
double precision :: p, q, delta

p = calP(prSig)

! delta = M_ext / M_comp
delta = (3.0d0 + sinPhi) / (3.0d0 - sinPhi)
q = prSig(1) + (delta - 1.0d0)*prSig(2) - delta * prSig(3)

f = (q/xM)**2 + p**2 - pc**2

end function calF_CAP_HS_MC_General

!--------------------------------------------------------------------------------------------------
double precision function calF_CAP_HS_MC_TriaxComp(xM, pc, prSig) result(f)
implicit none
double precision, intent(in):: xM, pc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: p, q

p = calP(prSig)

! this equation can be derived by setting both sig_2 and sig_3 to (sig_2 + sig_3) / 2 in the general formula
q = prSig(1) - 0.5d0*(prSig(2) + prSig(3))

f = (q/xM)**2 + p**2 - pc**2

end function calF_CAP_HS_MC_TriaxComp

!--------------------------------------------------------------------------------------------------
double precision function calF_CAP_HS_MC_TriaxExt(xM, pc, prSig, sinPhi) result(f)
implicit none
double precision, intent(in):: xM, pc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, intent(in):: sinPhi

! local variables:
double precision :: p, q, delta

p = calP(prSig)

! this equation can be derived by setting both sig_1 and sig_2 to (sig_1 + sig_2) / 2 in the general formula
! delta = M_ext / M_comp
delta = (3.0d0 + sinPhi) / (3.0d0 - sinPhi)
q = 0.5d0*delta*(prSig(1) + prSig(2)) - delta * prSig(3)

f = (q/xM)**2 + p**2 - pc**2

end function calF_CAP_HS_MC_TriaxExt

!--------------------------------------------------------------------------------------------------
! return df/dSigma of cap yield functions
subroutine caldFdS_CAP_HS(xM, prSig, dFdS, sinPhi, indexActivePartMC)
implicit none
double precision, intent(in)::  xM
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
double precision, optional, intent(in):: sinPhi
integer, optional, intent(in) :: indexActivePartMC

if (USE_MOHR_COULOMB_SHAPE_CAP .and. present(sinPhi)) then
  call caldFdS_CAP_HS_MC(xM, prSig, dFdS, sinPhi, indexActivePartMC)
else
  call caldFdS_CAP_HS_DP(xM, prSig, dFdS)
endif

end subroutine caldFdS_CAP_HS

!--------------------------------------------------------------------------------------------------
! return df/dSigma of cap yield functions (circlar)
subroutine caldFdS_CAP_HS_DP(xM, prSig, dFdS)
implicit none
double precision, intent(in)::  xM
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
double precision p, q
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dPdS, dQds

p = calP(prSig)
call getdPdSigma(dPdS)
q = calQ(prSig)
call getdQdSigma(p, q, prSig, dQds)
!f = (q/xM)**2 + p**2 - pc**2
dFdS = ((1.0d0/xM**2))*2.0d0*q*dQds + 2.0d0*p*dPdS

end subroutine caldFdS_CAP_HS_DP

!--------------------------------------------------------------------------------------------------
! return df/dSigma of cap yield functions
subroutine caldFdS_CAP_HS_MC(xM, prSig, dFdS, sinPhi, indexActivePartMC_I)
implicit none
double precision, intent(in)::  xM
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
double precision, intent(in):: sinPhi
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    ! Triaxial compression
    call caldFdS_CAP_HS_MC_TriaxComp(xM, prSig, dFdS)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    ! Triaxial extension
    call caldFdS_CAP_HS_MC_TriaxExt(xM, prSig, dFdS, sinPhi)
  case(INDEX_COULOMB_YIELD_GENERAL)
    ! between triaxial compression and extension
    call caldFdS_CAP_HS_MC_General(xM, prSig, dFdS, sinPhi)
  case default
    call SetError('unknown indexActivePartMC in caldFdS_CAP_HS_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  ! between triaxial compression and extension
  call caldFdS_CAP_HS_MC_General(xM, prSig, dFdS, sinPhi)
endif

end subroutine caldFdS_CAP_HS_MC

!--------------------------------------------------------------------------------------------------
! return df/dSigma of cap yield functions
subroutine caldFdS_CAP_HS_MC_General(xM, prSig, dFdS, sinPhi)
implicit none
double precision, intent(in)::  xM
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
double precision, intent(in):: sinPhi

! local variables:
double precision p, q, delta
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dPdS, dQds

p = calP(prSig)
call getdPdSigma(dPdS)

! delta = M_ext / M_comp
delta = (3.0d0 + sinPhi) / (3.0d0 - sinPhi)
q = prSig(1) + (delta - 1.0d0)*prSig(2) - delta * prSig(3)

dQds(1) =  1.0d0
dQds(2) =  delta - 1.0d0
dQds(3) = -delta

!f = (q/xM)**2 + p**2 - pc**2
dFdS = ((1.0d0/xM**2))*2.0d0*q*dQds + 2.0d0*p*dPdS

end subroutine caldFdS_CAP_HS_MC_General

!--------------------------------------------------------------------------------------------------
! return df/dSigma of cap yield functions
subroutine caldFdS_CAP_HS_MC_TriaxComp(xM, prSig, dFdS)
implicit none
double precision, intent(in)::  xM
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
double precision p, q
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dPdS, dQds

p = calP(prSig)
call getdPdSigma(dPdS)

q = prSig(1) - 0.5d0*prSig(2) - 0.5d0*prSig(3)

dQds(1) = 1.0d0
dQds(2) = -0.5d0
dQds(3) = -0.5d0

!f = (q/xM)**2 + p**2 - pc**2
dFdS = ((1.0d0/xM**2))*2.0d0*q*dQds + 2.0d0*p*dPdS

end subroutine caldFdS_CAP_HS_MC_TriaxComp

!--------------------------------------------------------------------------------------------------
! return df/dSigma of cap yield functions
subroutine caldFdS_CAP_HS_MC_TriaxExt(xM, prSig, dFdS, sinPhi)
implicit none
double precision, intent(in)::  xM
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
double precision, intent(in):: sinPhi

! local variables:
double precision:: p, q, delta
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dPdS, dQds

p = calP(prSig)
call getdPdSigma(dPdS)

! this equation can be derived by setting both sig_1 and sig_2 to (sig_1 + sig_2) / 2 in the general formula
! delta = M_ext / M_comp
delta = (3.0d0 + sinPhi) / (3.0d0 - sinPhi)
q = 0.5d0*delta*(prSig(1) + prSig(2)) - delta * prSig(3)

dQds(1) = 0.5d0 * delta
dQds(2) = 0.5d0 * delta
dQds(3) = -delta

!f = (q/xM)**2 + p**2 - pc**2
dFdS = ((1.0d0/xM**2))*2.0d0*q*dQds + 2.0d0*p*dPdS

end subroutine caldFdS_CAP_HS_MC_TriaxExt

!--------------------------------------------------------------------------------------------------
! derivative of f with respect to pc parameter for Cap model
double precision function calDfDh_CAP_HS(pc) result(dFdh)
implicit none
double precision, intent(in):: pc

dFdh = - 2.0*pc

end function calDfDh_CAP_HS

!--------------------------------------------------------------------------------------------------
! return qf (q at failure) of cap with Mohr-Coulomb criterion
double precision function calQf_MC(sinPhi, pApex, prSig, indexActivePartMC_I) result(qf)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    ! Triaxial compression
    qf = calQf_MC_TriaxialCom(sinPhi, pApex, prSig)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    ! Triaxial extension
    qf = calQf_MC_TriaxialExt(sinPhi, pApex, prSig)
  case(INDEX_COULOMB_YIELD_GENERAL)
    ! between triaxial compression and extension
    qf = calQf_MC_General(sinPhi, pApex, prSig)
  case default
    call SetError('unknown indexActivePartMC in calQf_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  ! between triaxial compression and extension
  qf = calQf_MC_General(sinPhi, pApex, prSig)
endif

end function calQf_MC

!--------------------------------------------------------------------------------------------------
! return qf (q at failure) of cap with Mohr-Coulomb criterion
double precision function calQf_MC_General(sinPhi, pApex, prSig) result(qf)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

qf = (2*sinPhi * (prSig(3) + pApex)) / (1.0d0 - sinPhi)

end function calQf_MC_General

!--------------------------------------------------------------------------------------------------
! return qf (q at failure) of cap with Mohr-Coulomb criterion
double precision function calQf_MC_TriaxialCom(sinPhi, pApex, prSig) result(qf)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

qf = (2*sinPhi * (0.5d0 * (prSig(2) + prSig(3)) + pApex)) / (1.0d0 - sinPhi)

end function calQf_MC_TriaxialCom

!--------------------------------------------------------------------------------------------------
! return qf (q at failure) of cap with Mohr-Coulomb criterion
double precision function calQf_MC_TriaxialExt(sinPhi, pApex, prSig) result(qf)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

qf = (2.0d0 * sinPhi * (prSig(3) + pApex)) / (1.0d0 - sinPhi)

end function calQf_MC_TriaxialExt

!--------------------------------------------------------------------------------------------------
! return qf (q at failure) of cap with Drucker-Prager criterion
double precision function calQf_DP(sinPhi, pApex, prSig) result(qf)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: p, xM

p = calP(prSig)
xM = calM_TriaxComp(sinPhi)

qf = xM * ( p +  pApex)

end function calQf_DP

!--------------------------------------------------------------------------------------------------
subroutine caldQfdS_MC(sinPhi, dQfdS, indexActivePartMC_I)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQfdS
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    ! Triaxial compression
    call caldQfdS_MC_TriaxialCom(sinPhi, dQfdS)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    ! Triaxial extension
    call caldQfdS_MC_TriaxialExt(sinPhi, dQfdS)
  case(INDEX_COULOMB_YIELD_GENERAL)
    ! between triaxial compression and extension
    call caldQfdS_MC_General(sinPhi, dQfdS)
  case default
    call SetError('unknown indexActivePartMC in caldQfdS_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  ! between triaxial compression and extension
  call caldQfdS_MC_General(sinPhi, dQfdS)
endif

end subroutine caldQfdS_MC
!--------------------------------------------------------------------------------------------------
subroutine caldQfdS_MC_General(sinPhi, dQfdS)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQfdS

dQfdS(1) = 0.0d0
dQfdS(2) = 0.0d0
dQfdS(3) = 2.0d0*sinPhi / (1.0d0 - sinPhi)

end subroutine caldQfdS_MC_General

!--------------------------------------------------------------------------------------------------
subroutine caldQfdS_MC_TriaxialCom(sinPhi, dQfdS)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQfdS

dQfdS(1) = 0.0d0
dQfdS(2) = sinPhi / (1.0d0 - sinPhi)
dQfdS(3) = sinPhi / (1.0d0 - sinPhi)

end subroutine caldQfdS_MC_TriaxialCom

!--------------------------------------------------------------------------------------------------
subroutine caldQfdS_MC_TriaxialExt(sinPhi, dQfdS)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQfdS

dQfdS(1) = 0.0d0
dQfdS(2) = 0.0d0
dQfdS(3) = 2.0d0*sinPhi / (1.0d0 - sinPhi)

end subroutine caldQfdS_MC_TriaxialExt

!--------------------------------------------------------------------------------------------------
! return f of hyperbolic yield function with Mohr-Coulomb criterion
double precision function calF_HS_MC(sinPhi, pApex, Rf, Ei, Eur, gammaP, prSig, indexActivePartMC) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex, Rf, Ei, Eur, gammaP
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
integer, optional, intent(in) :: indexActivePartMC

! local variables:
double precision:: q, qf, qa
double precision:: term1, term2, term3

q  = calQTX_MC(prSig, indexActivePartMC)
qf = calQf_MC(sinPhi, pApex, prSig, indexActivePartMC)
qa = qf/Rf

if (USE_NEW_F_HS_MC_FORMULA) then
  ! This reformulation does not suffer from devision by zero!
  term1 = (2.0d0 * q * qa) / Ei
  term2 = 1.0d0
  term3 = 2.0d0 * q * (qa - q) / Eur
  f = term1 * term2 - term3 - gammaP * (qa - q)

else
  if (q > ((1.0d0 - SMALL) * qa)) then
    f = calF_MC(sinPhi, pApex, prSig, indexActivePartMC)
  else
    term1 = (2.0d0 * q) / Ei
    term2 = qa / ((qa - q) + TINY)
    term3 = 2.0d0 * q / Eur

    f = term1 * term2 - term3 - gammaP
  endif
endif

end function calF_HS_MC

!--------------------------------------------------------------------------------------------------
! return df/dSigma of hyperbolic yield function with Mohr-Coulomb criterion
subroutine caldFdS_HS_MC(sinPhi, pApex, Rf, Ei, Eur, gammaP, prSig, dFdS, indexActivePartMC)
implicit none
double precision, intent(in)::  sinPhi, pApex, Rf, Ei, Eur, gammaP
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
integer, optional, intent(in) :: indexActivePartMC

! local variables:
double precision:: q, qf, qa
double precision:: term1, term2
double precision:: q2, qa2, denominator2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dQTxdS, dQfdSUnitVector

if (USE_NEW_F_HS_MC_FORMULA) then
  ! This reformulation does not suffer from devision by zero!
  call getdFdS_HS_MC_Reformulated(sinPhi, pApex, Rf, Ei, Eur, gammaP, prSig, dFdS, indexActivePartMC)
else
  dQTxdS          = reshape((/1.0d0,  0.0d0, -1.0d0/), shape(dQTxdS))
  dQfdSUnitVector = reshape((/0.0d0,  0.0d0,  1.0d0/), shape(dQfdSUnitVector))

  q  = calQTX_MC(prSig, indexActivePartMC)
  qf = calQf_MC(sinPhi, pApex, prSig, indexActivePartMC)
  qa = qf/Rf

  if (q > ((1.0d0 - SMALL) * qa)) then
    call caldFdS_MC(sinPhi, dFdS, indexActivePartMC)
  else
    denominator2 = 1.0d0 / ((qa - q) + TINY)
    denominator2 = denominator2 * denominator2

    q2  = q * q
    qa2 = qa * qa

    term1 =  (2.0d0 / Ei) * qa2 * denominator2 - (2.0d0 / Eur)
    term2 = -(2.0d0 / Ei) * q2  * denominator2 * (2.0d0*sinPhi / (Rf * (1.0d0 - sinPhi)))

    dFdS =  term1 * dQTxdS +  term2 * dQfdSUnitVector
  endif
endif

end subroutine caldFdS_HS_MC

!--------------------------------------------------------------------------------------------------
! return df/dSigma of hyperbolic yield function with Mohr-Coulomb criterion
subroutine getdFdS_HS_MC_Reformulated(sinPhi, pApex, Rf, Ei, Eur, gammaP, prSig, dFdS, indexActivePartMC)
implicit none
double precision, intent(in)::  sinPhi, pApex, Rf, Ei, Eur, gammaP
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
integer, optional, intent(in):: indexActivePartMC

! local variables:
double precision:: q, qf, qa
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dQdS, dQadS, dQfdS
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: qadQdS, qdQadS

q  = calQTX_MC(prSig, indexActivePartMC)
qf = calQf_MC(sinPhi, pApex, prSig, indexActivePartMC)
qa = qf/Rf

call caldQfdS_MC(sinPhi, dQfdS, indexActivePartMC)
dQadS = dQfdS / Rf

call caldQTXdS_MC(dQdS, indexActivePartMC)

qadQdS = qa * dQdS
qdQadS = q  * dQadS

dFdS =  (2.0d0 / Ei)  * (qadQdS + qdQadS) &
      - (2.0d0 / Eur) * (qadQdS + qdQadS) &
      + (4.0d0 / Eur) * q * dQdS          &
      - gammaP*(dQadS - dQdS)

end subroutine getdFdS_HS_MC_Reformulated

!--------------------------------------------------------------------------------------------------
! return f of hyperbolic yield function with Drucker-Prager criterion
double precision function calF_HS_DP(sinPhi, pApex, Rf, Ei, Eur, gammaP, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex, Rf, Ei, Eur, gammaP
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision:: q, qf, qa
double precision:: term1, term2, term3

q  = calQ(prSig)
qf = calQf_DP(sinPhi, pApex, prSig)
qa = qf/Rf

if (q > ((1.0d0 - RELATIVELY_SMALL) * qa)) then
  f = calF_DP(sinPhi, pApex, prSig)
else
  term1 = (2.0d0 * q) / Ei
  term2 = qa / ((qa - q) + TINY)
  term3 = 2.0d0 * q / Eur

  f = term1 * term2 - term3 - gammaP
endif

end function calF_HS_DP

!--------------------------------------------------------------------------------------------------
double precision function calQTX_MC(prSig, indexActivePartMC_I) result(q)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_YIELD_GENERAL)
    q = calqTX_General(prSig)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    q = calqTX_TriaxialCom(prSig)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    q = calqTX_TriaxialExt(prSig)
  case default
    call SetError('unknown indexActivePartMC in calQTX_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  q = calqTX_General(prSig)
endif

end function calQTX_MC

!--------------------------------------------------------------------------------------------------
subroutine caldQTXdS_MC(dQdS, indexActivePartMC_I)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQdS
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_YIELD_GENERAL)
    call caldQTXdS_MC_General(dQdS)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    call caldQTXdS_MC_TriaxialCom(dQdS)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    call caldQTXdS_MC_TriaxialExt(dQdS)
  case default
    call SetError('unknown indexActivePartMC in caldQTXdS_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  call caldQTXdS_MC_General(dQdS)
endif

end subroutine caldQTXdS_MC

!--------------------------------------------------------------------------------------------------
double precision function calqTX_General(prSig) result(q)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

q = prSig(1) - prSig(3)

end function calqTX_General

!--------------------------------------------------------------------------------------------------
subroutine caldQTXdS_MC_General(dQdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQdS

dQdS(1) =  1.0d0
dQdS(2) =  0.0d0
dQdS(3) = -1.0d0

end subroutine caldQTXdS_MC_General

!--------------------------------------------------------------------------------------------------
double precision function calqTX_TriaxialCom(prSig) result(q)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

q = prSig(1) - 0.5d0*(prSig(2) + prSig(3))

end function calqTX_TriaxialCom

!--------------------------------------------------------------------------------------------------
subroutine caldQTXdS_MC_TriaxialCom(dQdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQdS

dQdS(1) =  1.0d0
dQdS(2) = -0.5d0
dQdS(3) = -0.5d0

end subroutine caldQTXdS_MC_TriaxialCom

!--------------------------------------------------------------------------------------------------
double precision function calqTX_TriaxialExt(prSig) result(q)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

q = 0.5d0*(prSig(1) + prSig(2)) - prSig(3)

end function calqTX_TriaxialExt

!--------------------------------------------------------------------------------------------------
subroutine caldQTXdS_MC_TriaxialExt(dQdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dQdS

dQdS(1) =  0.5d0
dQdS(2) =  0.5d0
dQdS(3) = -1.0d0

end subroutine caldQTXdS_MC_TriaxialExt

!--------------------------------------------------------------------------------------------------
!--- Matsouka-Nakai-Convex model functions: -------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! return df/dSigma of Convex Matsouka-Nakai yield surfaces
subroutine caldFdS_MNC(sinPhi, apexHyperplolic, prSig, dFdS)
implicit none
double precision, intent(in) :: sinPhi
double precision, intent(in) :: apexHyperplolic
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
double precision :: J2, J3, theta, GTheta, dGdTheta, alpha, a, ksi, eta
double precision :: C1, C2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S, dSigmaBardSigma, dPdSigma

a = apexHyperplolic

call getPrincipalDeviatoricS(prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)
theta = getLodeAngle(J2,J3)

ksi    = calKsi_MNC(sinPhi)
GTheta = calGTheta_MNC(theta, ksi)
eta    = calEta_MNC(sinPhi)
dGdTheta = caldGdTheta_MNC(theta, ksi)

alpha = GTheta / sqrt(GTheta * GTheta + a*a * eta*eta)

C1 = - eta
call getdPdSigma(dPdSigma)

C2 = GTheta - tan(3.0d0*theta)*dGdTheta
call getdSigmaBardSigma(J2, S, dSigmaBardSigma)

dFdS = C1*dPdSigma + alpha*C2*dSigmaBardSigma

end subroutine caldFdS_MNC

!--------------------------------------------------------------------------------------------------
! returns f of Matsoua-Nakai-Convex yield surfaces
! paper: Pantephini and Lagioia (2014)
double precision function calF_MNC(sinPhi, apexHyperplolic, pApex, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi
double precision, intent(in):: apexHyperplolic
double precision, intent(in):: pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: J2, J3, theta, GTheta, p, a, ksi, eta
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

a = apexHyperplolic

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)

theta = getLodeAngle(J2,J3)

ksi    = calKsi_MNC(sinPhi)
GTheta = calGTheta_MNC(theta, ksi)
eta    = calEta_MNC(sinPhi)

f = -p*eta + sqrt( (J2 * GTheta*GTheta) + a*a * eta*eta) - pApex*eta

end function calF_MNC

!--------------------------------------------------------------------------------------------------
! derivative of f with respect to hardening parameter for Matsouka-Nakai-Convex model
! paper: Pantephini and Lagioia (2014)
double precision function calDfDh_MNC(sinPhi, apexHyperplolic, pApex, prSig) result(dfDh)
implicit none
double precision, intent(in):: sinPhi, apexHyperplolic, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: J2, J3, theta, GTheta, p, a, ksi, eta, sinPhi2, a2, sin3Theta
double precision :: dEtaDh, dKsiDh, dGThetaDKsi, dGThetaDh, denominator
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

sinPhi2 = sinPhi * sinPhi
a = apexHyperplolic
a2 = a*a

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)

theta = getLodeAngle(J2,J3)
sin3Theta = sin(3.0d0*theta)

ksi    = calKsi_MNC(sinPhi)
GTheta = calGTheta_MNC(theta, ksi)
eta    = calEta_MNC(sinPhi)

dEtaDh = 18.0d0 / sqrt((3.0d0 + sinPhi2)**3)
dKsiDh = 27.0d0*(1.0d0-sinPhi2) / sqrt((3.0d0 + sinPhi2)**5)

dGThetaDKsi = 2.0d0*sin3Theta * sin(acos(ksi*sin3Theta)/3.0d0)/(sqrt(3.0d0*(1-(sin3Theta*ksi)**2)))
dGThetaDh = dGThetaDKsi * dKsiDh
denominator = sqrt((J2 * GTheta*GTheta) + a2 * eta*eta)

dfDh = -p * dEtaDh + (a2*eta*dEtaDh) / denominator + J2*GTheta*dGThetaDh / denominator - pApex * dEtaDh

end function calDfDh_MNC

!--------------------------------------------------------------------------------------------------
! returns ksi in Matsouka-Nakai-Convex model
double precision function calKsi_MNC(sinPhi) result(ksi)
implicit none
double precision, intent(in):: sinPhi

! local variables:
double precision :: denominator

denominator = 3.0d0 + sinPhi * sinPhi
denominator = denominator * denominator * denominator

ksi = sinPhi * (9.0d0 - sinPhi * sinPhi) / sqrt(denominator)

end function calKsi_MNC

!--------------------------------------------------------------------------------------------------
! returns derivative of G with respect to theta in Matsouka-Nakai-Convex model
double precision function caldGdTheta_MNC(theta, ksi) result(dGdTheta)
implicit none
double precision, intent(in):: theta, ksi

! local variables:
double precision :: beta, cos3T, sin3T, factor

cos3T = cos(3.0d0 * theta)
sin3T = sin(3.0d0 * theta)
beta = acos(ksi * sin3T) / 3.0d0
factor = 2.0d0 * sqrt(3.0d0)

dGdTheta =  factor * cos3T * sin(beta) / sqrt(1.0d0 - ksi*ksi * sin3T*sin3T)

end function caldGdTheta_MNC

!--------------------------------------------------------------------------------------------------
! returns G(Theta) in Matsouka-Nakai-Convex model
double precision function calGTheta_MNC(theta, ksi) result(gTheta)
implicit none
double precision, intent(in):: theta, ksi

! local variables:
double precision :: beta

beta = acos(ksi * sin(3.0d0 * theta))
gTheta = 2.0d0 * sqrt(3.0d0) * cos( beta / 3.0d0)

end function calGTheta_MNC

!--------------------------------------------------------------------------------------------------
! returns eta in Matsouka-Nakai-Convex model
double precision function calEta_MNC(sinPhi) result(eta)
implicit none
double precision, intent(in):: sinPhi

eta = 6.0d0 * sinPhi / sqrt(3.0d0 + sinPhi * sinPhi)

end function calEta_MNC

!--------------------------------------------------------------------------------------------------
!--- Matsouka-Nakai model functions: --------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! returns f of Matsouka-Nakai yield surfaces
double precision function calF_MN(sinPhi, pApex, prSigI) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSigI

! local variables:
double precision :: xM, aI1, aI2, aI3
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: prSig

! appex is based on phi peak
prSig(:) = prSigI(:) + pApex

aI1 = getI1(prSig)
aI2 = getI2(prSig)
aI3 = abs(getI3(prSig))

xM = (9 - sinPhi*sinPhi)/((1.0 - sinPhi*sinPhi) + TINY)

if (IS_GRIFFITH_MATSUOKA_NAKAI) then
  ! compression is posititive
  f = (xM*(aI3/aI2) - aI1)
else
  ! compression is posititive
  f = (aI1*aI2 - xM*aI3)
endif

end function calF_MN

!--------------------------------------------------------------------------------------------------
! Derivative of Matsouka-Nakai yield surface w.r.t. principal stress
subroutine caldFdS_MN(sinPhi, pApex, prSigI, dFdS)
implicit none
double precision, intent(in)::  sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSigI
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
integer :: i
double precision :: xM, aI1, aI2, aI3
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: prSig, dI1dS, dI2dS, dI3dS

prSig(:) = prSigI(:) + pApex

aI1 = getI1(prSig)
aI2 = getI2(prSig)
aI3 = abs(getI3(prSig))
xM  = (9.0 - sinPhi*sinPhi)/(1.0 - sinPhi*sinPhi)

if (IS_GRIFFITH_MATSUOKA_NAKAI) then
  dFdS(1) = xM*((prSig(2)*prSig(3))/aI2)**2.0 - 1.0d0
  dFdS(2) = xM*((prSig(1)*prSig(3))/aI2)**2.0 - 1.0d0
  dFdS(3) = xM*((prSig(1)*prSig(2))/aI2)**2.0 - 1.0d0
else
  ! derivatives of stress invarients
  dI1dS = 1.0

  dI2dS(1) =             prSig(2) + prSig(3)
  dI2dS(2) = prSig(1)             + prSig(3)
  dI2dS(3) = prSig(1) + prSig(2)

  dI3dS(1) =             prSig(2) * prSig(3)
  dI3dS(2) = prSig(1)             * prSig(3)
  dI3dS(3) = prSig(1) * prSig(2)

  do i=1,N_PRINCIPAL_STRESS_VECTOR
    dFdS(i) = (dI1dS(i)*aI2 + aI1*dI2dS(i) - xM*dI3dS(i))
  enddo
endif

end subroutine caldFdS_MN

!--------------------------------------------------------------------------------------------------
! derivative of f with respect to hardening parameter for Matsouka-Nakai model
double precision function calDfDh_MN(sinPhiMob, prSig) result(dFdh)
implicit none
double precision, intent(in):: sinPhiMob
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision xI3, sin2PhiMob

xI3 = getI3(prSig)
sin2PhiMob = sinPhiMob * sinPhiMob

dFdh = - 16.0*sinPhiMob*xI3 / ( (1.0 - sin2PhiMob)*(1.0 - sin2PhiMob))

end function calDfDh_MN

!--------------------------------------------------------------------------------------------------
!--- Drucker-Prager model functions: --------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! returns f of Drucker-Prager yield surfaces (M at compression point): Default function
double precision function calF_DP(sinPhi, pApex, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision xM, p, q

q = calQ(prSig)
p = calP(prSig)

xM = calM_TriaxComp(sinPhi)

f = q - xM * (p + pApex)

end function calF_DP

!--------------------------------------------------------------------------------------------------
! returns f of Drucker-Prager yield surface (M at extension point)
double precision function calF_DPExtension(sinPhi, pApex, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: xM, p, q

q = calQ(prSig)
p = calP(prSig)

xM = calM_TriaxExt(sinPhi)

f = q - xM * (p + pApex)

end function calF_DPExtension

!--------------------------------------------------------------------------------------------------
! Derivative of Drucker-Prager yield surface w.r.t. principal stress (M at triaxial compression)
subroutine caldFdS_DP(sinPhi, prSig, dFdS)
implicit none
double precision, intent(in) :: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
integer :: i
double precision :: p, a, q

q = calQ(prSig)
p = calP(prSig)

! using Lode's angle at triaxial compression
a = 2.0d0*sinPhi/(3.0d0-sinPhi)

do i=1,N_PRINCIPAL_STRESS_VECTOR
  dFdS(i) = ((3.0d0*(prSig(i) - p)/(2.0d0*q)) - a)
enddo

end subroutine caldFdS_DP

!--------------------------------------------------------------------------------------------------
! Derivative of Drucker-Prager yield surface w.r.t. principal stress (M at triaxial extension)
subroutine caldFdS_DPExtension(sinPhi, prSig, dFdS)
implicit none
double precision, intent(in) :: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
integer :: i
double precision :: p, a, q

q = calQ(prSig)
p = calP(prSig)

! using Lode's angle at triaxial extension
a = 2.0d0*sinPhi/(3.0d0+sinPhi)

do i=1,N_PRINCIPAL_STRESS_VECTOR
  dFdS(i) = ((3.0d0*(prSig(i) - p)/(2.0d0*q)) - a)
enddo

end subroutine caldFdS_DPExtension

!--------------------------------------------------------------------------------------------------
! derivative of f with respect to hardening parameter for Drucker-Prager model
double precision function calDfDh_DP(sinPhiMob, pApex, prSig) result(dFdh)
implicit none
double precision, intent(in):: sinPhiMob, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: p

p = calP(prSig)

dFdh = - 18.0*(p + pApex) / ((3.0- sinPhiMob)*(3.0- sinPhiMob))

end function calDfDh_DP

!--------------------------------------------------------------------------------------------------
! derivative of f with respect to hardening parameter for Drucker-Prager model
double precision function calDfDh_DPExtension(sinPhiMob, pApex, prSig) result(dFdh)
implicit none
double precision, intent(in):: sinPhiMob, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: p

p = calP(prSig)

dFdh = - 18.0*(p + pApex) / ((3.0 + sinPhiMob)*(3.0 + sinPhiMob))

end function calDfDh_DPExtension

!--------------------------------------------------------------------------------------------------
double precision function calSinPhiMob_DP(pApex, prSig) result(sinPhiMob)
implicit none
double precision, intent(in):: pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision :: xM, p, q

q = calQ(prSig)
p = calP(prSig)

xM = q / max(p + pApex, TINY)

sinPhiMob = (3.0*xM) / (6.0 + xM)

end function calSinPhiMob_DP

!--------------------------------------------------------------------------------------------------
double precision function calSinPhiMob_DPExtension(pApex, prSig) result(sinPhiMob)
implicit none
double precision, intent(in):: pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

double precision xM, p, q

q = calQ(prSig)
p = calP(prSig)

xM = q / max(p + pApex, TINY)

sinPhiMob = getSinPhiFromMTriaxialExtension(xM)

end function calSinPhiMob_DPExtension

!--------------------------------------------------------------------------------------------------
!--- tensile yield functions: ---------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! returns f of tensile yield surfaces
double precision function calF_Ten(sigTen, prSig, indexActivePartMC_I) result(f)
implicit none
double precision, intent(in):: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_P_TENSION_FUNCTION) then
  f = calF_Ten_P(sigTen, prSig)
elseif (USE_ALWAYS_SIG3_TENSION_FUNCTION) then
  f = calF_Ten_Sig3_General(sigTen, prSig)
elseif (USE_ALWAYS_MIN_SIG_TENSION_FUNCTION) then
  f = calF_Ten_Min_Sig(sigTen, prSig)
else
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_YIELD_GENERAL)
    f = calF_Ten_Sig3_General(sigTen, prSig)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    f = calF_Ten_Sig3_TxComp(sigTen, prSig)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    f = calF_Ten_Sig3_TxExt(sigTen, prSig)
  case default
    call SetError('unknown indexActivePartMC in calF_Ten: ' // trim(string(indexActivePartMC)))
  end select
endif

end function calF_Ten

!--------------------------------------------------------------------------------------------------
subroutine caldFdS_Ten(prSig, dFdS, indexActivePartMC_I)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_P_TENSION_FUNCTION) then
  call caldFdS_Ten_P(dFdS)
elseif (USE_ALWAYS_SIG3_TENSION_FUNCTION) then
  call caldFdS_Ten_Sig3_General(dFdS)
elseif (USE_ALWAYS_MIN_SIG_TENSION_FUNCTION) then
  call caldFdS_Ten_Min_Sig(prSig, dFdS)
else
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_YIELD_GENERAL)
    call caldFdS_Ten_Sig3_General(dFdS)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    call caldFdS_Ten_Sig3_TxCom(dFdS)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    call caldFdS_Ten_Sig3_TxExt(dFdS)
  case default
    call SetError('unknown indexActivePartMC in caldFdS_Ten: ' // trim(string(indexActivePartMC)))
  end select
endif

end subroutine caldFdS_Ten

!--------------------------------------------------------------------------------------------------
subroutine caldFdS_Ten_P(dFdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

dFdS(1) = -1.0d0 / 3.0d0
dFdS(2) = -1.0d0 / 3.0d0
dFdS(3) = -1.0d0 / 3.0d0

end subroutine caldFdS_Ten_P

!--------------------------------------------------------------------------------------------------
subroutine caldFdS_Ten_Min_Sig(prSig, dFdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
integer :: indexMin

indexMin = minloc(prSig,1)
dFdS =  0.0d0
dFdS(indexMin) =  -1.0d0

end subroutine caldFdS_Ten_Min_Sig

!--------------------------------------------------------------------------------------------------
subroutine caldFdS_Ten_Sig3_General(dFdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

dFdS(1) =  0.0d0
dFdS(2) =  0.0d0
dFdS(3) = -1.0d0

end subroutine caldFdS_Ten_Sig3_General

!--------------------------------------------------------------------------------------------------
subroutine caldFdS_Ten_Sig3_TxCom(dFdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

dFdS(1) =  0.0d0
dFdS(2) = -0.5d0
dFdS(3) = -0.5d0

end subroutine caldFdS_Ten_Sig3_TxCom

!--------------------------------------------------------------------------------------------------
subroutine caldFdS_Ten_Sig3_TxExt(dFdS)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

dFdS(1) =  0.0d0
dFdS(2) =  0.0d0
dFdS(3) = -1.0d0

end subroutine caldFdS_Ten_Sig3_TxExt

!--------------------------------------------------------------------------------------------------
! returns f of tensile yield surface for triaxial compression condition
double precision function calF_Ten_Sig3_TxComp(sigTen, prSig) result(f)
implicit none
double precision, intent(in):: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f = sigTen - (prSig(2) + prSig(3)) * 0.5d0

end function calF_Ten_Sig3_TxComp

!--------------------------------------------------------------------------------------------------
! returns f of tensile yield surface for triaxial compression condition
double precision function calF_Ten_Sig3_TxExt(sigTen, prSig) result(f)
implicit none
double precision, intent(in):: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f = sigTen - prSig(3)

end function calF_Ten_Sig3_TxExt

!--------------------------------------------------------------------------------------------------
! returns f of tensile yield surfaces based on sig1
double precision function calF_Ten_Sig1(sigTen, prSig) result(f)
implicit none
double precision, intent(in):: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f = sigTen - prSig(1)

end function calF_Ten_Sig1

!--------------------------------------------------------------------------------------------------
! returns f of tensile yield surfaces based on sigMin
double precision function calF_Ten_Min_Sig(sigTen, prSig) result(f)
implicit none
double precision, intent(in):: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
integer :: indexMin

indexMin = minloc(prSig, 1)

f = sigTen - prSig(indexMin)

end function calF_Ten_Min_Sig

!--------------------------------------------------------------------------------------------------
! returns f of tensile yield surfaces based on sig3
double precision function calF_Ten_Sig3_General(sigTen, prSig) result(f)
implicit none
double precision, intent(in):: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f = sigTen - prSig(3)

end function calF_Ten_Sig3_General

!--------------------------------------------------------------------------------------------------
! returns f of tensile yield surfaces based on p
double precision function calF_Ten_P(sigTen, prSig) result(f)
implicit none
double precision, intent(in):: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f = sigTen - calP(prSig)

end function calF_Ten_P

!--------------------------------------------------------------------------------------------------
!--- Mohr-Coulomb model functions: ----------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
integer function FindActivePart_MC(prSig) result(res)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! By default use General condition
res = INDEX_COULOMB_YIELD_GENERAL

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  if (prSig(2) < prSig(3) .or. abs(prSig(2) - prSig(3)) < SMALL) then
    ! Triaxial compression
    if (.not. (abs(prSig(2)) < SMALL .and. abs(prSig(3)) < SMALL)) then
      res = INDEX_COULOMB_TRIAXIAL_COMP
    endif
  elseif (prSig(1) < prSig(2) .or. abs(prSig(1) - prSig(2)) < SMALL) then
    ! Triaxial extension
    if (.not. (abs(prSig(1)) < SMALL .and. abs(prSig(2)) < SMALL)) then
      res = INDEX_COULOMB_TRIAXIAL_EXTN
    endif
  else
    ! General condition
    res = INDEX_COULOMB_YIELD_GENERAL
  endif
endif

end function FindActivePart_MC

!--------------------------------------------------------------------------------------------------
! Yield funcion based on Mohr-Coulomb criterion
double precision function calF_MC(sinPhi, pApex, prSig, indexActivePartMC_I) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    ! Triaxial compression
    f = calF_MC_TxComp(sinPhi, pApex, prSig)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    ! Triaxial extension
    f = calF_MC_TxExt(sinPhi, pApex, prSig)
  case(INDEX_COULOMB_YIELD_GENERAL)
    ! between triaxial compression and extension
    f = calF_MC_General(sinPhi, pApex, prSig)
  case default
    call SetError('unknown indexActivePartMC in calF_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  ! between triaxial compression and extension
  f = calF_MC_General(sinPhi, pApex, prSig)
endif

end function calF_MC

!--------------------------------------------------------------------------------------------------
! Yield funcion based on Mohr-Coulomb criterion
! for the general part (between triaxial compression and extension)
double precision function calF_MC_General(sinPhi, pApex, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f =             0.5d0 * (prSig(1) - prSig(3)) &
    - sinPhi * (0.5d0 * (prSig(1) + prSig(3)) + pApex)

end function calF_MC_General

!--------------------------------------------------------------------------------------------------
! Yield funcion based on Mohr-Coulomb criterion
! for the triaxial compression part of the model
double precision function calF_MC_TxComp(sinPhi, pApex, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f =             0.25d0 * (2.0d0*prSig(1) - prSig(2) - prSig(3)) &
    - sinPhi * (0.25d0 * (2.0d0*prSig(1) + prSig(2) + prSig(3)) + pApex)

end function calF_MC_TxComp

!--------------------------------------------------------------------------------------------------
! Yield funcion based on Mohr-Coulomb criterion
! for the triaxial extension part of the model
double precision function calF_MC_TxExt(sinPhi, pApex, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

f =             0.25d0 * (prSig(1) + prSig(2) - 2.0d0*prSig(3)) &
    - sinPhi * (0.25d0 * (prSig(1) + prSig(2) + 2.0d0*prSig(3)) + pApex)

end function calF_MC_TxExt

!--------------------------------------------------------------------------------------------------
! derivative of Mohr-Coulomb criterion with respect to stress
subroutine caldFdS_MC(sinPhi, dFdS, indexActivePartMC_I)
implicit none
double precision, intent(in) :: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    ! Triaxial compression
    call cal_dFdS_MC_TxComp(sinPhi, dFdS)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    ! Triaxial extension
    call cal_dFdS_MC_TxExt(sinPhi, dFdS)
  case(INDEX_COULOMB_YIELD_GENERAL)
    ! between triaxial compression and extension
    call caldFdS_MC_General(sinPhi, dFdS)
  case default
    call SetError('unknown indexActivePartMC in caldFdS_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  ! between triaxial compression and extension
  call caldFdS_MC_General(sinPhi, dFdS)
endif

end subroutine caldFdS_MC

!--------------------------------------------------------------------------------------------------
! derivative of Mohr-Coulomb criterion with respect to stress
! for the triaxial compression part of the model
subroutine cal_dFdS_MC_TxComp(sinPhi, dFdS)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

dFdS(1) = 0.25d0*( 2.0 - 2.0d0*sinPhi)
dFdS(2) = 0.25d0*(-1.0 - 1.0d0*sinPhi)
dFdS(3) = 0.25d0*(-1.0 - 1.0d0*sinPhi)

end subroutine cal_dFdS_MC_TxComp

!--------------------------------------------------------------------------------------------------
! derivative of Mohr-Coulomb criterion with respect to stress
! for the general part (between triaxial compression and extension)
subroutine caldFdS_MC_General(sinPhi, dFdS)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

dFdS(1) = 0.5d0*( 1.0 - sinPhi)
dFdS(2) = 0.0
dFdS(3) = 0.5d0*(-1.0 - sinPhi)

end subroutine caldFdS_MC_General

!--------------------------------------------------------------------------------------------------
! derivative of Mohr-Coulomb criterion with respect to stress
! for the triaxial compression part of the model
subroutine cal_dFdS_MC_TxExt(sinPhi, dFdS)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

dFdS(1) = 0.25d0*( 1.0 - 1.0d0*sinPhi)
dFdS(2) = 0.25d0*( 1.0 - 1.0d0*sinPhi)
dFdS(3) = 0.25d0*(-2.0 - 2.0d0*sinPhi)

end subroutine cal_dFdS_MC_TxExt

!--------------------------------------------------------------------------------------------------
double precision function calSinPhiMob_MC(pApex, prSig, indexActivePartMC_I) result(sinPhiMob)
implicit none

double precision, intent(in) :: pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
integer, optional, intent(in) :: indexActivePartMC_I

! local variables:
integer :: indexActivePartMC

if (USE_SEPARATED_MOHR_COULOMB_ZONES) then
  indexActivePartMC = INDEX_COULOMB_YIELD_GENERAL
  if (present(indexActivePartMC_I)) indexActivePartMC = indexActivePartMC_I
  select case(indexActivePartMC)
  case(INDEX_COULOMB_TRIAXIAL_COMP)
    ! Triaxial compression
    sinPhiMob = calSinPhiMob_MC_TxComp(pApex, prSig)
  case(INDEX_COULOMB_TRIAXIAL_EXTN)
    ! Triaxial extension
    sinPhiMob = calSinPhiMob_MC_TxExt(pApex, prSig)
  case(INDEX_COULOMB_YIELD_GENERAL)
    ! between triaxial compression and extension
    sinPhiMob = calSinPhiMob_MC_General(pApex, prSig)
  case default
    call SetError('unknown indexActivePartMC in calSinPhiMob_MC: ' // trim(string(indexActivePartMC)))
  end select
else
  ! between triaxial compression and extension
  sinPhiMob = calSinPhiMob_MC_General(pApex, prSig)
endif

!sinPhiMob = max(sinPhiMob, SMALL)
!sinPhiMob = min(sinPhiMob, 1.0d0)

end function calSinPhiMob_MC

!--------------------------------------------------------------------------------------------------
double precision function calSinPhiMob_MC_General(pApex, prSig) result(sinPhiMob)
implicit none

double precision, intent(in) :: pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

! local variables:
double precision:: denom

denom = prSig(1) + prSig(3) + 2.0*pApex

if (abs(denom) > TINY) then
  sinPhiMob = (prSig(1) - prSig(3)) / (denom + TINY)
else
  sinPhiMob = 1.0d0
endif

sinPhiMob = max(sinPhiMob, SMALL)
sinPhiMob = min(sinPhiMob, 1.0d0)

end function calSinPhiMob_MC_General

!--------------------------------------------------------------------------------------------------
double precision function calSinPhiMob_MC_TxComp(pApex, prSig) result(sinPhiMob)
implicit none

double precision, intent(in) :: pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

! local variables:
double precision:: denom

denom = 0.25d0*(2.0d0*prSig(1) + prSig(2) + prSig(3)) + pApex

if (abs(denom) > TINY) then
  sinPhiMob = 0.25d0*(2.0d0*prSig(1) - prSig(2) - prSig(3)) / (denom + TINY)
else
  sinPhiMob = 1.0d0
endif

sinPhiMob = max(sinPhiMob, SMALL)
sinPhiMob = min(sinPhiMob, 1.0d0)

end function calSinPhiMob_MC_TxComp

!--------------------------------------------------------------------------------------------------
double precision function calSinPhiMob_MC_TxExt(pApex, prSig) result(sinPhiMob)
implicit none

double precision, intent(in) :: pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

! local variables:
double precision:: denom

denom = 0.25d0*(prSig(1) + prSig(2) + 2.0d0*prSig(3)) + pApex

if (abs(denom) > TINY) then
  sinPhiMob = 0.25d0*(prSig(1) + prSig(2) - 2.0d0*prSig(3)) / (denom + TINY)
else
  sinPhiMob = 1.0d0
endif

sinPhiMob = max(sinPhiMob, SMALL)
sinPhiMob = min(sinPhiMob, 1.0d0)

end function calSinPhiMob_MC_TxExt

!--------------------------------------------------------------------------------------------------
!--- Generalised Mohr-Coulomb model functions: ----------------------------------------------------
!--------------------------------------------------------------------------------------------------
! Grimstad, Ronningen & Nordal (2018): NUMGE2018
! Application of a generalized continious Mohr-coulomb criterion
double precision  function getSinPhi0GMC(sinPhi, A1, A2) result(sinPhi0)
implicit none
double precision, intent(in):: sinPhi, A1, A2
double precision, parameter :: ONE_THIRD = 1.0d0/3.0d0
double precision beta

beta = ONE_THIRD * asin(A1)

sinPhi0 = 2.0d0 * sqrt(3.0d0) * A2 * cos(beta) * sinPhi / (3.0d0*A2 - (A2-2.0d0*sin(beta))*sinPhi)

end function getSinPhi0GMC

!--------------------------------------------------------------------------------------------------
double precision function getSinPhi0MN(sinPhi) result(sinPhi0)
implicit none
double precision, intent(in):: sinPhi

sinPhi0 = 2.0d0 * sinPhi / sqrt(3.0d0 + sinPhi*sinPhi)

end function getSinPhi0MN

!--------------------------------------------------------------------------------------------------
double precision function getA1MN(sinPhi0) result(A1)
implicit none
double precision, intent(in):: sinPhi0

A1 = sinPhi0 * (3.0d0 - sinPhi0*sinPhi0) * 0.5d0

end function getA1MN

!--------------------------------------------------------------------------------------------------
!--- Cam-Clay model functions: --------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! return f of cap yield functions
double precision function calF_CamClay(xM, pc, pApex, prSig) result(f)
implicit none
double precision, intent(in):: xM, pc, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

! local variables:
double precision p, q

p = calP(prSig) + pApex
q = calQ(prSig)


if (USE_ALTERNATIVE_CAM_CLAY) then
  !f = (1/xM^2) * (q^2/p) + p - pc
  f = (1.0d0/xM**2) * (q**2/(p+TINY)) + p - pc
else
  !f = (q/xM)**2 + p**2 -p*pc
  f = (q/xM)**2 + p * (p - pc)
endif

end function calF_CamClay

!--------------------------------------------------------------------------------------------------
! return df/dSigma of Cam-Clay yield functions
subroutine caldFdS_CamClay(xM, pc, pApex, prSig, dfdSigma)
implicit none
double precision, intent(in)::  xM, pc, pApex
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out) :: dfdSigma

! local variables:
double precision:: p, q
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dPdS, dQds

p = calP(prSig)
q = calQ(prSig)

call getdPdSigma(dPdS)
call getdQdSigma(p, q, prSig, dQds)

p = p + pApex

if (USE_ALTERNATIVE_CAM_CLAY) then
  !f = (1/xM^2) * (q^2/p) + p - pc
  dfdSigma = (1.0d0/xM**2) * (2.0d0*q*dQds*p - dPdS * (q**2)) / (p**2 + TINY) + dPdS
else
  !f = (q/xM)**2 + p**2 -p*pc
  dfdSigma = (1.0d0/xM**2)*2.0d0*q*dQds + 2.0d0*p*dPdS - pc*dPdS
endif

end subroutine caldFdS_CamClay

!--------------------------------------------------------------------------------------------------
! derivative of f with respect to hardening parameter for Cam-Clay model
double precision function calDfDh_CamClay(prSig) result(dFdh)
implicit none
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

if (USE_ALTERNATIVE_CAM_CLAY) then
  !f = (1/xM^2) * (q^2/p) + p - pc
  dFdh = - 1.0d0
else
  !f = (q/xM)**2 + p**2 -p*pc
  dFdh = - calP(prSig)
endif

end function calDfDh_CamClay

!--------------------------------------------------------------------------------------------------
! returns f of generalized Mohr-Coulomb (GMC) yield functions
double precision function calGTheta_GMC(sinPhi, generalisedA1, generalisedA2, prSig) result(GTheta)
implicit none
double precision, intent(in):: sinPhi
double precision, intent(in):: generalisedA1, generalisedA2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig

! local variables:
double precision J2, J3, theta, KTheta, p
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)

theta = getLodeAngle(J2,J3)

theta = (1.0/3.0)*asin(generalisedA1*sin(3.0*theta))

KTheta   = getKTheta(theta, sinPhi, generalisedA2)

GTheta = KTheta / sinPhi

end function calGTheta_GMC

!-----------------------------------------------------------------
! returns f of generalized Mohr-Coulomb (GMC) yield functions
double precision function getF_GMC(sinPhi, apexHyperplolic, pApex, generalisedA1, generalisedA2, prSig) result(f)
implicit none
double precision, intent(in):: sinPhi, apexHyperplolic, pApex
double precision, intent(in):: generalisedA1, generalisedA2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision J2, J3, theta, KTheta, p, a
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

a = apexHyperplolic

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)

theta = getLodeAngle(J2,J3)

theta = (1.0/3.0)*asin(generalisedA1*sin(3.0*theta))

KTheta   = getKTheta(theta, sinPhi, generalisedA2)


!f = -p + sqrt( (J2 * KTheta*KTheta) / (sinPhi*sinPhi)  + a*a ) - pApex
f = -p*sinPhi + sqrt( (J2 * KTheta*KTheta) + a*a * (sinPhi*sinPhi)) - pApex*sinPhi

end function getF_GMC

!-----------------------------------------------------------------
! returns getGTheta_MNC of Convex Matsoua-Nakai yield surfaces
! paper: Pantephini and Lagioia (2014)
double precision function getGTheta_MNC(sinPhi, prSig) result(GTheta)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

! local variables:
double precision J2, J3, theta, p, ksi, eta
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)

theta = getLodeAngle(J2,J3)

ksi    = calKsi_MNC(sinPhi)
GTheta = calGTheta_MNC(theta, ksi)
eta    = calEta_MNC(sinPhi)

!f = -p*eta + sqrt( (J2 * GTheta*GTheta) + a*a * eta*eta) - pApex*eta
GTheta = GTheta / eta

end function getGTheta_MNC

!-----------------------------------------------------------------
! returns getGTheta_MC of Mohr-Coulomb yield surface
! paper: Pantephini and Lagioia & Panteghini (2016)
double precision function getGTheta_MC(sinPhi, prSig) result(GTheta)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig

double precision J2, J3, theta, p
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)

theta = getLodeAngle(J2,J3)

GTheta = (cos(theta) + (sin(theta)*sinPhi) / sqrt(3.0d0)) / (sinPhi + TINY)

end function getGTheta_MC

!-----------------------------------------------------------------
! returns getGTheta_DP of Mohr-Coulomb yield surface
! paper: Pantephini and Lagioia & Panteghini (2016)
double precision function getGTheta_DP(sinPhi) result(GTheta)
implicit none
double precision, intent(in)::sinPhi

! local variables:
double precision :: CompressionM

CompressionM = 6.0 * sinPhi / (3.0 - sinPhi)

GTheta = sqrt(3.0d0) / CompressionM

end function getGTheta_DP

!-----------------------------------------------------------------
subroutine getdFdSGMC(sinPhi, apexHyperplolic, generalisedA1, generalisedA2, prSig, dFdS)
implicit none
double precision, intent(in) :: sinPhi, apexHyperplolic
double precision, intent(in) :: generalisedA1, generalisedA2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
double precision :: J2, J3, theta, KTheta, dKdTheta, sigmaBar, alpha, a, sigmaK
double precision :: C1, C2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S, dSigmaBardSigma, dPdSigma

a = apexHyperplolic

call getPrincipalDeviatoricS(prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)
theta = getLodeAngle(J2,J3)

if (generalisedA1 > 0.0d0) then
  theta = (1.0d0/3.0d0)*asin(generalisedA1*sin(3.0d0*theta))
  KTheta   = getKTheta(theta, sinPhi, generalisedA2)
  dKdTheta = getdKdTheta(theta, sinPhi, generalisedA2)
else
  theta    = 0.0
  KTheta   = 1.0
  dKdTheta = 0.0
endif

sigmaBar = sqrt(J2)

sigmaK = sigmaBar * KTheta
alpha = sigmaK / sqrt(sigmaK * sigmaK + a*a * sinPhi*sinPhi)

C1 = - sinPhi
call getdPdSigma(dPdSigma)

C2 = KTheta - tan(3.0d0*theta)*dKdTheta
call getdSigmaBardSigma(J2, S, dSigmaBardSigma)

dFdS = C1*dPdSigma + alpha*C2*dSigmaBardSigma

end subroutine getdFdSGMC

!--------------------------------------------------------------------------------------------
double precision function getKTheta(theta, sinPhi, generalisedA2) result(KTheta)
implicit none
double precision, intent(in):: theta, sinPhi, generalisedA2

KTheta = cos(theta) + SIGN_FACTOR_GMC * (1.0d0/sqrt(3.0d0)) * sinPhi * (sin(theta)/generalisedA2)

end function getKTheta

!--------------------------------------------------------------------------------------------
double precision function getdKdTheta(theta, sinPhi, generalisedA2) result(dKdTheta)
implicit none
double precision, intent(in):: theta, sinPhi, generalisedA2

dKdTheta = -sin(theta) + SIGN_FACTOR_GMC * (1.0d0/sqrt(3.0d0)) * sinPhi * (cos(theta)/generalisedA2)

end function getdKdTheta

!-----------------------------------------------------------------
subroutine getdFdSDPLodeAngle(sinPhi, prSig, dFdS)
implicit none
double precision, intent(in) :: sinPhi
double precision, intent(In) :: prSig(N_PRINCIPAL_STRESS_VECTOR)
double precision, intent(out):: dFdS(N_PRINCIPAL_STRESS_VECTOR)

! local variables:
double precision p, q, theta, M, J2, J3
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dPdS, dQdS
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dSin3ThetadSigma
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dMdS

p = calP(prSig)
q = calQ(prSig)

call getdPdSigma(dPdS)
call getdQdSigma(p, q, prSig, dQds)

theta = getLodeAngle(prSig)
! using Lode's angle at triaxial compression
M = 6.0d0*sinPhi/(3.0d0 - sinPhi*sin(3.0d0*theta))

if (.false.) then
  call getPrincipalDeviatoricSWithP(p, prSig, S)
  J2 = getJ2(S)
  J3 = getJ3(S)
  call getdSin3ThetadSigma(S, J2, J3, dSin3ThetadSigma)
  dMdS = M*M*dSin3ThetadSigma
else
  dMdS = 0.0d0
endif

dFdS = dQds - dMdS * p - M * dPdS

end subroutine getdFdSDPLodeAngle

!-----------------------------------------------------------------
double precision function getFDPKinematic(xM, pApex, prSig, prBackStress) result(f)
implicit none
double precision, intent(in):: xM, pApex
double precision, intent(in):: prSig(N_PRINCIPAL_STRESS_VECTOR), prBackStress(N_PRINCIPAL_STRESS_VECTOR)

double precision f0, p, q
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

S = S - p*prBackStress

q = sqrt(3.0d0*getJ2(S))

f0 = pApex * xM

f = q - xM*p - f0

end function getFDPKinematic

!-----------------------------------------------------------------
subroutine getdFdSDPKinematic(xM, prSig, prBackStress, dFdS)
implicit none
double precision, intent(in):: xM
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)::  prSig, prBackStress
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS

! local variables:
double precision p, J2
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: S
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dPdSigma, dqdSigma
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dJ2dSigma

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

S  = S - p*prBackStress
J2 = getJ2(S)

call getdPdSigma(dPdSigma)
if (.true.) then
  call getdQSigmaKinematic(J2, S, prBackStress, dqdSigma)
  dFdS = dqdSigma - xM * dPdSigma
else
  call getdJ2dSigmaKinematic(S, prBackStress, dJ2dSigma)
  dFdS = 3.0d0*dJ2dSigma - 2.0d0 * xM * xM * p* dPdSigma
endif

end subroutine getdFdSDPKinematic

!-----------------------------------------------------------------
subroutine getdGtdS(DInverse, prSig, prSigApex, dGtdS)
implicit none
double precision, intent(in) :: DInverse(N_PRINCIPAL_STRESS_VECTOR,N_PRINCIPAL_STRESS_VECTOR)
double precision, intent(in) :: prSig(N_PRINCIPAL_STRESS_VECTOR), prSigApex(N_PRINCIPAL_STRESS_VECTOR)
double precision, intent(out):: dGtdS(N_PRINCIPAL_STRESS_VECTOR)

! local variables:
double precision:: vectorLength
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: dSig

dSig = prSig - prSigApex

! dg/dS = (D^-1)*(Sig_trial - Sig_Appex)
call MatrixToVectorProd(DInverse, dSig, N_PRINCIPAL_STRESS_VECTOR, dGtdS)

vectorLength = prStressLength(dGtdS)
dGtdS(:) = dGtdS(:) / (vectorLength + TINY)

end subroutine getdGtdS

!-----------------------------------------------------------------
! derivative of f with respect to hardening parameter for Generalised Mohr-Coulomb model
double precision function getDfDhGMC(sinPhi, apexHyperplolic, pApex, generalisedA1, generalisedA2, prSig) result(dfDh)
implicit none
double precision, intent(in):: sinPhi, apexHyperplolic, pApex
double precision, intent(in):: generalisedA1, generalisedA2
double precision, intent(in):: prSig(N_PRINCIPAL_STRESS_VECTOR)

double precision J2, J3, theta, KTheta, p, a, cosPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: S

cosPhi = cos(asin(sinPhi))
a = apexHyperplolic

p = calP(prSig)
call getPrincipalDeviatoricSWithP(p, prSig, S)

J2 = getJ2(S)
J3 = getJ3(S)

theta = getLodeAngle(J2,J3)
theta = (1.0/3.0)*asin(generalisedA1*sin(3.0*theta))

KTheta = getKTheta(theta, sinPhi, generalisedA2)

dfDh = -p + (a*a * sinPhi)/sqrt( (J2 * KTheta*KTheta) + a*a * (sinPhi*sinPhi)) - pApex

end function getDfDhGMC

!-----------------------------------------------------------------
! derivative of f with respect to hardening parameter for Mohr-Coulomb model
double precision function getDfDhMC(pApex, prSig) result(dFdh)
implicit none
double precision, intent(in):: pApex
double precision, intent(in):: prSig(N_PRINCIPAL_STRESS_VECTOR)

dFdh = - (0.5*(prSig(1) + prSig(3)) + pApex)

end function getDfDhMC

!-----------------------------------------------------------------
! complementary function (area which bring stress to apex)
double precision function calFComplementary(sinPhi, prSig, pApex) result(f)
implicit none
double precision, intent(in):: sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, intent(in):: pApex

double precision:: sinPhiComplementary, phi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: prStressNeg

if (sinPhi > 1.0 .or. sinPhi < -1.0) then
  call SetError('wrong sinPhi:' // trim(string(sinPhi)))
endif

if (sinPhi > 0.0d0 .and. calP(prSig) < pApex) then
  phi = asin(sinPhi)
  phi = PI*0.5d0 - phi
  sinPhiComplementary = sin(phi)
  prStressNeg = - prSig
  f = -calF_DP(sinPhiComplementary, pApex, prStressNeg)
else
  phi = asin(sinPhi)
  phi = phi + PI*0.5d0
  sinPhiComplementary = sin(phi)
  f = calF_DP(sinPhiComplementary, pApex, prSig)
endif

end function calFComplementary

end module ModYieldSurfaceLibrary