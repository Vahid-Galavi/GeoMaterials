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
module ElasticViscoPlastic_Module

use ModGlobalConstants
use ModMathLibrary
use ModDebug
use ModString
use ModYieldSurfaceLibrary
use ModSolversLibrary

implicit none

public:: ElasticViscoPlastic
public:: GetParamAndUnitElasticViscoPlastic
public:: GetParamCountElasticViscoPlastic
public:: GetStateVarCountElasticViscoPlastic
public:: GetModelNameElasticViscoPlastic

public:: NSTATS
public:: NPARAMS
public:: INDEX_PROP_E
public:: INDEX_PROP_POISSON_UR
public:: INDEX_PROP_COHESION
public:: INDEX_PROP_FRICTION_PEAK
public:: INDEX_PROP_DILATION_PEAK
public:: INDEX_PROP_SIG_T_CUT_OFF
public:: INDEX_PROP_YIELD_FUNCTION
public:: INDEX_PROP_UNDRAINED_NU
public:: STVAR_LOADING_TYPE

private

! flags
! CONSIDER_INITIAL_SUB_STEPPING should be true for MATSUOKA_NAKAI_YIELD. Needs to be checked!
logical, parameter:: CONSIDER_INITIAL_SUB_STEPPING      = .false.
logical, parameter:: CALCULATE_NSUB_BASED_ON_EPS        = .false.
logical, parameter:: CONSIDER_SUB_STEPPING              = .true.
logical, parameter:: REPEAT_THE_WHOLE_STEP              = .false.
logical, parameter:: USE_COUPLED_HARDENING_CROSS_POINTS = .true.

logical, parameter:: USE_DEFAULT_PLASTIC_POTENTIAL      = .false.
logical, parameter:: CHECK_ALL_YIELD_SURFACES_CONVERGED = .false.
logical, parameter:: ABORT_WHEN_NOT_CONVERGED           = .true.
logical, parameter:: WRITE_DEBUG_FILE_WHEN_NOT_CONVERGED= .false.

! limit constants
integer, parameter:: ITER_MAX = 100 ! MN models need more iterations. Needs to be checked!
integer, parameter:: N_SUB_MAX_MC = N_SUB_MAX

! indices constants
integer, parameter:: INDEX_MOHR_COULOMB_YIELD = 1
integer, parameter:: INDEX_TENSION_YIELD      = 2
integer, parameter:: INDEX_APEX_YIELD         = 3
integer, parameter:: N_YIELD_FUNCTIONS        = 3

! sub stepping constants:
double precision, parameter:: MIN_DEPS  = 1d-15

integer, parameter:: DEFAULT_PLASTIC_POTENTIAL           = DRUCKER_PRAGER_YIELD
integer, parameter:: DEFAULT_PLASTIC_POTENTIAL_EXTENSION = DRUCKER_PRAGER_YIELD_EXTENSION


! indices of property array
integer, parameter:: INDEX_PROP_E                  = 1
integer, parameter:: INDEX_PROP_POISSON_UR         = 2

integer, parameter:: INDEX_PROP_COHESION           = 3
integer, parameter:: INDEX_PROP_P_APEX             = 3

integer, parameter:: INDEX_PROP_FRICTION_PEAK      = 4
integer, parameter:: INDEX_PROP_SIN_FRICTION_PEAK  = 4

integer, parameter:: INDEX_PROP_DILATION_PEAK      = 5
integer, parameter:: INDEX_PROP_SIN_DILATION_PEAK  = 5

integer, parameter:: INDEX_PROP_SIG_T_CUT_OFF      = 6

integer, parameter:: INDEX_PROP_YIELD_FUNCTION     = 7
integer, parameter:: INDEX_PROP_UNDRAINED_NU       = 8

! viscosity
integer, parameter:: INDEX_PROP_VISCOSITY          = 9
integer, parameter:: INDEX_PROP_RATE_SENSITIVITY   = 10
integer, parameter:: INDEX_PROP_PERZYNA_F0         = 11

! total number of parameters
integer, parameter:: NPARAMS     = 11

!>    State variables
integer, parameter:: STVAR_EXCESS_PORE_PRESSURE = 1
integer, parameter:: STVAR_LOADING_TYPE         = 2
integer, parameter:: STVAR_LODE_ANGLE           = 3
integer, parameter:: STVAR_INITIALIZATION       = 4

integer, parameter:: NSTATS = STVAR_INITIALIZATION

contains

!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef IN_DEBUG
subroutine ElasticViscoPlastic( IDTask, isUndr,               &
                                    iStep, iTer, iEl, intPt,      &
                                    Time0, dTimeI,                &
                                    PropsI, SigI0, Swp0, stVarI0, &
                                    dEpsI, D, BulkW,              &
                                    Sig, Swp, stVar, ipl,         &
                                    nStat,                        &
                                    NonSym, iStrsDep, iTimeDep, iTang, &
                                    iAbort )

implicit none
integer, intent(in):: IDTask, isUndr, iStep, iTer, iEl, intPt
integer, intent(inout):: nStat
integer, intent(out)::  ipl, nonSym, iStrsDep, iTimeDep, iTang, iAbort
double precision, intent(in):: swp0, Time0, dTimeI
double precision, intent(inout):: swp, BulkW
double precision, intent(in):: propsI(NPARAMS), sigI0(*), dEpsI(*)
double precision, intent(inout):: sig(*), stVar(NSTATS), D(N_STRESS_VECTOR,N_STRESS_VECTOR), stVarI0(NSTATS)
double precision:: props(NPARAMS), stVar0(NSTATS), dEps(N_STRESS_VECTOR), sig0(N_STRESS_VECTOR)

props = (/100000.000000000e0      ,  0.300000000000000e0      ,  0.000000000000000E+000 ,   25.0000000000000e0      ,   25.0000000000000e0      ,  0.000000000000000E+000 ,   2.00000000000000      ,  0.000000000000000E+000/)

stVar0 = (/ 0.000000000000000E+000 ,   6.00000000000000      ,  -22.4109105310250      ,   1.00000000000000 /)

dEps = (/ 1.253547187174091E-002 ,  9.753484976574476E-003 ,  1.521177880296275E-004 ,  1.910702769386128E-002 ,  0.000000000000000E+000 ,  0.000000000000000E+000 /)

sig0 = (/ 2.515883959809619E-013 ,  2.031589549055022E-013 ,  0.000000000000000E+000 ,  1.688022977312425E-013 ,  0.000000000000000E+000 ,  0.000000000000000E+000 /)

call ElasticViscoPlastic_I(IDTask, isUndr,                    &
                               iStep, iTer, iEl, intPt,           &
                               Time0, dTimeI,                     &
                               Props, Sig0, Swp0, stVar0,         &
                               dEps, D, BulkW,                    &
                               Sig, Swp, stVar, ipl,              &
                               nStat,                             &
                               NonSym, iStrsDep, iTimeDep, iTang, &
                               iAbort )

end subroutine ElasticViscoPlastic_Temp

!--------------------------------------------------------------------------------------------------
subroutine ElasticViscoPlastic_I(IDTask, isUndr,                    &
                                     iStep, iTer, iEl, intPt,           &
                                     Time0, dTimeI,                     &
                                     PropsI, SigI0, Swp0, stVar0,       &
                                     dEps, D, BulkW,                    &
                                     Sig, Swp, stVar, ipl,              &
                                     nStat,                             &
                                     NonSym, iStrsDep, iTimeDep, iTang, &
                                     iAbort )

#else
subroutine ElasticViscoPlastic(IDTask, isUndr,                    &
                                   iStep, iTer, iEl, intPt,           &
                                   Time0, dTimeI,                     &
                                   PropsI, SigI0, Swp0, stVar0,       &
                                   dEps, D, BulkW,                    &
                                   Sig, Swp, stVar, ipl,              &
                                   nStat,                             &
                                   NonSym, iStrsDep, iTimeDep, iTang, &
                                   iAbort )
#endif

implicit none
!arguments:
! Task ID (see Plaxis manual or SolversLib.f90)
integer, intent(in):: IDTask
! undrained = 1, drained = 0
integer, intent(in):: isUndr
! global step number:
integer, intent(in):: iStep
! iteration number:
integer, intent(in):: iTer
! element number
integer, intent(in):: iEl
! local integration point number:
integer, intent(in):: intPt
! number of state variables:
integer, intent(inout):: nStat
! plasticity indicator (see Plaxis manual or SolversLib.f90)
integer, intent(out)  :: ipl
! is stiffness matrix non-symmetric? yes = 1
integer, intent(out)  :: nonSym
! is stiffness matrix stress-dependent? yes = 1
integer, intent(out)  :: iStrsDep
! is stiffness matrix time-dependent? yes = 1
integer, intent(out)  :: iTimeDep
! is stiffness matrix tangential? yes = 1
integer, intent(out)  :: iTang
! abort calculation, when set to 1
integer, intent(out)  :: iAbort
! excess pore pressure at the end of previous step:
double precision, intent(in) :: swp0
! time at the end of previous step:
double precision, intent(in) :: Time0
! increment of time:
double precision, intent(in) :: dTimeI
! current excess pore pressure for undrained calculations:
double precision, intent(inout):: swp
! current equivalent bulk modulus of water K/n for undrained calculations:
double precision, intent(inout):: BulkW
! vector of material properties:
double precision, dimension(NPARAMS), intent(in)   :: propsI
! stress at the end of previous step:
double precision, dimension(N_STRESS_VECTOR), intent(in)   :: sigI0
! current stress:
double precision, dimension(N_STRESS_VECTOR), intent(inout):: sig
! increment of strain:
double precision, dimension(N_STRESS_VECTOR), intent(in)   :: dEps
! state variables at the end of previous step:
double precision, dimension(NSTATS), intent(inout):: stVar0
! current state variables:
double precision, dimension(NSTATS), intent(inout):: stVar
! stiffness matrix:
double precision, dimension(N_STRESS_VECTOR, N_STRESS_VECTOR), intent(inout):: D

! Local variables
integer, parameter :: INT_DEBUG_LOCAL = -1
integer:: nSub, iSub
double precision:: GCurrent
double precision:: dTimeCurrent, dTime
double precision:: dSwp, dEpsV, dEpsQ, dTimeSub
double precision, dimension(N_STRESS_VECTOR,N_STRESS_VECTOR) :: DE
double precision, dimension(N_STRESS_VECTOR):: dSig, SigI, Sig0
double precision, dimension(N_STRESS_VECTOR):: dEpsSub
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: prSigElas, prSig, prSig0
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR):: xN1, xN2, xN3
double precision, dimension(NPARAMS):: props
double precision, dimension(NSTATS):: stVarSub0
double precision, dimension(N_YIELD_FUNCTIONS):: fYield, fYieldEnd
character(len=1023):: subroutineCalled
logical:: isConverged
logical:: IsDebug
logical:: restart

IsDebug = setIsDebug(iEl, intPt, INT_DEBUG_LOCAL)

iAbort = 0
ipl = ELASTIC_LOADING

call CheckAndCopyProps(propsI, props)
Sig0(1:N_STRESS_VECTOR) = SigI0(1:N_STRESS_VECTOR)

dTime = dTimeI
if (dTime < TINY) dTime = 1.0d0

select case (IDTask)
case (IDTASK_INITIALISATION)
  if (stVar0(STVAR_INITIALIZATION) /= 1.0) then
    call initialiseStateVariables(stVar)
    stVar0(1:NSTATS) = stVar(1:NSTATS)
  endif

case (IDTASK_STRESS)
  ! Calculate stresses
  if (stVar0(STVAR_INITIALIZATION) /= 1.0) then
    call initialiseStateVariables(stVar)
    stVar0(1:NSTATS) = stVar(1:NSTATS)
  endif

  ! Fill elastic material matrix
  GCurrent = calGRef(props)
  call DMatrixDrained(GCurrent, props(INDEX_PROP_POISSON_UR), DE)

  if (isUndr == 1) then
    call CalculateExcessPorePressure(GCurrent,                       &
                                     props(INDEX_PROP_POISSON_UR),   &
                                     props(INDEX_PROP_UNDRAINED_NU), &
                                     Swp0, getEpsV(dEps),            &
                                     BulkW, Swp)
    stVar(STVAR_EXCESS_PORE_PRESSURE) = Swp
  else
    if (props(INDEX_PROP_UNDRAINED_NU) > props(INDEX_PROP_POISSON_UR)) then
      call CalculateExcessPorePressure(GCurrent,                       &
                                       props(INDEX_PROP_POISSON_UR),   &
                                       props(INDEX_PROP_UNDRAINED_NU), &
                                       0.0d0, getEpsV(dEps),           &
                                       BulkW, Swp)

      ! calculate the effective stress
      Sig0(1:3) = Sig0(1:3) - stVar0(STVAR_EXCESS_PORE_PRESSURE)
    endif
    BulkW = 0
  endif

  ! sub-stepping
  nSub = 1
  if (ConsiderInitialSubStepping(props)) then
    if (CALCULATE_NSUB_BASED_ON_EPS) then
      nSub = max(getnSubStepsBasedOnEpsq(dEps), nSub)
    else
      call calElasticPrincipalStress(Sig0,           &
                                     dEps,           &
                                     Sig, prSigElas, &
                                     DE,             &
                                     isDebug)
      ! from mechanical sign to soil sign convention
      prSigElas = -prSigElas
      nSub = max(getnSubStepsBasedOnSigCone(props(INDEX_PROP_SIN_FRICTION_PEAK), &
                                            props(INDEX_PROP_P_APEX),            &
                                            prSig0,                              &
                                            prSigElas), nSub)
    endif
    nSub = min(nSub, N_SUB_MAX_INITIAL)
  endif

  if (IsDebug) then
    call WriteInDebugFile('nSub   : ' // trim(String(nSub)))
  endif

  stVarSub0(1:NSTATS) = stVar0(1:NSTATS)

  restart = .true.
  do while (restart)
    restart = .false.
    dTimeSub = dTime / float(nSub)
    dEpsSub(1:N_STRESS_VECTOR) = dEps(1:N_STRESS_VECTOR) / float(nSub)
    SigI(1:N_STRESS_VECTOR) = Sig0(1:N_STRESS_VECTOR)

    iSub = 0

    dTimeCurrent = dTimeSub
    stVar(1:NSTATS) = stVarSub0(1:NSTATS)

    do while (dTimeCurrent <= dTime)
    !do iSub= 1, nSub
      iSub = iSub + 1

      if (IsDebug) then
        call WriteInDebugFile('iSub   :' // trim(string(iSub)))
        call WriteInDebugFile('dEpsSub:' // trim(string(dEpsSub,1,N_STRESS_VECTOR)))
      endif

      call calPrincipalStresses(SigI, prSig0)

      ! from mechanical sign to soil sign convention
      prSig0 = -prSig0

      ! elastic stress sub_increment
      call MatrixToVectorProd(DE, dEpsSub, N_STRESS_VECTOR, dSig)
      ! elastic stress
      Sig(1:N_STRESS_VECTOR) = SigI(1:N_STRESS_VECTOR) + dSig(1:N_STRESS_VECTOR)
      ! calculate principal stresses and directions
      call calPrincipalStressesAndVectors(Sig, xN1, xN2, xN3, prSigElas)

      ! from mechanical sign to soil sign convention
      prSigElas = -prSigElas

      call calStress(iStep, iEl, intPt,  &
                     dTimeSub, GCurrent, &
                     prSigElas, prSig0,  &
                     props,              &
                     stVar,              &
                     prSig,              &
                     fYield,             &
                     isConverged,        &
                     subroutineCalled )

      ipl = stVar(STVAR_LOADING_TYPE)
      stVar(STVAR_LODE_ANGLE) = -getLodeAngle(prSig)*DEGREES
      if (CHECK_ALL_YIELD_SURFACES_CONVERGED) then
        call checkYieldSurfaces(props, prSig, IsDebug, isConverged)
      endif

      if (.not.isConverged) then
        ! to write debug data:
        if (AUTOMATIC_DEBUG_FILE) then
          IsDebug = .true.
        endif
        if (IsDebug) call WriteInDebugFile('GaussPoint:' // trim(string(intPt)))
        if (nSub < N_SUB_MAX_MC .and. ConsiderSubStepping(props)) then
          if (REPEAT_THE_WHOLE_STEP) then
            nSub = min(UP_SCALE_FACTOR * nSub, N_SUB_MAX_MC)
            dEpsSub(1:N_STRESS_VECTOR) = dEps(1:N_STRESS_VECTOR) / float(nSub)
            dEpsV = getEpsV(dEpsSub)
            dEpsQ = getEpsQ(dEpsSub)
            if (abs(dEpsV) >= MIN_DEPS .or. abs(dEpsQ) >= MIN_DEPS) then
              if (IsDebug) call WriteInDebugFile('Repeat the whole step with new nSub:' // trim(string(nSub)))
              stVarSub0(1:NSTATS) = stVar0(1:NSTATS)
              restart = .true.
              EXIT
            else
              if (IsDebug) call WriteInDebugFile('discontinue sup-stepping at step:' // trim(string(iSub)))
            endif
          else
            dTimeCurrent = dTimeCurrent - dTimeSub
            nSub = min(UP_SCALE_FACTOR * nSub, N_SUB_MAX_MC)
            dTimeSub = dTime / float(nSub)
            dEpsSub(1:N_STRESS_VECTOR) = dEps(1:N_STRESS_VECTOR) / float(nSub)
            dEpsV = getEpsV(dEpsSub)
            dEpsQ = getEpsQ(dEpsSub)
            if (abs(dEpsV) >= MIN_DEPS .or. abs(dEpsQ) >= MIN_DEPS) then
              if (IsDebug) call WriteInDebugFile('continue the rest of the step with new nSub:' // trim(string(nSub)))
              dTimeCurrent = dTimeCurrent + dTimeSub
              stVar(1:NSTATS) = stVarSub0(1:NSTATS)
              CYCLE
            else
              if (IsDebug) call WriteInDebugFile('discontinue sup-stepping at step:' // trim(string(iSub)))
            endif
          endif
        endif

        if (ABORT_WHEN_NOT_CONVERGED .or. WRITE_DEBUG_FILE_WHEN_NOT_CONVERGED) then
          call calAllYieldFs(props, prSig, fYieldEnd)
          call WriteInDebugFile('Mohr-Coulomb Model not converged',                                 ALWAYS_LEVEL)
          call WriteInDebugFile('subroutine :' // trim(subroutineCalled),                           ALWAYS_LEVEL)
          call WriteInDebugFile('props      :' // trim(string(propsI,1,NPARAMS)),                   ALWAYS_LEVEL)
          call WriteInDebugFile('GaussPoint :' // trim(string(intPt)),                              ALWAYS_LEVEL)
          call WriteInDebugFile('iEl        :' // trim(string(iEl)),                                ALWAYS_LEVEL)
          call WriteInDebugFile('dEpsV      :' // trim(string(dEpsV)),                              ALWAYS_LEVEL)
          call WriteInDebugFile('dEpsQ      :' // trim(string(dEpsQ)),                              ALWAYS_LEVEL)
          call WriteInDebugFile('stVar0     :' // trim(string(stVar0, 1,NSTATS)),                   ALWAYS_LEVEL)
          call WriteInDebugFile('stVar      :' // trim(string(stVar0, 1,NSTATS)),                   ALWAYS_LEVEL)
          call WriteInDebugFile('dEps       :' // trim(string(dEps,   1,N_STRESS_VECTOR)),          ALWAYS_LEVEL)
          call WriteInDebugFile('SigI0      :' // trim(string(SigI0,  1,N_STRESS_VECTOR)),          ALWAYS_LEVEL)
          call WriteInDebugFile('dEpsSub    :' // trim(string(dEpsSub,1,N_STRESS_VECTOR)),          ALWAYS_LEVEL)
          call WriteInDebugFile('new nSub   :' // trim(string(nSub)),                               ALWAYS_LEVEL)
          call WriteInDebugFile('prSig      :' // trim(string(-prSig,1,N_PRINCIPAL_STRESS_VECTOR)), ALWAYS_LEVEL)
          call WriteInDebugFile('fYield     :' // trim(string(fYield,1,N_YIELD_FUNCTIONS)),         ALWAYS_LEVEL)
          call WriteInDebugFile('fYieldEnd  :' // trim(string(fYieldEnd,1,N_YIELD_FUNCTIONS)),      ALWAYS_LEVEL)
          call WriteInDebugFile('lodeAngle  :' // trim(string(-getLodeAngle(prSig)*DEGREES)),       ALWAYS_LEVEL)
          call WriteInDebugFile('isConverged:' // trim(string(isConverged)),                        ALWAYS_LEVEL)
        endif

        if (ABORT_WHEN_NOT_CONVERGED) then
          iAbort =1
          RETURN
        else
          if (IsDebug) then
            call WriteInDebugFile('time  :' // trim(string(time0 + dTimeI)))
            call WriteInDebugFile('p     :' // trim(string(calP(prSig))))
            call WriteInDebugFile('q     :' // trim(string(calQ(prSig))))
          endif
        endif
      endif

      prSig0 = prSig

      ! from soil sign to mechanical sign convention
      prSig = -prSig
      ! back to Cartesian stresses
      call PrincipalStressToCartesian(prSig, xN1, xN2, xN3, sig)

      dTimeCurrent = dTimeCurrent + dTimeSub
      sigI(1:N_STRESS_VECTOR) = sig(1:N_STRESS_VECTOR)
      stVarSub0(1:NSTATS) = stVar(1:NSTATS)
    enddo
  enddo

  if (isUndr == 0) then
    if (props(INDEX_PROP_UNDRAINED_NU) > props(INDEX_PROP_POISSON_UR)) then
      stVar(STVAR_EXCESS_PORE_PRESSURE) = stVar0(STVAR_EXCESS_PORE_PRESSURE) + dSwp
      sig(1:3) = sig(1:3) + stVar(STVAR_EXCESS_PORE_PRESSURE)
    endif
  endif

  if (IsDebug) then
    call calPrincipalStresses(SigI0, prSig0)
    call calAllYieldFs(props, -prSig, fYieldEnd)
    call WriteInDebugFile('Mohr-Coulomb Model converged')
    call WriteInDebugFile('time       :' // trim(string(time0 + dTimeI)))
    call WriteInDebugFile('subroutine :' // trim(subroutineCalled))
    call WriteInDebugFile('GaussPoint :' // trim(string(intPt)))
    call WriteInDebugFile('iEl        :' // trim(string(iEl)))
    call WriteInDebugFile('stVar0     :' // trim(string(stVar0,  1, NSTATS)))
    call WriteInDebugFile('stVar      :' // trim(string(stVar0,  1, NSTATS)))
    call WriteInDebugFile('dEps       :' // trim(string(dEps,    1, N_STRESS_VECTOR)))
    call WriteInDebugFile('SigI0      :' // trim(string(SigI0,   1, N_STRESS_VECTOR)))
    call WriteInDebugFile('prSig0     :' // trim(string(-prSig0, 1, N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('prSig      :' // trim(string(-prSig,  1, N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('p          :' // trim(string(calP(-prSig))))
    call WriteInDebugFile('q          :' // trim(string(calQ(-prSig))))
    call WriteInDebugFile('lodeAngle  :' // trim(string(-getLodeAngle(-prSig)*DEGREES)))
    call WriteInDebugFile('fYield     :' // trim(string(fYield,   1, N_YIELD_FUNCTIONS)))
    call WriteInDebugFile('fYieldEnd  :' // trim(string(fYieldEnd,1, N_YIELD_FUNCTIONS)))
    call WriteInDebugFile('props      :' // trim(string(propsI,   1, NPARAMS)))
    call WriteInDebugFile('isConverged:' // trim(string(isConverged)))
  endif

case (IDTASK_DEP, IDTASK_DE)
  call calStiffnessMatrix(isUndr, BulkW, props, D)

case (IDTASK_NSTAT)
  ! Number of state parameters
  nStat = NSTATS

case (IDTASK_ATTRIBUTES)
  ! matrix type
  NonSym   = 0  ! 1 for non-symmetric D-matrix
  iStrsDep = 0  ! 1 for stress dependent D-matrix
  iTang    = 0  ! 1 for tangent D-matrix
  iTimeDep = 0  ! 1 for time dependent D-matrix
end select

end subroutine ElasticViscoPlastic

!--------------------------------------------------------------------------------------------------
logical function ConsiderInitialSubStepping(props) result(res)
implicit none
double precision, intent(in):: props(NPARAMS)
integer :: IDYieldFunc

IDYieldFunc = getYieldSurfaceFunction(props)
res = CONSIDER_INITIAL_SUB_STEPPING       .or. &
      IDYieldFunc == MATSUOKA_NAKAI_YIELD .or. &
      IDYieldFunc == MATSUOKA_NAKAI_CONVEX_YIELD

end function ConsiderInitialSubStepping

!--------------------------------------------------------------------------------------------------
logical function ConsiderSubStepping(props) result(res)
implicit none
double precision, intent(in):: props(NPARAMS)
integer :: IDYieldFunc

IDYieldFunc = getYieldSurfaceFunction(props)
res = CONSIDER_SUB_STEPPING               .or. &
      IDYieldFunc == MATSUOKA_NAKAI_YIELD .or. &
      IDYieldFunc == MATSUOKA_NAKAI_CONVEX_YIELD

end function ConsiderSubStepping

!--------------------------------------------------------------------------------------------------
subroutine calStiffnessMatrix(isUndr, BulkW, props, D)
implicit none
integer, intent(in):: isUndr
double precision, intent(inout):: BulkW
double precision, dimension(NPARAMS), intent(in):: props
double precision, dimension(N_STRESS_VECTOR, N_STRESS_VECTOR), intent(inout):: D

double precision:: GCurrent, undrainedPoissonRatio

GCurrent = calGRef(props)

call DMatrixDrained(GCurrent, props(INDEX_PROP_POISSON_UR), D)

if (isUndr == 1) then
  if (props(INDEX_PROP_UNDRAINED_NU) < TINY) then
    undrainedPoissonRatio = UNDRAINED_NU
  else
    undrainedPoissonRatio = props(INDEX_PROP_UNDRAINED_NU)
  endif
  ! bulk modulus of water
  BulkW = CalculateBulkModulusWater(GCurrent, props(INDEX_PROP_POISSON_UR), undrainedPoissonRatio)
else
  BulkW = 0
  if (props(INDEX_PROP_UNDRAINED_NU) > props(INDEX_PROP_POISSON_UR)) then
    ! add bulk modulus of water to D-Matrix
    call AddBulkWToDMatrixDrained(CalculateBulkModulusWater(GCurrent, &
                                                            props(INDEX_PROP_POISSON_UR), &
                                                            props(INDEX_PROP_UNDRAINED_NU)), &
                                                            D)
  endif
endif

end subroutine calStiffnessMatrix

!--------------------------------------------------------------------------------------------------
subroutine getPropsViscid(iOverSigFunc,         &
                          viscosity, rateVisco, &
                          perzynaF0,            &
                          props)
implicit none
integer, intent(out):: iOverSigFunc
double precision, intent (out):: viscosity, rateVisco, perzynaF0
double precision, dimension(NPARAMS), intent(in):: props

! viscosity parameters
viscosity = ETA_TO_MU_DP * props(INDEX_PROP_VISCOSITY)
rateVisco = props(INDEX_PROP_RATE_SENSITIVITY)
if (rateVisco > TINY) then
  iOverSigFunc = PERZYNA_MODEL
else
  iOverSigFunc = BINGHAM_MODEL
endif

perzynaF0 = props(INDEX_PROP_PERZYNA_F0)

end subroutine getPropsViscid

!--------------------------------------------------------------------------------------------------
subroutine getdGdS(IDYieldFunc,  &
                   prSig, props, &
                   dGdS,         &
                   indexActivePartMC)
implicit none
integer, intent(in):: IDYieldFunc
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(NPARAMS),                   intent(in) :: props
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dGdS
integer, intent(in):: indexActivePartMC

call cal_dFdS(IDYieldFunc,                        &
             props(INDEX_PROP_SIN_DILATION_PEAK), &
             prSig, props,                        &
             dGdS,                                &
             indexActivePartMC)

end subroutine getdGdS

!--------------------------------------------------------------------------------------------------
subroutine calStress(iStep, iEl, intPt, &
                     dTime, GCurrent,   &
                     prSigElas, prSig0, &
                     props,             &
                     stVar,             &
                     prSig,             &
                     fYield,            &
                     isConverged,       &
                     subroutineCalled)
implicit none

integer, intent(in):: iStep, iEl, intPt
double precision, intent(in):: dTime, GCurrent
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in)   :: prSigElas, prSig0
double precision, dimension(NPARAMS),                   intent(in)   :: props
double precision, dimension(NSTATS),                    intent(inout):: stVar
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out)  :: prSig
double precision, dimension(N_YIELD_FUNCTIONS),         intent(out)  :: fYield
logical,                                                intent(out)  :: isConverged
character(len=1023),                                    intent(out)  :: subroutineCalled

! local variables:
! debugging
integer, parameter:: INT_DEBUG_LOCAL = -1
! yield surfaces
integer :: IDYieldFunc
double precision prDE(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR)
logical :: isTensile = .false.
logical :: IsDebug = .false.
integer :: indexActivePartMC

! initialisation
IsDebug = setIsDebug(iEl, intPt, INT_DEBUG_LOCAL)

isConverged = .false.
isTensile   = .false.

subroutineCalled =''
prSig = prSigElas
indexActivePartMC = FindActivePart_MC(prSig)

! check different parts to select proper strategy for return algorithm
! 1. Tension part
fYield(INDEX_TENSION_YIELD) = calF(TENSION_YIELD, prSig, props, indexActivePartMC)

! 2. Primary part
! yield function type
IDYieldFunc     = getYieldSurfaceFunction(props)
fYield(INDEX_MOHR_COULOMB_YIELD) = calF(IDYieldFunc, prSig, props, indexActivePartMC)

! 3. Complementary part
fYield(INDEX_APEX_YIELD) = calFComplementary(props(INDEX_PROP_SIN_DILATION_PEAK), &
                                             prSig,                               &
                                             props(INDEX_PROP_P_APEX))

if (IsDebug) then
  call WriteInDebugFile('')
  call WriteInDebugFile('GaussPoint      : ' // trim(string(intPt)))
  call WriteInDebugFile('pApex           : ' // trim(string(props(INDEX_PROP_P_APEX))))
  call WriteInDebugFile('sinPhiPeak      : ' // trim(string(props(INDEX_PROP_SIN_FRICTION_PEAK))))
  call WriteInDebugFile('prSig0          : ' // trim(string(prSig0,1,N_PRINCIPAL_STRESS_VECTOR)))
  call WriteInDebugFile('prSigElas       : ' // trim(string(prSigElas,1,N_PRINCIPAL_STRESS_VECTOR)))
  call WriteInDebugFile('fComplementary  : ' // trim(string(fYield(INDEX_APEX_YIELD))))
  call WriteInDebugFile('fTension        : ' // trim(string(fYield(INDEX_TENSION_YIELD))))
  call WriteInDebugFile('fPrimary        : ' // trim(string(fYield(INDEX_MOHR_COULOMB_YIELD))))
endif

if (fYield(INDEX_APEX_YIELD)          <= getTolerance(COMPLEMENTARY_YIELD, prSig) .and. &
    fYield(INDEX_TENSION_YIELD)       <= getTolerance(TENSION_YIELD,       prSig) .and. &
    fYield(INDEX_MOHR_COULOMB_YIELD)  <= getTolerance(IDYieldFunc,         prSig) ) then
  isConverged = .true.
  RETURN
endif

! stiffness prDE matrices
call getPrincipalDMatrix(GCurrent, props(INDEX_PROP_POISSON_UR), prDE)

if (fYield(INDEX_MOHR_COULOMB_YIELD) > getTolerance(IDYieldFunc, prSig)) then
  ! cone
  if (IsDebug) then
    call WriteInDebugFile('before ConeYieldSurface 692')
  endif
  call setSubroutineCalled(IsDebug, subroutineCalled, PRIMARY)
  call ConeYieldSurface(iEl, intPt,                &
                        dTime,                   &
                        prDE,                    &
                        prSig0, prSigElas, &
                        props,                   &
                        stVar,                   &
                        prSig,                &
                        isConverged, isTensile)
  if (IsDebug) then
    call WriteInDebugFile('after ConeYieldSurface isConverged:' // trim(String(isConverged)))
    call WriteInDebugFile('after ConeYieldSurface isTensile  :' // trim(String(isTensile)))
  endif
endif

if (.not.isConverged .or. isTensile) then
  call calStressInTensionArea(iEl, intPt, &
                              dTime,      &
                              prSigElas,  &
                              props,      &
                              stVar,      &
                              prSig,      &
                              prDE,       &
                              isConverged,&
                              IsDebug,    &
                              subroutineCalled)
endif

end subroutine calStress

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine calStressInTensionArea(iEl, intPt, &
                                  dTime,      &
                                  prSigElas,  &
                                  props,      &
                                  stVar,      &
                                  prSig,      &
                                  prDE,       &
                                  isConverged,&
                                  IsDebug,    &
                                  subroutineCalled)
implicit none

integer,                                                                           intent(in)   :: iEl, intPt
double precision,                                                                  intent(in)   :: dTime
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR),                            intent(in)   :: prSigElas
double precision, dimension(NPARAMS),                                              intent(in)   :: props
double precision, dimension(NSTATS),                                               intent(inout):: stVar
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR),                            intent(out)  :: prSig
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(in)   :: prDE
logical,                                                                           intent(out)  :: isConverged
logical,                                                                           intent(in)   :: IsDebug
character(len=1023),                                                               intent(inout):: subroutineCalled

! local variables
integer :: IDYieldFunc
double precision :: frPrimary

isConverged = .false.

! Correct the stress based on the tension criterion
if (IsDebug) call WriteInDebugFile('before TensileYieldSurface 1773')
call setSubroutineCalled(IsDebug, subroutineCalled, TENSILE)
call TensileYieldSurface(iEl, intPt,       &
                         dTime,            &
                         prDE,             &
                         prSigElas, props, &
                         stVar,            &
                         prSig,            &
                         isConverged)
if (IsDebug) then
  call WriteInDebugFile('after TensileYieldSurface isConverged :' // trim(String(isConverged)))
  call WriteInDebugFile('prSig: ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
endif

!if (isConverged) then
  IDYieldFunc = getYieldSurfaceFunction(props)
  frPrimary = calF(IDYieldFunc, prSig, props, FindActivePart_MC(prSig))
  if (frPrimary > getTolerance(IDYieldFunc, prSig)) then
    if (IsDebug) then
      call WriteInDebugFile('after TensileYieldSurface frPrimary:' // trim(String(frPrimary)))
      call WriteInDebugFile('before : ConeWithTensionYieldSurfaces 900')
    endif
    call setSubroutineCalled(IsDebug, subroutineCalled, PRIMARY_TENSILE)
    call ConeWithTensionYieldSurfaces(iEl, intPt, &
                                      dTime,      &
                                      prDE,       &
                                      prSigElas,  &
                                      props,      &
                                      stVar,      &
                                      prSig,      &
                                      isConverged)
    if (IsDebug) then
      call WriteInDebugFile('after ConeWithTensionYieldSurfaces isConverged:' // trim(String(isConverged)))
      call WriteInDebugFile('prSig: ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
    endif
  endif
!endif

end subroutine calStressInTensionArea

!--------------------------------------------------------------------------------------------------
subroutine TensileYieldSurface(iEl, intPt,       &
                               dTime,            &
                               prDE,             &
                               prSigElas, props, &
                               stVar,            &
                               prSig,            &
                               isConverged)
implicit none

integer, intent(in):: iEl, intPt
double precision, intent(in):: dTime
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(in):: prDE
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSigElas
double precision, dimension(NPARAMS),   intent(in):: props
double precision, dimension(NSTATS), intent(inout):: stVar
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: prSig
logical, intent(out):: isConverged

! local variables
integer, parameter :: INT_DEBUG_LOCAL = -1
integer :: iOverSigFunc, IDYieldFunc, iter
double precision:: f, dLamda, dFdLamda, residual, alpha
double precision:: viscosity, rateVisco, perzynaF0
double precision:: phiF, OverStressDFactor
double precision:: sigTen
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: dGdS, DdGdS, dFdS
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: prSigElasMod

integer :: indexActivePartMC
logical :: isStressModified, isStressEverModified
logical :: IsDebug = .false.

IsDebug = setIsDebug(iEl, intPt, INT_DEBUG_LOCAL)

prSig = prSigElas
prSigElasMod = prSig
indexActivePartMC = FindActivePart_MC(prSig)

call getPropsViscid(iOverSigFunc, viscosity, rateVisco, perzynaF0, props)

! initialisation
isConverged = .false.

sigTen = props(INDEX_PROP_SIG_T_CUT_OFF)

IDYieldFunc = TENSION_YIELD
f = calF(IDYieldFunc, prSig, props, indexActivePartMC)

if (IsDebug) then
  call WriteInDebugFile('fTension    : ' // trim(string(f)))
  call WriteInDebugFile('prSigElas   : ' // trim(string(prSigElas,1,N_PRINCIPAL_STRESS_VECTOR)))
endif

if (f < getTolerance(IDYieldFunc, prSig)) then
  isConverged = .true.
  stVar(STVAR_LOADING_TYPE)= ELASTIC_LOADING
  RETURN
endif

! no viscosity for tension cut-off
viscosity = 0.0

! initialisation
iter = 0
dLamda = 0.0d0
isStressEverModified = .false.
isStressModified     = .false.

phiF = calOverStressPhi(iOverSigFunc, f, perzynaF0, rateVisco)
residual = calResidual(phiF, dLamda, dTime, viscosity)

do while (iter < ITER_MAX .and. abs(residual) > getTolerance(IDYieldFunc, prSig))
  ! loop until convergence
  if (iter > 0 .and. .not.isStressEverModified) then
    call CheckAndModifyPrincipalStressAtCornersMC(prSigElas, prSig, isStressModified, IsDebug)
    if (isStressModified) then
      dLamda = 0.0d0
      prSigElasMod = prSig
      isStressEverModified = .true.
    endif
  endif

  indexActivePartMC = FindActivePart_MC(prSig)

  ! dGdSigma
  call getdGdS(IDYieldFunc,  &
               prSig, props, &
               dGdS,         &
               indexActivePartMC)

  ! prDE*dGdSigma
  call MatrixToVectorProd(prDE, dGdS, N_PRINCIPAL_STRESS_VECTOR, DdGdS)

  ! derivative of the yield surface
  ! dFdSigma
  call cal_dFdS(IDYieldFunc,                        &
               props(INDEX_PROP_SIN_DILATION_PEAK), &
               prSig, props,                        &
               dFdS,                                &
               indexActivePartMC)

  ! NOT TESTED:
  OverStressDFactor = calOverStressDFactor(iOverSigFunc, f, perzynaF0, rateVisco)
  dFdS = OverStressDFactor * dFdS

  ! dfdSigma*prDE*dGdSigma
  dFdLamda = dot_product(DdGdS, dFdS)
  alpha = viscosity/dTime + dFdLamda
  dLamda = dLamda + residual / alpha

  ! prSig = prSigElas - dLamda * DdGdS
  call sumVectors(prSigElasMod, DdGdS, 1.0d0, -dLamda, N_PRINCIPAL_STRESS_VECTOR, prSig)

  if (IsDebug) then
    call WriteInDebugFile('iter       : ' // trim(string(iter)))
    call WriteInDebugFile('residualTen: ' // trim(string(residual)))
    call WriteInDebugFile('dGdS       : ' // trim(string(dGdS,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('dFdS       : ' // trim(string(dFdS,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('DdGdS      : ' // trim(string(DdGdS,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('dFdLamda   : ' // trim(string(dFdLamda)))
    call WriteInDebugFile('alpha      : ' // trim(string(alpha)))
    call WriteInDebugFile('dLamda     : ' // trim(string(dLamda)))
    call WriteInDebugFile('prSig      : ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
  endif

  if (isDiverged(prSig)) then
    isConverged = .false.
    RETURN
  endif

  !check yield function
  f = calF(IDYieldFunc, prSig, props, FindActivePart_MC(prSig))
  
  phiF = calOverStressPhi(iOverSigFunc, f, perzynaF0, rateVisco)
  residual = calResidual(phiF, dLamda, dTime, viscosity)

  if (IsDebug) then
    call WriteInDebugFile('residual : ' // trim(string(residual)))
  endif

  iter = iter + 1

  if (residual > DIVERGENCE_LIMIT) then
    ! diverging
    isConverged = .false.
    RETURN
  endif
enddo

if (abs(residual) < getTolerance(IDYieldFunc, prSig)) then
  isConverged = .true.
  stVar(STVAR_LOADING_TYPE) = TENSILE_LOADING
endif

end subroutine TensileYieldSurface

!--------------------------------------------------------------------------------------------------
subroutine ConeYieldSurface(iEl, intPt,        &
                            dTime,             &
                            prDE,              &
                            prSig0, prSigElas, &
                            props,             &
                            stVar,             &
                            prSig,             &
                            isConverged, isTensile)
implicit none

integer, intent(in):: iEl, intPt
double precision, intent(in):: dTime
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prDE
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig0, prSigElas
double precision, dimension(NPARAMS), intent(in):: props
double precision, dimension(NSTATS),  intent(inout):: stVar
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: prSig
logical, intent(out):: isConverged, isTensile

integer, parameter :: INT_DEBUG_LOCAL = -1

integer :: iOverSigFunc, IDYieldFunc, iter
integer :: IDPotentialFunc

double precision:: f,dLamda,dFdLamda,residual
double precision:: viscosity, rateVisco, perzynaF0
double precision:: phiF, OverStressDFactor

double precision:: fTension
double precision:: alpha

double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: dGdS, DdGdS, dFdS
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: dEpsP
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: prSigElasMod

logical :: isStressModified, isStressEverModified
logical :: IsDebug = .false.
integer :: indexActivePartMC

IsDebug = setIsDebug(iEl, intPt, INT_DEBUG_LOCAL)

!initialize
prSig       = prSigElas
prSigElasMod= prSigElas
indexActivePartMC = FindActivePart_MC(prSig)

dEpsP = 0.0

IDYieldFunc = getYieldSurfaceFunction(props)

call getPropsViscid(iOverSigFunc, viscosity, rateVisco, perzynaF0, props)

isConverged = .false.
isTensile   = .false.

! initialisation
iter   = 0
dLamda = 0.0d0
isStressEverModified = .false.
isStressModified     = .false.

f = calF(IDYieldFunc, prSig, props, indexActivePartMC)
phiF = calOverStressPhi(iOverSigFunc, f, perzynaF0, rateVisco)
residual = calResidual(phiF, dLamda, dTime, viscosity)

IDPotentialFunc = getPotentialSurfaceFunction(IDYieldFunc, prSig)

if (IsDebug) then
  call WriteInDebugFile('residual : ' // trim(string(residual)))
  call WriteInDebugFile('prSig0   : ' // trim(string(prSig0,1,N_PRINCIPAL_STRESS_VECTOR)))
  call WriteInDebugFile('prSigElas: ' // trim(string(prSigElas,1,N_PRINCIPAL_STRESS_VECTOR)))
  call WriteInDebugFile('prSig    : ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
  call WriteInDebugFile('lodeAngle: ' // trim(string(-getLodeAngle(prSigElas)*DEGREES)))
  call WriteInDebugFile('IDPotentialFunc: ' // trim(string(IDPotentialFunc)))
  if (IDYieldFunc == MOHR_COULOMB_YIELD) then
    call WriteInDebugFile('calF_MC_TxComp : ' // trim(string(calF_MC_TxComp(props(INDEX_PROP_SIN_FRICTION_PEAK), props(INDEX_PROP_P_APEX), prSig))))
    call WriteInDebugFile('calF_MC_General: ' // trim(string(calF_MC_General(props(INDEX_PROP_SIN_FRICTION_PEAK),  props(INDEX_PROP_P_APEX), prSig))))
    call WriteInDebugFile('calF_MC_TxExt  : ' // trim(string(calF_MC_TxExt(props(INDEX_PROP_SIN_FRICTION_PEAK), props(INDEX_PROP_P_APEX), prSig))))
  endif
endif

do while (iter < ITER_MAX .and. abs(residual) > getTolerance(IDYieldFunc, prSig))
  ! loop until convergence
  if (iter > 0 .and. .not.isStressEverModified) then
    call CheckAndModifyPrincipalStressAtCornersMC(prSigElas, prSig, isStressModified, IsDebug)
    if (isStressModified) then
      dLamda = 0.0d0
      prSigElasMod = prSig
      isStressEverModified = .true.
    endif
  endif

  indexActivePartMC = FindActivePart_MC(prSig)

  ! derivative of the plastic potential
  ! dGdSigma
  IDPotentialFunc = getPotentialSurfaceFunction(IDYieldFunc, prSig)
  call getdGdS(IDPotentialFunc, prSig, props, dGdS, indexActivePartMC)

  if (IsDebug) then
    call WriteInDebugFile('iter           : ' // trim(string(iter)))
    call WriteInDebugFile('f              : ' // trim(string(f)))
    call WriteInDebugFile('isStressModified: ' // trim(string(isStressModified)))
    call WriteInDebugFile('residual       : ' // trim(string(residual)))
    call WriteInDebugFile('prSig          : ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('prSigElasMod   : ' // trim(string(prSigElasMod,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('lodeAngle      : ' // trim(string(-getLodeAngle(prSig)*DEGREES)))
    call WriteInDebugFile('IDPotentialFunc: ' // trim(string(IDPotentialFunc)))
    call WriteInDebugFile('dGdS           : ' // trim(string(dGdS,1,N_PRINCIPAL_STRESS_VECTOR)))
  endif

  ! derivative of the yield surface
  ! dFdSigma
  call cal_dFdS(IDYieldFunc, props(INDEX_PROP_SIN_FRICTION_PEAK), prSig, props, dFdS, indexActivePartMC)

  OverStressDFactor = calOverStressDFactor(iOverSigFunc, f, perzynaF0, rateVisco)
  dFdS = OverStressDFactor * dFdS

  ! prDE*dGdSigma
  call MatrixToVectorProd(prDE, dGdS, N_PRINCIPAL_STRESS_VECTOR, DdGdS)

  ! - dfdSigma*prDE*dGdSigma
  dFdLamda = dot_product(DdGdS, dFdS)

  alpha = viscosity/dTime + dFdLamda
  dLamda = dLamda + residual / alpha

  ! prSig = prSigElas - dLamda * DdGdS
  call sumVectors(prSigElasMod, DdGdS, 1.0d0, -dLamda, N_PRINCIPAL_STRESS_VECTOR, prSig)

  if (IsDebug) then
    call WriteInDebugFile('DdGdS            : ' // trim(string(DdGdS,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('dFdLamda         : ' // trim(string(dFdLamda)))
    call WriteInDebugFile('alpha            : ' // trim(string(alpha)))
    call WriteInDebugFile('dLamda           : ' // trim(string(dLamda)))
    call WriteInDebugFile('residual / alpha : ' // trim(string(residual / alpha)))
  endif

  if (isDiverged(prSig)) then
    isConverged = .false.
    RETURN
  endif

  !check yield function
  f = calF(IDYieldFunc, prSig, props, FindActivePart_MC(prSig))
  phiF = calOverStressPhi(iOverSigFunc, f, perzynaF0, rateVisco)
  residual = calResidual(phiF, dLamda, dTime, viscosity)

  if (IsDebug) then
    call WriteInDebugFile('prSig   : ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('lodeAngle  : ' // trim(string(-getLodeAngle(prSig)*DEGREES)))
    call WriteInDebugFile('residual   : ' // trim(string(residual)))
    if (IDYieldFunc == MOHR_COULOMB_YIELD) then
      call WriteInDebugFile('calF_MC_TxComp: ' // trim(string(calF_MC_TxComp(props(INDEX_PROP_SIN_FRICTION_PEAK), props(INDEX_PROP_P_APEX), prSig))))
      call WriteInDebugFile('calF_MC_General : ' // trim(string(calF_MC_General(props(INDEX_PROP_SIN_FRICTION_PEAK),  props(INDEX_PROP_P_APEX), prSig))))
      call WriteInDebugFile('calF_MC_TxExt: ' // trim(string(calF_MC_TxExt(props(INDEX_PROP_SIN_FRICTION_PEAK), props(INDEX_PROP_P_APEX), prSig))))
    endif
  endif

  if (abs(residual) > DIVERGENCE_LIMIT) EXIT

  iter = iter + 1
enddo

if (IsDebug) then
  call WriteInDebugFile('iter     : ' // trim(string(iter)))
  call WriteInDebugFile('residual : ' // trim(string(residual)))
endif

if (abs(residual) < getTolerance(IDYieldFunc, prSig)) isConverged = .true.

if (isDiverged(prSig)) then
  isConverged = .false.
  RETURN
endif


if (isConverged) then
  fTension = calF(TENSION_YIELD, prSig, props, FindActivePart_MC(prSig))

  if (IsDebug) then
    call WriteInDebugFile('converged prSig end ConeYieldSurface: ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('converged fTension end ConeYieldSurface: ' // trim(string(fTension)))
  endif

  if (fTension > getTolerance(TENSION_YIELD, prSig)) then
    ! In tension zone, therefore it is not isConverged yet
    isConverged = .false.
    isTensile = .true.
  else
    stVar(STVAR_LOADING_TYPE) = FAILURE_LOADING
    isConverged = .true.
    isTensile = .false.
  endif
endif

end subroutine ConeYieldSurface

!--------------------------------------------------------------------------------------------------
subroutine ConeWithTensionYieldSurfaces(iEl, intPt, &
                                        dTime,      &
                                        prDE,       &
                                        prSigElas,  &
                                        props,      &
                                        stVar,      &
                                        prSig,      &
                                        isConverged)
implicit none

integer, intent(in):: iEl, intPt
double precision, intent(in):: dTime
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prDE
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSigElas
double precision, dimension(NPARAMS), intent(in):: props
double precision, dimension(NSTATS), intent(inout):: stVar
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: prSig
logical, intent(out):: isConverged

integer, parameter :: INT_DEBUG_LOCAL = -2018
integer, parameter :: FRICTION = 1
integer, parameter :: TENSION  = 2

integer :: iOverSigFunc, iter, i

double precision:: viscosity, rateVisco, perzynaF0
double precision:: OverStressDFactor

integer, dimension(N_CROSS_YIELD_SURFACES) :: IDYieldFunc, IDPotentialFunc

double precision, dimension(N_CROSS_YIELD_SURFACES) :: dFdh
double precision, dimension(N_CROSS_YIELD_SURFACES) :: dLamda, ddLamda
double precision, dimension(N_CROSS_YIELD_SURFACES) :: f, residual, phiF
double precision, dimension(N_CROSS_YIELD_SURFACES, N_CROSS_YIELD_SURFACES) :: dFdLamda, dFdLamdaInverse, Coupling

double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_CROSS_YIELD_SURFACES) :: dGdS, DdGdS, dFdS, dhdEpsP
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR, N_CROSS_YIELD_SURFACES) :: dEpsP
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR) :: prSigElasMod

integer :: indexActivePartMC
logical :: isStressModified, isStressEverModified
logical :: IsDebug = .false.

IsDebug = setIsDebug(iEl, intPt, INT_DEBUG_LOCAL)

!initialize

dEpsP    = 0.0
prSig = prSigElas
prSigElasMod = prSig
indexActivePartMC = FindActivePart_MC(prSig)

if (USE_COUPLED_HARDENING_CROSS_POINTS) then
  Coupling = 1.0d0
else
  Coupling = 0.0d0
  Coupling(1,1) = 1.0d0
  Coupling(2,2) = 1.0d0
endif

! cone
call getPropsViscid(iOverSigFunc, viscosity, rateVisco, perzynaF0, props)

IDYieldFunc(FRICTION)     = getYieldSurfaceFunction(props)
IDPotentialFunc(FRICTION) = getPotentialSurfaceFunction(IDYieldFunc(FRICTION), prSig)
IDYieldFunc(TENSION)      = TENSION_YIELD
IDPotentialFunc(TENSION)  = TENSION_YIELD

isConverged = .false.

! initialisation
iter          = 0
dLamda        = 0.0d0
isStressEverModified = .false.
isStressModified     = .false.
dFdh    = 0.0d0
dhdEpsP = 0.0d0

do i=1,N_CROSS_YIELD_SURFACES
  f(i) = calF(IDYieldFunc(i), prSig, props, indexActivePartMC)
  phiF(i) = calOverStressPhi(iOverSigFunc, f(i), perzynaF0, rateVisco)
  residual(i)= calResidual(phiF(i), dLamda(i), dTime, viscosity)
enddo

if (IsDebug) then
  call WriteInDebugFile('f(i)    : ' // trim(string(f,1,2)))
  call WriteInDebugFile('prSig   : ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
endif

do while (iter < ITER_MAX .and. &
          (abs(residual(FRICTION)) > getTolerance(IDYieldFunc(FRICTION), prSig) .or. &
           abs(residual(TENSION) ) > getTolerance(IDYieldFunc(TENSION),  prSig)      ))
  ! loop until convergence
  if (iter > 0 .and. .not.isStressEverModified) then
    call CheckAndModifyPrincipalStressAtCornersMC(prSigElas, prSig, isStressModified, IsDebug)
    if (isStressModified) then
      dLamda = 0.0d0
      prSigElasMod = prSig
      isStressEverModified = .false.
    endif
  endif

  indexActivePartMC = FindActivePart_MC(prSig)

  IDPotentialFunc(FRICTION) = getPotentialSurfaceFunction(IDYieldFunc(FRICTION), prSig)

  ! derivative of the plastic potential
  ! dGdSigma
  do i=1,N_CROSS_YIELD_SURFACES
    call getdGdS(IDPotentialFunc(i), prSig, props, dGdS(:, i), indexActivePartMC)
  enddo

  if (IsDebug) then
    call WriteInDebugFile('iter          : ' // trim(string(iter)))
    call WriteInDebugFile('dGdS(Friction): ' // trim(string(dGdS(:,1),1,N_PRINCIPAL_STRESS_VECTOR)))
    call WriteInDebugFile('dGdS(Tension) : ' // trim(string(dGdS(:,2),1,N_PRINCIPAL_STRESS_VECTOR)))
  endif

  ! derivative of the yield surface
  ! dFdSigma
  do i=1,N_CROSS_YIELD_SURFACES
    call cal_dFdS(IDYieldFunc(i), props(INDEX_PROP_SIN_FRICTION_PEAK), prSig, props, dFdS(:,i), indexActivePartMC)
    OverStressDFactor = calOverStressDFactor(iOverSigFunc, f(i), perzynaF0, rateVisco)
    dFdS(:,i) = OverStressDFactor * dFdS(:,i)
  enddo

  ! prDE*dGdSigma
  do i=1,N_CROSS_YIELD_SURFACES
    call MatrixToVectorProd(prDE, dGdS(:,i), N_PRINCIPAL_STRESS_VECTOR, DdGdS(:,i))
  enddo

  ! - dfdSigma*prDE*dGdSigma + dfdh*dhdEpsP*dGdSigma
  call getdFdLamdaCrossPoint(dFdh, Coupling, dGdS, DdGdS, dFdS, dhdEpsP, dFdLamda)
  dFdLamda = dFdLamda + viscosity/dTime

  call calInverseMatrix(dFdLamda, dFdLamdaInverse, 2)

  ! AInverse*residual
  call MatrixToVectorProd(dFdLamdaInverse, residual, 2, ddLamda)
  dLamda = dLamda + ddLamda

  do i=1, N_CROSS_YIELD_SURFACES
    dEpsP(:,i) = dLamda(i) * dGdS(:,i)
  enddo

  ! prSig = prSigElas - dLamdaFriction * DdGdS - dLamdaTension * DdGtdS
  call sumVectors(prSigElasMod, DdGdS(:,FRICTION), 1.0d0, -dLamda(FRICTION), N_PRINCIPAL_STRESS_VECTOR, prSig)
  call sumVectors(prSig,        DdGdS(:,TENSION),  1.0d0, -dLamda(TENSION),  N_PRINCIPAL_STRESS_VECTOR, prSig)

  if (isDiverged(prSig)) then
    isConverged = .false.
    RETURN
  endif

  if (IsDebug) then
    call WriteInDebugFile('prSig       : ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
  endif

  !check yield function

  do i=1,N_CROSS_YIELD_SURFACES
    f(i) = calF(IDYieldFunc(i), prSig, props, FindActivePart_MC(prSig))
    phiF(i)  = calOverStressPhi(iOverSigFunc, f(i), perzynaF0, rateVisco)
    residual(i) = calResidual(phiF(i), dLamda(i), dTime, viscosity)

    !if (IsDebug) call WriteInDebugFile('f(i)        : ' // trim(string(f(i))))
    if (IsDebug) then
      call WriteInDebugFile('residual(i) : ' // trim(string(residual(i))))
    endif

    if (abs(residual(i)) > DIVERGENCE_LIMIT) then
      isConverged = .false.
      RETURN
    endif
  enddo

  iter = iter + 1
enddo

if (IsDebug) then
  call WriteInDebugFile('iter     : ' // trim(string(iter)))
  call WriteInDebugFile('residual : ' // trim(string(residual,1,2)))
endif

if (abs(residual(FRICTION)) < getTolerance(IDYieldFunc(FRICTION), prSig) .and. &
    abs(residual(TENSION))  < getTolerance(IDYieldFunc(TENSION),  prSig)      ) isConverged = .true.

if (isDiverged(prSig)) then
  isConverged = .false.
  RETURN
endif

if (isConverged) then
  stVar(STVAR_LOADING_TYPE) = PRIMARY_TENSILE_LOADING
  isConverged = .true.
endif

end subroutine ConeWithTensionYieldSurfaces

!---------------------------------------------------------------------------------------------------------
subroutine initialiseStateVariables(stVar)
implicit none
double precision, intent(out):: stVar(NSTATS)

stVar(1 : NSTATS) = 0.0

stVar(STVAR_INITIALIZATION) = 1.0

end subroutine initialiseStateVariables

!--------------------------------------------------------------------------------------------------
double precision function calGRef(props) result(GRef)
implicit none
double precision, intent(in) :: props(NPARAMS)

double precision :: EurRef, nu

nu     = props(INDEX_PROP_POISSON_UR)
EurRef = props(INDEX_PROP_E)

GRef = 0.5d0*EurRef / (1.0d0 + nu)

end function calGRef

!--------------------------------------------------------------------------------------------------
double precision function calF(IDYieldFunc, prSig, props, indexActivePartMC) result(f)
implicit none
integer, intent(in):: IDYieldFunc
double precision, intent(in):: prSig(N_PRINCIPAL_STRESS_VECTOR)
double precision, intent(in):: props(NPARAMS)
integer, intent(in) :: indexActivePartMC

! local variables
integer:: IDYieldConeFunc
double precision:: pApex, sinPhi

pApex = props(INDEX_PROP_P_APEX)
sinPhi = props(INDEX_PROP_SIN_FRICTION_PEAK)

select case(IDYieldFunc)
case(COMPLEMENTARY_YIELD)
  f = calFComplementary(sinPhi, prSig, pApex)

case(DRUCKER_PRAGER_YIELD)
  f = calF_DP(sinPhi, pApex, prSig)

case(DRUCKER_PRAGER_YIELD_EXTENSION)
  f = calF_DPExtension(sinPhi, pApex, prSig)

case(MATSUOKA_NAKAI_YIELD)
  f = calF_MN(sinPhi, pApex, prSig)

case(MOHR_COULOMB_YIELD)
  f = calF_MC(sinPhi, pApex, prSig, indexActivePartMC)

case(MATSUOKA_NAKAI_CONVEX_YIELD)
  f = calF_MNC(sinPhi, 0.0d0, pApex, prSig)

case(TENSION_YIELD)
  IDYieldConeFunc = getYieldSurfaceFunction(props)
  if (IDYieldConeFunc == MATSUOKA_NAKAI_YIELD        .or. &
      IDYieldConeFunc == MATSUOKA_NAKAI_CONVEX_YIELD .or. &
      IDYieldConeFunc == DRUCKER_PRAGER_YIELD) then
    f = calF_Ten_P(props(INDEX_PROP_SIG_T_CUT_OFF), prSig)
  else
    f = calF_Ten(props(INDEX_PROP_SIG_T_CUT_OFF), prSig, indexActivePartMC)
  endif

case default
  call SetError('unknown yield function: '  // trim(string(IDYieldFunc)))
end select

end function calF

!--------------------------------------------------------------------------------------------------
subroutine cal_dFdS(IDYieldFunc, &
                    sinPhi,       &
                    prSig, props, &
                    dFdS,         &
                    indexActivePartMC)
implicit none
integer, intent(in):: IDYieldFunc
double precision, intent(in)::  sinPhi
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in) :: prSig
double precision, dimension(NPARAMS),                   intent(in) :: props
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(out):: dFdS
integer, intent(in):: indexActivePartMC

! local variables:
integer:: IDYieldConeFunc

select case(IDYieldFunc)
case(DRUCKER_PRAGER_YIELD)
  call caldFdS_DP(sinPhi, prSig, dFdS)

case(DRUCKER_PRAGER_YIELD_EXTENSION)
  call caldFdS_DPExtension(sinPhi, prSig, dFdS)

case(MATSUOKA_NAKAI_YIELD)
  call caldFdS_MN(sinPhi, props(INDEX_PROP_P_APEX), prSig, dFdS)

case(MOHR_COULOMB_YIELD)
  call caldFdS_MC(sinPhi, dFdS, indexActivePartMC)

case(TENSION_YIELD)
  IDYieldConeFunc = getYieldSurfaceFunction(props)
  if (IDYieldConeFunc == MATSUOKA_NAKAI_YIELD        .or. &
      IDYieldConeFunc == MATSUOKA_NAKAI_CONVEX_YIELD .or. &
      IDYieldConeFunc == DRUCKER_PRAGER_YIELD) then
    call caldFdS_Ten_P(dFdS)
  else
    call caldFdS_Ten(dFdS, indexActivePartMC)
  endif

case(MATSUOKA_NAKAI_CONVEX_YIELD)
  call caldFdS_MNC(sinPhi, 0.0d0, prSig, dFdS )

case default
  call SetError('unknown yield function in cal_dFdS: '  // trim(string(IDYieldFunc)))
end select

end subroutine cal_dFdS

!---------------------------------------------------------------------------------------------------------
subroutine checkYieldSurfaces(props, prSig, IsDebug, isConverged)
implicit none
double precision, intent(in):: props(NPARAMS)
double precision, intent(in):: prSig(N_PRINCIPAL_STRESS_VECTOR)
logical, intent(in) :: IsDebug
logical, intent(out) :: isConverged

integer :: IDYieldFunc, indexActivePartMC
double precision:: fPrimary, fTension

isConverged = .true.

IDYieldFunc = getYieldSurfaceFunction(props)
indexActivePartMC = FindActivePart_MC(prSig)

! check different parts to select proper strategy for return algorithm
! 1. Tension part
fTension = calF(TENSION_YIELD, prSig, props, indexActivePartMC)

! 2. Primary part
fPrimary = calF(IDYieldFunc, prSig, props, indexActivePartMC)

if (IsDebug) then
  call WriteInDebugFile('checkYieldSurfaces...')
  call WriteInDebugFile('prSig       : ' // trim(string(prSig,1,N_PRINCIPAL_STRESS_VECTOR)))
  call WriteInDebugFile('IDYieldFunc : ' // trim(string(IDYieldFunc)))
  call WriteInDebugFile('p           : ' // trim(string(calP(prSig))))
  call WriteInDebugFile('q           : ' // trim(string(calQ(prSig))))
  call WriteInDebugFile('fTension0   : ' // trim(string(fTension)))
  call WriteInDebugFile('fPrimary0   : ' // trim(string(fPrimary)))
  call WriteInDebugFile('isConverged : ' // trim(string(isConverged)))
endif

if (fPrimary  > getTolerance(IDYieldFunc, prSig)) then
  if (IsDebug) call WriteInDebugFile('high fPrimary')
  isConverged = .false.
endif

if (fTension > getTolerance(TENSION_YIELD, prSig)) then
  if (IsDebug) call WriteInDebugFile('high fTension')
  isConverged = .false.
endif

if (isDiverged(prSig)) then
  isConverged = .false.
  if (IsDebug) call WriteInDebugFile('isDiverged')
endif

end subroutine checkYieldSurfaces

!--------------------------------------------------------------------------------------------------
subroutine CheckAndCopyProps(propsI, props)
implicit none
double precision, intent(in):: propsI(NPARAMS)
double precision, intent(out):: props(NPARAMS)

props = propsI

props(INDEX_PROP_P_APEX)        = abs(props(INDEX_PROP_P_APEX))
props(INDEX_PROP_COHESION)      = max(props(INDEX_PROP_COHESION), 0.0d0)
props(INDEX_PROP_FRICTION_PEAK) = max(props(INDEX_PROP_FRICTION_PEAK), SMALL)

! last operations because same indeces are used
props(INDEX_PROP_P_APEX) = props(INDEX_PROP_COHESION) * (1.0d0/tan(props(INDEX_PROP_FRICTION_PEAK) * RADIANS))
props(INDEX_PROP_SIG_T_CUT_OFF)  = max(-abs(props(INDEX_PROP_SIG_T_CUT_OFF)), -props(INDEX_PROP_P_APEX))

props(INDEX_PROP_SIN_FRICTION_PEAK) = sin(props(INDEX_PROP_FRICTION_PEAK) * RADIANS)
props(INDEX_PROP_SIN_DILATION_PEAK) = sin(props(INDEX_PROP_DILATION_PEAK) * RADIANS)

end subroutine CheckAndCopyProps

!--------------------------------------------------------------------------------------------------
integer function GetParamCountElasticViscoPlastic() result(res)
implicit none

res = NPARAMS

end function GetParamCountElasticViscoPlastic

!--------------------------------------------------------------------------------------------------
integer function GetStateVarCountElasticViscoPlastic() result(res)
implicit none

res = NSTATS

end function GetStateVarCountElasticViscoPlastic

!--------------------------------------------------------------------------------------------------
subroutine calAllYieldFs(props, prSig, fYield)
implicit none
double precision, intent(in):: props(NPARAMS)
double precision, dimension(N_PRINCIPAL_STRESS_VECTOR), intent(in):: prSig
double precision, dimension(N_YIELD_FUNCTIONS), intent(out):: fYield

integer :: IDYieldFunc, indexActivePartMC

! yield function type
IDYieldFunc = getYieldSurfaceFunction(props)

indexActivePartMC = FindActivePart_MC(prSig)

! 1. Tension part
fYield(INDEX_TENSION_YIELD) = calF(TENSION_YIELD, prSig, props, indexActivePartMC)

! 2. Primary part
fYield(INDEX_MOHR_COULOMB_YIELD) = calF(IDYieldFunc, prSig, props, indexActivePartMC)

! 3. Complementary part
fYield(INDEX_APEX_YIELD) = calFComplementary(props(INDEX_PROP_SIN_DILATION_PEAK), &
                                             prSig,                               &
                                             props(INDEX_PROP_P_APEX))
 
end subroutine calAllYieldFs

!---------------------------------------------------------------------------------------------------------
integer function getYieldSurfaceFunction(props) result(IDYieldFuncion)
implicit none
double precision, dimension(NPARAMS), intent(in):: props

IDYieldFuncion = max(1, nint(props(INDEX_PROP_YIELD_FUNCTION)))

end function getYieldSurfaceFunction

!---------------------------------------------------------------------------------------------------------
integer function getPotentialSurfaceFunction(IDYieldFunc, prSig) result(IDPotentialFunc)
implicit none
integer, intent(in):: IDYieldFunc
double precision, intent(in):: prSig(N_PRINCIPAL_STRESS_VECTOR)

double precision :: lodeAngle

if (USE_DEFAULT_PLASTIC_POTENTIAL) then
  IDPotentialFunc = DEFAULT_PLASTIC_POTENTIAL
else
  if (IDYieldFunc == DRUCKER_PRAGER_YIELD .or. IDYieldFunc == MOHR_COULOMB_YIELD) then
    IDPotentialFunc = IDYieldFunc
  else
    lodeAngle = -getLodeAngle(prSig)
    if (lodeAngle > LODE_ANGLE_TRIAXIAL_LIMIT) then
      IDPotentialFunc = DEFAULT_PLASTIC_POTENTIAL
    elseif (lodeAngle < -LODE_ANGLE_TRIAXIAL_LIMIT) then
      IDPotentialFunc = DEFAULT_PLASTIC_POTENTIAL_EXTENSION
    else
      select case(IDYieldFunc)
      case(MATSUOKA_NAKAI_CONVEX_YIELD, MATSUOKA_NAKAI_YIELD)
        IDPotentialFunc = MOHR_COULOMB_YIELD
      case default
        IDPotentialFunc = IDYieldFunc
      end select
    endif
  endif
endif

end function getPotentialSurfaceFunction

!--------------------------------------------------------------------------------------------------
subroutine GetParamAndUnitElasticViscoPlastic(iParam, ParamName, Units)
!
! return the parameters name and units of the different models
!
! Units: use F for force unit
!            L for length unit
!            T for time unit
!
implicit none
integer, intent(in):: iParam
character (Len=255) ParamName, Units

select case (iParam)
  ! stiffness parameters
  case(INDEX_PROP_E)
    ParamName = 'E'                 ; Units = 'F/L^2#'

  case (INDEX_PROP_FRICTION_PEAK)
    ParamName = '@f#_peak#'         ; Units = 'deg'

  case (INDEX_PROP_DILATION_PEAK)
    ParamName = '@y#_peak#'         ; Units = 'deg'

  case (INDEX_PROP_COHESION)
    ParamName = 'c'''               ; Units = 'F/L^2#'

  case (INDEX_PROP_POISSON_UR)
    ParamName = '@n#_ur'            ; Units = '-'

  case (INDEX_PROP_SIG_T_CUT_OFF)
    ParamName = '|@s#_t, cut-off#|' ; Units = 'F/L^2#'

  case (INDEX_PROP_UNDRAINED_NU)
    ParamName = '@n#_un# (UMAT)'    ; Units = '-'

  case (INDEX_PROP_YIELD_FUNCTION)
    ParamName = 'yield function (MC=1; DP=2; MNC=3; MN=4)' ; Units   = '-'

  ! viscosity parameters
  case (INDEX_PROP_VISCOSITY)
    ParamName = '@h#'                   ; Units = 'FT/L^2#'
  case (INDEX_PROP_RATE_SENSITIVITY)
    ParamName = 'N'                     ; Units = '-'
  case (INDEX_PROP_PERZYNA_F0)
    ParamName = 'f_0'                   ; Units = '-'

  case default
    ! parameter not in model
    ParamName = ' N/A '        ; Units     = ' N/A '

end select ! iParam

ParamName = trim(string(iParam)) //'. '// trim(ParamName)

end subroutine GetParamAndUnitElasticViscoPlastic

!--------------------------------------------------------------------------------------------------
subroutine GetModelNameElasticViscoPlastic(ModelName)
use ModGlobalConstants
implicit none
character (Len= * ) ModelName

ModelName= 'ElasticViscoPlastic'

end subroutine GetModelNameElasticViscoPlastic

end module ElasticViscoPlastic_Module

