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
subroutine User_Mod( IDTask, iMod, IsUndr,                    &
                     iStep, iTer, iEl, intPt,                 &
                     X, Y, Z,                                 &
                     Time0, dTime,                            &
                     Props, Sig0, Swp0, StVar0,               &
                     dEps, D, BulkW,                          &
                     Sig, Swp, StVar, ipl,                    &
                     nStat, NonSym, iStrsDep, iTimeDep,iTang, &
                     iPrjDir, iPrjLen, iAbort )
!
! Purpose: User supplied soil model
!
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : return number of state variables
!                       5 : inquire matrix properties
!                           return switch for non-symmetric D-matrix
!                           stress/time dependent matrix
!                       6 : calculate elastic material stiffness matrix
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iMod     I   I    : model number (1..10)
!  IsUndr   I   I    : =1 for undrained, 0 otherwise
!  iStep    I   I    : Global step number
!  iter     I   I    : Global iteration number
!  iel      I   I    : Global element number
!  intPt    I   I    : Global integration point number
!  X        I   R    : X-Position of integration point
!  Y        I   R    : Y-Position of integration point
!  Z        I   R    : Z-Position of integration point
!  Time0    I   R    : Time at start of step
!  dTime    I   R    : Time increment
!  Props    I   R()  : List with model parameters
!  Sig0     I   R()  : Stresses at start of step
!  Swp0     I   R    : Excess pore pressure start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulkmodulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : Non-Symmetric D-matrix ?
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iTang    O   I    : =1 for tangent matrix
!  iAbort   O   I    : =1 to force stopping of calculation
!
use ModMathLibrary
use ModString
use ModDebug
use ModGlobalConstants
use ElasticPerfectlyPlastic_Module
use ElasticViscoPlastic_Module

implicit none

!arguments:
! model number:
integer, intent(in):: iMod
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
! size of project directory string
integer, intent(in):: iPrjLen
! string of project directory
integer, dimension(*), intent(in):: iPrjDir
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
! coordiantes of integration point
double precision, intent(in) :: X, Y, Z
! excess pore pressure at the end of previous step:
double precision, intent(in) :: swp0
! time at the end of previous step:
double precision, intent(in) :: Time0
! increment of time:
double precision, intent(in) :: dTime
! current excess pore pressure for undrained calculations:
double precision, intent(inout):: swp
! current equivalent bulk modulus of water K/n for undrained calculations:
double precision, intent(inout):: BulkW
! vector of material properties:
double precision, dimension(*), intent(in)   :: props
! stress at the end of previous step:
double precision, dimension(N_STRESS_VECTOR), intent(in)   :: sig0
! current stress:
double precision, dimension(N_STRESS_VECTOR), intent(inout):: sig
! increment of strain:
double precision, dimension(N_STRESS_VECTOR), intent(in)   :: dEps
! state variables at the end of previous step:
double precision, dimension(*), intent(inout):: stVar0
! current state variables:
double precision, dimension(*), intent(inout):: stVar
! stiffness matrix:
double precision, dimension(N_STRESS_VECTOR, N_STRESS_VECTOR), intent(inout):: D


!DEC$ ATTRIBUTES DLLExport,StdCall,reference :: User_Mod

call OpenDebugFileAndWriteParameters(iPrjDir, iPrjLen, Props, iMod)

select case (iMod)
case(IMOD_ELASTIC_PERFECTLY_PLASTIC)
  call ElasticPerfectlyPlastic(IDTask, isUndr,                    &
                               iStep, iTer, iEl, intPt,           &
                               Time0, dTime,                      &
                               Props, Sig0, Swp0, stVar0,         &
                               dEps, D, BulkW,                    &
                               Sig, Swp, stVar, ipl,              &
                               nStat,                             &
                               NonSym, iStrsDep, iTimeDep, iTang, &
                               iAbort )
case(IMOD_ELASTIC_VISCO_PLASTIC)
  call ElasticViscoPlastic(IDTask, isUndr,                    &
                           iStep, iTer, iEl, intPt,           &
                           Time0, dTime,                      &
                           Props, Sig0, Swp0, stVar0,         &
                           dEps, D, BulkW,                    &
                           Sig, Swp, stVar, ipl,              &
                           nStat,                             &
                           NonSym, iStrsDep, iTimeDep, iTang, &
                           iAbort )
case Default
  call SetError('invalid model number in UsrMod:'//trim(String(iMod)))
  call SetError('IDTask: '//trim(String(IDTask)))
end select

if (CatchError()) iAbort = 1

end subroutine User_Mod

