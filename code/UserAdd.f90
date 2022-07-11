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
! Subroutines in this file:
! Exported:
!  subroutine GetModelCount( nMod )
!  subroutine GetModelName ( iMod , ModelName )
!  subroutine GetParamCount( iMod , nParam )
!  subroutine GetParamName ( iMod , iParam, ParamName )
!  subroutine GetParamUnit ( iMod , iParam, Units )
!  subroutine GetStateVarCount( iMod , nVar )
!  subroutine GetStateVarName ( iMod , iVar, Name )
!  subroutine GetStateVarUnit ( iMod , iVar, Unit )
!
! Local:
!  subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
!  subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
!
subroutine GetModelCount(nMod)
use ModGlobalConstants
!
! return the maximum model number (iMod) in this DLL
!
implicit none
integer (Kind=4) nMod

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetModelCount

nMod = N_MODELS ! Maximum model number (iMod) in current DLL

end subroutine GetModelCount

!--------------------------------------------------------------------------------------------------
subroutine GetModelName( iMod , ModelName )
use ModGlobalConstants
use ElasticPerfectlyPlastic_Module
use ElasticViscoPlastic_Module
!
! return the name of the different models
!
implicit none
integer  iMod
character (Len= * ) ModelName
character (Len=255) tName
integer lt

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetModelName

select case (iMod)
  case(IMOD_ELASTIC_PERFECTLY_PLASTIC)
    call GetModelNameElasticPerfectlyPlastic(tName)
  case(IMOD_ELASTIC_VISCO_PLASTIC)
    call GetModelNameElasticViscoPlastic(tName)
  case default
    tName = 'N/A'
end select
LT = len_trim(tName)
ModelName= char(lt) // tName(1:Lt)

end subroutine GetModelName

!--------------------------------------------------------------------------------------------------
subroutine GetParamCountInternal( iMod , nParam )
use ModGlobalConstants
!
! return the number of parameters of the different models
!
use ElasticPerfectlyPlastic_Module,   only: GetParamCountElasticPerfectlyPlastic
use ElasticViscoPlastic_Module, only: GetParamCountElasticViscoPlastic


implicit none
integer, intent(in):: iMod
integer, intent(out) :: nParam

select case (iMod)
  case(IMOD_ELASTIC_PERFECTLY_PLASTIC)
    nParam = GetParamCountElasticPerfectlyPlastic()
  case(IMOD_ELASTIC_VISCO_PLASTIC)
    nParam = GetParamCountElasticViscoPlastic()
  case default
    nParam = 0
end select

end subroutine GetParamCountInternal

!--------------------------------------------------------------------------------------------------
subroutine GetParamCount( iMod , nParam )
use ModGlobalConstants
!
! return the number of parameters of the different models
!
implicit none
integer, intent(in):: iMod
integer, intent(out) :: nParam
external GetParamCountInternal

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetParamCount
call GetParamCountInternal( iMod , nParam )

end subroutine GetParamCount

!--------------------------------------------------------------------------------------------------
subroutine GetParamName( iMod , iParam, ParamName )
!
! return the parameters name of the different models
!
implicit none
integer, intent(in):: iMod, iParam

character (Len=255) ParamName, Units
external GetParamAndUnit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetParamName
call GetParamAndUnit(iMod,iParam,ParamName,Units)

end subroutine GetParamName

!--------------------------------------------------------------------------------------------------
subroutine GetParamUnit( iMod , iParam, Units )
!
! return the units of the different parameters of the different models
!
implicit none
integer, intent(in):: iMod, iParam
character (Len=255) ParamName, Units
external GetParamAndUnit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetParamUnit

call GetParamAndUnit(iMod,iParam,ParamName,Units)

end subroutine GetParamUnit

!--------------------------------------------------------------------------------------------------
subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
use ModGlobalConstants
!
! return the parameters name and units of the different models
!
! Units: use F for force unit
!            L for length unit
!            T for time unit
!
use ElasticPerfectlyPlastic_Module
use ElasticViscoPlastic_Module

implicit none
integer, intent(in):: iMod, iParam
character (Len=255) ParamName, Units, tName
integer :: lt

select case (iMod)
  case(IMOD_ELASTIC_PERFECTLY_PLASTIC)
    call GetParamAndUnitElasticPerfectlyPlastic(iParam, ParamName, Units)
  case(IMOD_ELASTIC_VISCO_PLASTIC)
    call GetParamAndUnitElasticViscoPlastic(iParam, ParamName, Units)
  case default
    ! model not in DLL
    ParamName = ' N/A '        ; Units     = ' N/A '
end select
tName    = ParamName
LT       = len_trim(tName)
ParamName= char(lt) // tName(1:Lt)

tName = Units
LT    = len_trim(tName)
Units = char(lt) // tName(1:Lt)

end subroutine GetParamAndUnit

!--------------------------------------------------------------------------------------------------
subroutine GetStateVarCount( iMod , nVar )
use ModGlobalConstants
!
! return the number of state variables of the different models
!
use ElasticPerfectlyPlastic_Module, &
  only: GetStateVarCountElasticPerfectlyPlastic
use ElasticViscoPlastic_Module, &
  only: GetStateVarCountElasticViscoPlastic

implicit none
integer, intent(in):: iMod
integer, intent(out):: nVar

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetStateVarCount

select case (iMod)
case(IMOD_ELASTIC_PERFECTLY_PLASTIC)
  nVar = GetStateVarCountElasticPerfectlyPlastic()
case(IMOD_ELASTIC_VISCO_PLASTIC)
  nVar = GetStateVarCountElasticViscoPlastic()
case default
  nVar = 0
end select

end subroutine GetStateVarCount

!--------------------------------------------------------------------------------------------------
subroutine GetStateVarName( iMod , iVar, Name )
!
! return the name of the different state variables
! of the different models
!
use ModGlobalConstants
implicit none
integer, intent(in):: iMod, iVar

character (Len=255) Name, Unit
external GetStateVarNameAndUnit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetStateVarName

select case(iMod)
!case(IMOD_ABC_MODEL)
!  call GetStateVarNameABCModel( iMod , iVar, Name )
case default
  call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
end select

end subroutine GetStateVarName

!--------------------------------------------------------------------------------------------------
subroutine GetStateVarUnit( iMod , iVar, Unit )
!
! return the units of the different state variables of the different models
!
implicit none
integer, intent(in):: iMod, iVar
character (Len=255) Name, Unit
external GetStateVarNameAndUnit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetStateVarUnit

call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )

end subroutine GetStateVarUnit

!--------------------------------------------------------------------------------------------------
subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
use ModGlobalConstants

!
! return the name and unit of the different state variables of the different models
!
implicit none
integer, intent(in):: iMod, iVar
character (Len=255) Name, Unit, tName
integer :: lt
select case (iMod)
!case (IMOD_DELTA_SAND_VERSION_6)
!  call GetStateVarNameAndUnitDeltaSandVersion_6(iVar, Name, Unit)
case default
  Name='N/A'                   ; Unit = '?'
end select

tName = Unit
LT    = len_trim(tName)
Unit  = char(lt) // tName(1:Lt)

end subroutine GetStateVarNameAndUnit
