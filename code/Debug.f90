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
module ModDebug
use ModString

implicit none
private

public :: WriteInDebugFile
public :: OpenDebugFile
public :: OpenDebugFileAndWriteParameters
public :: FlushInFile
public :: setIsDebug
public :: RELEASE_LEVEL, DEBUG_LEVEL, ALWAYS_LEVEL
public :: IS_MINOR_DEBUG, AUTOMATIC_DEBUG_FILE
public :: CatchError
public :: SetError

! debugging
logical, parameter :: OPEN_DEBUG_FILE = .false.
logical, parameter :: IS_MINOR_DEBUG  = .false.
logical, parameter :: AUTOMATIC_DEBUG_FILE = .false.
integer, parameter :: INT_DEBUG_ALL_DATA = -1
integer, parameter :: IEL_DEBUG_ALL_DATA = -1

!> Feedback level
integer, parameter :: RELEASE_LEVEL  = 1
integer, parameter :: DEBUG_LEVEL    = 2
integer, parameter :: ALWAYS_LEVEL   = -huge(0)
integer, parameter :: FEEDBACK_LEVEL = DEBUG_LEVEL

! file units
integer, parameter :: DebugUnit = 11

! error
logical:: IsError = .false.

contains

!--------------------------------------------------------------------------------------------------
logical function setIsDebug(iEl, intPt, INT_DEBUG_LOCAL) result(IsDebug)
implicit none
integer, intent(in):: iEl, intPt, INT_DEBUG_LOCAL

if (iEl == IEL_DEBUG_ALL_DATA .and. (intPt == INT_DEBUG_ALL_DATA .or. intPt == INT_DEBUG_LOCAL)) then
  IsDebug = .true.
else
  IsDebug = .false.
endif

end function setIsDebug

!--------------------------------------------------------------------------------------------------
subroutine WriteInDebugFile(message, FeedbackLevel)
implicit none
! arguments:
character*(*), intent(in) :: message
integer, optional, intent(in):: FeedbackLevel
! local variables
integer level
integer, external :: GetFeedbackLevel

if (OPEN_DEBUG_FILE) then
  if (present(FeedbackLevel)) then
    level = FeedbackLevel
  else
    level = DEBUG_LEVEL
  endif

  if (level <= FEEDBACK_LEVEL) then
    !$OMP CRITICAL
    write (DebugUnit, fmt='(a)') trim(message)
    !$OMP end CRITICAL
  endif
endif

end subroutine WriteInDebugFile

!--------------------------------------------------------------------------------------------------
subroutine FlushInFile()
implicit none

flush(DebugUnit)

end subroutine FlushInFile

!--------------------------------------------------------------------------------------------------
subroutine OpenDebugFile(iPrjDir, iPrjLen, FileName)
implicit none
integer, intent(in) :: iPrjLen, iPrjDir(*)
character*255, optional, intent(in):: FileName

character*100 BaseName
character*255 PrjDir, debugName
integer i, nErr, ios

if (OPEN_DEBUG_FILE) then
  if (.not.present(FileName)) then
    BaseName = 'UDSMDebug'

    PrjDir=' '
    do i=1,iPrjLen
      PrjDir(i:i) = char(iPrjDir(i))
    end do

    debugName=PrjDir(:iPrjLen)//'data.'//trim(BaseName)//'.rr0'
    nErr=0
    do while(nErr < 10)
      !$OMP CRITICAL
      open(Unit= DebugUnit, File= debugName, iostat=ios)
      if (ios == 0) close(Unit=DebugUnit, Status='delete',iostat=ios)
      !$OMP end CRITICAL

      if (ios /= 0) then
        nErr=nErr+1
        debugName = PrjDir(:iPrjLen)//'data.'// trim(BaseName)//char(48+nErr)//'.rr0'
      else
        exit
      end if
    enddo
  else
    debugName=trim(FileName)
  endif

  open(unit=DebugUnit, file= debugName)
endif

end subroutine OpenDebugFile

!--------------------------------------------------------------------------------------------------
subroutine OpenDebugFileAndWriteParameters(iPrjDir, iPrjLen, Props, iMod)
implicit none
integer, intent(in) :: iPrjLen
integer, intent(in) :: iMod
integer,          dimension(*), intent(in):: iPrjDir
double precision, dimension(*), intent(in):: Props

external GetParamCountInternal

! local variables
integer :: nParam = 50
logical, save :: isOpen = .false.

if (OPEN_DEBUG_FILE .and. .not.isOpen) then
  ! open a file for debugging purposes
  call openDebugFile(iPrjDir, iPrjLen)

  call WriteInDebugFile('File 1 opened')
  call GetParamCountInternal( iMod , nParam )
  call WriteInDebugFile('props:' // trim(String(props, 1, nParam)))
  isOpen = .true.
endif

end subroutine OpenDebugFileAndWriteParameters

!--------------------------------------------------------------------------------------------------
logical function CatchError() result(res)
implicit none

res = IsError

end function CatchError

!--------------------------------------------------------------------------------------------------
subroutine SetError(message)
implicit none
character*(*), intent(in) :: message

IsError = .true.
call WriteInDebugFile(message, ALWAYS_LEVEL)

end subroutine SetError

end module ModDebug

