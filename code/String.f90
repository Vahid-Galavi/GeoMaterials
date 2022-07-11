!
!    Main authors:    Vahid Galavi
!
module ModString
!**********************************************************************
!
!    function:  returns string out of any types of variable.
!
!**********************************************************************
implicit none

! String functions
public:: String         ! makes string out of different types of variables
public:: Copy           !  (generic) char array to string or string to char array
public:: Clen           !  returns same as len      unless last non-blank char = null
public:: Clen_trim      !  returns same as len_trim    "              "
public:: Ctrim          !  returns same as TRIM        "              "
public:: CountItems     !  in string that are blank or comma separated
public:: ReduceBlanks   !  in string to 1 blank between items, last char not blank
public:: ReplaceText    !  in all occurances in string with replacement string
public:: Spack          !  pack string's chars == extract string's chars
public:: Tally          !  occurances in string of text arg
public:: Translate      !  text arg via indexed code table
public:: Upper          !  case the text arg
public:: Lower          !  case the text arg
public:: ContainLetters ! returns true if at least one charcter is letter
public:: GetItem        ! returns a given item in a string
public:: GetFirstLocationNextItem ! returns location of the first character of the next item
public:: GetLastLocationNextItem  ! returns location of the last character of the next item
public:: HasKey

!> return integer value from a string
public:: StringToInteger


interface String
module procedure StringReal
module procedure StringInt
module procedure StringInt8
module procedure StringLogical
module procedure StringString
module procedure StringRealArray
module procedure StringIntArray
module procedure StringRealMatrix
end interface String

interface Copy    ! generic
module procedure CopyArrayToString, CopyStringToArray
end interface Copy

private

contains

! -----------------------------------------------------------------------------------------
character(len = 64) function StringReal(Variable, fmt) result(str)
implicit none
double precision, intent(in):: Variable
character(*), intent(in), optional :: fmt

if (present(fmt)) then
  write(str, fmt) Variable
else
  !write(str, '(E16.6E3)') Variable
  write(str, *) Variable
endif

end function StringReal

! -----------------------------------------------------------------------------------------
character(len = 16) function StringInt(Variable, fmt) result(str)
implicit none
integer, intent(in):: Variable
character(*), intent(in), optional :: fmt

!$OMP CRITICAL
if (present(fmt)) then
  write(str, fmt) Variable
else
  if (Variable < 1000000) then
    write(str, '(I6)') Variable
  else
    write(str, '(I12)') Variable
  endif
endif
!$OMP end CRITICAL

end function StringInt

! -----------------------------------------------------------------------------------------
character(len = 16) function StringLogical(Variable, fmt) result(str)
implicit none
logical, intent(in):: Variable
character(*), intent(in), optional :: fmt
!$OMP CRITICAL
if (present(fmt)) then
  write(str, fmt) Variable
else
  write(str, *) Variable
endif
!$OMP end CRITICAL

end function StringLogical

! -----------------------------------------------------------------------------------------
character(len = 24) function StringInt8(Variable, fmt) result(str)
implicit none
integer(8), intent(in):: Variable
character(*), intent(in), optional :: fmt
!$OMP CRITICAL

if (present(fmt)) then
  write(str, fmt) Variable
else
  write(str, '(I14)') Variable
endif
!$OMP end CRITICAL

end function StringInt8

! -----------------------------------------------------------------------------------------
character(len = 4096) function StringString(Variable, fmt) result(str)
implicit none
character(*), intent(in):: Variable
character(*), intent(in), optional :: fmt
!$OMP CRITICAL

if (present(fmt)) then
  write(str, fmt) Variable
else
  write(str, *) Variable
endif
!$OMP end CRITICAL

end function StringString

! -----------------------------------------------------------------------------------------
character(len = 4096) function StringRealArray(Variable, first, last) result(str)
implicit none
double precision, intent(in):: Variable(*)
integer, intent(in):: first, last
character(len = 256) :: strTemp
integer :: i
!$OMP CRITICAL

str = ''
do i=first, last-1
  write(strTemp, *) Variable(i), ', '
  str = trim(str) // trim(strTemp)
enddo
write(strTemp, *) Variable(last)
str = trim(str) // trim(strTemp)
!$OMP end CRITICAL

end function StringRealArray

! -----------------------------------------------------------------------------------------
character(len = 4096) function StringRealMatrix(Variable, lastRow, lastCol) result(str)
implicit none
double precision, intent(in):: Variable(:,:)
integer, intent(in):: lastRow, lastCol
character(len = 256) :: strTemp
integer :: i, j

!$OMP CRITICAL
str = '['
do i=1, lastRow
  str = trim(str) // '('
  do j=1,lastCol-1
    write(strTemp, *) Variable(i,j), ', '
    str = trim(str) // trim(strTemp)
  enddo
  write(strTemp, *) Variable(i,lastCol)
  if (i==lastRow) then
    str = trim(str) // trim(strTemp) // ')]'
  else
    str = trim(str) // trim(strTemp) // ')' // NEW_LINE('A')
  endif
enddo
!$OMP end CRITICAL

end function StringRealMatrix

! -----------------------------------------------------------------------------------------
character(len = 4096) function StringIntArray(Variable, first, last) result(str)
implicit none
integer, intent(in):: Variable(*)
integer, intent(in):: first, last
!$OMP CRITICAL

write(str, *) Variable(first:last)
!$OMP end CRITICAL

end function StringIntArray

! -----------------------------------------------------------------------------------------
logical function ContainLetters(s) result(res)
implicit none
character(*),intent(IN) :: s
integer :: i

res = .false.
do i=65,90
  ! Capital letters
  res = scan(s, char(i)) > 0
  if (res) RETURN
enddo

do i=97,122
  ! small letters
  res = scan(s, char(i)) > 0
  if (res) RETURN
enddo

end function ContainLetters

! -----------------------------------------------------------------------------------------
pure function CopyArrayToString(a) result(s)    ! copy char array to string
implicit none
character,intent(IN) :: a(:)
character(SIZE(a)) :: s
integer :: i
do i = 1,SIZE(a)
  s(i:i) = a(i)
end do

end function CopyArrayToString

! -----------------------------------------------------------------------------------------
pure function CopyStringToArray(s) result(a)   ! copy s(1:Clen(s)) to char array
implicit none
character(*),intent(IN) :: s
character :: a(len(s))
integer :: i
do i = 1,len(s)
  a(i) = s(i:i)
end do
end function CopyStringToArray

! -----------------------------------------------------------------------------------------
pure integer function Clen(s)      ! returns same result as len unless:
implicit none
character(*),intent(IN) :: s       ! last non-blank char is null
integer :: i
Clen = len(s)
i = len_trim(s)
if (s(i:i) == char(0)) Clen = i-1  ! len of C string

end function Clen

! -----------------------------------------------------------------------------------------
pure integer function Clen_trim(s) ! returns same result as len_trim unless:
implicit none
character(*),intent(IN) :: s       ! last char non-blank is null, if true:
integer :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
i = len_trim(s) ; Clen_trim = i
if (s(i:i) == char(0)) Clen_trim = Clen(s)   ! len of C string

end function Clen_trim

! -----------------------------------------------------------------------------------------
function Ctrim(s1) result(s2)     ! returns same result as TRIM unless:
implicit none
character(*),intent(IN)  :: s1     ! last non-blank char is null in which
character(Clen_trim(s1)) :: s2     ! case trailing blanks prior to null
s2 = s1                            ! are output

end function Ctrim

! -----------------------------------------------------------------------------------------
integer function CountItems(s)  ! in string or C string that are blank or comma separated
implicit none
character(*), intent(in) :: s
!character(Clen(s1)) :: s
integer :: i, k

k = 0
do i = 1,len_trim(s)
  if (.not.IsDelimiter(s(i:i)) .and. IsDelimiter(s(i+1:i+1))) k = k+1
end do
CountItems = k

end function CountItems

! -----------------------------------------------------------------------------------------
subroutine GetItem(s, item, outs)
implicit none
character(*), intent(in):: s
character(*), intent(out):: outs
integer, intent(in) :: item
integer itemCounter, length, iFirstOfWord, iLastOfWord

outs = ' '

! find location of given item:
iFirstOfWord = GetFirstLocationNextItem(s)
iLastOfWord  = GetLastLocationNextItem(s)

do itemCounter=1,item
  if (itemCounter == item) then
    length = iLastOfWord - iFirstOfWord + 1
    outs(1:length) = s(iFirstOfWord:iLastOfWord)
    RETURN
  else
    iFirstOfWord = GetFirstLocationNextItem(s,iLastOfWord)
    iLastOfWord  = GetLastLocationNextItem(s,iLastOfWord)
  endif
enddo

end subroutine GetItem

! -----------------------------------------------------------------------------------------
function ReduceBlanks(s) result(outs)
implicit none
character(*)      :: s
character(len_trim(s)) :: outs
integer           :: i, k, n

n = 0  ; k = len_trim(s)          ! k=index last non-blank (may be null)
do i = 1,k-1                      ! dont process last char yet
  n = n+1 ; outs(n:n) = s(i:i)
  if (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
end do
n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
if (n < k) outs(n+1:) = ' '       ! pad trailing blanks

end function ReduceBlanks

! -----------------------------------------------------------------------------------------
function ReplaceText(s,text,rep) result(outs)
implicit none
character(*)        :: s,text,rep
character(len(s)+100) :: outs     ! provide outs with extra 100 char len
integer             :: i, nt, nr

outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
do
  i = index(outs,text(:nt)) ; if (i == 0) EXIT
  outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
end do

end function ReplaceText

! -----------------------------------------------------------------------------------------
function Spack(s,ex) result(outs)
implicit none
character(*) :: s,ex
character(len(s)) :: outs
character :: aex(len(ex))   ! array of ex chars to extract
integer   :: i, n

n = 0  ;  aex = Copy(ex)
do i = 1,len(s)
  if (.NOT.any(s(i:i) == aex)) CYCLE   ! dont pack char
  n = n+1 ; outs(n:n) = s(i:i)
end do
outs(n+1:) = ' '     ! pad with trailing blanks

end function Spack

! -----------------------------------------------------------------------------------------
integer function Tally(s,text)
implicit none
character(*) :: s, text
integer :: i, nt

Tally = 0 ; nt = len_trim(text)
do i = 1,len(s)-nt+1
  if (s(i:i+nt-1) == text(:nt)) Tally = Tally+1
end do

end function Tally

! -----------------------------------------------------------------------------------------
function Translate(s1,codes) result(s2)
implicit none
character(*)       :: s1, codes(2)
character(len(s1)) :: s2
character          :: ch
integer            :: i, j

do i = 1,len(s1)
  ch = s1(i:i)
  j = index(codes(1),ch) ; if (j > 0) ch = codes(2)(j:j)
  s2(i:i) = ch
end do

end function Translate

! -----------------------------------------------------------------------------------------
function Upper(s1) result(s2)
implicit none
character(*)       :: s1
character(len(s1)) :: s2
character          :: ch
integer,parameter  :: DUC = ichar('A') - ichar('a')
integer            :: i

do i = 1,len(s1)
  ch = s1(i:i)
  if (ch >= 'a'.and.ch <= 'z') ch = char(ichar(ch) + DUC)
  s2(i:i) = ch
end do

end function Upper

! -----------------------------------------------------------------------------------------
function Lower(s1) result (s2)
implicit none
character(*)       :: s1
character(len(s1)) :: s2
character          :: ch
integer, parameter :: DUC = ichar('A') - ichar('a')
integer            :: i

do i = 1,len(s1)
  ch = s1(i:i)
  if (ch >= 'A'.and.ch <= 'Z') ch = char(ichar(ch) - DUC)
  s2(i:i) = ch
end do
end function Lower
! -----------------------------------------------------------------------------------------
integer function GetFirstLocationNextItem(s, LocLastCharacter) result(location)
implicit none
character(*), intent(in):: s
integer, intent(in), optional :: LocLastCharacter
integer i, locLast
        
  if (present(LocLastCharacter)) then
  locLast = LocLastCharacter
else
  locLast = 0
endif

! find location of first character:
location = len_trim(s)
do i = locLast+1,len_trim(s)
  if (.not.IsDelimiter(s(i:i))) then
    location = i
    RETURN
  endif
enddo

end function GetFirstLocationNextItem

! -----------------------------------------------------------------------------------------
integer function GetLastLocationNextItem(s, LocLastCharacter) result(location)
implicit none
character(*), intent(in):: s
integer, intent(in), optional :: LocLastCharacter
integer i, locLast,iFirst

if (present(LocLastCharacter)) then
  locLast = LocLastCharacter
else
  locLast = 0
endif

! find location of first character:
location = len_trim(s)
do i = locLast+1,len_trim(s)
  if (.not.IsDelimiter(s(i:i))) then
    iFirst = i
    EXIT
  endif
enddo

do i = iFirst,len_trim(s)
  if (.not.IsDelimiter(s(i:i)) .and. IsDelimiter(s(i+1:i+1)) ) then
    location = i
    RETURN
  endif
enddo

end function GetLastLocationNextItem

! -----------------------------------------------------------------------------------------
logical function IsDelimiter(c) result(res)
implicit none
character, intent(in):: c

res =      c == ' '  &
      .or. c == ','  &
      .or. c == '('  &
      .or. c == ')'  &
      .or. c == '['  &
      .or. c == ']'  &
      .or. c == '{'  &
      .or. c == '}'  &
      .or. c == '''' &
      .or. c == '"'  &
      .or. c == ';'  &
      .or. c == ':'

end function IsDelimiter

! -----------------------------------------------------------------------------------------
logical function HasKey(s,key) result(res)
implicit none
character(*), intent(in):: s
character(*), intent(in):: key
character(len_trim(s)) :: itemString
integer :: item

res = .false.
do item=1,CountItems(s)
  call GetItem(s, item, itemString)
  if (trim(lower(key)) == trim(lower(itemString))) then
    res = .true.
    RETURN
  endif
enddo
end function HasKey

! -----------------------------------------------------------------------------------------
integer function StringToInteger(str, stat) result(int)
implicit none
character(*), intent(in):: str
integer,intent(out)     :: stat

read(str,*,iostat=stat)  int

end function StringToInteger

end module ModString
