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
module ModGlobalConstants
implicit none

!> model IDs
integer, parameter:: IMOD_ELASTIC_PERFECTLY_PLASTIC = 1
integer, parameter:: IMOD_ELASTIC_VISCO_PLASTIC     = 2

!> model IDs
integer, parameter:: N_MODELS = 2

! Mathematical constants
double precision, parameter :: PI = 3.141592653589793238462643383279502884197169399d0
double precision, parameter :: RADIANS = PI/180.0d0
double precision, parameter :: DEGREES = 180.0d0 / PI

! stress-strain constants
integer, parameter :: N_STRESS_VECTOR = 6
integer, parameter :: N_PRINCIPAL_STRESS_VECTOR = 3

! limit constants
double precision, parameter :: RELATIVELY_SMALL = 1d-4
double precision, parameter :: SMALL = 1d-15
double precision, parameter :: TINY  = 1d-60
double precision, parameter :: LARGE = 1D10
double precision, parameter :: RELATIVELY_LARGE = 1D6

end module ModGlobalConstants

