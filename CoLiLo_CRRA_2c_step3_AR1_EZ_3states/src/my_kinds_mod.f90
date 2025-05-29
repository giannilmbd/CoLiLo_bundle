module my_kinds_mod
use iso_fortran_env, only: wp => real64 !real32!
use iso_fortran_env, only: sp => real32
use iso_fortran_env, only: dp => real64
use iso_fortran_env, only: i4 => int32
use iso_fortran_env, only: i8 => int64
implicit none
public wp,dp, sp, i4, i8
end module 