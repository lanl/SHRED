module types
    use, intrinsic :: iso_c_binding
implicit none
private
public :: sp, dp, hp, op

integer, parameter :: dp=kind(0.d0)            ! double precision
integer, parameter :: hp=selected_real_kind(15)    ! high precision
integer, parameter :: sp = kind(0.)                ! single precision

integer, parameter :: op=kind(0.d0)                ! double precision

contains

end module