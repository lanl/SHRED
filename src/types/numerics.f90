module numerics_type
    use types, only : dp

    implicit none
    
    private

    public  :: numerics_struct, pulay_numerics, diis_numerics, pulay_pointers, pulay_real_pointers

type lobpcg_numerics
    integer :: inner_steps
    real(dp) ::  soft_lock_thresh
    integer :: n_init, PC_type
end type

type cheby_numerics
    integer :: inner_steps
end type

type pulay_pointers
    complex(dp), pointer :: ppointer(:)
endtype
type pulay_real_pointers
    real(dp), pointer :: ppointer(:)
endtype

type pulay_numerics
    integer :: n, k, max_iter, n_skip
    real(dp) ::  L2_eps, eig_eps, L2
    type(pulay_pointers), allocatable :: R_i(:,:), F_i(:,:)
end type

type diis_numerics
    integer :: n_steps, PC_type
    real(dp) :: min_resnorm
end type
type kerker_numerics
    logical :: plasma
    real(dp) :: A, B, C, Amin
end type

type numerics_struct
    integer :: SCF_type, n_cg_init
    logical :: precondition_density
    type(lobpcg_numerics) :: lobpcg
    type(cheby_numerics) :: cheby
    type(pulay_numerics)  :: pulay
    type(diis_numerics) :: diis
    type(kerker_numerics) :: kerker
end type






end module