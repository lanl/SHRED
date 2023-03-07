module stopping_type
    use types, only : dp
    use atom_type, only : atom_struct
    use element_type, only : element_struct

    implicit none
    
    private 

    public stopping_struct
    
    type stopping_struct
        logical  :: do_stopping=.false., remove_t0=.false.
        type(atom_struct) :: atom(1)
        type(element_struct) :: element(1)
        integer :: output_unit
        integer :: stop_at, stop_el
        real(dp) :: power, power_average, t0, R0(3)=-1.0_dp
        real(dp), allocatable ::  density_t(:)
    end type

end module