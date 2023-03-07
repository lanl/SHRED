module Density_of_states_type
   use types, only : dp

   implicit none
   
   public Density_of_states_struct
   
   type:: Density_of_states_struct
      real(dp) :: a
      real(dp) dw
      integer  :: nw, output_unit
      real(dp), allocatable :: dos(:), fdos(:)
      real(dp), allocatable :: det_dos(:), det_fdos(:)
      real(dp), allocatable :: stoc_dos(:), stoc_fdos(:)
   end type Density_of_states_struct

end module