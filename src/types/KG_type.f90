module KG_type
   use types, only : dp
   use odp_type, only : field_struct

   
   implicit none
   
   public KG_struct
   

   integer :: SIL=1
   integer :: Cheb=2

   type:: KG_struct
      real(dp) delta, err_allowed, dw
      integer  :: nw, output_unit, output_unit_time
      real(dp) energy_gap_min
      complex(dp), allocatable :: Onsager(:,:,:)
      integer :: nt, SIL_rank, Prop
      real(dp) :: dt
      type(field_struct), allocatable :: Q(:,:,:,:)
      logical :: filter_stoc, project_out_J
   end type KG_struct

end module