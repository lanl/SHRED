module td_type
   use types, only : dp
   use odp_type, only : field_struct   
   implicit none
   
   public td_struct, TD_BOMD, TD_RT, TD_None, IsoKinetic, IsoEnergy

   integer :: TD_BOMD=1
   integer :: TD_RT=2
   integer :: TD_None=3

   integer :: IsoKinetic=1,IsoEnergy=2

   
   type:: td_struct
      real(dp) :: time, last_time, dt, err_allowed
      real(dp) :: total_time
      real(dp) :: E0(3), A(3), Ex(3)
      real(dp) :: t_max_field, w_field, tw, t0, phase_field
      integer :: nt, SIL_rank, Nc
      integer :: type, ETRS_steps
      integer :: thermostat, it, p_print_mod, e_print_mod
      type(field_struct), allocatable :: Q(:,:,:,:)
      real(dp), allocatable :: Ebar(:,:), deltaE(:,:)
      complex(dp), allocatable :: Cheby_Mn(:,:,:), Cheby_Cn(:,:,:)
   end type td_struct

end module