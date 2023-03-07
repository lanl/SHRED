module tf_type
   use types, only : dp

   
   implicit none
   
   public tf_struct
   
   type:: tf_struct
      real(dp) lambda
      real(dp) gamma
      integer max_iter
      real(dp) energy_eps
      real(dp) brent_eps
      integer update_type
      logical :: dynamic_kedp=.false.
   end type tf_struct

end module