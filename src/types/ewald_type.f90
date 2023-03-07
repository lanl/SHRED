module ewald_type
   use types, only : dp

   
   implicit none
   
   public ewald_struct, ntable
  
   integer, parameter :: ntable=5000

   
   type:: ewald_struct
      real(dp) :: energy, alpha, preuk, unitk(3), kmaxsq, delrsq
      integer :: kmax(3) 
      real(dp) ::  delta(3,3), stress(3,3)
      real(dp), allocatable :: urtable(:), frtable(:)
      real(dp), allocatable :: sinkr(:,:,:), coskr(:,:,:), ukarray(:,:)
      real(dp), allocatable :: cosat(:)
      real(dp), allocatable :: sinat(:)
      real(dp), allocatable :: ctmp(:), stmp(:), v(:)
      real(dp), allocatable :: E_field(:,:)

   end type ewald_struct

end module