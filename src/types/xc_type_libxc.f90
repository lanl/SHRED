module xc_type
    use xc_f03_lib_m
    use iso_c_binding
   
   implicit none
   
   public xc_struct
   
   type:: xc_struct
      TYPE(xc_f03_func_t) :: func
      TYPE(xc_f03_func_info_t) :: info
      integer :: id, vmajor, vminor, vmicro
      integer(c_int) :: family, xc_kind
   end type xc_struct

end module