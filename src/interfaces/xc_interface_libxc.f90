module xc_interface

!    use xc_f03_types_m
    use xc_f03_lib_m
    use iso_c_binding
   
    use types, only: dp
    use xc_type
   implicit none

   public update_xc, init_xc


   interface update_xc
    module procedure update_xc_3D
   end interface

   contains
   subroutine update_xc_3D(rho, xc, sigma, orbitals, exc, vxc, sxc)
    use constants, only: i_
    real(dp), intent(in) :: rho(:,:,:,:) ! charge density
    type(xc_struct), intent(in) :: xc(:) !1-correlation/exchnage_correlation 2-exchange
    real(dp), intent(in), optional :: orbitals(:,:,:), sigma(:,:,:,:)! the orbitals, sigma
    real(dp), intent(out), optional :: exc(:,:,:), vxc(:,:,:,:), sxc(:,:,:,:)! XC energy density, potential, stress tensor
    real(dp), dimension(:,:,:,:), allocatable :: tmp_v, tmp_vs
    real(dp), dimension(:,:,:), allocatable ::  tmp_e
    integer(c_size_t) :: size_rho
    integer :: i

    if(present(orbitals)) then
        print *, 'No Meta GGA support'
    endif

    if(present(exc)) exc=0
    if(present(vxc)) vxc=0
    if(present(sxc)) sxc=0

    if(size(xc).gt.1) then
        if(present(exc)) allocate(tmp_e(size(exc, 1),size(exc, 2),size(exc, 3)))
        if(present(vxc)) allocate(tmp_v(size(vxc, 1),size(vxc, 2),size(vxc, 3),size(vxc, 4)))
        if( any(xc(:)%family.eq.XC_FAMILY_GGA) .and. present(vxc)) &
             allocate(tmp_vs(size(sigma, 1),size(sigma, 2),size(sigma, 3),size(sigma, 4)))
    endif
   size_rho=size(rho)

   do i=1,size(xc)
    select case (xc(i)%family)
    case (XC_FAMILY_LDA);
        if(present(vxc).and.present(exc)) then
                if(size(xc).gt.1) then
                 call xc_f03_lda_exc_vxc(xc(i)%func, size_rho, rho, tmp_e, tmp_v)
                 exc=exc+tmp_e
                 vxc=vxc+tmp_v
                else
                 call xc_f03_lda_exc_vxc(xc(i)%func, size_rho, rho, exc, vxc)
                endif
        else if(present(vxc)) then
                if(size(xc).gt.1) then
                 call xc_f03_lda_vxc(xc(i)%func, size_rho, rho, tmp_v)
                 vxc=vxc+tmp_v
                 else
                 call xc_f03_lda_vxc(xc(i)%func, size_rho, rho, vxc)
                 endif
                else if(present(exc)) then
                    if(size(xc).gt.1) then
                        call xc_f03_lda_exc(xc(i)%func, size_rho, rho, tmp_e)
                        exc=exc+tmp_e
                     else
                     call xc_f03_lda_exc(xc(i)%func, size_rho, rho, exc)
                     endif
        else
        endif
    case (XC_FAMILY_GGA);
        if(present(vxc).and.present(exc)) then
            if(size(xc).gt.1) then
            call xc_f03_gga_exc_vxc(xc(i)%func, size_rho, rho, sigma, &
                                    tmp_e, tmp_v, tmp_vs)
                exc=exc+tmp_e
                vxc=vxc+tmp_v
                sxc=sxc+tmp_vs
            else
                call xc_f03_gga_exc_vxc(xc(i)%func, size_rho, rho, sigma, exc, vxc, sxc)
            endif
        else if(present(vxc)) then
            if(size(xc).gt.1) then
                call xc_f03_gga_vxc(xc(i)%func, size_rho, rho, sigma, tmp_v, tmp_vs)
                vxc=vxc+tmp_v
                sxc=sxc+tmp_vs
            else
                call xc_f03_gga_vxc(xc(i)%func, size_rho, rho, sigma, vxc, sxc)
            endif
        else if(present(exc)) then
            if(size(xc).gt.1) then
                call xc_f03_gga_exc(xc(i)%func, size_rho, rho, sigma, tmp_e)
                exc=exc+tmp_e
            else
                call xc_f03_gga_exc(xc(i)%func, size_rho, rho, sigma, exc)
            endif
        else
        endif
    case (XC_FAMILY_HYB_GGA);
    case (XC_FAMILY_MGGA);
    case (XC_FAMILY_HYB_MGGA);
    case default;   
    end select

enddo

end subroutine

subroutine init_xc(xc, functional, T)
   ! use xc_f03_types_m
    use xc_f03_lib_m
      character(len=*),intent(in) :: functional
      type(xc_struct), intent(inout) :: xc
      real(dp), optional, intent(in) :: T
      real(dp) :: params(20)
   
      call xc_f03_version(xc%vmajor, xc%vminor, xc%vmicro)
   !   write(*,'("Libxc version: ",I1,".",I1,".",I1)') xc%vmajor, xc%vminor, xc%vmicro
      xc%id=xc_f03_functional_get_number(trim(functional))
   !   call xc_f03_func_init(xc%func, xc%info, xc%id, XC_UNPOLARIZED)
      call xc_f03_func_init(xc%func, xc%id, XC_UNPOLARIZED)
   
      xc%info=xc_f03_func_get_info(xc%func)
      xc%family=xc_f03_func_info_get_family(xc%info)
      xc%xc_kind=xc_f03_func_info_get_kind(xc%info)
   
      if(functional.eq.'LDA_XC_KSDT') then
           params(1)=T
           call xc_f03_func_set_ext_params(xc%func, params(1))
      endif
   
   end subroutine

end module

   
   
