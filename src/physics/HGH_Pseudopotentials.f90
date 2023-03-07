module HGH_Pseudopotentials
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: Calculate_Projector_HGH_RL, Calculate_deriv_Projector_HGH_RL, HGH_calc_dVen0GdG2

    contains

    subroutine Calculate_Projector_HGH_RL(Rvec, element, ll, mm, j, h_or_k, projector)
        use  element_type, only : element_struct
        use Spherical_Harmonic_mod, only : real_spherical_harmonic2

        type(element_struct), intent(inout) :: element
        real(dp), intent(in) :: Rvec(3)
        real(dp),intent(inout) :: projector
        integer, intent(in) :: ll, mm, j
        character, intent(in) :: h_or_k
        integer ::  i
        real(dp) :: R, tmp

        R=sqrt(Rvec(1)**2+Rvec(2)**2+Rvec(3)**2)
        projector=0.0
        tmp=0.0
        if(h_or_k.eq.'h') then
            do i=1,element%HGH%nlh(ll+1) 
            tmp= tmp + element%HGH%U_HE(i,j,ll+1)*HGH_Real_Radial(R, i, ll, element)
            enddo 
        else
            do i=1,element%HGH%nlk(ll+1) 
                tmp= tmp + element%HGH%U_HK(i,j,ll+1)*HGH_Real_Radial(R, i, ll, element)
            enddo 
        endif
        projector = real_spherical_harmonic2(ll, mm, Rvec)
        projector = tmp*projector

    end subroutine

    subroutine Calculate_Projector_HGH_GNL(element, ll, mm, j, h_or_k, projector, G2, G)
        use  element_type, only : element_struct
        use Spherical_Harmonic_mod, only : real_spherical_harmonic2
        use  grids_type, only : grid_struct

        type(element_struct), intent(inout) :: element
        real(dp),intent(inout) :: projector(:,:,:)
        integer, intent(in) :: ll, mm, j
        character, intent(in) :: h_or_k
        real(dp), intent(in) :: G2(:,:,:), G(:,:,:,:)
        integer ::  i, x,y,z
        real(dp) :: tmp(size(G2,1),size(G2,2),size(G2,3))
        
        projector=0.0
        tmp=0.0
        if(h_or_k.eq.'h') then
            do i=1,element%HGH%nlh(ll+1) 
                tmp= tmp + element%HGH%U_HE(i,j,ll+1)*HGH_projector_GNL(G2, i, ll, element)
            enddo 
        else
            do i=1,element%HGH%nlk(ll+1) 
                tmp= tmp + element%HGH%U_HK(i,j,ll+1)*HGH_projector_GNL(G2, i, ll, element)
            enddo 
        endif

        do z=1,size(G2,3); do y=1,size(G2,2); do x=1,size(G2,1)
            projector = real_spherical_harmonic2(ll, mm, G(x,y,z,:) )
        enddo;enddo;enddo

        projector = tmp*projector

    end subroutine

    subroutine Calculate_deriv_Projector_HGH_RL(Rvec, element, ll, mm, j, h_or_k, projector)
        use  element_type, only : element_struct
        use Spherical_Harmonic_mod, only : real_spherical_harmonic_and_derivatives

        type(element_struct), intent(inout) :: element
        real(dp), intent(in) :: Rvec(3)
        real(dp),intent(inout) :: projector(:)
        integer, intent(in) :: ll, mm, j
        character, intent(in) :: h_or_k
        
        integer ::  i
        real(dp) :: R2, R, P
        real(dp) :: derP(3)
        real(dp) :: Ylm, dYlm(3)

        R2=Rvec(1)**2+Rvec(2)**2+Rvec(3)**2
        if(R2.lt.tiny(1.0_dp)) then
            projector(:)=0.0_dp
            return
        endif
        
        R=sqrt(R2)
        call real_spherical_harmonic_and_derivatives(ll, mm, Rvec, Ylm, dYlm)

        !The minus here is because it is derivative wrt R_atom and Rvec = R-R_atom
        projector(:)=0.0_dp
        if(h_or_k.eq.'h') then
            do i=1,element%HGH%nlh(ll+1)
                P = HGH_Real_Radial(R, i, ll, element)
                derP(:) = -dYlm(:)*P - Rvec(:)*((-2.0_dp+2*i+ll)/R2-1.0_dp/element%HGH%rl(ll+1)**2)*P*Ylm
                projector(:)= projector(:) + element%HGH%U_HE(i,j,ll+1)*derP(:)
            enddo
        else
            do i=1,element%HGH%nlk(ll+1) 
                P = HGH_Real_Radial(R, i, ll, element)
                derP(:) = -dYlm(:)*P -Rvec(:)*((-2.0_dp+2*i+ll)/R2-1.0/element%HGH%rl(ll+1)**2)*P*Ylm
                projector(:)= projector(:) + element%HGH%U_HK(i,j,ll+1)*derP(:)
            enddo 
        endif

    end subroutine

    real(dp)    function HGH_Real_Radial(R, j, l, element)  result(P_of_r)
        use  element_type, only : element_struct
        type(element_struct), intent(in) :: element
        real(dp), intent(in) :: R
        integer, intent(in) :: j,l

        real(dp) :: tmp

        tmp=l*1.0_dp+(4*j-1)*0.5_dp
        if(R.lt.tiny(1.0_dp)) then
                if((j.eq.1).and.(l.eq.0)) then
                        P_of_r=2.0_dp/element%HGH%rl(l+1)**(1.5_dp)/pi**(0.25_dp) 
                else
                        P_of_r=0.0
                endif
                return
        endif  
        P_of_r=sqrt(2.0_dp)*R**(l+2*(j-1))*exp(-0.5_dp*R**2/element%HGH%rl(l+1)**2)&
                                 / element%HGH%rl(l+1)**(tmp) / sqrt(gamma(tmp))

        
    end function

    real(dp) elemental function HGH_projector_GNL(G2, j, l, element)
    use constants, only: pi54
    use  element_type, only : element_struct
    type(element_struct), intent(in) :: element
    
    integer,intent(in) ::  l, j
    real(dp),intent(in) :: G2
     
    real(dp) :: w, factor, r2, rl
        
        rl=element%HGH%rl(l+1)
        r2=rl**2
        factor=pi54*exp(-0.5_dp*G2*r2)
        
        HGH_projector_GNL=0.0_dp
        if(l.eq.0) then
                if(j.eq.1) then
                    HGH_projector_GNL=Sqrt(32.0_dp*rl**3)
                else if(j.eq.2) then 
                    HGH_projector_GNL=Sqrt(128.0_dp*rl**3/15.0_dp)*(3.0_dp-G2*r2)
                else if(j.eq.3) then
                    HGH_projector_GNL=Sqrt(512.0_dp*rl**3/945.0_dp)*(15.0_dp-10.0_dp*G2*r2+(G2*r2)**2)
                else 
                endif
        else if(l.eq.1) then 
                w=sqrt(G2)
                if(j.eq.1) then
                    HGH_projector_GNL=Sqrt(64.0_dp*rl**5/3.0_dp)*w
                else if(j.eq.2) then
                    HGH_projector_GNL=Sqrt(256.0_dp*rl**5/105.0_dp)*w*(5.0_dp-G2*r2)
                else if(j.eq.3) then
                    HGH_projector_GNL=Sqrt(1024.0_dp*rl**5/10395.0_dp)*w*(35.0_dp-14.0_dp*G2*r2+(G2*r2)**2)
                else 
                endif 
        else if(l.eq.2) then
                if(j.eq.1) then
                    HGH_projector_GNL=Sqrt(128.0_dp*rl**7/15.0_dp)*G2
                else if(j.eq.2) then
                    HGH_projector_GNL=Sqrt(512.0_dp*rl**7/945.0_dp)*G2*(7.0_dp-G2*r2)
                else 
                endif 
        else if(l.eq.3) then
                w=sqrt(G2)
                if(j.eq.1) then
                    HGH_projector_GNL=Sqrt(256.0_dp*rl**9/105.0_dp)*G2*w
                else 
                endif 
        else 
        endif
        HGH_projector_GNL=HGH_projector_GNL*factor
 end function

    subroutine HGH_calc_dVen0GdG2(dVen0GdG2, Ven0G,Zion, rloc, C, G2)
        use constants, only: pi
          real(dp), intent(in) :: rloc, C(4), Zion
          real(dp), intent(in) :: G2(:,:,:)
          real(dp), intent(in) :: Ven0G(:,:,:)
          real(dp), intent(out) :: dVen0GdG2(:,:,:)
          integer :: i,j,k
          real(dp) :: w2
            do k = 1, size(G2,3)
            do j = 1, size(G2,2)
            do i = 1, size(G2,1)
                w2 = G2(i,j,k)
                if (w2 < tiny(1._dp)) then
                    dVen0GdG2(i,j,k) = 0.0_dp
                else
                    dVen0GdG2(i,j,k) = -rloc**2/2.0_dp*Ven0G(i,j,k) &
                                       +4.0_dp*pi*Zion/w2**2*exp(-w2*rloc**2/2.0_dp) &
                                   + sqrt(8.0_dp*pi**3)*rloc**3*exp(-w2*rloc**2/2.0_dp) & 
                              *(C(2)*(-rloc**2) &
                               +C(3)*(-10.0_dp*rloc**2 + 2.0*w2*rloc**4) &
                               +C(4)*(-105.0_dp*rloc**2 + 42.0_dp*w2*rloc**4 -3.0_dp*w2**2*rloc**6) &
                               )
                endif
            enddo
            enddo
            enddo
    end subroutine
        
end module
