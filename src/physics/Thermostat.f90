module Thermostat_mod
    use types, only:dp
    use constants, only: pi
    implicit none  
    public IsoKinetic_Thermostat
    contains
       
       subroutine IsoKinetic_Thermostat(atoms, atoms_P_out, elements, t_bostep)  
        use atom_type, only : atom_struct 
        use system_type, only : system_struct
        use element_type, only : element_struct  
  
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout), target :: elements(:)
        real(dp), intent(inout) :: atoms_P_out(:,:)
        real(dp), intent(in) :: t_bostep
        type(element_struct), pointer :: element
        integer :: i, e

        elements(:)%isokin_ec=0.0_dp
        elements(:)%isokin_a=0.0_dp
        elements(:)%isokin_b=0.0_dp
        
        do i=1,size(atoms)
            if(atoms(i)%frozen_R) cycle
            if(atoms(i)%frozen_V) cycle
            e=atoms(i)%element
            element=>elements(e)
            element%isokin_ec=element%isokin_ec+sum(atoms(i)%P(:)*atoms(i)%P(:))/element%M
            element%isokin_a=element%isokin_a+sum(atoms(i)%F(:)*atoms(i)%P(:))/element%M
            element%isokin_b=element%isokin_b+sum(atoms(i)%F(:)*atoms(i)%F(:))/element%M
        enddo
        do e=1,size(elements)
            element=>elements(e)
            if(element%all_frozen_R) cycle
            if(element%all_frozen_V) cycle
            if(element%isokin_ec.gt.tiny(1.0_dp)) then
                element%isokin_a=element%isokin_a/element%isokin_ec
                element%isokin_b=element%isokin_b/element%isokin_ec
                element%isokin_s=element%isokin_a/element%isokin_b*(cosh(sqrt(element%isokin_b)*t_bostep)-1.0_dp) &
                                +sinh(sqrt(element%isokin_b)*t_bostep)/sqrt(element%isokin_b)
                element%isokin_sdot=element%isokin_a/sqrt(element%isokin_b)* &
                                    sinh(sqrt(element%isokin_b)*t_bostep)+cosh(sqrt(element%isokin_b)*t_bostep) 
            else
                element%isokin_sdot=1.0_dp
                element%isokin_s=1.0_dp
            endif
        enddo
        do i=1,size(atoms)
            if(atoms(i)%frozen_V) cycle
            e=atoms(i)%element
            element=>elements(e)
            atoms_P_out(:,i)=(atoms(i)%P(:)+element%isokin_s*atoms(i)%F(:))/element%isokin_sdot
        enddo

       end subroutine

    end module