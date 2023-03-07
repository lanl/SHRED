module Entropy
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: Calc_Entropy

    contains

    subroutine Calc_Entropy(sim)
        use simulation_type, only : simulation_struct
        use Thomas_Fermi, only: Thomas_Fermi_Entropy
        use Stochastic_Mod, only : Stoc_Entropy
        use odp_type, only : orbital_struct

        type(simulation_struct), intent(inout), target :: sim

        if(sim%system%orbital_free) then
            sim%entropy=Thomas_Fermi_Entropy(sim%density, sim%grids, sim%parallel, sim%system) 
        else

            sim%entropy=KS_Entropy_Deterministic(sim%orbitals, sim%parallel, sim%system)
            sim%entropy=sim%entropy+Stoc_Entropy(sim)
        endif
        if(sim%parallel%myid.eq.0) print *, '-T*Entropy(a.u.): ', sim%entropy
    end subroutine

    function KS_Entropy_Deterministic(psi, parallel, system) result (Entrop)
        use odp_type, only: orbital_struct
        use system_type, only: system_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task
        use odp_type, only : deterministic
        type(parallel_struct), intent(in) :: parallel
        type(system_struct), intent(in) :: system
        type(orbital_struct), intent(in) :: psi(:,:,:)

        real(dp):: Entrop, tmp1, tmp2
        integer :: i,j,k,s
        
        Entrop=0.0_dp
        !Note that i and s both relate account for spin,
        !but only one of the two should be able to be > 1 at a time
        do i=1, size(psi,3); do j=1, size(psi,2); do k=1, size(psi,1)
            if(psi(k,j,i)%type.eq.deterministic) then
            do s=1,psi(k,j,i)%n_spinor
                if(abs(psi(k,j,i)%filter(s)).lt.1.0E-12) cycle
                if(abs(psi(k,j,i)%occ(s)/psi(k,j,i)%degeneracy/psi(k,j,i)%filter(s)).lt.1.0E-12 .or. &
                   abs(1.0_dp-psi(k,j,i)%occ(s)/psi(k,j,i)%degeneracy/psi(k,j,i)%filter(s)).lt.1.0E-12 ) cycle   
                tmp1=psi(k,j,i)%occ(s)/psi(k,j,i)%degeneracy*psi(k,j,i)%filter(s)
                tmp2=1.0_dp-tmp1
                Entrop = Entrop + psi(k,j,i)%weight(s)*psi(k,j,i)%degeneracy*psi(k,j,i)%filter(s)*( &
                    tmp1*log(tmp1) + tmp2*log(tmp2))
            end do
            endif
        end do;end do;end do

        call parallel_task('sum', Entrop, parallel, 'space')

        Entrop=system%temperature*Entrop
    end function
    

end module