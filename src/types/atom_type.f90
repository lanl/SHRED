module atom_type
    use types, only : dp
    use shared_type, only: shared_struct
    use m_paw_an, only : paw_an_type
    use m_paw_ij, only : paw_ij_type
    use m_pawfgrtab,only : pawfgrtab_type
    use m_pawrhoij, only :pawrhoij_type
    use m_pawcprj, only : pawcprj_type
    use m_pawrad, only: pawrad_type
    use m_pawtab, only: pawtab_type

    implicit none
    
    private 

    public atom_struct

    type atom_Real_Local
        real(dp), pointer :: projector(:,:), deriv_projector(:,:,:)
        type(shared_struct) :: proj_shared, deriv_proj_shared
        integer, allocatable :: map(:,:)
        real(dp), allocatable :: Rs(:,:)
        integer :: ns
    end type

    type atom_PAW
        type(paw_an_type), pointer :: an
        type(paw_ij_type), pointer  :: ij
        type(pawfgrtab_type), pointer  :: fgrtab
        type(pawrhoij_type), pointer  :: rhoij
        type(pawcprj_type), pointer :: cprj(:)
        type(pawrad_type), pointer :: rad
        type(pawtab_type), pointer :: tab
        integer, pointer :: dimcprj, typat
        real(dp), pointer :: projector_ort(:,:), deriv_projector_ort(:,:,:)
        type(shared_struct) :: proj_shared_ort, deriv_proj_shared_ort
        real(dp), allocatable :: Ei(:), obar_i(:)
        real(dp), allocatable :: Dij(:,:,:), Uij(:,:), Sij(:,:)
        complex(dp), allocatable :: rho_ij(:,:,:), rho_ij_save(:,:,:)
        integer :: index
        real(dp), allocatable  :: Sinv_block(:,:),  Rs_fine(:,:)
        integer :: ns_fine
        integer, allocatable :: map_fine(:,:)
    end type
    
    type atom_struct
        real(dp) :: R(3)
        real(dp) :: P(3)
        real(dp) :: F_nn(3), F_loc(3), F_nonloc(3), F(3), F_comp(3),F_xc_core(3)
        integer :: element
        integer :: atom_of_element
        logical :: frozen_R, frozen_V, update_me, update_me_force, mine
        type(atom_Real_Local), allocatable :: RL
        real(dp), allocatable :: R_mix(:,:)
        type(atom_PAW), allocatable :: PAW
    end type
    
end module