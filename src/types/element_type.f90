module element_type
    use types, only : dp
    use shared_type, only: shared_struct
    use m_pawrad,    only : pawrad_type
    use m_pawtab,    only : pawtab_type 
    use m_pawang,    only : pawang_type

    use m_libpaw_defs, only: fnlen
    use m_pawxmlps, only: paw_setup_t

    implicit none
    
    private 

    public element_struct, Real_local, Euler_Spline, Reciprocal, None, HGH, PAW, local_only

    integer :: local_only=0
    integer :: Real_local=1
    integer :: Euler_Spline=2
    integer :: Reciprocal=3


    integer :: None=1
    integer :: HGH=2
    integer :: PAW=3
  

    type element_HGH
        real(dp) ::  epsatm
        real(dp) :: rloc, C(4)
        integer  :: nloc, nlh(4), nlk(4), n_proj
        real(dp) :: H(3,3,4), K(3,3,4), U_HE(3,3,4), EKB(3,4), U_HK(3,3,4), KKB(3,4), rl(4)
        integer, allocatable :: l(:), m(:), i(:)
        character, allocatable :: h_or_k(:)
    end type

    type element_None
    end type

    type element_PAW
        type(pawang_type), pointer :: ang
        type(pawtab_type), pointer :: tab
        type(pawrad_type), pointer :: rad
        type(paw_setup_t), pointer :: setup


        logical :: xml_file=.true.
        character(len=fnlen) :: filepsp  ! name of the psp file
        integer  :: icoulomb=0
        integer  :: ixc, nspden, lnmax, spnorb
        integer  :: xclevel
        integer  :: pawxcdev=1,usewvl=0,usexcnhat=-1
        real(dp) :: xc_denpos=1.0E-14
        integer :: mqgrid_ff , mqgrid_vl
        real(dp) :: epsatm,xcccrc
        real(dp), allocatable :: qgrid_ff(:),qgrid_vl(:)
        real(dp), allocatable :: ffspl(:,:,:), ffspl_Ort(:,:,:)
        real(dp), allocatable :: vlspl(:,:)
        real(dp), allocatable :: sij(:,:), eijkl(:,:,:,:), dij0(:,:)
        real(dp), allocatable :: tproj_RO(:,:) , o_i(:)
        integer :: index, n_proj

        integer :: max_sphere_fine(3)
        real(dp), allocatable :: fine_grid(:,:)
        real(dp) :: radius
    end type

    type element_Real_Local
        integer :: max_sphere(3),  max_sphere_dense(3), OH_factor
        !Ono-Hirose Transform splines
        real(dp), pointer :: Beta_x(:,:),Beta_y(:,:),Beta_z(:,:)
        type(shared_struct) :: beta_x_shared, beta_y_shared, beta_z_shared
        real(dp), allocatable :: coarse_grid(:,:), dense_grid(:,:)
        real(dp) :: radius, dense_dR(3)
    end type

    type element_Reciprocal_Non_local
        real(dp), allocatable :: projector_G(:,:,:,:)
        real(dp), allocatable :: projector_G_ort(:,:,:,:)
    end type

    type element_struct
        character(len=256) :: PP_file
        integer  :: PP_type
        integer  :: PA_type
        integer  :: n_atoms_of_element
        integer :: my_natom, my_atom_min, my_atom_max
        integer  :: n_proj
        integer, allocatable :: counts(:), displs(:)
        real(dp) :: M
        real(dp) :: Znuc, Zion, vcore
        type(element_None)  :: None
        type(element_HGH), allocatable  :: HGH
        type(element_PAW), allocatable  :: PAW
        type(element_Real_Local), allocatable  :: RL
        type(element_Reciprocal_Non_local), allocatable  :: GNL(:)
        real(dp), allocatable :: Ven0G(:,:,:), dVen0GdG2(:,:,:), Den0G(:,:,:), Den1G(:,:,:), dDen0GdG2(:,:,:)
        real(dp) :: mass_tot, isokin_ec, isokin_a, isokin_b, &
                   isokin_s, isokin_sdot
        logical :: all_frozen_R, all_frozen_V
        logical :: all_thawed_R, all_thawed_V
        logical :: recip_apply
    end type

end module
