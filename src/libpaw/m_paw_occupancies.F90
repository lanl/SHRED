!!****m* m_paw_occupancies/m_paw_occupancies
!! NAME
!!  m_paw_occupancies
!!
!! FUNCTION
!!  This module contains routines related to the computation of PAW on-site occupancies (rhoij).
!!
!! COPYRIGHT
!! Copyright (C) 2018-2021 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "libpaw.h"

MODULE m_paw_occupancies

  USE_DEFS
  USE_MSG_HANDLING
  USE_MPI_WRAPPERS
  USE_MEMORY_PROFILING

 use m_libpaw_defs
 
 use m_pawtab,     only : pawtab_type
 use m_pawrhoij,   only : pawrhoij_type,pawrhoij_init_unpacked,pawrhoij_mpisum_unpacked, &
&                         pawrhoij_alloc,pawrhoij_free,pawrhoij_inquire_dim
 use m_pawcprj,    only : pawcprj_type,pawcprj_alloc,pawcprj_get, &
&                         pawcprj_gather_spin, pawcprj_free, pawcprj_mpi_send, &
&                         pawcprj_mpi_recv, pawcprj_copy, pawcprj_unpack, pawcprj_pack
 use m_paw_io,     only : pawio_print_ij
 use m_paral_atom, only : get_my_atmtab,free_my_atmtab

 implicit none

 private

!public procedures.
 public :: initrhoij   ! Initialize the PAW occupancies rhoij from atomic data

CONTAINS  !========================================================================================

!----------------------------------------------------------------------

!!****f* m_paw_occupancies/initrhoij
!! NAME
!! initrhoij
!!
!! FUNCTION
!! Initialize PAW rhoij occupancies (in packed storage)
!! from atomic ones
!!
!! INPUTS
!!  cpxocc=1 if rhoij are real, 2 if they are complex
!!  lexexch(ntypat)=l on which local exact-exchange is applied for a given type of atom
!!  lpawu(ntypat)=l on which U is applied for a given type of atom (PAW+U)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of atom types
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated in PAW augmentation regions
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!                                     (containing initial rhoij)
!!  qphase=2 if rhoij have a exp(iqR) phase, 1 if not (typical use: response function at q<>0)
!!  spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!!  typat(natom)=type of each atom
!!  === Optional arguments
!!    ngrhoij=number of gradients to be allocated (OPTIONAL, default=0)
!!    nlmnmix=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!!    use_rhoij_=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!!    use_rhoijres=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)

!!
!! OUTPUT
!!  pawrhoij(natom) <type(pawrhoij_type)>=rhoij quantities for each atom
!!                                        in packed storage
!!
!! PARENTS
!!      m_gstate,m_nonlinear,m_positron,m_respfn_driver
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawrhoij_alloc,pawrhoij_inquire_dim
!!
!! SOURCE

subroutine initrhoij(cpxocc,lexexch,lpawu,my_natom,natom,nspden,nspinor,nsppol,&
&                    ntypat,pawrhoij,pawspnorb,pawtab,qphase,spinat,typat,&
&                    ngrhoij,nlmnmix,use_rhoij_,use_rhoijres,& ! optional arguments
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cpxocc,my_natom,natom,nspden,nspinor,nsppol,ntypat,pawspnorb,qphase
 integer,intent(in),optional :: comm_atom,ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
 character(len=500) :: message
!arrays
 integer,intent(in) :: lexexch(ntypat),lpawu(ntypat)
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: spinat(3,natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!Arrays
!scalars
 integer :: cplex_rhoij,iatom,iatom_rhoij,ilmn,ispden,itypat,j0lmn,jl,jlmn,jspden
 integer :: klmn,klmn1,ln,lnspinat0,my_comm_atom
 integer :: ngrhoij0,nlmnmix0,nselect,nselect1,nspden_rhoij,qphase_rhoij
 integer :: use_rhoij_0,use_rhoijres0
 real(dp) :: ratio,ro,roshift,zratio,zz
 logical :: my_atmtab_allocated,paral_atom,spinat_zero,test_exexch,test_pawu,test_lnspinat
!arrays
 integer,pointer :: my_atmtab(:),lnspinat(:)
 real(dp),allocatable :: occ(:)
!************************************************************************

!PAW+U and local exact-exchange restriction
 do itypat=1,ntypat
   if (lpawu(itypat)/=lexexch(itypat).and. lpawu(itypat)/=-1.and.lexexch(itypat)/=-1) then
     message = ' lpawu must be equal to lexexch !'
     LIBPAW_ERROR(message)
   end if
 end do

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                          nspden=nspden,spnorb=pawspnorb,cpxocc=cpxocc,cplex=qphase)

 ratio=one;if (nspden_rhoij==2) ratio=half
 spinat_zero=all(abs(spinat(:,:))<tol10)

 if (my_natom>0) then
   ngrhoij0=0;if (present(ngrhoij)) ngrhoij0=ngrhoij
   nlmnmix0=0;if (present(nlmnmix)) nlmnmix0=nlmnmix
   use_rhoij_0=0;if (present(use_rhoij_)) use_rhoij_0=use_rhoij_
   use_rhoijres0=0;if (present(use_rhoijres)) use_rhoijres0=use_rhoijres
   if (paral_atom) then
     call pawrhoij_alloc(pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,nsppol,typat,&
&     ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0,&
&     qphase=qphase_rhoij,pawtab=pawtab,comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
   else
     call pawrhoij_alloc(pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,nsppol,typat,qphase=qphase_rhoij,&
&     pawtab=pawtab,ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0)
   end if
 end if

 do iatom_rhoij=1,my_natom
   iatom=iatom_rhoij;if (paral_atom) iatom=my_atmtab(iatom_rhoij)
   itypat=typat(iatom)
   nselect=0
   LIBPAW_ALLOCATE(lnspinat,(pawtab(itypat)%basis_size))
   lnspinat=-1
! Determine occupancies of each orbital
   if (nspden_rhoij==2) then
     LIBPAW_ALLOCATE(occ,(pawtab(itypat)%basis_size))
     occ=zero
     do jlmn=1,pawtab(itypat)%lmn_size
       ln=pawtab(itypat)%indlmn(5,jlmn)
       klmn=jlmn*(jlmn+1)/2
       occ(ln)=occ(ln)+pawtab(itypat)%rhoij0(klmn)
     end do
     do ln=1,pawtab(itypat)%basis_size
       if(pawtab(itypat)%orbitals(ln)==0.and.occ(ln)==1) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==1.and.(occ(ln)>=1.and.occ(ln)<=5)) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==2.and.(occ(ln)>=1.and.occ(ln)<=9)) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==3.and.(occ(ln)>=1.and.occ(ln)<=13)) lnspinat(ln)=ln
     end do
     LIBPAW_DEALLOCATE(occ)
   end if
   lnspinat0=maxval(lnspinat)
   lnspinat0=-1

!  Determine Z (trace of rhoij0 or part of it)
   zz=zero
   do jlmn=1,pawtab(itypat)%lmn_size
     jl=pawtab(itypat)%indlmn(1,jlmn)
     ln=pawtab(itypat)%indlmn(5,jlmn)
     j0lmn=jlmn*(jlmn-1)/2
     test_lnspinat=(lnspinat0==-1.or.lnspinat(ln)==ln)
     test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
     test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       if ((ilmn==jlmn).and.test_pawu.and.test_exexch.and.test_lnspinat) &
&       zz=zz+pawtab(itypat)%rhoij0(klmn)
     end do
   end do

!  Compute rhoij from tabulated value and magnetization
   do ispden=1,nspden_rhoij

     zratio=zero
     roshift=one
     ratio=one
     if (nspden_rhoij==2) then
       ratio=half
       if ((spinat(3,iatom)>zero.and.ispden==1).or.&
&       (spinat(3,iatom)<zero.and.ispden==2)) then
         if(abs(zz)>tol12)then
           zratio=two*abs(spinat(3,iatom))/zz
         else
           zratio=zero
         end if
       end if
     else if (nspden_rhoij==4.and.ispden>=2) then
       roshift=zero
       if(abs(zz)>tol12)then
         zratio=spinat(ispden-1,iatom)/zz
       else
         zratio=zero
       end if
     end if

     nselect=0;nselect1=1-cpxocc
     do jlmn=1,pawtab(itypat)%lmn_size
       jl=pawtab(itypat)%indlmn(1,jlmn)
       ln=pawtab(itypat)%indlmn(5,jlmn)
       j0lmn=jlmn*(jlmn-1)/2
       test_lnspinat=(lnspinat0==-1.or.lnspinat(ln)==ln)
       test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
       test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ro=pawtab(itypat)%rhoij0(klmn)
         if ((ilmn==jlmn).and.test_pawu.and.test_exexch.and.test_lnspinat) then
           ro=ro*ratio*(roshift+zratio)
         else
           ro=ro*ratio*roshift
         end if

         klmn1=cpxocc*(klmn-1)+1
         if (abs(ro)>tol10) then
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=ro
         else
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=zero
         end if

         if (ispden==nspden_rhoij) then
           if (any(abs(pawrhoij(iatom_rhoij)%rhoijp(klmn1,:))>tol10)) then
             nselect=nselect+1;nselect1=nselect1+cpxocc
             pawrhoij(iatom_rhoij)%rhoijselect(nselect)=klmn
             do jspden=1,nspden_rhoij
               pawrhoij(iatom_rhoij)%rhoijp(nselect1,jspden)=pawrhoij(iatom_rhoij)%rhoijp(klmn1,jspden)
             end do
           end if
         end if

       end do
     end do

   end do
   pawrhoij(iatom_rhoij)%nrhoijsel=nselect
   if (nselect<pawrhoij(iatom_rhoij)%lmn2_size) &
&    pawrhoij(iatom_rhoij)%rhoijselect(nselect+1:pawrhoij(iatom_rhoij)%lmn2_size)=0

!  Non-collinear magnetism: avoid zero magnetization, because it produces numerical instabilities
!    Add a small real to the magnetization ; not yet activated => must be tested.
!   if (pawrhoij(iatom_rhoij)%nspden==4.and.spinat_zero) then
!     pawrhoij(iatom_rhoij)%rhoijp(:,4)=pawrhoij(iatom_rhoij)%rhoijp(:,4)+tol10
!   end if
   LIBPAW_DEALLOCATE(lnspinat)
 end do ! iatom_rhoij

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)


end subroutine initrhoij
!!***

!----------------------------------------------------------------------

END MODULE m_paw_occupancies
!!***
