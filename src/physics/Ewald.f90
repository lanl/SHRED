module Ewald_mod
    use types, only: dp
    use constants, only: pi
    implicit none  
    private outer_prod, erfc2, ewald_in
    public Ewald_Calc, init_ewald, free_ewald
    contains
       
        subroutine outer_prod(ab,a,b)
          real(dp), dimension(3), intent(in) :: a,b
          real(dp), dimension(3,3), intent(out) :: ab

               ab(1:3,1)= a(1:3)*b(1)
               ab(1:3,2)= a(1:3)*b(2)
               ab(1:3,3)= a(1:3)*b(3)

        end subroutine outer_prod
     
        function erfc2(arg)

          real(dp) :: t,tmp,arg,erfc2
    
          t=1.0_dp/(1.0_dp+.3275911_dp*arg)
          tmp=t*(.254829592_dp+t*(-.284496736_dp+t*(1.421413741_dp+t*(-1.453152027_dp+&
               t*1.061405429_dp))))
          erfc2=tmp*exp(-arg**2)
    
        end function erfc2

       subroutine init_ewald(ewald, Box_Length, atoms)
          use ewald_type, only : ewald_struct, ntable
          use atom_type, only : atom_struct
          type(atom_struct), intent(in) :: atoms(:)
          type(ewald_struct), intent(inout) :: ewald
          real(dp) :: Box_Length(3), R

          integer :: N, i

          N=size(atoms)

          ewald%alpha=7.0_dp/minval(Box_Length)  
          ewald%kmax=int(1.0_dp + 3.7_dp*ewald%alpha*Box_Length/pi) 
          ewald%preuk=8.0_dp*pi/product(Box_Length) 
          ewald%unitk=2.0_dp*pi/Box_Length 
          ewald%kmaxsq=maxval((ewald%unitk*ewald%kmax)**2.0_dp)
          ewald%delrsq=sum((Box_Length/2.0_dp)**2.0_dp)/((ntable-1)*1.0_dp)      

          ! delta function 3D
          ewald%delta(:,:)=0.0_dp
          ewald%delta(1,1)=1.0_dp
          ewald%delta(2,2)=1.0_dp
          ewald%delta(3,3)=1.0_dp

          allocate(ewald%urtable(ntable))
          allocate(ewald%frtable(ntable))
          !Allocate working space
          allocate(ewald%E_field(3,N))
          allocate(ewald%cosat(N))
          allocate(ewald%sinat(N))
          allocate(ewald%v(N))
          allocate(ewald%ctmp(N))
          allocate(ewald%stmp(N))

          allocate(ewald%coskr(0:maxval(ewald%kmax),3,N))
          allocate(ewald%sinkr(0:maxval(ewald%kmax),3,N))
          allocate(ewald%ukarray(0:maxval(ewald%kmax),3))
      
          do i=0,maxval(ewald%kmax)
               ewald%ukarray(i,:)=exp(-(i*ewald%unitk/(2.0_dp*ewald%alpha))**2.0_dp) ! \exp(-\frac{k^{2}}{4 \alpha^{2}})
          enddo
 
          do i=1,ntable
               R=sqrt(i*ewald%delrsq)
               ewald%urtable(i)=erfc2(ewald%alpha*R)/R ! urtable(r)=erfc(ar)/r
               ewald%frtable(i)=(ewald%urtable(i)+ &
                    2.0_dp/(sqrt(pi))*ewald%alpha*exp(-(ewald%alpha*R)**2.0_dp))/R**2.0_dp ! frtable(r)=-1/r d(urtable(r))/dr
          enddo
          
       end subroutine
       
       subroutine free_ewald(ewald)
          use ewald_type, only : ewald_struct
          type(ewald_struct), intent(inout) :: ewald

          deallocate(ewald%cosat)
          deallocate(ewald%sinat)
          deallocate(ewald%coskr)
          deallocate(ewald%sinkr)
          deallocate(ewald%v)
          deallocate(ewald%ctmp)
          deallocate(ewald%stmp)
          deallocate(ewald%urtable)
          deallocate(ewald%frtable)
          deallocate(ewald%ukarray)
          deallocate(ewald%E_field)
       end subroutine

       subroutine ewald_in(ewald, atoms, Zion, Box_Length, calc_force, calc_stress)
          use ewald_type, only : ewald_struct, ntable
          use atom_type, only : atom_struct
          type(atom_struct), intent(in) :: atoms(:)
          type(ewald_struct), intent(inout) :: ewald
          real(dp), intent(in) :: Zion(:)
          real(dp), intent(in) :: Box_Length(3)
          logical, intent(in) :: calc_stress, calc_force

          ! Local variables.
          real(dp) :: r_tensor(3,3),k_tensor(3,3), dr(3)
          real(dp) :: k_vector(3)
          real(dp) :: uk, fk, k_square
          real(dp) :: Zdotcos, Zdotsin
          integer :: dir, i, j, kx, ky, kz, itab, perm
          real(dp) :: dr2, itab_real, vtab, etab
             
          !Energy
          ewald%v(:)=-pi/(ewald%alpha**2.0_dp*product(Box_Length))* &
                         sum(Zion)-2.0_dp*ewald%alpha/sqrt(pi)*Zion(:)
          !Force
          ewald%E_field(:,:)=0.0_dp
          !Stress
          ewald%stress(:,:)=-pi/(2.0_dp*ewald%alpha**2.0_dp*product(Box_Length)) &
                              *sum(Zion)**2.0_dp*ewald%delta(:,:)
          !kx/y/z=0
          ewald%coskr(0,:,:)=1.0_dp
          ewald%sinkr(0,:,:)=0.0_dp
          !kx/y/z=1
          do dir=1,3
               ewald%coskr(1,dir,:)=cos(ewald%unitk(dir)*atoms(:)%R(dir))
               ewald%sinkr(1,dir,:)=sin(ewald%unitk(dir)*atoms(:)%R(dir))
          enddo

          do kx=2,maxval(ewald%kmax)
             do dir=1,3
                    ewald%coskr(kx,dir,:)=ewald%coskr(kx-1,dir,:)*ewald%coskr(1,dir,:)-ewald%sinkr(kx-1,dir,:)*ewald%sinkr(1,dir,:) 
                    ewald%sinkr(kx,dir,:)=ewald%sinkr(kx-1,dir,:)*ewald%coskr(1,dir,:)+ewald%coskr(kx-1,dir,:)*ewald%sinkr(1,dir,:)
             enddo
          enddo

          do kx=0, ewald%kmax(1)
             k_vector(1)=ewald%unitk(1)*kx
             do ky=0, ewald%kmax(2)
                k_vector(2)=ewald%unitk(2)*ky
                do kz=0, ewald%kmax(3)      
                   k_vector(3)=ewald%unitk(3)*kz
                   k_square=sum(k_vector**2.0_dp)
                   if(k_square.lt.tiny(1.0_dp)) cycle !skip the 0, 0, 0 
                   if(k_square<ewald%kmaxsq) then
                         ! Permutations
                         ! k->-k (sin(k)->-sin(k))
                         ! 1=(kx,ky,kz)
                         ! 2=(kx,-ky,kz)
                         ! 3=(kx,-ky,-kz)
                         ! 4=(kx,ky,-kz)
                         ! 5= reset and back to 1
                         do perm=1,5
                              !if 1 proceed as is
                              if(perm.gt.1) then
                                   !If 2 components are zero, dont need symmetric negative value
                                   if(ky.eq.0.and.kx.eq.0) cycle !No - contribution for 0
                                   if(kz.eq.0.and.kx.eq.0) cycle !No - contribution for 0
                                   if(ky.eq.0.and.kz.eq.0) cycle !No - contribution for 0
                              endif
                              if(perm.eq.2) then
                                   ewald%sinkr(ky,2,:)=-ewald%sinkr(ky,2,:)
                                   k_vector(2)=-k_vector(2)
                                   if(ky.eq.0) cycle !No - contribution for 0
                              endif
                              if(perm.gt.2.and.kx.eq.0) then
                                   !If kx component is zero, dont need symmetric negative value
                                   ewald%sinkr(ky,2,:)=-ewald%sinkr(ky,2,:)
                                   k_vector(2)=-k_vector(2)
                                   cycle
                              endif
                                   
                              if(perm.eq.3) then
                                   ewald%sinkr(kz,3,:)=-ewald%sinkr(kz,3,:)
                                   k_vector(3)=-k_vector(3)
                                   if(ky.eq.0) cycle
                                   if(kz.eq.0) cycle !No - contribution for 0
                              endif
                              if(perm.eq.4) then
                                   ewald%sinkr(ky,2,:)=-ewald%sinkr(ky,2,:)
                                   k_vector(2)=-k_vector(2)
                                   if(kz.eq.0) cycle !No - contribution for 0
                              endif
                              if(perm.eq.5) then
                                   ewald%sinkr(kz,3,:)=-ewald%sinkr(kz,3,:)
                                   k_vector(3)=-k_vector(3)
                                   cycle
                              endif
                              ewald%ctmp(:)= ewald%coskr(ky,2,:)*ewald%coskr(kz,3,:) - ewald%sinkr(ky,2,:)*ewald%sinkr(kz,3,:)
                              ewald%stmp(:)= ewald%sinkr(ky,2,:)*ewald%coskr(kz,3,:) + ewald%coskr(ky,2,:)*ewald%sinkr(kz,3,:)

                              ewald%cosat(:)= ewald%coskr(kx,1,:)*ewald%ctmp - ewald%sinkr(kx,1,:)*ewald%stmp     
                              ewald%sinat(:)= ewald%sinkr(kx,1,:)*ewald%ctmp + ewald%coskr(kx,1,:)*ewald%stmp

                              Zdotcos=sum(Zion(:)*ewald%cosat(:))
                              Zdotsin=sum(Zion(:)*ewald%sinat(:))

                              uk=ewald%preuk*ewald%ukarray(kx,1)*ewald%ukarray(ky,2)*ewald%ukarray(kz,3)/k_square     

                              ewald%v(:)=ewald%v(:) + &
                                        uk*(Zdotcos*ewald%cosat(:)+Zdotsin*ewald%sinat(:))

                              !print *, kx, ky, kz, perm, uk*(Zdotcos*ewald%cosat(:)+Zdotsin*ewald%sinat(:))

                              if(calc_force) then
                                   fk=k_vector(1)*uk
                                   ewald%E_field(1,:)=ewald%E_field(1,:) + &
                                        fk*(Zdotcos*ewald%sinat(:)-Zdotsin*ewald%cosat(:)) 
                                   fk=k_vector(2)*uk
                                   ewald%E_field(2,:)=ewald%E_field(2,:) + &
                                        fk*(Zdotcos*ewald%sinat(:)-Zdotsin*ewald%cosat(:))
                                   fk=k_vector(3)*uk
                                   ewald%E_field(3,:)=ewald%E_field(3,:)+&
                                        fk*(Zdotcos*ewald%sinat(:)-Zdotsin*ewald%cosat(:))
                              endif
                              if(calc_stress) then
                                   call outer_prod(k_tensor(:,:),k_vector,k_vector)
                                   
                                   ewald%stress(:,:)=ewald%stress(:,:)+&
                                        (ewald%delta(:,:)-(0.5_dp/(ewald%alpha**2.0_dp)+2.0_dp/k_square)*k_tensor(:,:)) * &
                                        uk*sum(0.5_dp*Zion(:)*(Zdotcos*ewald%cosat(:)+Zdotsin*ewald%sinat(:)))

                              endif
                         enddo ! sign permutations
                    endif
                enddo
             enddo
          enddo

          ! Real space part
          do i=1,size(atoms)
             do j=i+1,size(atoms)
                dr(:)=atoms(i)%R(:)-atoms(j)%R(:)
                dr(:)=dr(:)-Box_Length(:)*idnint(dr(:)/Box_Length(:))
                dr2=sum(dr(:)**2)
    
                itab_real=dr2/ewald%delrsq
                ! Make sure itab falls on the table
                itab=idnint(itab_real)
                itab=min(itab,ntable-1)
                itab=max(itab,1)

                !table extrapolation/interpolation
                itab_real=itab_real-1.0_dp*itab

                !Energy
                vtab=itab_real*ewald%urtable(itab+1)+(1.0_dp-itab_real)*ewald%urtable(itab) 
                ewald%v(i)=ewald%v(i)+Zion(j)*vtab 
                ewald%v(j)=ewald%v(j)+Zion(i)*vtab
                !Force
                if(calc_force.or.calc_stress)  &
                    etab=itab_real*ewald%frtable(itab+1)+(1.0_dp-itab_real)*ewald%frtable(itab)
                if(calc_force) then
                    print *, 'Recip space part:', i, j,  ewald%E_field(:,i), ewald%E_field(:,j)

                    ewald%E_field(:,i)=ewald%E_field(:,i)+Zion(j)*dr(:)*etab
                    ewald%E_field(:,j)=ewald%E_field(:,j)-Zion(i)*dr(:)*etab

                    print *, 'Real space part:', i, j,  Zion(j)*dr(:)*etab, -Zion(i)*dr(:)*etab
                endif
                !Stress
                if(calc_stress) then
                    r_tensor(:,:)=0.0_dp
                    r_tensor(:,1)=dr(:)*dr(1)
                    r_tensor(:,2)=dr(:)*dr(2)
                    r_tensor(:,3)=dr(:)*dr(3)
                    ewald%stress(:,:)=ewald%stress(:,:)+Zion(i)*Zion(j)*etab*r_tensor(:,:)
                endif
             enddo
          enddo
    
          ewald%energy=0.5_dp*sum(Zion(:)*ewald%v(:))
          ewald%stress(2,1)=ewald%stress(1,2)
          ewald%stress(3,1)=ewald%stress(1,3)
          ewald%stress(3,2)=ewald%stress(2,3)
          ewald%stress=-ewald%stress

          return
        end subroutine ewald_in
    
        subroutine Ewald_Calc(atoms, elements, grids, ewald, energies, calc_force, stress)
            use grids_type, only : grid_struct
            use simulation_type, only : stress_struct, energies_struct
            use element_type, only : element_struct
            use atom_type, only : atom_struct
            use ewald_type, only : ewald_struct

            type(atom_struct), intent(inout) :: atoms(:)
            type(element_struct), intent(inout) :: elements(:)
            type(grid_struct), intent(inout) :: grids(:)
            type(stress_struct), intent(inout), optional :: stress
            type(energies_struct), intent(inout) :: energies
            type(ewald_struct), intent(inout) :: ewald
            logical, intent(in) :: calc_force
            integer :: i
            real(dp) :: gmet(3,3), rmet(3,3), rprimd(3,3), grewtn(3,size(atoms)), xred(3,size(atoms)), ab_stress(6)
            
            if(.true.) then
               gmet=0
               gmet(1,1)=grids(1)%Box_Length(1)
               gmet(2,2)=grids(1)%Box_Length(2)
               gmet(3,3)=grids(1)%Box_Length(3)
               gmet=gmet/product(grids(1)%Box_Length(:))

               rmet=0
               rmet(1,1)=grids(1)%Box_Length(1)**2
               rmet(2,2)=grids(1)%Box_Length(2)**2
               rmet(3,3)=grids(1)%Box_Length(3)**2

               xred(1,:)=atoms(:)%R(1)/grids(1)%Box_Length(1)
               xred(2,:)=atoms(:)%R(2)/grids(1)%Box_Length(2)
               xred(3,:)=atoms(:)%R(3)/grids(1)%Box_Length(3)

               call abinit_ewald(ewald%energy,gmet,grewtn,size(atoms),size(elements), &
                    rmet, atoms(:)%element, product(grids(1)%Box_Length(:)), xred, elements(:)%Zion)
               
               energies%nuclear_pot=ewald%energy

                    if(calc_force) then
                         do i=1, size(atoms)
                              atoms(i)%F_nn(:)=-grewtn(:,i)/grids(1)%Box_Length(:)
                         enddo
                    endif
               if(present(stress)) then
                    rprimd=0
                    rprimd(1,1)=grids(1)%Box_Length(1)
                    rprimd(2,2)=grids(1)%Box_Length(2)
                    rprimd(3,3)=grids(1)%Box_Length(3)
                    call abinit_ewald2(gmet,size(atoms),size(elements),rmet,rprimd,ab_stress, &
                         atoms(:)%element,product(grids(1)%Box_Length),xred,elements(:)%Zion)
                    stress%nuclear_pot(1,1)=ab_stress(1)
                    stress%nuclear_pot(2,2)=ab_stress(2)
                    stress%nuclear_pot(3,3)=ab_stress(3)
                    stress%nuclear_pot(3,2)=ab_stress(4)
                    stress%nuclear_pot(3,1)=ab_stress(5)
                    stress%nuclear_pot(2,1)=ab_stress(6)
                    stress%nuclear_pot(2,3)=stress%nuclear_pot(3,2)
                    stress%nuclear_pot(1,3)=stress%nuclear_pot(3,1)
                    stress%nuclear_pot(1,2)=stress%nuclear_pot(2,1)
                    stress%nuclear_pot=stress%nuclear_pot*product(grids(1)%Box_Length)
               endif

            else
               if(present(stress)) then
                call ewald_in(ewald, atoms, &
                     elements(atoms(:)%element)%Zion,grids(1)%Box_Length, calc_force, .true.)
                     stress%nuclear_pot=ewald%stress
               else
                 call ewald_in(ewald, atoms, &
                     elements(atoms(:)%element)%Zion,grids(1)%Box_Length, calc_force, .false.)
               endif

                energies%nuclear_pot=ewald%energy
               if(calc_force) then
                do i=1, size(atoms)
                         print *, "Total", i, elements(atoms(i)%element)%Zion, ewald%E_field(:,i),  &
                              elements(atoms(i)%element)%Zion*ewald%E_field(:,i)
                        atoms(i)%F_nn(:)=elements(atoms(i)%element)%Zion*ewald%E_field(:,i)
                enddo 
               endif
            endif
          end subroutine   

!!****f* m_ewald/ewald
!!
!! NAME
!! ewald
!!
!! FUNCTION
!! Compute Ewald energy and derivatives with respect to dimensionless
!! reduced atom coordinates xred.
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in unit cell
!! ntypat=numbe of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! eew=final ewald energy in hartrees
!! grewtn(3,natom)=grads of eew wrt xred(3,natom), hartrees.
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!      matr3inv,spline
!!
!! SOURCE

subroutine abinit_ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)

      implicit none
     
     !Arguments ------------------------------------
     !scalars
      integer,intent(in) :: natom,ntypat
      real(dp),intent(in) :: ucvol
      real(dp),intent(out) :: eew
     !arrays
      integer,intent(in) :: typat(natom)
      real(dp),intent(in) :: gmet(3,3),rmet(3,3),xred(3,natom),zion(ntypat)
      real(dp),intent(out) :: grewtn(3,natom)
     
     !Local variables-------------------------------
     !scalars
      integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
      real(dp) :: arg,c1i,ch,chsq,derfc_arg,direct,drdta1,drdta2,drdta3,eta,fac
      real(dp) :: fraca1,fraca2,fraca3,fracb1,fracb2,fracb3,gsq,gsum,phi,phr,r1
      real(dp) :: minexparg
      real(dp) :: r1a1d,r2,r2a2d,r3,r3a3d,recip,reta,rmagn,rsq,sumg,summi,summr,sumr
      real(dp) :: t1,term
      !character(len=500) :: message
     
     ! *************************************************************************
     
     
     !This is the minimum argument of an exponential, with some safety
      minexparg=log(tiny(0._dp))+5.0_dp
     
     !Add up total charge and sum of $charge^2$ in cell
     
      chsq=0._dp
      ch=0._dp
      do ia=1,natom
        ch=ch+zion(typat(ia))
        chsq=chsq+zion(typat(ia))**2
      end do
     
     !Compute eta, the Ewald summation convergence parameter,
     !for approximately optimized summations:
      direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
     & rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
      recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
     & gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
     !A bias is introduced, because G-space summation scales
     !better than r space summation ! Note : debugging is the most
     !easier at fixed eta.
      eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)
     
     !Conduct reciprocal space summations
      fac=pi**2/eta
      gsum=0._dp
      grewtn(:,:)=0.0_dp
     
     !Sum over G space, done shell after shell until all
     !contributions are too small.
      ng=0
      do
        ng=ng+1
        newg=0
     !   if (ng > 20 .and. mod(ng,10)==0) then
     !      write (message,'(3a,I10)') "Very large box of G neighbors in ewald: you probably do not want to do this.", ch10,&
     !&       " If you have a metal consider setting dipdip 0.  ng = ", ng
     !      MSG_WARNING(message)
     !   end if
     
        do ig3=-ng,ng
          do ig2=-ng,ng
            do ig1=-ng,ng
     
     !        Exclude shells previously summed over
              if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng .or. ng==1 ) then
     
     !          gsq is G dot G = |G|^2
                gsq=gmet(1,1)*dble(ig1*ig1)+gmet(2,2)*dble(ig2*ig2)+&
     &           gmet(3,3)*dble(ig3*ig3)+2._dp*(gmet(2,1)*dble(ig1*ig2)+&
     &           gmet(3,1)*dble(ig1*ig3)+gmet(3,2)*dble(ig3*ig2))
     
     !          Skip g=0:
                if (gsq>1.0d-20) then
                  arg=fac*gsq
     
     !            Larger arg gives 0 contribution because of exp(-arg)
                  if (arg <= -minexparg ) then
     !              When any term contributes then include next shell
                    newg=1
                    term=exp(-arg)/gsq
                    summr = 0.0_dp
                    summi = 0.0_dp
     !              Note that if reduced atomic coordinates xred drift outside
     !              of unit cell (outside [0,1)) it is irrelevant in the following
     !              term, which only computes a phase.
                    do ia=1,natom
                      arg=2.0_dp*pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
     !                Sum real and imaginary parts (avoid complex variables)
                      summr=summr+zion(typat(ia))*cos(arg)
                      summi=summi+zion(typat(ia))*sin(arg)
                    end do
     
     !              The following two checks avoid an annoying underflow error message
                    if (abs(summr)<1.d-16) summr=0.0_dp
                    if (abs(summi)<1.d-16) summi=0.0_dp
     
     !              The product of term and summr**2 or summi**2 below
     !              can underflow if not for checks above
                    t1=term*(summr*summr+summi*summi)
                    gsum=gsum+t1
     
                    do ia=1,natom
     !                Again only phase is computed so xred may fall outside [0,1).
                      arg=2.0_dp*pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
                      phr= cos(arg)
                      phi=-sin(arg)
     !                (note: do not need real part, commented out)
     !                c1r=(phr*summr-phi*summi)*(term*zion(typat(ia)))
                      c1i=(phi*summr+phr*summi)*(term*zion(typat(ia)))
     !                compute coordinate gradients
                      grewtn(1,ia)=grewtn(1,ia)-c1i*ig1
                      grewtn(2,ia)=grewtn(2,ia)-c1i*ig2
                      grewtn(3,ia)=grewtn(3,ia)-c1i*ig3
                    end do
     
                  end if ! End condition of not larger than -minexparg
                end if ! End skip g=0
              end if ! End triple loop over G s and associated new shell condition
     
            end do
          end do
        end do
     
     !  Check if new shell must be calculated
        if (newg==0) exit
     
      end do !  End the loop on ng (new shells). Note that there is one exit from this loop.
     
      sumg=gsum/(2.0_dp*pi*ucvol)
     
     !Stress tensor is now computed elsewhere (ewald2) hence do not need
     !length scale gradients (used to compute them here).
     
     !normalize coordinate gradients by unit cell volume ucvol
      term=-2._dp/ucvol
      grewtn(:,:)=grewtn(:,:)*term
     !call DSCAL(3*natom,term,grewtn,1)
     
     !Conduct real space summations
      reta=sqrt(eta)
      fac=2._dp*sqrt(eta/pi)
      sumr=0.0_dp
     
     !In the following a summation is being conducted over all
     !unit cells (ir1, ir2, ir3) so it is appropriate to map all
     !reduced coordinates xred back into [0,1).
     !
     !Loop on shells in r-space as was done in g-space
      nr=0
      do
        nr=nr+1
        newr=0
     !   if (nr > 20 .and. mod(nr,10)==0) then
     !      write (message,'(3a,I10)') "Very large box of R neighbors in ewald: you probably do not want to do this.", ch10,&
     !&       " If you have a metal consider setting dipdip 0.  nr = ", nr
     !      MSG_WARNING(message)
     !   end if
     !
        do ir3=-nr,nr
          do ir2=-nr,nr
            do ir1=-nr,nr
              if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr&
     &         .or. nr==1 )then
     
                do ia=1,natom
     !            Map reduced coordinate xred(mu,ia) into [0,1)
                  fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
                  fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
                  fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
                  drdta1=0.0_dp
                  drdta2=0.0_dp
                  drdta3=0.0_dp
     
                  do ib=1,natom
                    fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
                    fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
                    fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
                    r1=dble(ir1)+fracb1-fraca1
                    r2=dble(ir2)+fracb2-fraca2
                    r3=dble(ir3)+fracb3-fraca3
                    rsq=rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3+&
     &               2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r1*r3)
     
     !              Avoid zero denominators in 'term':
                    if (rsq>=1.0d-24) then
     
     !                Note: erfc(8) is about 1.1e-29, so do not bother with larger arg.
     !                Also: exp(-64) is about 1.6e-28, so do not bother with larger arg**2 in exp.
                      term=0._dp
                      if (eta*rsq<64.0_dp) then
                        newr=1
                        rmagn=sqrt(rsq)
                        arg=reta*rmagn
     !                  derfc is the real(dp) complementary error function
                        derfc_arg = erfc(arg)
                        term=derfc_arg/rmagn
                        sumr=sumr+zion(typat(ia))*zion(typat(ib))*term
                        term=zion(typat(ia))*zion(typat(ib))*&
     &                   (term+fac*exp(-eta*rsq))/rsq
     !                  Length scale grads now handled with stress tensor in ewald2
                        r1a1d=rmet(1,1)*r1+rmet(1,2)*r2+rmet(1,3)*r3
                        r2a2d=rmet(2,1)*r1+rmet(2,2)*r2+rmet(2,3)*r3
                        r3a3d=rmet(3,1)*r1+rmet(3,2)*r2+rmet(3,3)*r3
     !                  Compute terms related to coordinate gradients
                        drdta1=drdta1+term*r1a1d
                        drdta2=drdta2+term*r2a2d
                        drdta3=drdta3+term*r3a3d
                      end if
                    end if ! End avoid zero denominators in'term'
                  end do ! end loop over ib:
                  !   print *, 'long(3) and short(3) ewald force terms', &
                 !            ia, nr, newr, grewtn(:,ia), drdta1,drdta2,drdta3
                  grewtn(1,ia)=grewtn(1,ia)+drdta1
                  grewtn(2,ia)=grewtn(2,ia)+drdta2
                  grewtn(3,ia)=grewtn(3,ia)+drdta3
                end do ! end loop over ia:
              end if
            end do ! end triple loop over real space points and associated condition of new shell
          end do
        end do
     
     !  Check if new shell must be calculated
        if(newr==0) exit
      end do ! End loop on nr (new shells). Note that there is an exit within the loop
     !
      sumr=0.5_dp*sumr
      fac=pi*ch**2/(2.0_dp*eta*ucvol)
     
     !Finally assemble Ewald energy, eew
      eew=sumg+sumr-chsq*reta/sqrt(pi)-fac
     
     !DEBUG
     !write(std_out,*)'eew=sumg+sumr-chsq*reta/sqrt(pi)-fac'
     !write(std_out,*)eew,sumg,sumr,chsq*reta/sqrt(pi),fac
     !ENDDEBUG
     
     !Length scale grads handled with stress tensor, ewald2
     
     !Output the final values of ng and nr
     ! write(message, '(a,a,i4,a,i4)' )ch10,' ewald : nr and ng are ',nr,' and ',ng
     ! call wrtout(std_out,message,'COLL')
     
     end subroutine abinit_ewald

     !!****f* m_ewald/ewald2
!!
!! NAME
!! ewald2
!!
!! FUNCTION
!! Compute the part of the stress tensor coming from the Ewald energy
!! which is calculated by derivating the Ewald energy with respect to strain.
!! See Nielsen and Martin, Phys. Rev. B 32, 3792 (1985).
!! Definition of stress tensor is $(1/ucvol)*d(Etot)/d(strain(a,b))$.
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in umit cell
!! ntypat=number of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2) (inverse transpose of gmet)
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! $stress(6)=(1/ucvol)*gradient$ of Ewald energy with respect to strain,
!!      in hartrees/bohr^3
!! Cartesian components of stress are provided for this symmetric
!! tensor in the order 11 22 33 32 31 21.
!!
!! PARENTS
!!      stress
!!
!! CHILDREN
!!      matr3inv,spline
!!
!! SOURCE

subroutine abinit_ewald2(gmet,natom,ntypat,rmet,rprimd,stress,typat,ucvol,xred,zion)
     
      implicit none
     
     !Arguments ------------------------------------
     !scalars
      integer,intent(in) :: natom,ntypat
      real(dp),intent(in) :: ucvol
     !arrays
      integer,intent(in) :: typat(natom)
      real(dp),intent(in) :: gmet(3,3),rmet(3,3),rprimd(3,3),xred(3,natom)
      real(dp),intent(in) :: zion(ntypat)
      real(dp),intent(out) :: stress(6)
     
     !Local variables-------------------------------
     !scalars
      integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
      real(dp) :: arg1,arg2,arg3,ch,dderfc,derfc_arg,direct,eta,fac,fraca1
      real(dp) :: fraca2,fraca3,fracb1,fracb2,fracb3,g1,g2,g3,gsq,r1,r1c,r2,r2c
      real(dp) :: minexparg
      real(dp) :: r3,r3c,recip,reta,rmagn,rsq,summi,summr,t1,t2,t3,t4,t5,t6,term1
      real(dp) :: term2,term3,term4
     !arrays
      real(dp) :: gprimd(3,3),strg(6),strr(6)
     
     ! *************************************************************************
     
     !Define dimensional reciprocal space primitive translations gprimd
     !(inverse transpose of rprimd)

     ! call matr3inv(rprimd,gprimd)
      !SHRED USES ONLY ORTHOGONAL COORDINATES
      gprimd=0
      gprimd(1,1)=1.0_dp/rprimd(1,1)
      gprimd(2,2)=1.0_dp/rprimd(2,2)
      gprimd(3,3)=1.0_dp/rprimd(3,3)

     !This is the minimum argument of an exponential, with some safety
      minexparg=log(tiny(0._dp))+5.0_dp
     
     !Add up total charge and sum of charge^2 in cell
      ch=0._dp
      do ia=1,natom
        ch=ch+zion(typat(ia))
      end do
     
     !Compute eta, the Ewald summation convergence parameter,
     !for approximately optimized summations:
      direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
     & rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
      recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
     & gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
     !Here, a bias is introduced, because G-space summation scales
     !better than r space summation !
      eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)
     
      fac=pi**2/eta
     
     !Conduct reciprocal space summations
      strg(1:6)=0.0_dp
     
     !Sum over G space, done shell after shell until all
     !contributions are too small
      ng=0
      do
        ng=ng+1
        newg=0
     
        do ig3=-ng,ng
          do ig2=-ng,ng
            do ig1=-ng,ng
     
     !        Exclude shells previously summed over
              if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng .or. ng==1 ) then
     
     !          Compute Cartesian components of each G
                g1=gprimd(1,1)*ig1+gprimd(1,2)*ig2+gprimd(1,3)*ig3
                g2=gprimd(2,1)*ig1+gprimd(2,2)*ig2+gprimd(2,3)*ig3
                g3=gprimd(3,1)*ig1+gprimd(3,2)*ig2+gprimd(3,3)*ig3
     !          Compute |G|^2 (no pi factors)
                gsq=(g1**2+g2**2+g3**2)
     
     !          skip g=0:
                if (gsq>1.0d-20) then
                  arg1=fac*gsq
     
     !            larger arg1 gives 0 contribution because of exp(-arg1)
                  if (arg1<= -minexparg) then
     !              When any term contributes then include next shell
                    newg=1
                    term1=exp(-arg1)/arg1
                    summr = 0.0_dp
                    summi = 0.0_dp
                    do ia=1,natom
                      arg2=2.0_dp*pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
     !                Sum real and imaginary parts (avoid complex variables)
                      summr=summr+zion(typat(ia))*cos(arg2)
                      summi=summi+zion(typat(ia))*sin(arg2)
                    end do
     
     !              Avoid underflow error messages
                    if (abs(summr)<1.d-16) summr=0.0_dp
                    if (abs(summi)<1.d-16) summi=0.0_dp
     
                    term2=(2._dp/gsq)*(1._dp+arg1)
                    t1=term2*g1*g1-1._dp
                    t2=term2*g2*g2-1._dp
                    t3=term2*g3*g3-1._dp
                    t4=term2*g2*g3
                    t5=term2*g1*g3
                    t6=term2*g1*g2
                    term3=term1*(summr*summr+summi*summi)
                    strg(1)=strg(1)+t1*term3
                    strg(2)=strg(2)+t2*term3
                    strg(3)=strg(3)+t3*term3
                    strg(4)=strg(4)+t4*term3
                    strg(5)=strg(5)+t5*term3
                    strg(6)=strg(6)+t6*term3
     
                  end if ! End condition not being larger than -minexparg
                end if ! End skip g=0
     
              end if ! End triple loop and condition of new shell
            end do
          end do
        end do
     
     !  Check if new shell must be calculated
        if (newg==0) exit
      end do ! End loop on new shell. Note that there is an "exit" instruction within the loop
     
     
     !Conduct real space summations
      reta=sqrt(eta)
      strr(1:6)=0.0_dp
     
     !Loop on shells in r-space as was done in g-space
      nr=0
      do
        nr=nr+1
        newr=0
     
        do ir3=-nr,nr
          do ir2=-nr,nr
            do ir1=-nr,nr
              if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1 )then
     
                do ia=1,natom
     !            Convert reduced atomic coordinates to [0,1)
                  fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
                  fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
                  fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
                  do ib=1,natom
                    fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
                    fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
                    fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
                    r1=ir1+fracb1-fraca1
                    r2=ir2+fracb2-fraca2
                    r3=ir3+fracb3-fraca3
     !              Convert from reduced to cartesian coordinates
                    r1c=rprimd(1,1)*r1+rprimd(1,2)*r2+rprimd(1,3)*r3
                    r2c=rprimd(2,1)*r1+rprimd(2,2)*r2+rprimd(2,3)*r3
                    r3c=rprimd(3,1)*r1+rprimd(3,2)*r2+rprimd(3,3)*r3
     !              Compute |r|^2
                    rsq=r1c**2+r2c**2+r3c**2
                    rmagn=sqrt(rsq)
     
     !              Avoid zero denominators in 'term':
                    if (rmagn>=1.0d-12) then
     
     !                Note: erfc(8) is about 1.1e-29, so do not bother with larger arg.
     !                Also: exp(-64) is about 1.6e-28, so do not bother with larger arg**2 in exp.
                      arg3=reta*rmagn
                      if (arg3<8.0_dp) then
                        newr=1
     !                  derfc computes the complementary error function
     !                  dderfc is the derivative of the complementary error function
                        dderfc=(-2/sqrt(pi))*exp(-eta*rsq)
                        derfc_arg = erfc(arg3)
                        term3=dderfc-derfc_arg/arg3
                        term4=zion(typat(ia))*zion(typat(ib))*term3
                        strr(1)=strr(1)+term4*r1c*r1c/rsq
                        strr(2)=strr(2)+term4*r2c*r2c/rsq
                        strr(3)=strr(3)+term4*r3c*r3c/rsq
                        strr(4)=strr(4)+term4*r2c*r3c/rsq
                        strr(5)=strr(5)+term4*r1c*r3c/rsq
                        strr(6)=strr(6)+term4*r1c*r2c/rsq
                      end if ! End the condition of not being to large
                    end if ! End avoid zero denominator
     
                  end do ! End loop over ib:
                end do  ! End loop over ia:
     
              end if ! End triple loop overs real space points, and associated new shell condition
            end do
          end do
        end do
     
     !  Check if new shell must be calculated
        if(newr==0) exit
      end do ! End loop on new shells
     
     !Finally assemble stress tensor coming from Ewald energy, stress
     !(note division by unit cell volume in accordance with definition
     !found in Nielsen and Martin, Phys. Rev. B 32, 3792 (1985).)
     
      fac = pi/(2._dp*ucvol*eta)
      stress(1)=(0.5_dp*reta*strr(1)+fac*(strg(1)+(ch**2)))/ucvol
      stress(2)=(0.5_dp*reta*strr(2)+fac*(strg(2)+(ch**2)))/ucvol
      stress(3)=(0.5_dp*reta*strr(3)+fac*(strg(3)+(ch**2)))/ucvol
      stress(4)=(0.5_dp*reta*strr(4)+fac*strg(4))/ucvol
      stress(5)=(0.5_dp*reta*strr(5)+fac*strg(5))/ucvol
      stress(6)=(0.5_dp*reta*strr(6)+fac*strg(6))/ucvol
     
     end subroutine abinit_ewald2

    end module