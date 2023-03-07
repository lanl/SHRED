module Spherical_Harmonic_mod
    use types, only: dp
    use constants, only : i_, pi
    implicit none

    public :: real_spherical_harmonic, real_spherical_harmonic_and_derivatives, &
    real_Dx_spherical_harmonic, real_Dy_spherical_harmonic, real_Dz_spherical_harmonic, &
    complex_spherical_harmonic, complex_spherical_harmonic_and_derivatives, &
    complex_Dx_spherical_harmonic, complex_Dy_spherical_harmonic, complex_Dz_spherical_harmonic

contains

subroutine complex_spherical_harmonic_and_derivatives(l,m,vec, Ylm, dYlm)
        integer, intent(in) :: l, m
        real(dp), intent(in) :: vec(3)
        complex(dp), intent(out) :: Ylm, dYlm(3)

        Ylm=complex_spherical_harmonic(l, m,vec)
        dYlm(1)=complex_Dx_spherical_harmonic(l,m,vec(:),Ylm)
        dYlm(2)=complex_Dy_spherical_harmonic(l,m,vec(:),Ylm)
        dYlm(3)=complex_Dz_spherical_harmonic(l,m,vec(:),Ylm)

end subroutine

subroutine real_spherical_harmonic_and_derivatives(l,m,vec, Ylm, dYlm)
        integer, intent(in) :: l, m
        real(dp), intent(in) :: vec(3)
        real(dp), intent(out) :: Ylm, dYlm(3)
        complex(dp) :: Ylm_c, dYlm_c(3)

        call complex_spherical_harmonic_and_derivatives(l,m,vec, Ylm_c, dYlm_c)
        if(m.lt.0) then
                Ylm=sqrt(2._dp)*(-1.0_dp)**m*aimag(Ylm_c)
                dYlm=sqrt(2._dp)*(-1.0_dp)**m*aimag(dYlm_c)
        else if(m.gt.0) then
                Ylm=sqrt(2._dp)*(-1.0_dp)**m*real(Ylm_c)
                dYlm=sqrt(2._dp)*(-1.0_dp)**m*real(dYlm_c)
        else
                Ylm=real(Ylm_c)
                dYlm=real(dYlm_c)
        endif
end subroutine

real(dp) function real_spherical_harmonic2(l,m,vec) result(Ylm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
complex(dp) :: Ylm_c

Ylm_c=complex_spherical_harmonic(l,m,vec)

if(m.lt.0) then
        Ylm=sqrt(2._dp)*(-1.0_dp)**m*aimag(Ylm_c)
else if(m.gt.0) then
        Ylm=sqrt(2._dp)*(-1.0_dp)**m*real(Ylm_c)
else
        Ylm=real(Ylm_c)
endif

end function

real(dp) function real_Dx_spherical_harmonic(l,m,vec) result(dYlm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
complex(dp) :: Ylm_c, dYlm_c

Ylm_c=complex_spherical_harmonic(l,m,vec)
dYlm_c=complex_Dx_spherical_harmonic(l,m,vec,Ylm_c)

if(m.lt.0) then
        dYlm=sqrt(2._dp)*(-1.0_dp)**m*aimag(dYlm_c)
else if(m.gt.0) then
        dYlm=sqrt(2._dp)*(-1.0_dp)**m*real(dYlm_c)
else
        dYlm=real(Ylm_c)
endif

end function

real(dp) function real_Dy_spherical_harmonic(l,m,vec) result(dYlm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
complex(dp) :: Ylm_c, dYlm_c

Ylm_c=complex_spherical_harmonic(l,m,vec)
dYlm_c=complex_Dy_spherical_harmonic(l,m,vec,Ylm_c)

if(m.lt.0) then
        dYlm=sqrt(2._dp)*(-1.0_dp)**m*aimag(dYlm_c)
else if(m.gt.0) then
        dYlm=sqrt(2._dp)*(-1.0_dp)**m*real(dYlm_c)
else
        dYlm=real(Ylm_c)
endif

end function

real(dp) function real_Dz_spherical_harmonic(l,m,vec) result(dYlm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
complex(dp) :: Ylm_c, dYlm_c

Ylm_c=complex_spherical_harmonic(l,m,vec)
dYlm_c=complex_Dz_spherical_harmonic(l,m,vec,Ylm_c)

if(m.lt.0) then
        dYlm=sqrt(2._dp)*(-1.0_dp)**m*aimag(dYlm_c)
else if(m.gt.0) then
        dYlm=sqrt(2._dp)*(-1.0_dp)**m*real(dYlm_c)
else
        dYlm=real(Ylm_c)
endif

end function

real(dp) function real_spherical_harmonic(l,m,vec) result(Ylm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
real(dp) :: x,y,z, r2, r

x=vec(1)
y=vec(2)
z=vec(3)

if((abs(m).gt.l).or.(l.lt.0)) then
        print *, "check your Spherical Harmonic definition"
        stop
else if (l.gt.3) then
        print*, "Spherical Harmonics maxed at l=3"
        stop
endif

Ylm=0.0_dp
r2=x**2+y**2+z**2
if(r2.lt.tiny(1.0_dp)) then
        if((m.eq.0).and.(l.eq.0)) then 
                Ylm=Sqrt(0.25_dp/pi)
        endif
        return
endif

if(l.eq.0) then
        Ylm=Sqrt(0.25_dp/pi)
else if(l.eq.1) then
        r=sqrt(r2)
             if(m.eq.-1) then
                Ylm=Sqrt(0.75_dp/pi)*y/r
        else if(m.eq.0) then
                Ylm=Sqrt(0.75_dp/pi)*z/r
        else if(m.eq.1) then
                Ylm=Sqrt(0.75_dp/pi)*x/r
        endif
else if(l.eq.2) then
             if(m.eq.-2) then
                Ylm=Sqrt(3.75_dp/pi)*x*y/r2
        else if(m.eq.-1) then
                Ylm=Sqrt(3.75_dp/pi)*y*z/r2
        else if(m.eq.0) then
                Ylm=Sqrt(1.25_dp/pi)*(2*z**2-x**2-y**2)/r2
        else if(m.eq.1) then
                Ylm=Sqrt(3.75_dp/pi)*x*z/r2
        else if(m.eq.2) then
                Ylm=Sqrt(0.9375/pi)*(x**2-y**2)/r2
        endif
else if(l.eq.3) then
        r=sqrt(r2)
             if(m.eq.-3) then
                Ylm=Sqrt(1.09375_dp/pi)*(3*x**2-y**2)*y/r**3
        else if(m.eq.-2) then
                Ylm=Sqrt(26.25_dp/pi)*x*y*z/r**3
        else if(m.eq.-1) then
                Ylm=Sqrt(0.65625_dp/pi)*y*(4*z**2-x**2-y**2)/r**3
        else if(m.eq.0) then
                Ylm=Sqrt(0.4375_dp/pi)*z*(2*z**2-3*x**2-3*y**2)/r**3
        else if(m.eq.1) then
                Ylm=Sqrt(0.65625_dp/pi)*x*(4*z**2-x**2-y**2)/r**3
        else if(m.eq.2) then
                Ylm=Sqrt(6.5625_dp/pi)*(x*2-y**2)*z/r**3
        else if(m.eq.3) then
                Ylm=Sqrt(1.09375_dp/pi)*(x**2-3*y**2)*x/r**3
        endif
endif

end function

complex(dp) function complex_spherical_harmonic(l,m,vec) result(Ylm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
real(dp) :: x,y,z, r2, r

x=vec(1)
y=vec(2)
z=vec(3)

if((abs(m).gt.l).or.(l.lt.0)) then
        print *, "check your Spherical Harmonic definition"
        stop
else if (l.gt.3) then
        print *, "Spherical Harmonics maxed at l=3"
        stop
endif

Ylm=0.0_dp
r2=x**2+y**2+z**2
if(r2.lt.tiny(1.0_dp)) then
        if((m.eq.0).and.(l.eq.0)) then 
                Ylm=Sqrt(0.25_dp/pi)
        endif
        return
endif

if(l.eq.0) then
        Ylm=Sqrt(0.25_dp/pi)
else if(l.eq.1) then
        r=sqrt(r2)
             if(m.eq.-1) then
                Ylm=Sqrt(0.375_dp/pi)*(x-i_*y)/r
        else if(m.eq.0) then
                Ylm=Sqrt(0.75_dp/pi)*z/r
        else if(m.eq.1) then
                Ylm=-Sqrt(0.375_dp/pi)*(x+i_*y)/r
        endif
else if(l.eq.2) then
             if(m.eq.-2) then
                Ylm=Sqrt(0.46875_dp/pi)*(x-i_*y)**2/r2
        else if(m.eq.-1) then
                Ylm=Sqrt(1.875_dp/pi)*(x-i_*y)*z/r2
        else if(m.eq.0) then
                Ylm=Sqrt(0.3125_dp/pi)*(2*z**2-x**2-y**2)/r2
        else if(m.eq.1) then
                Ylm=-Sqrt(1.875_dp/pi)*(x+i_*y)*z/r2
        else if(m.eq.2) then
                Ylm=Sqrt(0.46875_dp/pi)*(x+i_*y)**2/r2
        endif
else if(l.eq.3) then
        r=sqrt(r2)
             if(m.eq.-3) then
                Ylm=Sqrt(0.546875_dp/pi)*((x-i_*y)/r)**3
        else if(m.eq.-2) then
                Ylm=Sqrt(3.28125_dp/pi)*(x-i_*y)**2*z/r**3
        else if(m.eq.-1) then
                Ylm=Sqrt(0.328125_dp/pi)*(x-i_*y)*(4*z**2-x**2-y**2)/r**3
        else if(m.eq.0) then
                Ylm=Sqrt(0.4375_dp/pi)*z*(2*z**2-3*x**2-3*y**2)/r**3
        else if(m.eq.1) then
                Ylm=-Sqrt(0.328125_dp/pi)*(x+i_*y)*(4*z**2-x**2-y**2)/r**3
        else if(m.eq.2) then
                Ylm=Sqrt(3.28125_dp/pi)*(x+i_*y)**2*z/r**3
        else if(m.eq.3) then
                Ylm=-Sqrt(0.546875_dp/pi)*((x+i_*y)/r)**3
        endif
endif

end function

complex(dp) function complex_Dx_spherical_harmonic(l,m,vec,Ylm) result(dYlm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
complex(dp),intent(in) :: Ylm
real(dp) :: x,y,z, r2, r, dr_inv
complex(dp) :: tmp

x=vec(1)
y=vec(2)
z=vec(3)

if((abs(m).gt.l).or.(l.lt.0)) then
        print *, "check your Spherical Harmonic definition"
        stop
else if (l.gt.3) then
        print *, "Spherical Harmonics maxed at l=3"
        stop
endif

r2=x**2+y**2+z**2
if(r2.lt.tiny(1.0_dp).or.(l.eq.0)) then
        dYlm=0.0_dp
        return
endif
r=sqrt(r2)
dr_inv=-l*x/r2 !missing 1/r^l comes from Y
dYlm=dr_inv*Ylm !Ylm is result of spherical_harmonic(l,m,vec)
tmp=0.0
if(l.eq.1) then
        if(m.eq.-1) then
                tmp=Sqrt(0.375_dp/pi)/r
        else if(m.eq.1) then
                tmp=-Sqrt(0.375_dp/pi)/r
        endif
else if(l.eq.2) then
        if(m.eq.-2) then
                tmp=Sqrt(0.46875_dp/pi)*2.0_dp*(x-i_*y)/r2
        else if(m.eq.-1) then
                tmp=Sqrt(1.875_dp/pi)*z/r2
        else if(m.eq.0) then
                tmp=Sqrt(0.3125_dp/pi)*(-2.0_dp*x)/r2
        else if(m.eq.1) then
                tmp=-Sqrt(1.875_dp/pi)*z/r2
        else if(m.eq.2) then
                tmp=Sqrt(0.46875_dp/pi)*2.0_dp*(x+i_*y)/r2
        endif
else if(l.eq.3) then
        if(m.eq.-3) then
                tmp=Sqrt(0.546875_dp/pi)*3.0_dp*(x-i_*y)**2/r**3
        else if(m.eq.-2) then
                tmp=Sqrt(3.28125_dp/pi)*2.0_dp*(x-i_*y)*z/r**3
        else if(m.eq.-1) then
                tmp=Sqrt(0.328125_dp/pi)*((5.0_dp*z**2-r2)-(x-i_*y)*(2.0_dp*x))/r**3
        else if(m.eq.0) then
                tmp=Sqrt(0.4375_dp/pi)*z*(-6.0_dp*x)/r**3
        else if(m.eq.1) then
                tmp=-Sqrt(0.328125_dp/pi)*((5.0_dp*z**2-r2)-(x+i_*y)*(2.0_dp*x))/r**3
        else if(m.eq.2) then
                tmp=Sqrt(3.28125_dp/pi)*2.0_dp*(x+i_*y)*z/r**3
        else if(m.eq.3) then
                tmp=-Sqrt(0.546875_dp/pi)*3.0_dp*(x+i_*y)**2/r**3
        endif
endif

dYlm=dYlm+tmp

end function
complex(dp) function complex_Dy_spherical_harmonic(l,m,vec,Ylm) result(dYlm)
integer, intent(in) :: l, m
complex(dp),intent(in) :: Ylm
real(dp), intent(in) :: vec(3)
real(dp) :: x,y,z, r2, r, dr_inv
complex(dp) :: tmp

x=vec(1)
y=vec(2)
z=vec(3)

if((abs(m).gt.l).or.(l.lt.0)) then
        print *, "check your Spherical Harmonic definition"
        stop
else if (l.gt.3) then
        print *, "Spherical Harmonics maxed at l=3"
        stop
endif

r2=x**2+y**2+z**2
if(r2.lt.tiny(1.0_dp).or.(l.eq.0)) then
        dYlm=0.0_dp
        return
endif

r=sqrt(r2)
dr_inv=-l*y/r2 !missing 1/r^l comes from Y
dYlm=dr_inv*Ylm !Ylm is result of spherical_harmonic(l,m,vec)
tmp=0.0
if(l.eq.1) then
        if(m.eq.-1) then
                tmp=Sqrt(0.375_dp/pi)*(-i_)/r
        else if(m.eq.1) then
                tmp=-Sqrt(0.375_dp/pi)*(i_)/r
        endif
else if(l.eq.2) then
        if(m.eq.-2) then
                tmp=Sqrt(0.46875_dp/pi)*(-i_)*2.0_dp*(x-i_*y)/r2
        else if(m.eq.-1) then
                tmp=Sqrt(1.875_dp/pi)*(-i_)*z/r2
        else if(m.eq.0) then
                tmp=Sqrt(0.3125_dp/pi)*(-2.0_dp*y)/r2
        else if(m.eq.1) then
                tmp=-Sqrt(1.875_dp/pi)*i_*z/r2
        else if(m.eq.2) then
                tmp=Sqrt(0.46875_dp/pi)*i_*2.0_dp*(x+i_*y)/r2
        endif
else if(l.eq.3) then
        if(m.eq.-3) then
                tmp=Sqrt(0.546875_dp/pi)*(-i_)*3.0_dp*(x-i_*y)**2/r**3
        else if(m.eq.-2) then
                tmp=Sqrt(3.28125_dp/pi)*(-i_)*2.0_dp*(x-i_*y)*z/r**3
        else if(m.eq.-1) then
                tmp=Sqrt(0.328125_dp/pi)*((-i_)*(5.0_dp*z**2-r2)-(x-i_*y)*(2.0_dp*y))/r**3
        else if(m.eq.0) then
                tmp=Sqrt(0.4375_dp/pi)*z*(-6.0_dp*y)/r**3
        else if(m.eq.1) then
                tmp=-Sqrt(0.328125_dp/pi)*(i_*(5.0*z**2-r2)-(x+i_*y)*(2.0_dp*y))/r**3
        else if(m.eq.2) then
                tmp=Sqrt(3.28125_dp/pi)*(i_*2.0_dp)*(x+i_*y)*z/r**3
        else if(m.eq.3) then
                tmp=-Sqrt(0.546875_dp/pi)*i_*3.0_dp*(x+i_*y)**2/r**3
        endif

endif
dYlm=dYlm+tmp

end function

complex(dp) function complex_Dz_spherical_harmonic(l,m,vec,Ylm) result(dYlm)
integer, intent(in) :: l, m
real(dp), intent(in) :: vec(3)
complex(dp),intent(in) :: Ylm
real(dp) :: x,y,z, r2, r, dr_inv
complex(dp) :: tmp

x=vec(1)
y=vec(2)
z=vec(3)

if((abs(m).gt.l).or.(l.lt.0)) then
        print *, "check your Spherical Harmonic definition"
        stop
else if (l.gt.3) then
        print *, "Spherical Harmonics maxed at l=3"
        stop
endif

r2=x**2+y**2+z**2
if(r2.lt.tiny(1.0_dp).or.(l.eq.0)) then
        dYlm=0.0_dp
        return
endif

r=sqrt(r2)
dr_inv=-l*z/r2 !missing 1/r^l comes from Y
dYlm=dr_inv*Ylm !Ylm is result of spherical_harmonic(l,m,vec)

tmp=0.0
if(l.eq.1) then
        if(m.eq.0) then
                tmp=Sqrt(0.75_dp/pi)/r
        endif
else if(l.eq.2) then
        if(m.eq.-1) then
                tmp=Sqrt(1.875_dp/pi)*(x-i_*y)/r2
        else if(m.eq.0) then
                tmp=Sqrt(0.3125_dp/pi)*4.0_dp*z/r2
        else if(m.eq.1) then
                tmp=-Sqrt(1.875_dp/pi)*(x+i_*y)/r2
        endif
else if(l.eq.3) then
        if(m.eq.-2) then
                tmp=Sqrt(3.28125_dp/pi)*(x-i_*y)**2/r**3
        else if(m.eq.-1) then
                tmp=Sqrt(0.328125_dp/pi)*(x-i_*y)*(8.0_dp*z)/r**3
        else if(m.eq.0) then
                tmp=Sqrt(0.4375_dp/pi)*(4.0_dp*z**2)/r**3
        else if(m.eq.1) then
                tmp=-Sqrt(0.328125_dp/pi)*(x+i_*y)*(8.0_dp*z)/r**3
        else if(m.eq.2) then
                tmp=Sqrt(3.28125_dp/pi)*(x+i_*y)**2/r**3
        endif
endif

dYlm=dYlm+tmp

end function

    
end module