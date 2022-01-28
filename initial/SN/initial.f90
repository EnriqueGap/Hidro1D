module initial
use globals
implicit none
contains
  !=======================================================================
  ! generates initial condition
  subroutine initflow(time, tprint, itprint)
    use globals
    implicit none
    real, intent(out) :: time, tprint
    integer, intent (out) :: itprint
    !internal variables
    real, parameter :: n= 100., T=1000., u1=0, u0=9.5e8
    real :: x
    integer :: i
    real :: R, E1, rho1, E0, rho0
    R = 20*Km !20Km
    E1 = n*T*boltz
    rho1 = n*mh
    E0 = 3.*10d51/(4.*pi*R**3) ! E = 10^51 ergs / V
    rho0 = 3*SM/(4.*pi*R**3) ! rho = 1SM / V
    !  fill the vector u
    do i=0, nx+1
      x=real(i)*dx   ! obtain the position $x_i$
      if( x <= R )  then
        ! rho
        u(1,i)=rho0
        ! v
        u(2,i)=u0*rho0
        ! pressure/(gamma-1)
        u(neq,i)=E0
      else
        u(1,i)=rho1
        u(2,i)=u1*rho1
        u(neq,i)=E1
      end if
  
!      if( (x-0.5*dx <= R).and.(x+0.5*dx >= R) ) then
!        u(1,i)=(rho0+rho1)/2.
!        u(2,i)=(u0+u1)/2.
!        u(neq,i)=(E0+E1)
!      end if
    end do
    ! reset the counters and time to 0
    time=0.
    tprint=0.
    itprint=0
    return
  end subroutine initflow
end module initial
