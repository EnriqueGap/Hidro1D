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
    real, parameter :: u0=0.0, u1=0., p0=0.15, p1=1.
    real, parameter :: rho0=0.1, rho1=1.,x1=0.5
    real :: x
    integer :: i
  
    !  fill the vector u
    do i=0, nx+1
      x=real(i)*dx   ! obtain the position $x_i$
      if( x < x1 )  then
        ! rho
        u(1,i)=rho1
        ! v
        u(2,i)=u1
        ! pressure/(gamma-1)
        u(neq,i)=p1/(gamma-1.)
      else
        u(1,i)=rho0
        u(2,i)=u0
        u(neq,i)=p0/(gamma-1.)
      end if
  
      if( (x-0.5*dx <= x1).and.(x+0.5*dx >= x1) ) then
        u(1,i)=(rho0+rho1)/2.
        u(2,i)=(u0+u1)/2.
        u(neq,i)=(p0+p1)/(2.*(gamma-1.))
      end if
    end do
    ! reset the counters and time to 0
    time=0.
    tprint=0.
    itprint=0
    return
  end subroutine initflow
end module initial
