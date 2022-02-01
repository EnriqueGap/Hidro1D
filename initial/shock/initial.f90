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
    real, parameter :: ux0=0.0, uy0=0., ux1=0.0, uy1=0.                   !(vx,vy)=v0 for x<xm & (vx,vy)=v1 for x>xm
    real, parameter :: rho0=1., rho1=0.1, p0=1., p1=0.15, xm=0.5
    real :: x
    integer :: i
    !  fill the vector u
    do i=0, nx+1
      ! obtain the position $x_i$
      x=real(i)*dx
      if( x < xm )  then
        u(1,i,:)=rho0
        u(2,i,:)=ux0
        u(3,i,:)=uy0
        u(neq,i,:)=p0/(gamma-1.)
      else
        u(1,i,:)=rho1
        u(2,i,:)=ux1
        u(3,i,:)=uy1
        u(neq,i,:)=p1/(gamma-1.)
      end if
      ! middle
      if( (x-0.5*dx <= xm).and.(x+0.5*dx >= xm) ) then
        u(1,i,:)=(rho0+rho1)/2.
        u(2,i,:)=(ux0+ux1)/2.
        u(3,i,:)=(uy0+uy1)/2.
        u(neq,i,:)=(p0+p1)/(2.*(gamma-1.))
      end if
    end do
    ! reset the counters and time to 0
    time=0.
    tprint=0.
    itprint=0
    return
  end subroutine initflow
end module initial
