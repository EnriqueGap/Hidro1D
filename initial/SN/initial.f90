module initial
use globals
use mesh
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
    real, parameter :: n= 100., T=1000., u1=0
    real :: x, R, R0, E1, rho1, E0, rho0, u0
    integer :: i
    R = 100*Km				 ! 100Km SNR initial radius
    u0 = 9500.*Km			 ! SNR velocity 9500Km/s
    E1 = n*T*boltz/(gamma-1.)		 ! PV=nKT -> P/(gamma -1) = Et = Etot ISM
    E0 = 3.*1d51/(4.*pi*R**3) 	 ! E = 10^51 ergs / Vesf
    rho1 = n*mh			 ! ISM mass density
    R0 = 0.5*R
    rho0 = 3*0.25*SM/(4.*pi*R0**3.)	 ! mass density
    !  fill the vector u
    do i=0, nx+1
      x=real(i)*dx   ! obtain the position $x_i$
      if(x <= R0)  then
        ! rho
        u(1,i)=rho0
        ! momentum
        u(2,i)=rho0*u0*x/R
        ! Energy density
        u(neq,i)=E0
      else if(x <= R .and. x>R0) then
        u(1,i)=0.75*SM/(4.*pi*(R-R0)*R**2.)*(R/x)**2
        u(2,i)=u0*(x/R)*u(1,i)
        u(neq,i)=E0
      else
        u(1,i)=rho1
        u(2,i)=u1*rho1
        u(neq,i)=E1
      end if
    end do
    ! reset the counters and time to 0
    time=0.
    tprint=0.
    itprint=0
    return
  end subroutine initflow
end module initial
