module method
use globals
use general
implicit none
contains
!=======================================================================
  ! integration from t to t+dt with the method of Lax
  subroutine tstep(dt,time)
    use globals
    implicit none
    real,            intent(in) :: dt, time
    !internal variables
    real, dimension(neq,0:nx+1) :: up, upp
    real                        :: dtx, eta
    integer                     :: i
  
    !  obtain the fluxes
    eta=0.09
    !
    call fluxes(u,f)
  
    !   Here is the Lax method, notice that the values at the extremes can
    !   not be calculated, we need to enter then as boundary conditions
    dtx=dt/dx
  !
    do i=1,nx
      up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))
    end do 
  !
  !   Boundary conditions to the U^n+1
    call boundaries(up)
  
    ! copy the up to the u
    u(:,:)=up(:,:)
  
    return
  end subroutine tstep
end module method
