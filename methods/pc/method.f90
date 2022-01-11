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
  !  do i=1,nx
  !    up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))
  !  end do 
    !	MacCormack Method
    !	Predictor Step 
    do i=1,nx
      upp(:,i) = u(:,i) - dtx*(f(:,i+1)-f(:,i)) + eta*(u(:,i-1) + u(:,i+1) - 2*u(:,i))
    end do
    call boundaries(upp)
    call fluxes(upp,f)
  !
    !	Corrector Step
    do i=1,nx
      up(:,i)=0.5*(u(:,i)+upp(:,i) - dtx*(f(:,i)-f(:,i-1))) + eta*(upp(:,i-1) + upp(:,i+1) - 2*upp(:,i))
    end do
  !   Boundary conditions to the U^n+1
    call boundaries(up)
  
    ! copy the up to the u
    u(:,:)=up(:,:)
  
    return
  end subroutine tstep
end module method