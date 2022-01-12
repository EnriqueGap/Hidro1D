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
    real, dimension(neq,0:nx+1,0:ny+1) :: up
    real                        :: dtx,dty
    integer                     :: i,j
  
    !  obtain the fluxes
    !
    call fluxes(u,f,g)
  
    !   Here is the Lax method, notice that the values at the extremes can
    !   not be calculated, we need to enter then as boundary conditions
    dtx=dt/dx
    dty=dt/dy
  !
    do i=1,nx
      do j=1,ny
        up(:,i,j)=0.5*(u(:,i-1,j)+u(:,i+1,j) - dtx*(f(:,i+1,j)-f(:,i-1,j)) - dty*(g(:,i,j+1)-g(:,i,j-1)))
      end do
    end do
  !
  !   Boundary conditions to the U^n+1
    call boundaries(up)
  
    ! copy the up to the u
    u(:,:,:)=up(:,:,:)
  
    return
  end subroutine tstep
  end module method
