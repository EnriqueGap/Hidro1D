  module method
  use globals
  use general
  implicit none
  contains
  !=======================================================================
  ! integration from t to t+dt with the method of MacCormack
  subroutine tstep(dt,time)
    use globals
    implicit none
    real,            intent(in) :: dt, time
    !internal variables
    real, dimension(neq,0:nx+1,0:ny+1) :: up,upp
    real                        :: dtx,dty,eta
    integer                     :: i,j
  
    eta=0.09
    !  obtain the fluxes
    !
    call fluxes(u,f,g)

    dtx=dt/dx
    dty=dt/dy
  ! Predictor step
    do i=1,nx
      do j=1,ny
        upp(:,i,j) = u(:,i,j) - dtx*(f(:,i+1,j)-f(:,i,j)) - dty*(g(:,i,j+1)-g(:,i,j)) &
        + eta*(u(:,i-1,j) + u(:,i+1,j) + u(:,i,j-1) + u(:,i,j+1) - 4*u(:,i,j))
      end do
    end do
    
    call boundaries(upp)
    call fluxes(upp,f,g)
    !Corrector Step
    do i=1,nx
      do j=1,ny
      	up(:,i,j)=0.5*(u(:,i,j)+upp(:,i,j) &
      	- dtx*(f(:,i,j)-f(:,i-1,j)) &
      	- dty*(g(:,i,j)-g(:,i,j-1))) &
      	+ eta*(upp(:,i-1,j) + upp(:,i+1,j) + upp(:,i,j-1) + upp(:,i,j+1) - 4*upp(:,i,j))
      end do
    end do  
    call boundaries(up)
    ! copy the up to the u
    u(:,:,:)=up(:,:,:)
  
    return
  end subroutine tstep
  end module method
