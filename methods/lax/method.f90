module method
use globals
use physics
implicit none
contains
!=======================================================================
  ! integration from t to t+dt with the method of Lax
  subroutine tstep(dt,time)
    use globals
    use physics
    implicit none
    real,            intent(in) :: dt, time
    !internal variables
    ! auxiliar
    real, dimension(neq,0:nx+1) :: up
    ! sources
    real, dimension(neq) :: ss
    ! delta t/ delta x
    real                        :: dtx
    integer                     :: i
    !  obtain the fluxes
    call fluxes(u,f)
    !   Here is the Lax method, notice that the values at the extremes can
    !   not be calculated, we need to enter then as boundary conditions
    dtx=dt/dx
    do i=1,nx
      up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i)) 
      call sources(u,i,ss)
!      up(:,i)=up(:,i)-dt*ss(1:neq)
    end do 
    !  Boundary conditions to the U^n+1
    call boundaries(up)
    ! copy the up to the u
    u(:,:)=up(:,:)  
    return
  end subroutine tstep
end module method
