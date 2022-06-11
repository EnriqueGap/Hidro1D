module method
use globals
use physics
use mesh
implicit none
contains
!=======================================================================
  ! integration from t to t+dt with the method of MacCormack
  subroutine tstep(dt,time)
    use globals
    use physics
    implicit none
    real,            intent(in) :: dt, time
    ! internal variables
    ! auxiliar
    real, dimension(neq,0:nx+1) :: up, upp
    !sources
    real, dimension(neq) :: ss
    ! delta t/delta x & artificial viscosity
    real                        :: dtx, eta
    integer                     :: i
    ! eta in (0 , 0.5)
    eta=0.16
    ! obtain the fluxes
    call fluxes(u,f)
    ! MacCormack Method is a two step method. The predictor step pushes forward
    ! We save the information in upp
    dtx=dt/dx
    !	Predictor Step 
    do i=1,nx
      upp(:,i) = u(:,i) - dtx*(f(:,i+1)-f(:,i)) + eta*(u(:,i-1) + u(:,i+1) - 2*u(:,i))
      call sources(u,i,ss)
      upp(:,i)=upp(:,i)-dt*ss(1:neq)
    end do
    call boundaries(upp)
    call fluxes(upp,f)
    ! The corrector step corrects the previous step
    ! We save the information in up
    !	Corrector Step
    do i=1,nx
      up(:,i)=0.5*(u(:,i)+upp(:,i) - dtx*(f(:,i)-f(:,i-1))) + eta*(upp(:,i-1) + upp(:,i+1) - 2*upp(:,i))
      call sources(upp,i,ss)
      up(:,i)=up(:,i)-dt*ss(1:neq)
    end do
    ! Boundary conditions to the U^n+1
    call boundaries(up)
    ! copy the up to the u
    u(:,:)=up(:,:)
    return
  end subroutine tstep
end module method
