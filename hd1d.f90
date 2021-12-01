!=======================================================================
!   This program solves the hidrodinamics equation 
!=======================================================================
!   This module contains global variables
module globals
  implicit none
  !
  !   This is the number of points used to discretize X
  !
  integer, parameter :: nx=200
  integer, parameter :: ny=200
  !   This is the number of equation
  integer, parameter :: neq=4
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=1.
  real, parameter :: ymax=1.
  real, parameter :: dx=xmax/real(nx)
  real, parameter :: dy=ymax/real(ny)
  ! The simulation times
  real, parameter :: tmax= 1.             ! maximumn integration time
  real, parameter :: dtprint=0.1          ! interval between outputs
  ! Courant number
  real, parameter :: Co=0.5

  ! simulation constants
  real, parameter :: gamma=5./3.
  real, parameter :: mu=1.4
  ! universal constants
  real, parameter :: boltz=1.38e-16
  real, parameter :: mh=1.67e-24
  !   This is a vector that contains u(x)
  real,dimension(neq,0:nx+1,0:ny+1) :: u,f !U(x,y) & F(x,y)

end module globals
!=======================================================================
!   main program
program hd_1d
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real            :: time, dt             !  t, $\Delta t$
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output

! This subroutine generates the initial conditions
  call initflow(time, tprint, itprint)

  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)

  ! output at tprint intervals
    if(time.ge.tprint) then
      print*,'itprint=',itprint, time,tmax,dt
      call output(itprint)
      tprint=tprint+dtprint
      itprint=itprint+1
    end if

    ! Obtain the $\Delta t$ allowed by the CFL criterium
    call timestep(dt)
    !
    ! Integrate u fom t to t+dt
    call tstep(dt,time)
    ! time counter increases
    time=time+dt

  end do

  stop
end program hd_1d
!
!=======================================================================
! generates initial condition
subroutine initflow(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  !internal variables
  real, parameter :: ux0=0.0, uy0=0., ux1=0.0, uy1=0. 
  real, parameter :: rho0=1., rho1=0.1, p0=1., p1=0.15, xm=0.5
  real :: x
  integer :: i

  !  fill the vector u
  do i=0, nx+1
    x=real(i)*dx   ! obtain the position $x_i$
    if( x < xm )  then
      u(1,i,:)=rho0
      u(2,i,:)=ux0
      u(3,i,:)=uy0
      u(neq,i,:)=p0
    else
      u(1,i,:)=rho1
      u(2,i,:)=ux1
      u(3,i,:)=uy1
      u(neq,i,:)=p1
    end if

    if( (x-0.5*dx <= xm).and.(x+0.5*dx >= xm) ) then
      u(1,i,:)=(rho0+rho1)/2.
      u(2,i,:)=(ux0+ux1)/2.
      u(3,i,:)=(uy0+uy1)/2.
      u(neq,i,:)=(p0+p1)/2.
    end if
  end do

  print*,u(1,1) !??????????????????????????????????????????????

  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0

  return
end subroutine initflow

!=======================================================================
! output to file
subroutine output(itprint)
  use globals
  implicit none
  integer, intent(in) :: itprint
  !internal variables
  character (len=20) file1
  real                :: temp
  real,dimension(neq) :: prim
  integer :: i

  ! open output file
  write(file1,'(a,i2.2,a)') 'hd-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')

  ! writes x and u
  do i=1,nx
    call uprim(u(:,i),prim,temp)                        !Revisar la subrutina uprim primero
    write(10,*) real(i)*dx,prim(1),prim(2),prim(3)
  end do

  ! closes output file
  close(10)

  return
end subroutine output

!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine timestep(dt)
  use globals
  implicit none
  real, intent(out) ::dt
  !internal variables
  real :: temp,cs,csound,del
  real,dimension(neq) :: prim
  integer :: i
  !
  del=1.e+30
  do i=1,nx
    do j=1,ny
      call uprim(u(:,i,j),prim,temp)
      cs=csound(prim(1),prim(3))
      del=min(del, dx/abs(prim(2)+cs), dy/abs(prim(3)+cs))
    enddo
  enddo
  dt=Co*del
  return
end subroutine timestep
!
function csound(n_g,press)
  use globals, only : gamma
  implicit none
  real, intent(in) :: n_g,press
! function declaration
  real             :: csound
  !
  csound=sqrt(gamma*press/n_g)
  !
  return
end function
!
subroutine uprim(uu,prim,temp)
  use globals
  implicit none
  real, dimension(neq), intent(in)  :: uu
  real, dimension(neq), intent(out) :: prim
  !internal variables
  real                 :: ek, et
  real, intent(out)    :: temp
  !

  prim(1)=uu(1)
  prim(2)=uu(2)/prim(1)
  prim(3)=uu(3)/prim(1)
  ek=0.5*prim(1)*(prim(2)+prim(3))**2.
  et=uu(neq)-ek
  prim(neq)=et/(gamma-1.)
  temp=prim(neq)/(prim(1)*boltz/(mu*mh))
  !
  return
end subroutine uprim
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  implicit none
  real,            intent(in) :: dt, time
  !internal variables
  real, dimension(neq,0:nx+1) :: up
  real                        :: dtx
  integer                     :: i

  !  obtain the fluxes
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

!=======================================================================
! Obtain the fluxes F
subroutine fluxes(u,f)
  use globals, only :neq,nx,ny,gamma
  implicit none
  real,dimension(neq,0:nx+1,0:ny+1),intent(in) :: u
  real,dimension(neq,0:nx+1,0:ny+1),intent(out) :: f
  !internal variables
  real, dimension(neq) :: prim
  integer :: i
  real :: temp,etot

   do i=0,nx+1
    do j=0,ny+1
      call uprim(u(:,i,j),prim,temp)
      Etot=0.5*prim(1)*(prim(2)+prim(3))**2.+prim(neq)/(gamma-1.)
      f(1,i,j)=prim(1)*prim(2)
      !f(2,i,j)=prim(1)*prim(2)**2.+prim(3)
      f(2,i,j)=prim(1)*prim(2)**2.+prim(1)*prim(2)*prim(3)
      f(3,i,j)=prim(1)*prim(3)**2.+prim(1)*prim(2)*prim(3)
      f(neq,i,j)=(prim(2)+prim(3))*(etot+prim(3))
    enddo
!    print*,u(1,i),prim(1),f(1,i),prim(2)
!    print*,prim(1),f(2,i),prim(2),prim(3)
  enddo

  return
end subroutine fluxes

!=======================================================================
! Set boundary conditions
subroutine boundaries(u)
  use globals, only : nx,ny,neq
  implicit none
  real,dimension(neq,0:nx+1,0:ny+1), intent(inout) :: u
  ! free outflow (salida libre)
  !x
  u(:,0,:)=u(:,1,:)
  u(:,nx+1,:)=u(:,nx,:)
  !y
  u(:,:,0)=u(:,:,1)
  u(:,:,ny+1)=u(:,:,ny)
  return
end subroutine boundaries
!=======================================================================
