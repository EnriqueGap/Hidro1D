module general
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
        u(1,i)=rho1
        u(2,i)=u1
        u(3,i)=p1/(gamma-1.)
      else
        u(1,i)=rho0
        u(2,i)=u0
        u(3,i)=p0/(gamma-1.)
      end if
  
      if( (x-0.5*dx <= x1).and.(x+0.5*dx >= x1) ) then
        u(1,i)=(rho0+rho1)/2.
        u(2,i)=(u0+u1)/2.
        u(3,i)=(p0+p1)/(2.*(gamma-1.))
      end if
    end do
  
    print*,u(1,1)
  
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
    write(file1,'(a,i2.2,a)') 'vishd-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
  
    ! writes x and u
    do i=1,nx
      call uprim(u(:,i),prim,temp)
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
      call uprim(u(:,i),prim,temp)
      cs=sqrt(gamma*prim(3)/prim(1))    
      del=min(del,dx/abs(prim(2)+cs))
    enddo
    dt=Co*del
    return
  end subroutine timestep
  !========================================================================
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
    ek=0.5*prim(1)*prim(2)**2.
    et=uu(3)-ek
    prim(3)=et*(gamma-1.)
    temp=prim(3)/(prim(1)*boltz/(mu*mh))
    !
    return
  end subroutine uprim
  !=======================================================================
  ! Obtain the fluxes F
  subroutine fluxes(u,f)
    use globals, only :neq,nx,gamma
    implicit none
    real,dimension(neq,0:nx+1),intent(in) :: u
    real,dimension(neq,0:nx+1),intent(out) :: f
    !internal variables
    real, dimension(neq) :: prim
    integer :: i
    real :: temp,etot
  
     do i=0,nx+1
      call uprim(u(:,i),prim,temp)
      Etot=0.5*prim(1)*prim(2)**2.+prim(3)/(gamma-1.)
      f(1,i)=prim(1)*prim(2)
      f(2,i)=prim(1)*prim(2)**2.+prim(3)
      f(3,i)=prim(2)*(etot+prim(3))
  !    print*,u(1,i),prim(1),f(1,i),prim(2)
  !    print*,prim(1),f(2,i),prim(2),prim(3)
    enddo
  
    return
  end subroutine fluxes
  
  !=======================================================================
  ! Set boundary conditions
  subroutine boundaries(u)
    use globals, only : nx,neq
    implicit none
    real,dimension(neq,0:nx+1), intent(inout) :: u
    ! free outflow (salida libre)
    u(:,0)=u(:,1)
    u(:,nx+1)=u(:,nx)
    return
  end subroutine boundaries
  !=======================================================================
  end module general
