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
    real, parameter :: ux0=0.0, uy0=0., ux1=0.0, uy1=0.                   !(vx,vy)=v0 for x<xm & (vx,vy)=v1 for x>xm
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
        u(neq,i,:)=p0/(gamma-1.)
      else
        u(1,i,:)=rho1
        u(2,i,:)=ux1
        u(3,i,:)=uy1
        u(neq,i,:)=p1/(gamma-1.)
      end if
  
      if( (x-0.5*dx <= xm).and.(x+0.5*dx >= xm) ) then
        u(1,i,:)=(rho0+rho1)/2.
        u(2,i,:)=(ux0+ux1)/2.
        u(3,i,:)=(uy0+uy1)/2.
        u(neq,i,:)=(p0+p1)/(2.*(gamma-1.))
      end if
    end do
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
    real,dimension(1:nx,1:ny) :: rho, vx, vy, pressure
    integer :: i,j
  
    ! open output file
    write(file1,'(a,i2.2,a)') 'rho-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
  
    do i=1,nx
      do j=1,ny
        call uprim(u(:,i,j),prim,temp)
      rho(i,j)=prim(1)
      end do
    end do
    do j=1,ny
      write(10,*) (rho(i,j),i=1,nx)
    end do
  
    ! closes output file
    close(10)
  
  ! open output file
    write(file1,'(a,i2.2,a)') 'vx-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
  
    do i=1,nx
      do j=1,ny
        call uprim(u(:,i,j),prim,temp)
      vx(i,j)=prim(2)
      end do
    end do
    do j=1,ny
      write(10,*) (vx(i,j),i=1,nx)
    end do
  
    ! closes output file
    close(10)
  
  ! open output file
    write(file1,'(a,i2.2,a)') 'vy-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
  
    do i=1,nx
      do j=1,ny
        call uprim(u(:,i,j),prim,temp)
      vy(i,j)=prim(3)
      end do
    end do
    do j=1,ny
      write(10,*) (vy(i,j),i=1,nx)
    end do
    ! closes output file
    close(10)
  
  ! open output file
    write(file1,'(a,i2.2,a)') 'pressure-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
  
    do i=1,nx
      do j=1,ny
        call uprim(u(:,i,j),prim,temp)
      pressure(i,j)=prim(neq)
      end do
    end do
    do j=1,ny
      write(10,*) (pressure(i,j),i=1,nx)
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
    integer :: i,j
    !
    del=1.e+30
    do i=1,nx
      do j=1,ny
        call uprim(u(:,i,j),prim,temp)
        cs=sqrt(gamma*prim(neq)/prim(1))    
        del=min(del, dx/abs(prim(2)+cs))
        del=min(del, dy/abs(prim(3)+cs))
      enddo
    enddo
    dt=Co*del
    return
  end subroutine timestep
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
    ek=0.5*prim(1)*(prim(2)**2. +prim(3)**2.)
    et=uu(neq)-ek
    prim(neq)=et*(gamma-1.)
    temp=prim(neq)/(prim(1)*boltz/(mu*mh))
    !
    return
  end subroutine uprim
  !=======================================================================
  ! Obtain the fluxes F
  subroutine fluxes(u,f,g)
    use globals, only :neq,nx,ny,gamma
    implicit none
    real,dimension(neq,0:nx+1,0:ny+1),intent(in) :: u
    real,dimension(neq,0:nx+1,0:ny+1),intent(out) :: f
    real,dimension(neq,0:nx+1,0:ny+1),intent(out) :: g
    !internal variables
    real, dimension(neq) :: prim
    integer :: i,j
    real :: temp,etot
  
    do i=0,nx+1
      do j=0,ny+1
        call uprim(u(:,i,j),prim,temp)
        Etot=0.5*prim(1)*(prim(2)**2. + prim(3)**2.) + prim(neq)/(gamma-1.)
        f(1,i,j)=prim(1)*prim(2)
        f(2,i,j)=prim(1)*prim(2)**2.+prim(neq)
        f(3,i,j)=prim(1)*prim(2)*prim(3)
        f(neq,i,j)=prim(2)*(etot+prim(3))
      enddo
    enddo
  
    do i=0,nx+1
      do j=0,ny+1
        call uprim(u(:,i,j),prim,temp)
        Etot=0.5*prim(1)*(prim(2)**2. + prim(3)**2.) + prim(neq)/(gamma-1.)
        g(1,i,j)=prim(1)*prim(3)
        g(3,i,j)=prim(1)*prim(3)**2.+prim(neq)
        g(2,i,j)=prim(1)*prim(2)*prim(3)
        g(neq,i,j)=prim(3)*(etot+prim(3))
      enddo
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
  end module general
