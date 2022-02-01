module physics
use globals
implicit none
contains
!========================================================
  subroutine uprim(uu,prim,temp)
    use globals
    implicit none
    real, dimension(neq), intent(in)  :: uu
    real, dimension(neq), intent(out) :: prim
    !internal variables
    real                 :: ek, et
    real, intent(out)    :: temp
    ! prim1 = rho = u1
    prim(1)=uu(1)
    ! prim2 = vx = rho*vx/rho = u(2)/u(1)
    prim(2)=uu(2)/prim(1)
    ! prim3 = vy = rho*vy/rho = u(3)/u(1)
    prim(3)=uu(3)/prim(1)
    ! Kinetic energy = 1/2 rho * (vx^2 + vy^2)
    ek=0.5*prim(1)*(prim(2)**2. +prim(3)**2.)
    ! Et = Energy - Ek
    et=uu(neq)-ek
    ! Pressurre = Et*(gamma-1)
    prim(neq)=et*(gamma-1.)
    ! P =nKbT -> T = P/nKb & n= rho/mu*mh numerical density
    temp=prim(neq)/(prim(1)*boltz/(mu*mh))
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
    ! prim1 = rho, prim2,3,4 = v, prim_neq = P
    do i=0,nx+1
      do j=0,ny+1
        call uprim(u(:,i,j),prim,temp)
        ! Etot = Ek+ Pressurre/(gamma-1)
        Etot=0.5*prim(1)*(prim(2)**2. + prim(3)**2.) + prim(neq)/(gamma-1.)
        ! Fluxes with vx
        ! rho*vx
        f(1,i,j)=prim(1)*prim(2)
        ! rho*vx^2 + P
        f(2,i,j)=prim(1)*prim(2)**2.+prim(neq)
        ! rho*vx*vy
        f(3,i,j)=prim(1)*prim(2)*prim(3)
        ! vx*(Etot + P)
        f(neq,i,j)=prim(2)*(etot+prim(neq))
        ! Fluxes with vy
        ! rho*vy
        g(1,i,j)=prim(1)*prim(3)
        ! rho*vx*vy
        g(2,i,j)=prim(1)*prim(2)*prim(3)
        ! rho*vy^2 + P
        g(3,i,j)=prim(1)*prim(3)**2.+prim(neq)
        ! vy*(Etot + P)
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
end module physics
