module physics
use globals
use mesh
implicit none
contains
  !========================================================================
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
    ! prim2 = v = rho*v/rho = u(2)/u(1)
    prim(2)=uu(2)/prim(1)
    ! Kinetic energy = 1/2 rho*v²
    ek=0.5*prim(1)*prim(2)**2.
    ! Energy = uu3 = ek + et = ek + P/(gamma -1) 
    ! et = P/(gamma-1) = Energy - ek = uu3 - ek
    et=uu(neq)-ek
    ! prim_neq = P = et*(gamma-1)
    prim(neq)=et*(gamma-1.)
    ! P = n*Kb*T -> T = P/n*Kb & n = rho/mu*mh numerical density
    temp=prim(neq)/(prim(1)*boltz/(mu*mh))
    ! prim1 = rho, prim2,3,4 = v, prim_neq = P
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
    ! prim1 = rho, prim2,3,4 = v, prim_neq = P
     do i=0,nx+1
      call uprim(u(:,i),prim,temp)
      ! Etot= Ek + P/(gamma-1)
      Etot=0.5*prim(1)*prim(2)**2.+prim(neq)/(gamma-1.)
      ! flux1 = rho*v
      f(1,i)=prim(1)*prim(2)
      ! flux2,3,4 = rho*v*v + P
      f(2,i)=prim(1)*prim(2)**2.+prim(neq)
      ! flux_neq = v*(Etot + P)
      f(neq,i)=prim(2)*(etot+prim(neq))
    enddo
  
    return
  end subroutine fluxes
  !=======================================================================
  ! Obtain the sources
  subroutine sources(u,i,ss)
    use globals, only :neq,dx,gamma
    implicit none
    real,dimension(neq,0:nx+1),intent(in) :: u
    integer, intent(in) :: i
    real,dimension(neq),intent(out) :: ss
    ! Internal variables
    real :: R, alfa, term
    real, dimension(neq) :: prim
    real :: temp,etot
    alfa=2.
    R=float(i)*dx
    term=alfa/R
    call uprim(u(:,i),prim,temp)
    ! prim1 = rho, prim2,3,4 = v, prim_neq = P
    ! Etot= Ek + P/(gamma-1)
    Etot=0.5*prim(1)*prim(2)**2.+prim(neq)/(gamma-1.)
    ! source1 = eta*rho*v/x
    ss(1)=term*prim(1)*prim(2)
    ! source2 = eta*rho*v²/x
    ss(2)=term*prim(1)*prim(2)**2.
    ! source_neq = eta*v*(Etot+P)/x
    ss(neq)=term*prim(2)*(etot+prim(neq))
    return
  end subroutine sources
  !=======================================================================
  ! Set boundary conditions
  subroutine boundaries(u)
    use globals, only : nx,neq
    implicit none
    real,dimension(neq,0:nx+1), intent(inout) :: u
    u(:,0)=u(:,1)
    u(2,0)=-u(2,1)
    u(:,nx+1)=u(:,nx)
    return
  end subroutine boundaries
  !=======================================================================
end module physics
