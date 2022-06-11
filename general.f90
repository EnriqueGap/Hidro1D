module general
use globals
use mesh
use physics
implicit none
contains
  !=======================================================================
  ! output to file
  subroutine output(itprint)
    use globals
    use physics
    implicit none
    integer, intent(in) :: itprint
    !internal variables
    character (len=20) file1
    real                :: temp, dist
    real,dimension(neq) :: prim
    integer :: i
    ! open output file
    write(file1,'(a,i2.2,a)') 'hd-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
    ! writes x and u
    do i=1,nx
      call uprim(u(:,i),prim,temp)
      ! x(Km), rho, v, T
      dist = real(i)*dx/Km
      write(10,*) dist,prim(1),prim(2),temp
    end do
    ! closes output file
    close(10)
    return
  end subroutine output
  !=======================================================================
  ! computes the timestep allowed by the CFL criterium
  subroutine timestep(dt)
    use globals
    use physics
    implicit none
    real, intent(out) ::dt
    !internal variables
    real :: temp,cs,csound,del
    real,dimension(neq) :: prim
    integer :: i
    ! Starts here
    del=1.e+30
    do i=1,nx
      call uprim(u(:,i),prim,temp)
      !prim1 = rho, prim2,3,4 = v, prim_neq = P
      !sound = sqrt(gamma*P/rho)
      cs=sqrt(gamma*prim(neq)/prim(1))    
      del=min(del,dx/abs(prim(2)+cs))
    enddo
    dt=Co*del
    return
  end subroutine timestep
end module general
