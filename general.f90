module general
use globals
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
    real                :: temp
    real,dimension(neq) :: prim
    real,dimension(1:nx,1:ny) :: rho, vx, vy, pressure
    integer :: i,j
    ! first we calculate the correct data to export
    do i=1,nx
      do j=1,ny
        call uprim(u(:,i,j),prim,temp)
        rho(i,j)=prim(1) ! mass density
        vx(i,j)=prim(2)  ! velocity in x
        vy(i,j)=prim(3)  ! velocity in y
        pressure(i,j)=prim(neq) ! pressure
      end do
    end do
    ! write density
    write(file1,'(a,i2.2,a)') 'rho-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
    do j=1,ny
      write(10,*) (rho(i,j),i=1,nx)
    end do
    ! closes output file
    close(10)
    ! write vx
    write(file1,'(a,i2.2,a)') 'vx-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
    do j=1,ny
      write(10,*) (vx(i,j),i=1,nx)
    end do
    ! closes output file
    close(10)
    ! write vy
    write(file1,'(a,i2.2,a)') 'vy-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
    do j=1,ny
      write(10,*) (vy(i,j),i=1,nx)
    end do
    ! closes output file
    close(10)
    ! write pressure
    write(file1,'(a,i2.2,a)') 'pressure-',itprint,'.dat'
    open(unit=10,file=file1,status='unknown')
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
    use physics
    implicit none
    real, intent(out) ::dt
    !internal variables
    real :: temp,cs,del
    real,dimension(neq) :: prim
    integer :: i,j
    del=1.e+30
    do i=1,nx
      do j=1,ny
        call uprim(u(:,i,j),prim,temp)
        !prim1 = rho, prim2,3,4 = v, prim_neq = P
        !sound = sqrt(gamma*P/rho)
        cs=sqrt(gamma*prim(neq)/prim(1))    
        del=min(del, dx/abs(prim(2)+cs))
        del=min(del, dy/abs(prim(3)+cs))
      enddo
    enddo
    dt=Co*del
    return
  end subroutine timestep
end module general
