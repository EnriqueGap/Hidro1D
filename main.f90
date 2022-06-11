!=======================================================================
!   This program solves the hidrodinamics equation 
!=======================================================================
!   main program
program hd_1d
  use globals
  use initial
  use mesh
  use general
  use method
  use physics
  implicit none
  ! declaration of some variables needed by the main program
  real            :: time, dt             !  t, $Delta t$
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output

! This subroutine generates the initial conditions
  call initflow(time, tprint, itprint)

  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)				!time+dt<=tmax

  ! output at tprint intervals
    if(time.ge.tprint) then				!time<=tprint
      call output(itprint)
      tprint=tprint+dtprint
      itprint=itprint+1
    end if

    ! Obtain the $Delta t$ allowed by the CFL criterium
    call timestep(dt)
    print*,'itprint=',itprint, time,tmax,dt
    ! Integrate u fom t to t+dt
    call tstep(dt,time)
    ! time counter increases
    time=time+dt
    !print*,'ifout=',itprint, time,tmax,dt

  end do

  stop
end program hd_1d
