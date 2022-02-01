!   This module contains global variables
module globals
    implicit none
    ! This is the number of points used to discretize X
    integer, parameter :: nx=200
    integer, parameter :: ny=200
    ! This is the number of equation
    integer, parameter :: neq=4
    ! Here we set the extent of X and calculate $\Delta x$
    real, parameter :: xmax=1.
    real, parameter :: ymax=1.
    real, parameter :: dx=xmax/real(nx)
    real, parameter :: dy=ymax/real(ny)
    ! The simulation times
    real, parameter :: tmax= 1.             ! maximumn integration time
    real, parameter :: dtprint=0.1          ! interval between outputs
    ! Courant number
    real, parameter :: Co=0.8
    ! simulation constants
    real, parameter :: gamma=5./3.
    real, parameter :: mu=1.4
    ! universal constants
    !boltzman in cgs
    real, parameter :: boltz=1.38e-16
    ! hydrogen mass
    real, parameter :: mh=1.67e-24
    ! pi
    real, parameter :: pi=3.141592653589
    ! Solar Mass
    real, parameter :: SM=1.988435e33
    ! Kilometers
    real, parameter :: Km=1e5
    ! This is a vector that contains u(x)
    real,dimension(neq,0:nx+1,0:ny+1) :: u,f,g !U(x,y) & F(x,y)
  end module globals