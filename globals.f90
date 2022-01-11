!   This module contains global variables
module globals
    implicit none
    !
    !   This is the number of points used to discretize X
    !
    integer, parameter :: nx=200
    !   This is the number of equation
    integer, parameter :: neq=3
    !   Here we set the extent of X and calculate $\Delta x$
    real, parameter :: xmax=1.
    real, parameter :: dx=xmax/real(nx)
    ! The simulation times
    real, parameter :: tmax= 1.             ! maximumn integration time
    real, parameter :: dtprint=0.1          ! interval between outputs
    ! Courant number
    real, parameter :: Co=0.8
  
    ! simulation constants
    real, parameter :: gamma=5./3.
    real, parameter :: mu=1.4
    ! universal constants
    real, parameter :: boltz=1.38e-16
    real, parameter :: mh=1.67e-24
    !   This is a vector that contains u(x)
    real,dimension(neq,0:nx+1) :: u,f
  
  end module globals