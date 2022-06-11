!   This module contains global variables
module globals
    use mesh
    implicit none
    !   This is the number of equations
    integer, parameter :: neq=3
    ! The simulation times
    real, parameter :: tmax= 1             ! maximumn integration time
    real, parameter :: dtprint=0.1!0.005          ! interval between outputs
    ! Courant number
    real, parameter :: Co=0.3
    ! simulation constants
    ! adiabatic index
    real, parameter :: gamma=5./3.
    ! mean weight
    real, parameter :: eta=0.
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
    !   This is a vector that contains u(x)
    real,dimension(neq,0:nx+1) :: u,f
  end module globals
