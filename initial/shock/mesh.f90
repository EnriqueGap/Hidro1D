module mesh
implicit none
!   This is the number of points used to discretize X
integer, parameter :: nx=200
!   Here we set the extent of X and calculate $\Delta x$
real, parameter :: xmax=1
real, parameter :: dx=xmax/real(nx)
end module mesh
