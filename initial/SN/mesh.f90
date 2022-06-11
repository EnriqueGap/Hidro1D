module mesh
implicit none
!   This is the number of points used to discretize X
integer, parameter :: nx=1000
!   Here we set the extent of X and calculate $\Delta x$
real, parameter :: xmax=5e7 !1000 Km
real, parameter :: dx=xmax/real(nx)
end module mesh