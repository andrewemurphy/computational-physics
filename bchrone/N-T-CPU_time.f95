!Andrew Murphy, 2020

!uses subroutines in variation.f95 to measure the CPU time of allocation, initialization, sinusoidal variation, and action calculation,
!while the numbner of points in our path ('N') increses.

program main
use variation
implicit none
integer :: N, m =1
real :: start, finish
REAL(16):: T, lambda 
REAL(16), dimension(2)::end_point
REAL(16), dimension(:,:), allocatable :: P


lambda = 0.401
end_point = (/1.0_16, 0.0_16/)

do N = 1000, 100000, 1000
	call cpu_time(start)
	allocate(P(N,2))
	call initialize(end_point,N,P)
	call sine_variation(end_point,N,P,m,lambda)
	call calc_action(end_point,N,P,T)
	deallocate(P)
	call cpu_time(finish)
	print *, N, T, finish - start
	
end do
end program
