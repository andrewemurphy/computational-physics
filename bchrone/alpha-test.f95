program main
use variation
implicit none
integer :: N, m, j, counter=0
real :: start, finish
real(16) :: alpha, alpha_ex, lambda, T
real(16), dimension(2) :: end_point
real(16), dimension(:,:), allocatable :: P

end_point = (/ 1.0_16, 0.0_16 /)
N = 100
lambda = 0.401

allocate(P(N,2))

do alpha_ex = 0.0_16, 10.0_16, 0.1_16
	call cpu_time(start)
	alpha = exp(-alpha_ex)
	call initialize(end_point,N,P)
	call sine_variation(end_point,N,P,m,lambda)
	counter = 0
	do j=1,1000
		call random_variation(end_point,N,P,alpha, counter)
	end do
	call calc_action(end_point,N,P,T)
	call cpu_time(finish)
	print *, alpha, T, counter, finish - start
end do

end program