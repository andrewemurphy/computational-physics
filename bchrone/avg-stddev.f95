subroutine avg_stddev(counter_vals,avg,stddev)
integer, dimension(100) :: counter_vals
integer:: i, total=0 , square_total=0
real :: avg, stddev
do i=1,100
	total = total + counter_vals(i)
	square_total = square_total + counter_vals(i)*counter_vals(i)
end do
avg = real(total)/100.0
stddev = sqrt(real(square_total)/100.0 - avg*avg)
end subroutine


program main
use variation
implicit none
integer :: N, m=1 ,i, j, counter=0
integer, dimension(100):: counter_vals
real :: start, finish, avg, stddev
real(16) :: alpha, alpha_ex, lambda, T
real(16), dimension(2) :: end_point
real(16), dimension(:,:), allocatable :: P

end_point = (/ 1.0_16, 0.0_16 /)
N = 100
lambda = 0.4
alpha= .1
allocate(P(N,2))
do i=1,100
	call cpu_time(start)
	call initialize(end_point,N,P)
	call sine_variation(end_point,N,P,m,lambda)
	counter = 0
	do j=1,10000
		call random_variation(end_point,N,P,alpha, counter)
	end do
	counter_vals(i) = counter
	call calc_action(end_point,N,P,T)
	call cpu_time(finish)
	print *, alpha, T, counter, finish - start
end do

call avg_stddev(counter_vals, avg, stddev)
print *, avg, stddev
end program
