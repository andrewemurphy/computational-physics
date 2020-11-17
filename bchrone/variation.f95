!Andrew Murphy, 2020

!Module with subroutines pertaining to the Brachistochrone Problem
!The objective is to reduce the variational problem of minimizing the integral sqrt((1+y'**2)/y)dx over some path,
!where the paths are limited to straight line segments between a set of N - (x,y) cooridnates ('P') 
!along with the two endpoints, the origin and 'end_point'. The action integral, is thereby reduced to a sum of the
! 'differential action' calculated along each of the straight line segments. 

module variation
 implicit none
 
contains
!---------------------------------------------------------------------------------
! 1. Subroutines pertaining to the generation of paths:
!---------------------------------------------------------------------------------

subroutine initialize(end_point,N,P)
! creates a starting path of points, spaced equidistant along the straight line from the origin to end_point
! a, b are the start and end points (respecitvely), P is the array of points between the endpoints
implicit none
integer :: i, N ! 'i' = dummy index, and 'N' = number of intermediate points
REAL(16), dimension(2):: end_point ! 'end_point' = (x,y) coordinate of terminal point
REAL(16), dimension(N,2) :: P ! 'P' = Nx2 array of intermediate points

do i=1,N 
	P(i,1) = i*(end_point(1)/(N+1))
	P(i,2) = i*(end_point(2)/(N+1))
end do

end

!----------------------------------------------------------------------------------

subroutine sine_variation(end_point,N,P,m,lambda)
!introduces sinusoidal variation of order m to an existing path
implicit none
integer :: i , N , m ! 'i'= dummy variable, 'N' = number of points, 'm' = order of sinulsoidal variation
REAL(16) :: L, lambda, pi=4*atan(1.0_16) ! 'L' = charateristic length, 'lambda'= coefficent of variation
REAL(16), dimension(2) :: end_point !end_point of path
REAL(16), dimension(N,2) :: P ! set of path points
L = end_point(1)

do i=1,N
	P(i,1) = P(i,1)
	P(i,2) = P(i,2) + lambda*sin((m*pi/L)*P(i,1))
end do
end subroutine

!-----------------------------------------------------------------------------------

subroutine cycloid(R,t_max, N, P)
!Produces a cycloid path of N points ('P'), parametrically, given 'R' and 't_max'
implicit none
integer :: i, N !'i' =dummy variable, 'N' = number of intermediate points
Real(16)::R, t_max , t  ! 'R' = radius of generated cricle, 
real(16), dimension(N,2) :: P

do i = 1
	t = real(i)/real(N + 1)*t_max
	P(i,1) = R*(t - sin(t))
	P(i,2) = R*(1 - cos(t))
end do
end

!-----------------------------------------------------------------------------------
! 2. Subroutines pertaining to the calculation of the action along a path:
!------------------------------------------------------------------------------------ 

subroutine differential_action(a,b,dT)
!calculates the differential action (dT) over a straight-line segment from a to b, given a start from the origin, 
! takes inputs 'a' , 'b' as the endpoints of the straight line segment, and returns 'dT' the differential action along the segment 
REAL(16) ::  dT, M ! 'dT' = differential action, 'M' = slope of the line segment
REAL(16), dimension(2) :: a, b !two end points
M = (b(2) - a(2))/(b(1)-a(1)) !calculate slope
if (M == 0) then
	dT = abs((b(1) -a(1))/sqrt(a(2))) !case if slope =0
else
	dT = sqrt(1.0 + (1.0/M)**2)*abs(sqrt(a(2)) - sqrt(b(2))) !see notes for derivation of this formula
end if

end subroutine

!----------------------------------------------------------------------------------------

subroutine calc_action(end_point,N, P, T)
!Calculates action over entire path, starting at (0,0) and ending at end_point, with intermediate points 'P'
implicit none
integer :: i, N
REAL(16) :: T, dT
REAL(16), dimension(2) :: end_point, p1, p2 
REAL(16), dimension(N,2) :: P

T = 0.0 !set T to 0

!Calculate action at starting endpoint
p1 = (/ 0, 0 /)
p2 = (/ P(1,1), P(1,2) /)
call differential_action(p1,p2,dT)
T = T + dT

!Calculate action at final endpoint
p1 = (/ P(N,1) , P(N,2) /)
p2 = end_point
call differential_action(p1,p2,dT)
T = T + dT

!calculate action over intermediate points
do i =2, N
	p1 = (/ P(i-1,1) , P(i-1,2) /)
	p2 = (/ P(i,1) , P(i,2) /)
	call differential_action(p1,p2,dT)
	T = T + dT
end do
end subroutine

!------------------------------------------------------------------------------------------------
! 3. Subroutines pertaining to random variation of the path points
!------------------------------------------------------------------------------------------------

subroutine random_variation(end_point,N,P, alpha, counter)
implicit none
integer :: N, random_index, counter
real(16):: alpha, r, T , Tnew
real(16), dimension(2) :: dr , end_point
real(16), dimension(N,2) :: P, Pnew
Pnew = P
do
	call random_seed()
	call random_number(r)
	random_index = 1 + floor(N*r)
	call random_seed()
	call random_number(dr)
	dr = alpha*2*(dr - (/0.5 , 0.5/) )
	Pnew(random_index,1) = Pnew(random_index,1) + dr(1)
	Pnew(random_index,1) = Pnew(random_index,2) + dr(2)
	if (Pnew(random_index,2) > 0) then
		exit
	end if
end do
call calc_action(end_point,N,P,T)
call calc_action(end_point,N,Pnew,Tnew)
if (Tnew < T) then
	P = Pnew
	counter = counter + 1
end if
end subroutine

end module