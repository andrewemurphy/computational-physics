!------------------------------------------------------------
! INTRODUCTION
!------------------------------------------------------------

! This program simulates Rutherford Scattering in 3 Dimensions.

!Written by Andrew Murphy, 2020



!----------------------------------------------------------
!	SUBROUTINES
!----------------------------------------------------------
!This subroutine represents Hamilton's Equations for a repulsive Kepler-Force.
!'pc' is the 6-dimesional phase-space coordinate (x, y, z, px, py, pz) of the particle
!'pq-dot' contains  the time-derivatives of the phase-space coordinates, given by Hamilton's Equations (dH/dx = -dpx/dt; dH/dpx = dx/dt)
subroutine Hamiltons_eqs(pc, pq_dot, k, m)
implicit none
real(16) :: k, m
real(16), dimension(6):: pc, pq_dot
pq_dot(1) = pc(4)/m
pq_dot(2) = pc(5)/m
pq_dot(3) = pc(6)/m
pq_dot(4) = k*pc(1)/((pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))**1.5_16) 
pq_dot(5) = k*pc(2)/((pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))**1.5_16) 
pq_dot(6) = k*pc(3)/((pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))**1.5_16) 
end subroutine Hamiltons_eqs

subroutine scattering_process(pc_0, k, m, , N_steps, dt, theta)
implicit none
integer :: N_steps, i
real(16):: k, m, dt, theta
real(16), dimension(6):: pc, pq_dot, pc_0
pc = pc_0
do i =1, N_steps
	call Hamiltons_eqs(pc, pq_dot, k, m)
	pc = pc + pq_dot*dt

!---------------------------------------------------
!	MAIN CODE
!---------------------------------------------------

program rscattering
implicit none
real(16), parameter :: k=100.0, m=1.0, s_max =1.0, z_0 = -100.0_16, p_0=100.0_16, pi=4*atan(1.0_16)
integer, parameter:: N_particles=1000 ,N_steps=1000
!our parameters include: 'k' = force-constant (ZZ'e^2 for nuclear scattering), 'm' = mass
!'s_max' = beam-radius (maximum impact parameter), 'z_0' = initial z-coordinate
!'p_0' = initial z-momentum, and 'pi' . Note: since pi is only used to generate random initial coordinate, a precise value is not too neccessary
!'N_particles' = number of particles to be simulated, 'N_steps' = number of time-steps in our Euler-Method Algortihm

real(16):: s_0, phi_0, H, dt
real(16), dimension(6):: pc_0, pc, pq_dot
real(16), dimension(N_particles) :: s_vals , theta_vals
integer :: i, j

dt = abs((4*m*p_0/(z_0))/N_steps) !our maximum time is abs(4m*p_0/z_0), which is ""sufficeiently large"" to capture the scattering dynamics 

do i=1,N_particles
	call random_seed()
	call random_number(s_0)
	call random_number(phi_0)
	s_0 = s_max * s_0										!Generate random impact parameter
	phi_0 = 2*pi *phi_0											!Generate random azimuthal angle 
	pc_0 = (/ s_0*cos(phi_0), s_0*sin(phi_0), z_0, 0.0_16, 0.0_16, p_0 /)	!create initial phase-space coordinate
	pc = pc_0											!record initial impact-parameter
	
	do j=1,N_steps								!Begin Euler-Method Algortihm for approximating the trajectory 
		call Hamiltons_eqs(pc, pq_dot, k, m)	!Call subroutine to generate phase-space differentials
		pc = pc + pq_dot*dt						!Increment the phase-space coordinate of the particle
	end do
	
	theta_vals(i) = acos(pc(3)/sqrt(pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3)))
	s_vals(i) = s_0

end do	

open(10, file="data\theta_vals.txt", status='replace')	!list of scattering angles
open(11, file="data\s_vals.txt", status='replace')		!list of initial impact parameters


do i=1, N_particles				!Write data to files
	write(10,*) theta_vals(i)
	write(11,*) s_vals(i)


end do
end program rscattering