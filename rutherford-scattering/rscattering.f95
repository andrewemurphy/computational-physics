module rscattering
implicit none

contains
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
end subroutine

subroutine initialize(s_max, z_0, p_0, pc_0)
real(16) :: s_max, z_0, p_0, s, phi_0, pi = 4*atan(1.0_16)
real(16), dimension(6) :: pc_0
call random_seed()
call random_number(s)
call random_number(phi)
pc_0(1) = s*s_0*cos(2*pi*phi)
pc_0(2) = s*s_0*sin(2*pi*phi)
pc_0(3) = z_0
pc_0(4) = 0.0_16
pc_0(5) = 0.0_16
pc_0(6) = p_0
end subroutine

subroutine scattering_process(pc_0, k, m, , N_steps, dt, theta)
implicit none
integer :: N_steps, i
real(16):: k, m, dt, theta
real(16), dimension(6):: pc, pq_dot, pc_0
pc = pc_0
do i =1, N_steps
	call Hamiltons_eqs(pc, pq_dot, k, m)
	pc = pc + pq_dot*dt
end do
theta = acos(pc(3)/sqrt(pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3)))
end subroutine

end module
