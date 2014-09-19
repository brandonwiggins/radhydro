program hydro
 

 use initial
 use global_data
 use advection

 integer :: i, tstep_num, j
 !real :: rho(30)
 real :: vel(30)
 real :: C
 N = 30
 length = 1
 tstep_num = 3
 C = 0.7

 call initial_conditions(N,rho, vel)
 ! do i = 1, N
 !   write(*,*) 'rho = ', rho(i)
 !end do
 
 do i = 1, tstep_num
    call enforce_boundary_conditions(N, rho)
    call get_timestep(length, C, N, dt, vel)
    call interface_states(rho, N, length, dt, vel)
    do j = 1, N
     write(*,*) 'rho(out) = ', rho(j)
    end do
 end do
 


end program hydro
