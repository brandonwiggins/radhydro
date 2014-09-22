program hydro
 

 use initial
 use global_data
 use advection

 integer :: i, tstep_num, j, ad_type

 ad_type = 1
 num = 99
 n_ghost = 2
 N = num+2*n_ghost
 
 index_lo = n_ghost + 1
 index_hi = num+n_ghost

 allocate(rho(N))
 allocate(vel(N))
 length = 1
 tstep_num = 3
 C = 0.7

 call initial_conditions(num,rho, vel, index_lo)
 
 
 do i = 1, tstep_num
    call enforce_boundary_conditions(rho, n_ghost, index_lo, index_hi)
    call get_timestep(length, C, num, dt, vel)
    if (ad_type.eq.1) call interface_states(rho, num, length, dt, vel, index_lo, index_hi)
    if (ad_type.eq.2) call ppm_states(rho, num, length, dt, vel, index_lo, index_hi, n_ghost)
 end do
 
deallocate(rho)
deallocate(vel)


end program hydro
