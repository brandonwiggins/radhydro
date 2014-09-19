module advection
  implicit none

contains
  subroutine enforce_boundary_conditions(N, rho)
    integer :: N
    real, dimension(:) :: rho
    real :: rho_lo1, rho_lo2, rho_high1, rho_high2

    !rho(1) = rho(N) !Enforcing periodic boundary conditions

    !Enforcing periodic boundary conditions with 2 ghost cells
    
    rho_lo1 = rho(N)
    rho_lo2 = rho(N-1)

    rho_high1 = rho(1)
    rho_high2 = rho(2)
    
    
  end subroutine enforce_boundary_conditions

!================================================
  subroutine get_timestep(length, C, N, dt, vel)
    real :: length, C, dt
    real, dimension(:) :: vel
    integer :: N
    !dt = 0.0001
    dt = C*length/N*1/maxval(vel)
    write(*,*) 'dt ', dt

  end subroutine get_timestep
!================================================
  subroutine interface_states(rho, N, length, dt, vel)
    real :: length, dt, dadx1, dadx2, rho_lo1, rho_lo2, rho_high1, rho_high2, rho_old
    real, dimension(:) :: rho, vel
    real :: a_l(N+1), a_r(N+1), a(N+1)
    integer :: i, N, j
    
    dadx1 = 0.5*(rho(1) - rho_lo2)*length/N
    a_l(1) = rho(1) + length/N*0.5*(1-dt*N/length*vel(1))*dadx1

    dadx1 = 0.5*(rho(2) - rho_lo1)*length/N
    a_l(2) = rho(2) + length/N*0.5*(1-dt*N/length*vel(2))*dadx1

    dadx2 = 0.5*(rho(3)- rho(1))*length/N
    a_r(1) = rho(2) - length/N*0.5*(1+dt*N/length*vel(1))*dadx2

    dadx2 = 0.5*(rho(4) - rho(2))*length/N
    a_r(2) = rho(3) - length/N*0.5*(1+dt*N/length*vel(2))*dadx2


!!Ghost cell action with right boundary

    dadx1 = 0.5*(rho(1) - rho_lo2)*length/N
    a_l(N+1) = rho(1) + length/N*0.5*(1-dt*N/length*vel(1))*dadx1

    dadx1 = 0.5*(rho(2) - rho_lo1)*length/N
    a_l(N) = rho(2) + length/N*0.5*(1-dt*N/length*vel(2))*dadx1

    dadx2 = 0.5*(rho_high2- rho(N))*length/N
    a_r(N+1) = rho_high1 - length/N*0.5*(1+dt*N/length*vel(N))*dadx2

    dadx2 = 0.5*(rho_high1 - rho(N-1))*length/N
    a_r(N) = rho(N) - length/N*0.5*(1+dt*N/length*vel(N-1))*dadx2
    

    !! Calculating Interface States for
    !dadx2 = 0.5*(rho(3) - rho(1))*length/N
    !a(1) =  rho(2) - length/N*0.5*(1+dt*N/length*vel(1))*dadx2

    !dadx1 = 0.5*(rho(N) - rho(N-2))*length/N
    !a(N-1) = rho(N-1) + length/N*0.5*(1-dt*N/length*vel(N-1))*dadx1
    
    do i = 3,N-1
       dadx1 = 0.5*(rho(i+1) - rho(i-1))*length/N
       a_l(i) = rho(i) + length/N*0.5*(1-dt*N/length*vel(i))*dadx1
       !write(*,*) 'a_l(i)' ,i, a_l(i)
       dadx2 = 0.5*(rho(i+2) - rho(i))*length/N
       a_r(i) = rho(i+1) - length/N*0.5*(1+dt*N/length*vel(i+1))*dadx2
       !write(*,*) 'a_r(i)', i, a_r(i)

       !! Solve the Riemann Problem


    end do


    !do i = 1, N
    !   write(*,*) 'a(i) ',i, a(i)
    !end do
    
    do i = 1, N
       if (vel(i) > 0) then
           a(i) = a_l(i)
           !write(*,*) 'a_l(i)',i,a_l(i)
           !write(*,*) 'a(i)', i, a(i)
       else
           a(i) = a_r(i)
       end if
       
    end do
    !! Update quantities
    
    !rho(1) = rho(1) - a(1)*N/length*dt
    !rho(N) = rho(N) + a(N-1)*N/length*dt
    do i = 1, N
       rho_old = rho(i)
       !write(*,*) 'rho(i)' , i , rho(i)
       rho(i) = rho_old - (a(i+1) - a(i))*N*dt/length
       !write(*,*) 'a - a, last_term, rho(i), rho(next) ' , i, a(i+1) - a(i), (a(i+1) - a(i))*N*dt/length, rho(i)
       !write(*,*) 'last_term' ,i, (a(i+1) - a(i))*N*dt/length
       !write(*,*) 'rho is' , i, rho(i)
    end do
     !do j = 1, N
     !write(*,*) 'rho = ', rho(i)
    !end do
  end subroutine interface_states
end module advection
