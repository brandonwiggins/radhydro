module initial
  implicit none

contains
  subroutine initial_conditions(N, rho, vel, index_lo)
    integer :: N, i, index_lo
    real, dimension(:) :: vel, rho
    write(*,*) 'index_lo ', index_lo
    write(*,*) 'N/3+index_lo', N/3+index_lo
!!Top Hat initial conditions
 do i = index_lo, N/3+index_lo
    rho(i) = 0.0001
    vel(i) = 1
    !write(*,*) 'rho(i) ', i, rho(i), vel(i)
 end do
 do i = N/3+1+index_lo, 2*N/3+index_lo
    rho(i) = 10
    vel(i) = 1
 end do
 do i = 2*N/3+1+index_lo, N+index_lo
    rho(i) = 0.0001
    vel(i) = 1
 end do
 do i = index_lo, N+index_lo
    write(*,*) 'rho = ', rho(i), ' vel ', vel(i)
 end do
    write(*,*) 'done setting up ICs'
 end subroutine initial_conditions
end module initial
