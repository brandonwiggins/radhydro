module initial
  implicit none

contains
  subroutine initial_conditions(N, rho, vel)
    integer :: N, i
    real :: rho(N)
    real, dimension(:) :: vel


!!Top Hat initial conditions
 do i = 1, N/3
    rho(i) = 0.0001
    vel(i) = 1
 end do
 do i = N/3+1, 2*N/3
    rho(i) = 10
    vel(i) = 1
 end do
 do i = 2*N/3+1, N
    rho(i) = 0.0001
    vel(i) = 1
 end do
 do i = 1, N
    write(*,*) 'rho = ', rho(i), ' vel ', vel(i)
 end do

 end subroutine initial_conditions
end module initial
