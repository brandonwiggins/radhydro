module global_data
  implicit none

  integer :: N, n_ghost, num, index_lo, index_hi
  real :: length, dt, C
  real, dimension(:), allocatable :: rho, vel
end module global_data
