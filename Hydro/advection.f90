module advection
  implicit none

contains
  subroutine enforce_boundary_conditions(rho, n_ghost, index_lo, index_hi)
    integer :: N, n_ghost, index_lo, index_hi, m
    real, dimension(:) :: rho

    !Enforcing periodic boundary conditions with 2 ghost cells
   
    do m = 1, n_ghost
       rho(m) = rho(index_hi-n_ghost + m)
       rho(index_hi + m) = rho(index_lo + m -1)
    end do
    
  end subroutine enforce_boundary_conditions

!================================================
  subroutine get_timestep(length, C, N, dt, vel)
    real :: length, C, dt
    real, dimension(:) :: vel
    integer :: N

    dt = C*length/N*1/maxval(vel)
    write(*,*) 'dt ', dt

  end subroutine get_timestep
!================================================








!================================================
  subroutine interface_states(rho, N, length, dt, vel, index_lo, index_hi)
    real :: length, dt, dadx1, dadx2, rho_old
    real, dimension(:) :: rho, vel
    real :: a_l(N+1), a_r(N+1), a(N+1)
    integer :: i, N, j, index_lo, index_hi
  
    do i = index_lo-1, index_hi
       dadx1 = 0.5*(rho(i+1) - rho(i-1))*length/N
       a_l(i) = rho(i) + length/N*0.5*(1-dt*N/length*vel(i))*dadx1
       dadx2 = 0.5*(rho(i+2) - rho(i))*length/N
       a_r(i) = rho(i+1) - length/N*0.5*(1+dt*N/length*vel(i+1))*dadx2
    end do

    !! Riemann Problem

    do i = index_lo-1, index_hi
       if (vel(i) > 0) then
           a(i) = a_l(i)
       else
           a(i) = a_r(i)
       end if
    end do

    !! Update quantities
 
    do i = index_lo, index_hi
       rho_old = rho(i)
       rho(i) = rho_old - (a(i) - a(i-1))*N*dt/length
    
    end do
    
    do j = index_lo, index_hi
     write(*,*) 'rho(inside) = ', j, rho(j)
    end do


  end subroutine interface_states

!==============================================================
  subroutine ppm_states(rho, N, length, dt, vel, index_lo, index_hi, n_ghost)
    real :: length, dt, dadx1, dadx2, rho_old, delta_rho1, delta_rho2, dmaj, flux_l, flux_r, thing1, thing2
    real :: al_sub, ar_sub
    real, dimension(:) :: rho, vel
    real :: a_l(N+1), a_r(N+1), a(N+1), a_from_plus(N+1), a_from_minus(N+1)
    integer :: i, N, j, index_lo, index_hi, n_ghost
  

    !C&W eq. 1.7 for a uniform grid
    do i = index_lo - n_ghost +1, index_hi + n_ghost - 1
       delta_rho1 = 0.5*(rho(i+1) - rho(i))
       delta_rho2 = 0.5*(rho(i) - rho(i-1))

       if ((rho(i+1) - rho(i))*(rho(i) - rho(i-1)) < 0) then
          delta_rho1 = sign(1.0, delta_rho1)*min(abs(delta_rho1), 2*abs(rho(i+1) - rho(i)), 2*abs(rho(i) - rho(i-1)))
       else 
          delta_rho1 = 0
       end if

       if ((rho(i+1) - rho(i))*(rho(i) - rho(i-1)) < 0) then
          delta_rho2 = sign(1.0, delta_rho2)*min(abs(delta_rho2), 2*abs(rho(i+1) - rho(i)), 2*abs(rho(i) - rho(i-1)))
       else 
          delta_rho2  = 0
       end if
       !! here we fit a quartic (cubic?) spline using a 5 cell stensil to get a_{j+1/2}. We use the same stencil to predict the a_{j+3/2} value as well which is stored in the a_from_minus value.

       a_from_plus(i) = 0.5*(rho(i) + rho(i+1)) - (1/6)*(delta_rho2-delta_rho1)
       a_from_minus(i+1) = a_from_plus(i)
    end do
    
    do j = index_lo - n_ghost + 2, index_hi + n_ghost

          thing1 = (a_from_plus(j) - a_from_minus(j))**2/6
          thing2 = (a_from_plus(j) - a_from_minus(j) )* (rho(j) - 0.5*(a_from_minus(j) + a_from_plus(j)))
          if ((a_from_plus(j)-rho(j))*(rho(j) - a_from_minus(j)).le.0) then
             a_from_minus(j) = rho(j)
             a_from_plus(j) = rho(j)
          else if ((a_from_plus(j) - a_from_minus(j) )* (rho(j) - 0.5*(a_from_minus(j) + a_from_plus(j))).gt.thing1) then 
             a_from_minus(j) = 3*rho(j) - 2*a_from_plus(j)
          else if ( -(a_from_plus(j) - a_from_minus(j))**2/6.gt.thing2) then
             a_from_plus(j) = 3*rho(j) - 2*a_from_minus(j)
          end if

    end do

    do j = index_lo, index_hi+2
       al_sub = (1.0-2.0/3.0*vel(j-1)*(dt*N/length))*6.0*(rho(j-1)-0.5*(a_from_minus(j-1)+a_from_plus(j-1)))
       ar_sub = (1.0-2.0/3.0*vel(j)  *(dt*N/length))*6.0*(rho(j) - 0.5*(a_from_minus(j) + a_from_plus(j)))
       a_l(j) = a_from_plus(j-1) - 0.5*vel(j-1)*(dt*N/length)*((a_from_plus(j-1) - a_from_minus(j-1))- al_sub)
       a_r(j) = a_from_minus(j) -  0.5*vel(j)  *(dt*N/length)*((a_from_plus(j)   - a_from_minus(j)) + ar_sub)
    end do
    
    !!Then do Riemann Problem like before

    do i = index_lo-1, index_hi
       if (vel(i) > 0) then
           a(i) = vel(i)*a_l(i)
       else
           a(i) = vel(i)*a_r(i)
       end if
    end do

    !! Update quantities
 
    do i = index_lo, index_hi
       rho_old = rho(i)
       rho(i) = rho_old - (a(i-1) - a(i))*N*dt/length
    
    end do
    
    do j = index_lo, index_hi
     write(*,*) 'rho(inside) = ', j, rho(j)
    end do




 end subroutine ppm_states

end module advection
