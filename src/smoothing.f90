module smoothing
  use boundaries
  implicit none

contains

  !Laplace smoothing of algebraic grid
  subroutine smooth_grid(x, y, x_phys, eta, NI, NJ, max_iter, tol)
    integer, intent(in)    :: NI, NJ, max_iter
    real(8), intent(in)    :: tol
    real(8), intent(in)    :: x_phys(NI), eta(NJ)
    real(8), intent(inout) :: x(NI,NJ), y(NI,NJ)

    real(8) :: x_old(NI,NJ), y_old(NI,NJ)
    real(8) :: x_xi, x_eta, y_xi, y_eta
    real(8) :: alpha, beta_, gamma_
    real(8) :: x_c, y_c, denom
    real(8) :: Tx, Ty, Tmag, Nx, Ny, ds
    real(8) :: t, xw, fw, dfw
    real(8) :: res
    integer :: i, j, iter, k, i_min
          
    do iter = 1, max_iter

       x_old = x
       y_old = y

       !Interior Laplace update
       do i = 2, NI-1
          do j = 2, NJ-1

             x_xi  = (x_old(i+1,j) - x_old(i-1,j)) / 2.0d0
             x_eta = (x_old(i,j+1) - x_old(i,j-1)) / 2.0d0
             y_xi  = (y_old(i+1,j) - y_old(i-1,j)) / 2.0d0
             y_eta = (y_old(i,j+1) - y_old(i,j-1)) / 2.0d0

             !Coefficients, lagged
             alpha  = x_eta**2 + y_eta**2
             beta_  = x_xi*x_eta + y_xi*y_eta 
             gamma_ = x_xi**2 + y_xi**2

             x_c = x_old(i+1,j+1) - x_old(i+1,j-1) - x_old(i-1,j+1) + x_old(i-1,j-1)
             y_c = y_old(i+1,j+1) - y_old(i+1,j-1) - y_old(i-1,j+1) + y_old(i-1,j-1)

             denom = 2.0d0*(alpha + gamma_)

             x(i,j) = (alpha*(x(i+1,j)+x(i-1,j)) - 0.5d0*beta_*x_c + &
                  gamma_*(x(i,j+1)+x(i,j-1))) / denom

             y(i,j) = (alpha*(y(i+1,j)+y(i-1,j)) - 0.5d0*beta_*y_c + &
                  gamma_*(y(i,j+1)+y(i,j-1))) / denom

          end do
       end do

       !Project j=1 onto bottom wall along normal from j=2
       do i = 2, NI-1

          !Tangent
          Tx = (x(i+1,2) - x(i-1,2)) / 2.0d0
          Ty = (y(i+1,2) - y(i-1,2)) / 2.0d0
          Tmag = sqrt(Tx*Tx + Ty*Ty)
          Nx = -Ty/Tmag;  Ny = Tx/Tmag
          if (Ny < 0.0d0) then; Nx = -Nx; Ny = -Ny; end if

          !Newton, distance to bottom
          t = sqrt((x(i,2)-x(i,1))**2 + (y(i,2)-y(i,1))**2)
          do k = 1, 20
             xw  = x(i,2) - t*Nx
             fw  = y(i,2) - t*Ny - y_bottom(xw)
             dfw = -Ny + Nx*dy_bottom(xw)
             if (abs(dfw) < 1.0d-14) exit
             t   = t - fw/dfw
          end do

          x(i,1) = x(i,2) - t*Nx
          y(i,1) = y_bottom(x(i,1))

       end do

       !Project j=NJ onto top wall along normal from j=NJ-1
       do i = 2, NI-1

          !Tangent
          Tx = (x(i+1,NJ-1) - x(i-1,NJ-1)) / 2.0d0
          Ty = (y(i+1,NJ-1) - y(i-1,NJ-1)) / 2.0d0
          Tmag = sqrt(Tx*Tx + Ty*Ty)
          Nx = -Ty/Tmag;  Ny = Tx/Tmag
          if (Ny > 0.0d0) then; Nx = -Nx; Ny = -Ny; end if

          !Newton, distance to top
          t = sqrt((x(i,NJ-1)-x(i,NJ))**2 + (y(i,NJ-1)-y(i,NJ))**2)
          do k = 1, 20
             xw  = x(i,NJ-1) - t*Nx
             fw  = y(i,NJ-1) - t*Ny - y_top(xw)
             dfw = -Ny + Nx*dy_top(xw)
             if (abs(dfw) < 1.0d-14) exit
             t   = t - fw/dfw
          end do

          x(i,NJ) = x(i,NJ-1) - t*Nx
          y(i,NJ) = y_top(x(i,NJ))

       end do

       !Snap closest bottom wall point to x=2 and x=3
       i_min = 2
       do i = 3, NI-1
          if (abs(x(i,1)-2.0d0) < abs(x(i_min,1)-2.0d0)) i_min = i
       end do
       x(i_min,1) = 2.0d0;  y(i_min,1) = y_bottom(2.0d0)

       i_min = 2
       do i = 3, NI-1
          if (abs(x(i,1)-3.0d0) < abs(x(i_min,1)-3.0d0)) i_min = i
       end do
       x(i_min,1) = 3.0d0;  y(i_min,1) = y_bottom(3.0d0)

       !Snap closest top wall point to x=2 and x=3
       i_min = 2
       do i = 3, NI-1
          if (abs(x(i,NJ)-2.0d0) < abs(x(i_min,NJ)-2.0d0)) i_min = i
       end do
       x(i_min,NJ) = 2.0d0;  y(i_min,NJ) = y_top(2.0d0)

       i_min = 2
       do i = 3, NI-1
          if (abs(x(i,NJ)-3.0d0) < abs(x(i_min,NJ)-3.0d0)) i_min = i
       end do
       x(i_min,NJ) = 3.0d0;  y(i_min,NJ) = y_top(3.0d0)
       
       !Inlet / outlet, keep lines horizontal
       do j = 2, NJ-1
          ds = sqrt((x(2,j)-x(1,j))**2 + (y(2,j)-y(1,j))**2)
          x(2,j) = x(1,j) + ds;  y(2,j) = y(1,j)

          ds = sqrt((x(NI-1,j)-x(NI,j))**2 + (y(NI-1,j)-y(NI,j))**2)
          x(NI-1,j) = x(NI,j) - ds;  y(NI-1,j) = y(NI,j)
       end do

       !Fix inlet/outlet columns
       do j = 1, NJ
          x(1,j)  = 0.0d0;  y(1,j)  = eta(j)
          x(NI,j) = 5.0d0;  y(NI,j) = eta(j)
       end do      

       res = maxval(abs(x - x_old)) + maxval(abs(y - y_old))
       if (res < tol) then
          print *, "Converged at iteration", iter, " res =", res
          exit
       end if

    end do

  end subroutine smooth_grid
end module smoothing
