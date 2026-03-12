program grid
  use boundaries
  use grid_init
  use smoothing
  use io_grid
  implicit none

  integer, parameter :: NI       = 51
  integer, parameter :: NJ       = 11
  integer, parameter :: max_iter = 20000
  real(8), parameter :: tol      = 1.0d-7

  real(8) :: x(NI,NJ), y(NI,NJ)
  real(8) :: x_phys(NI), eta(NJ)
  integer :: i, j

  !Build coordinate arrays
  do i = 1, NI
     x_phys(i) = (i-1) * 5.0d0 / (NI-1)
  end do
  do j = 1, NJ
     eta(j) = (j-1) * 1.0d0 / (NJ-1)
  end do

  !Algebraic initialization
  call initialize_grid(x, y, x_phys, eta, NI, NJ)

  !Laplace smoothing
  call smooth_grid(x, y, x_phys, eta, NI, NJ, max_iter, tol)

  !Write output
  call write_grid(x, y, NI, NJ, "grid.dat")

end program grid
