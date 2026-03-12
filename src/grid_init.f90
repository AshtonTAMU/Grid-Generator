module grid_init
  use boundaries
  implicit none

contains

  !Build algebraic grid
  subroutine initialize_grid(x, y, x_phys, eta, NI, NJ)
    integer,  intent(in)  :: NI, NJ
    real(8),  intent(in)  :: x_phys(NI), eta(NJ)
    real(8),  intent(out) :: x(NI,NJ), y(NI,NJ)

    integer :: i, j
    real(8) :: xp, yb, yt

    do i = 1, NI
       xp = x_phys(i)
       yb = y_bottom(xp)
       yt = y_top(xp)
       do j = 1, NJ
          x(i,j) = xp
          y(i,j) = yb + eta(j)*(yt - yb)
       end do
    end do

  end subroutine initialize_grid

end module grid_init
