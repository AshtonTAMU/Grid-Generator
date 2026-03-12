module io_grid
  implicit none

contains

  !Write grid data to dat file
  subroutine write_grid(x, y, NI, NJ, filename)
    integer,          intent(in) :: NI, NJ
    real(8),          intent(in) :: x(NI,NJ), y(NI,NJ)
    character(len=*), intent(in) :: filename

    integer :: i, j

    open(unit=10, file=filename, status="replace")
    write(10,*) 'VARIABLES = "X","Y"'
    write(10,*) 'ZONE I=', NI, ', J=', NJ, ', F=POINT'

    do j = 1, NJ
       do i = 1, NI
          write(10,'(2E20.12)') x(i,j), y(i,j)
       end do
    end do

    close(10)
    print *, "Grid written to ", filename

  end subroutine write_grid

end module io_grid
