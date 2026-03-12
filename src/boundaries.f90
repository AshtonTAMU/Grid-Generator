module boundaries
  implicit none
  real(8), parameter :: pi = 4.0d0*atan(1.0d0)

contains

  !Bottom sine curve, floor
  real(8) function y_bottom(xp)
    real(8), intent(in) :: xp
    if (xp > 2.0d0 .and. xp < 3.0d0) then
       y_bottom = 0.17d0*sin((xp-2.0d0)*pi)
    else
       y_bottom = 0.0d0
    end if
  end function y_bottom

  !Top sine curve, roof
  real(8) function y_top(xp)
    real(8), intent(in) :: xp
    if (xp > 2.0d0 .and. xp < 3.0d0) then
       y_top = 1.0d0 - 0.17d0*sin((xp-2.0d0)*pi)
    else
       y_top = 1.0d0
    end if
  end function y_top

  !Derivative of floor, for normal
  real(8) function dy_bottom(xp)
    real(8), intent(in) :: xp
    if (xp >= 2.0d0 .and. xp <= 3.0d0) then
       dy_bottom = 0.17d0*pi*cos((xp-2.0d0)*pi)
    else
       dy_bottom = 0.0d0
    end if
  end function dy_bottom

  !Derivative of roof, for normal
  real(8) function dy_top(xp)
    real(8), intent(in) :: xp
    if (xp >= 2.0d0 .and. xp <= 3.0d0) then
       dy_top = -0.17d0*pi*cos((xp-2.0d0)*pi)
    else
       dy_top = 0.0d0
    end if
  end function dy_top

end module boundaries
