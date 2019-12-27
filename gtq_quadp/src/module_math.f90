module module_math
  implicit none

  private
  public factorial

contains

  recursive real(8) function factorial(n) result(out)
    integer,intent(in) :: n

    if(n.lt.0) then
       stop "aho"
    elseif(n.eq.1.or.n.eq.0) then
       out=1.d0
    else
       out=dble(n)*factorial(n-1)
    end if
  end function factorial
  
end module module_math

