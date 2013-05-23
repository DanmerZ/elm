     subroutine SolveEq(start,finish,imax,x_,ReQ_,ImQ_,RSol_,ISol_)
      use Profiles
      implicit none
      integer, intent(in) :: imax	  
	  real(8), intent(in) :: start  
	  real(8), intent(in) :: finish   
	  real(8), dimension(imax), intent(in) :: ReQ_, ImQ_
	  real(8), dimension(imax), intent(in) :: x_ 
	  real(8), dimension(imax), intent(inout) :: RSol_,ISol_  
	  	  
	  real(8) :: u1p,u1n,u2p,u2n, v1p,v1n,v2p,v2n  !u - real y,  v - imaginary y
	  real(8) :: sh1, sh2  !shooting method's start values
	  integer :: i
	  real(8) :: dx, right
	  dx = -(x_(2) - x_(1))
	  
	  !Initial conditions: RSol(0)=u1=0, Rsol'(0)=u2=sh1
	  !                    ISol(0)=v1=0, ISol'(0)=v2=sh2
	  sh1 = 1.0d+0; sh2 = -1.2d-2
	  
	  u1p = 1.0d0; u2p = sh1
	  v1p = 1.0d0; v2p = sh2
	  
	  do i = imax,1,-1
	      RSol_(i) = u1p; ISol_(i) = v1p
	      right = 1. !-3.*u2p/x_(i) + m*u1p*(m-1.)/(x_(i))**2 
	      u2n = u2p + dx*right
	      u1n = u1p + dx*u2p
	      right = -3.*v2p/x_(i) + m*v1p*(m-1.)/(x_(i))**2 
	      v2n = v2p + dx*right
	      v1n = v1p + dx*v2p
	      
	      u1p = u1n; u2p = u2n
	      v1p = v2n; v2p = v1n	      
	  end do
	  	  
	  contains
	  
	  function euler(f,y,dx)
	   implicit none
	   real(8), intent(in) :: f,y,dx
	   real(8) :: euler
	   euler = y + f*dx	   	   
	  end function euler
      
     end subroutine SolveEq
