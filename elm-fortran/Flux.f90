﻿subroutine CalculateFlux 
 use Profiles
 use ArDef
 implicit none

 integer i
 real(8) x, dx, s, x1, psi1
		
 dx = (finish-start)/1.0d0/imax
		
 psi1 = 0.0d0
 if (psi1 == 0.0d0) then
	s = 0.0d0
	x1 = 0.0d0
	do while (x1 < 1.0d0)
		s = s + dx*x1/q(x1)
		x1 = x1 + dx			
	enddo
		psi1 = s		 
 end if	
		 
 s = 0.0d0
 x1 = start
 i = 1
 do while ( i <= imax) 
	s = s + dx*x1/q(x1)
	x_(i) = x1
	NFlux_(i) = s/psi1			
	i = i + 1
	x1 = x1 + dx			
 end do
		
 x1 = start
 i = 1
 DFlux_(1) = (NFlux_(2) - NFlux_(1))/dx
		
 do while ( i <= imax) 
	if (i > 1) then
		DFlux_(i) = (NFlux_(i) - NFlux_(i-1))/dx
				
	end if
	x1 = x1 + dx
	i = i + 1
 end do	


end subroutine CalculateFlux
