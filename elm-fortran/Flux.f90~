subroutine CalculateFlux
 use Profiles
 use ArDef
 implicit none
 

 integer i
 real(8) x, dx, s, x1, psi1, psi2
 real(8) s1
		
 dx = (finish-start)/1.0d0/imax
		
 psi1 = 0.0d0; psi2 = 0.0
 !if (psi1 == 0.0d0) then
 s = 0.0d0; s1 = 0.0d0;
 x1 = start
 do while (x1 < 1.0d0)
	s = s + dx*x1/q(x1)
	s1 = s1 + dx*x1/q1(x1)
	x1 = x1 + dx			
 enddo
 psi1 = s; psi2 = s1		 
 !end if	
		 
 s = 0.0d0; s1 = 0.0d0
 x1 = start
 
 do i = 1,imax 
	s = s + dx*x1/q(x1)
	s1 = s1 + dx*x1/q1(x1)
	x_(i) = x1
	NFlux_(i) = s/psi1	
	NFlux1_(i) = s1/psi2		
	x1 = x1 + dx			
 end do
		
 x1 = start
 
 DFlux_(1) = (NFlux_(2) - NFlux_(1))/dx
		
 do i = 2,imax 	
	DFlux_(i) = (NFlux_(i) - NFlux_(i-1))/dx
	x1 = x1 + dx	
 end do	


end subroutine CalculateFlux
