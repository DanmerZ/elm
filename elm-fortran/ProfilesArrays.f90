﻿	subroutine ProfilesArrays(imax,start,finish,x_,NFlux_,ne_,Te_,Ti_,Pe_,Pi_,E_,DPe_,DPi_)
		use Profiles
		implicit none

		integer, intent(in) :: imax
		real(8), intent(in) :: start, finish
		real(8), dimension(imax), intent(in) :: NFlux_, x_
		real(8), dimension(imax), intent(out) :: ne_,Te_,Ti_,Pe_,Pi_,E_,DPe_,DPi_

		integer i
		real(8) x, dx
		real(8) psi1, psi2	
		
		100	format(4E14.5)
		dx = (finish-start)/1.0d0/imax 
		
		DPe_(1) = (Pe(NFlux_(2)) - Pe(NFlux_(1))) / dx 
		DPi_(1) = (Pi(NFlux_(2)) - Pi(NFlux_(1))) / dx

		do  i = 1, imax
					
			x = x_(i)
			
			ne_(i) = ne(NFlux_(i))
			Te_(i) = Te(NFlux_(i))
			Ti_(i) = Ti(NFlux_(i))
			Pe_(i) = Pe(NFlux_(i))
			Pi_(i) = Pi(NFlux_(i))
			if (i > 1) then
				psi1 = NFlux_(i); psi2 = NFlux_(i-1)
				DPe_(i) = (Pe(psi1)-Pe(psi2))/dx;
				DPi_(i) = (Pi(psi1)-Pi(psi2))/dx;				
			end if
			E_(i) = Ef(NFlux_(i))
			
		end do

	end subroutine ProfilesArrays