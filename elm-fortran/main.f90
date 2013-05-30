program elm
 use Profiles
 implicit none

 integer, parameter :: imax = 10000	  
 real(8), parameter :: start = 0.00d0  ! integration interval
 real(8), parameter :: finish = 1.0d0  ! ....................
	  
 real(8), dimension(imax) :: x_ 
 real(8), dimension(imax) :: NFlux_,DFlux_,ne_,Te_,Ti_,Pe_,Pi_,E_   ! '_' => arrays 
 real(8), dimension(imax) :: DPe_, DPi_
 real(8), dimension(imax) :: ReQ_, ImQ_, ReQ1_,ImQ1_
 !real(8), dimension(imax) :: RePressurePert_, ImPressurePert_, RePressurePert1_, ImPressurePert1_
 real(8), dimension(imax) :: znam_, E3_
 real(8), dimension(imax) :: RSol_,ISol_ 
			  		  

 call CalculateFlux(start, finish, imax, x_, NFlux_, DFlux_ )
 call ProfilesArrays(imax,start,finish,x_,NFlux_,ne_,Te_,Ti_,Pe_,Pi_,E_,DPe_,DPi_)
 call CalculateQ(imax,x_,NFlux_,ne_,Te_,Ti_,Pe_,Pi_,E_,DPe_,DPi_,ReQ_,ImQ_,ReQ1_,ImQ1_)		  
 !call SolveEq(start,finish,imax,x_,ReQ_,ImQ_,RSol_,ISol_)
 call ProfilesToFiles(imax, NFlux_, DFlux_,x_,ne_,Te_,Ti_,Pe_,Pi_,E_,DPe_,DPi_,ReQ_,ImQ_,RSol_,ISol_,ReQ1_,ImQ1_)

		 
end program elm
