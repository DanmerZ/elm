module ArDef
 public
 integer, parameter :: imax = 10000	  
 real(8), parameter :: start = 0.00d0  ! calculation interval
 real(8), parameter :: finish = 1.050d0  ! ....................
	  
 real(8), dimension(imax) :: x_ 
 real(8), dimension(imax) :: NFlux_,NFlux1_,DFlux_,ne_,Te_,Ti_,Pe_,Pi_,E_   ! '_' => arrays 
 real(8), dimension(imax) :: DPe_, DPi_, DTe_, DTi_, Dne_
 real(8), dimension(imax) :: ReQ_, ImQ_, ReQ1_,ImQ1_
 real(8), dimension(imax) :: ReP_,ImP_
 !real(8), dimension(imax) :: RSol_,ISol_ 
 
 integer, parameter ::  i1=940, i2=1000
 
end module ArDef
