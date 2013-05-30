module ArDef
 public
 integer, parameter :: imax = 10000	  
 real(8), parameter :: start = 0.00d0  ! calculation interval
 real(8), parameter :: finish = 1.0d0  ! ....................
	  
 real(8), dimension(imax) :: x_ 
 real(8), dimension(imax) :: NFlux_,DFlux_,ne_,Te_,Ti_,Pe_,Pi_,E_   ! '_' => arrays 
 real(8), dimension(imax) :: NFlux1_
 real(8), dimension(imax) :: DPe_, DPi_, DTe_, DTi_, Dne_
 real(8), dimension(imax) :: ReQ_, ImQ_, ReQ1_,ImQ1_
 !real(8), dimension(imax) :: RePressurePert_, ImPressurePert_, RePressurePert1_, ImPressurePert1_
 real(8), dimension(imax) :: RSol_,ISol_ 
end module ArDef
