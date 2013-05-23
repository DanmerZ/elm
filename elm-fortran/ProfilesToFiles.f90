	subroutine ProfilesToFiles(imax, NFlux_, DFlux_, x_,ne_,Te_,Ti_,Pe_,Pi_,E_,DPe_,DPi_,ReQ_,ImQ_,RSol_,ISol_)
		use Profiles
		implicit none

		integer, intent(in) :: imax
		real(8), dimension(imax), intent(in) :: NFlux_,DFlux_,x_,ne_,Te_,Ti_,Pe_,Pi_,E_,DPe_,DPi_,ReQ_,ImQ_,RSol_,ISol_

		integer :: i

		100	format(4E14.5)	
		open(11, file="data/Flux.dat")
		open(12, file="data/ne.dat")		
		open(13, file="data/T.dat")		
		open(14, file="data/P.dat")
		open(15, file="data/DP.dat")
		open(16, file="data/E.dat")
		open(17, file="data/Q.dat")
		open(18, file="data/Sol.dat")

		do i = 1, imax
			write(11,100) x_(i), NFlux_(i), DFlux_(i)
			write(12,100) x_(i), ne_(i)
			write(13,100) x_(i), Te_(i), Ti_(i)
			write(14,100) x_(i), Pe_(i), Pi_(i), Pe_(i) + Pi_(i)
			write(15,100) x_(i), DPe_(i), DPi_(i), DPe_(i) + DPi_(i)
			write(16,100) x_(i), E_(i)
			write(17,100) x_(i), ReQ_(i), ImQ_(i)
			write(18,100) x_(i), RSol_(i), 0.5*x_(i)*(x_(i)+1.) !ISol_(i)
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)
		close(17)
		close(18)
	

	end subroutine ProfilesToFiles
