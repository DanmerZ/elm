subroutine ProfilesToFiles 
 use Profiles
 use ArDef
 implicit none



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
 open(19, file="data/Q1.dat")

do i = 1, imax
    write(11,100) x_(i), NFlux_(i), DFlux_(i)
    write(12,100) x_(i), ne_(i), Fm(x_(i))
    write(13,100) x_(i), Te_(i), DTe_(i) !Te_(i), Ti_(i)
    write(14,100) x_(i), Pe_(i), Pi_(i), Pe_(i) + Pi_(i)
    write(15,100) x_(i), DPe_(i), DPi_(i), DPe_(i) + DPi_(i)
    write(16,100) x_(i), E_(i)
    write(17,100) x_(i), ReQ_(i), ImQ_(i)
    !write(18,100) x_(i), RSol_(i), 0.5*x_(i) !ISol_(i)
    write(19,100) x_(i), ReQ1_(i),ImQ1_(i)
end do

 close(11)
 close(12)
 close(13)
 close(14)
 close(15)
 close(16)
 close(17)
 close(18)
 close(19)
	
end subroutine ProfilesToFiles
