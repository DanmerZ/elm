subroutine CalculateQ 
 use Profiles
 use ArDef
 implicit none
 

 real(8), dimension(imax) :: Km_, A_
 real(8), dimension(imax) :: znam_
 integer i
 real(8) V01, d1, znam, x 
 

 V01 = 0.0d0 !4500000.0d0
 ReP_(1) = 0.0d0; ImP_(1) = 0.0d0;
 ReQ_(1) = 0.0d0; ImQ_(1) = 0.0d0;
 x = x_(1)
 Km_(1) = (DPi_(1)*Ti_(1)/Pi_(1)/(ap/100.0) - E_(1))/B0/30000.0d0 + Fm(x)*x*ap*V01/(m*aR*c)	
 A_(1) = 80.0*PiNum*m*m*x* (DPe_(1)+DPi_(1))*(1.0d0/q(x)/q(x) - 9.0) /B0/B0

 do i = 2,imax
     x = x_(i)
     Km_(i) = (DPi_(i)*Ti_(i)/Pi_(i)/(ap/100.0) - E_(i))/B0/30000.0d0 + Fm(x)*x*ap*V01/(m*aR*c)	
     A_(i) = 80.0*PiNum*m*m*x* (DPe_(i)+DPi_(i))*(1.0d0/q(x)/q(x) - 9.0) /B0/B0  !80.0 <-- from Pa to erg/cm^-3
     d1 = c/(4.0*PiNum*sigma(NFlux_(i))*x*ap)
     znam = 1.0d0/((m*Km_(i)*(Fm(x))**2)**2 + (A_(i)*d1)**2)
     znam_(i) = znam 
     ReQ_(i) =  znam*Km_(i)*Km_(i)*A_(i)*Fm(x)*Fm(x)*m  
     ImQ_(i) =  znam*Km_(i)*A_(i)*A_(i)*d1  
     ReQ1_(i) = znam*Km_(i)*Km_(i)*Fm(x)*Fm(x)*Fm(x)*m  
     ImQ1_(i) = znam*Km_(i)*Fm(x)*A_(i)*d1 

    ReP_(i) = -x*(DPe_(i)+DPi_(i))*(aR/ap)*m*Km_(i)*Fm(x)*A_(i)*znam*d1
	ImP_(i) = x*(DPe_(i)+DPi_(i))*(aR/ap)*m*m*Km_(i)*Km_(i)*Fm(x)*Fm(x)*Fm(x)*znam
 end do
 
 200 format(A20,E10.2)
 write(*,200) 'an', i1/1./i2
 write(*,*) '*********************************************'
 write(*,200) 'ne, cm-3', 1.d+13*ne_(i1*imax/i2)
 write(*,200) 'Te, erg', 1.6022d-12*Te_(i1*imax/i2)
 write(*,200) 'Ti, erg', 1.6022d-12*Ti_(i1*imax/i2)
 write(*,200) 'P, dyna/cm^2', 10.*(Pe_(i1*imax/i2) + Pi_(i1*imax/i2))
 write(*,200) 'dP, dyna/cm^2', 10.*(DPe_(i1*imax/i2)+DPi_(i1*imax/i2))
 write(*,200) 'ReQ', ReQ_(i1*imax/i2)
 write(*,200) 'ImQ', ImQ_(i1*imax/i2)
 write(*,200) 'A', A_(i1*imax/i2)
 write(*,200) 'Fm', Fm(x_(i1*imax/i2))
 write(*,200) 'Km', Km_(i1*imax/i2)
 write(*,200) 'ReQ1', ReQ1_(i1*imax/i2)
 write(*,200) 'ImQ1', ImQ1_(i1*imax/i2)
 write(*,*) '**********************************************'
 
 call A_Km_ToFiles
		
 contains
		
 subroutine A_Km_ToFiles
  implicit none
  integer i

  open(11, file="data/Km_A_.dat")
  open(12, file="data/znam_.dat")

  do i = 1, imax
	write(11,*) x_(i), Km_(i), A_(i)
	write(12,*) x_(i), ReP_(i),ImP_(i) !znam_(i)
  end do

  close(11)
  close(12)
		
 end subroutine A_Km_ToFiles		
		

end subroutine CalculateQ

	
