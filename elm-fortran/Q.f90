subroutine CalculateQ 
 use Profiles
 use ArDef
 implicit none
 

 real(8), dimension(imax) :: Km_, A_
 real(8), dimension(imax) :: znam_
 integer i
 real(8) V01, d1, znam, x

 V01 = 0.0d0 !4500000.0d0
 ReQ_(1) = 0.0d0; ImQ_(1) = 0.0d0;
 x = x_(1)
 Km_(1) = (DPi_(1)*Ti_(1)/Pi_(1)/(ap/100.0) - E_(1))/B0/30000.0d0 + Fm(x)*x*ap*V01/(m*aR*c)	
 A_(1) = 80.0*PiNum*m*m*x* (DPe_(1)+DPi_(1))*(1.0d0/q1(x)/q1(x) - 8.0) /B0/B0

 do i = 2,imax
     x = x_(i)
     Km_(i) = (DPi_(i)*Ti_(i)/Pi_(i)/(ap/100.0) - E_(i))/B0/30000.0d0 + Fm(x)*x*ap*V01/(m*aR*c)	
     A_(i) = 80.0*PiNum*m*m*x* (DPe_(i)+DPi_(i))*(1.0d0/q1(x)/q1(x) - 8.0) /B0/B0
     d1 = c/(4.0*PiNum*sigma(NFlux_(i))*x*ap)
     znam = (1.0d0/((m*Km_(i)*(Fm(x))**2)**2 + (A_(i)*d1)**2))
     znam_(i) = znam 
     ReQ_(i) = znam*m*Km_(i)*Km_(i)*A_(i)*Fm(x)*Fm(x)  
     ImQ_(i) = znam*Km_(i)*A_(i)*A_(i)*d1  
     ReQ1_(i) = znam*m*Km_(i)*Km_(i)*Fm(x)*Fm(x)*Fm(x)  
     ImQ1_(i) = znam*Km_(i)*Fm(x)*A_(i)*d1 
 end do

 call A_Km_ToFiles
		
 contains
		
 subroutine A_Km_ToFiles
  implicit none
  integer i

  open(11, file="data/Km_A_.dat")
  open(12, file="data/znam_.dat")

  do i = 1, imax
	write(11,*) x_(i), Km_(i), A_(i)
	write(12,*) x_(i), znam_(i)
  end do

  close(11)
  close(12)
		
 end subroutine A_Km_ToFiles		
		

end subroutine CalculateQ
