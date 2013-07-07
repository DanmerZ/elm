 subroutine Evans60
  use Profiles
  use ArDef
  implicit none
  integer i
  real(8) :: dx,x1,x,d1,znam,sigma1
  real(8), dimension(imax) :: ne60_,Te60_,Ti60_,Pe60_,Pi60_,E60_,DPe60_,DPi60_
  real(8), dimension(imax) :: Km60_, A60_, znam60_
  real(8), dimension(imax) :: ReQ,ImQ,ReP,ImP,ReQ1,ImQ1
  dx = (finish-start)/1.0d0/imax
  x1=start
  100	format(4E14.5)
  open(11,file="data/evans60/ne60.dat")
  open(12,file="data/evans60/T60.dat")  
  open(13,file="data/evans60/P60.dat")
  open(14,file="data/evans60/E60.dat")

  do i=1,imax
      x_(i)=x1
	  ne60_(i)=ne60(x1)
	  Te60_(i)=Te60(x1)
	  Ti60_(i)=Ti60(x1)
	  Pe60_(i)=1.6*ne60_(i)*Te60_(i)
	  Pi60_(i)=1.6*ne60_(i)*Ti60_(i)
      E60_(i)=Ef60(x1)
	  
	  write(11,100) x_(i),ne60_(i) 
	  write(12,100) x_(i),Te60_(i),Ti60_(i)	  
	  write(13,100) x_(i),Pe60_(i),Pi60_(i),Pe60_(i)+Pi60_(i)
	  write(14,100) x_(i),E60_(i)
	   
	  x1=x1+dx
  end do
  close(11)
  close(12)
  close(13)
  close(14)


  DPe60_(1) = (Pe60_(2) - Pe60_(1)) / dx 
  DPi60_(1) = (Pi60_(2) - Pi60_(1)) / dx   

  open(15,file="data/evans60/dP60.dat")   
  open(16,file="data/evans60/K60.dat")
  do i = 2,imax    
    DPe60_(i) = (Pe60_(i) - Pe60_(i-1)) / dx 
    DPi60_(i) = (Pi60_(i) - Pi60_(i-1)) / dx
	Km60_(i) = (DPi60_(i)*Ti60_(i)/Pi60_(i)/(ap/100.0) - E60_(i))/B0/30000.0d0 !+ Fm(x)*x*ap*V01/(m*aR*c) 
    
	write(15,100) x_(i),DPe60_(i),DPi60_(i),DPe60_(i)+DPi60_(i)
    write(16,100) x_(i),Km60_(i)
  end do
  close(15)
  close(16)

  ReQ(1) = 0.0d0; ImQ(1) = 0.0d0;
  ReQ1(1) = 0.0d0; ImQ1(1) = 0.0d0;
  
  open(17,file="data/evans60/ppert.dat")
  open(18,file="data/evans60/Q.dat")
  open(19,file="data/evans60/A.dat")
  open(20,file="data/evans60/Q1.dat")
  do i=2,imax
    x=x_(i)
    A60_(i) = 80.0*PiNum*m*m*x*(DPe60_(i)+DPi60_(i))*(1.0d0/q(x)/q(x) - 9.0)/B0/B0
    sigma1 = 1.2d+17 * (Te60(x)/500.0d0)**1.5
    d1 = c/(4.0*PiNum*sigma1*x*ap)
    znam=1.0d0/((m*Km60_(i)*(Fm(x))**2)**2 + (A60_(i)*d1)**2)
    ReQ(i) = znam*m*Km60_(i)*Km60_(i)*A60_(i)*Fm(x)*Fm(x)  
    ImQ(i) = znam*Km60_(i)*A60_(i)*A60_(i)*d1

    ReQ1(i) = znam*m*Km60_(i)*Km60_(i)*Fm(x)*Fm(x)*Fm(x)
    ImQ1(i) = znam*Km60_(i)*Fm(x)*A60_(i)*d1

    ReP(i) = -x*(DPe60_(i)+DPi60_(i))*(aR/ap)*m*Km60_(i)*Fm(x)*A60_(i)*znam*d1
	ImP(i) = x*(DPe60_(i)+DPi60_(i))*(aR/ap)*m*m*Km60_(i)*Km60_(i)*Fm(x)*Fm(x)*Fm(x)*znam

	write(17,100) x,ReP(i),ImP(i)
	write(18,100) x,ReQ(i),ImQ(i)
	write(19,100) x, A60_(i)
	write(20,100) x,ReQ1(i),ImQ1(i)

  end do
  close(17)
  close(18)
  close(19)
  close(20)
  

  

 end subroutine Evans60