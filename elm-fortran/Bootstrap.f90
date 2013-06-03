subroutine Bootstrap
 use Profiles 
 use ArDef
 implicit none
 
 integer :: i
 real(8), dimension(imax) :: jb, djb, jbr, jbi
 real(8) :: te1,ti1,ne1,dte1,dti1,dne1
 
 real(8) :: dx, j1, Ig
 
 dx = (finish - start)/imax
 
 do i = 1,imax
     ne1 = 1.d+13*ne_(i)    !cm^-3
     te1 = 1.6022d-12*Te_(i) !erg
     ti1 = 1.6022d-12*Ti_(i) !erg
     dne1 = 1.d+13*Dne_(i)
     dte1 = 1.6022d-12*DTe_(i)
     dti1 = 1.6022d-12*DTi_(i)
     jb(i) = ne1*q(x_(i))*(2.44*(te1+ti1)*dne1/ne1 + 0.69*dte1 - 0.42*dti1)
     
 end do 
 !print*, 'Te', Te_(9500), 'Ti',Ti_(9500), 'ne', ne_(9500)
 !print*, 'dTe', DTe_(9500), 'dTi', DTi_(9500), 'dne', Dne_(9500)
 
 djb(1) = (jb(2) - jb(1)) / dx
 
 do i = 2,imax
     djb(i) = (jb(i) - jb(i-1)) / dx
 end do
 
 do i = 1,imax
     j1=((aR/ap)**1.5)*(4.*PiNum*m*x_(i)**2)*djb(i)/B0/B0
     jbr(i) = j1*ReQ1_(i)
     jbi(i)=j1*ImQ1_(i)
 end do
 
 open(11, file="data/Bootstrap.dat")
 do i = 1, imax
     write(11,*) x_(i), jbr(i), jbi(i)
     !write(11,*) x_(i), -jb(i)*c*((ap/aR)**(-0.5)) / B0
 end do
 close(11)
 
 !full bootstrap current Ig
 !Ig = 0.0d0
 !do i = 1, imax
 !    Ig = Ig + dx*x_(i)*jb(i)*c*((ap/aR)**(-0.5)) / B0
 !end do
 !print *, 2.*PiNum*Ig*3.d-9
 
end subroutine Bootstrap
