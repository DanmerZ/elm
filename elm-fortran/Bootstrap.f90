subroutine Bootstrap
 use Profiles 
 use ArDef
 implicit none
 
 integer :: i
 real(8), dimension(imax) :: jb, djb, jbr, jbi
 real(8) :: te1,ti1,ne1,dte1,dti1,dne1
 real(8), dimension(imax) :: current, dcurrent, QGr, QGi
 
 real(8) :: dx, j1, Ig
 
 dx = (finish - start)/1.0d0/imax
 
 do i = 1,imax
     ne1 = 1.d+13*ne_(i)    !cm^-3
     te1 = 1.6022d-12*Te_(i) !erg
     ti1 = 1.6022d-12*Ti_(i) !erg
     dne1 = 1.d+13*Dne_(i)
     dte1 = 1.6022d-12*DTe_(i)
     dti1 = 1.6022d-12*DTi_(i)
     jb(i) = ne1*q(x_(i))*(2.44*(te1+ti1)*dne1/ne1 + 0.69*dte1 - 0.42*dti1)
     current(i)=-jb(i)*((ap/aR)**(-0.5))*c/B0
 end do 
 
 djb(1) = (jb(2) - jb(1)) / dx
 dcurrent(1) = (current(2) - current(1)) / dx
 
 do i = 2,imax
     djb(i) = (jb(i) - jb(i-1)) / dx
     dcurrent(i) = (current(i) - current(i-1)) / dx
 end do
 
 do i = 1,imax
     j1=((aR/ap)**1.5)*(4.*PiNum*m*(x_(i))**2)*djb(i)/B0/B0
     jbr(i) = j1*ReQ1_(i)
     jbi(i)=j1*ImQ1_(i)
     
     !QGr(i)=ReQ_(i)+4.*PiNum*m*(aR/ap)*x_(i)*x_(i)*dcurrent(i)*ReQ1_(i)/(c*B0)
     !QGi(i)=ImQ_(i)+4.*PiNum*m*(aR/ap)*x_(i)*x_(i)*dcurrent(i)*ImQ1_(i)/(c*B0)
     
     QGr(i)=4.*PiNum*m*(aR/ap)*x_(i)*x_(i)*dcurrent(i)*ReQ1_(i)/(c*B0)
     QGi(i)=4.*PiNum*m*(aR/ap)*x_(i)*x_(i)*dcurrent(i)*ImQ1_(i)/(c*B0)
 end do
 
 open(11, file="data/Bootstrap.dat")
 open(12, file="data/QG.dat")
 do i = 1, imax
     write(11,*) x_(i),jbr(i),jbi(i)
     write(12,*) x_(i),QGr(i),QGi(i)
     !write(11,*) x_(i), -jb(i)*c*((ap/aR)**(-0.5)) / B0
 end do
 close(11) 
 close(12)
 
 200 format(A20,E10.2)
 write(*,200) 'ne, cm-3', 1.d+13*ne_(i1*imax/i2)
 write(*,200) 'Te, erg', 1.6022d-12*Te_(i1*imax/i2)
 write(*,200) 'Ti, erg', 1.6022d-12*Ti_(i1*imax/i2)
  write(*,200) 'dne, cm-3', 1.d+13*Dne_(i1*imax/i2)
 write(*,200) 'dTe, erg', 1.6022d-12*DTe_(i1*imax/i2)
 write(*,200) 'dTi, erg', 1.6022d-12*DTi_(i1*imax/i2)
 write(*,200) 'P, dyne/cm^2', 1.d+13*ne_(i1*imax/i2)*1.6022d-12* (Te_(i1*imax/i2) + Ti_(i1*imax/i2)) 
 write(*,200) 'ne*q*(2.44*......',jb(i1*imax/i2) 
 write(*,200) 'd[n*q*(2.44...)]/dan = ', djb(i1*imax/i2)
 write(*,*) 'Finally calculate:'
write(*,200) 'Re', jbr(i1*imax/i2)
 write(*,200) 'Im', jbi(i1*imax/i2)
 
 
 !full bootstrap current Ig
 !Ig = 0.0d0
 !do i = 1, imax
 !    Ig = Ig + dx*x_(i)*jb(i)*c*((ap/aR)**(-0.5)) / B0
 !end do
 !print *, 2.*PiNum*Ig*3.d-9
 
end subroutine Bootstrap
