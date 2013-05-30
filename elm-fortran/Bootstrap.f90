subroutine Bootstrap
 use Profiles 
 use ArDef
 implicit none
 
 integer :: i
 real(8), dimension(imax) :: jb, djb, jbr, jbi
 real(8) :: te1,ti1,ne1,dte1,dti1,dne1
 
 real(8) :: dx, j1
 
 dx = (finish - start)/imax
 
 do i = 1,imax
     ne1 = 1.d+13*ne_(i)    !cm^-3
     te1 = 1.1604d+4*Te_(i) !K
     ti1 = 1.1604d+4*Ti_(i)  !K
     dne1 = 1.d+13*Dne_(i)
     dte1 = 1.1604d+4*DTe_(i)
     dti1 = 1.1604d+4*DTi_(i)
     jb(i) = ne1*q1(x_(i))*(2.44*(te1+ti1)*dne1/ne1 + 0.69*dte1 - 0.42*dti1)
     
 end do 
 
 djb(1) = (jb(2) - jb(1)) / dx
 
 do i = 2,imax
     djb(i) = (jb(i) - jb(i-1)) / dx
 end do
 
 do i = 1,imax
     j1=((aR/ap)**1.5)*(4.*PiNum*m*x_(i)**2)*djb(i)/B0
     jbr(i) = j1*ReQ1_(i)
     jbi(i)=j1*ImQ1_(i)
 end do
 
 
end subroutine Bootstrap
