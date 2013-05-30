module Profiles
 implicit none

 integer, parameter, public :: m =-11      
 integer, parameter, public :: n = 3	
 real(8), parameter, public :: PiNum = 3.141592654 
 real(8), parameter, public :: c = 2.9979245800d+10   !light speed cm/sec
 real(8), parameter, public :: ap = 61		      !minor radius, cm
 real(8), parameter, public :: aR = 168		      !major r, cm
 real(8), parameter, public :: B0 = 20000.0d0	      ! = 2 Tesla 	
				
 real(8), private :: deltaX =  0.0435		 !general shift for ne, Te, Ti
 real(8), parameter, private :: cq1 = 3.6	 ! q = 1 + cq1*x**cq2
 real(8), parameter, private :: cq2 = 5.6	 !6.60032	  
		

contains
!=====================  q safety factor =====================================
function q(x)
 implicit none
 real(8) q
 real(8), intent(in) :: x
 real(8) :: b
 b = 1.00d0  !change q here!!!		 
 q = b + cq1*x**cq2 
 !q = 1.0 + cq1*x**cq2
 !q = 1.0 + cq1*x**cq2 + 0.3*exp(-50.0*x**2) 
end function q

!this q for calculation profiles, attention, not change!!
function q1(x)			
 implicit none
 real(8) q1
 real(8), intent(in) :: x					 
 q1 = 1.0d0 + cq1*x**cq2 
end function q1	
!========  Fm(an)  ===========================
function Fm(x)  !Fm(an)
 implicit none
 real(8), intent(in) :: x
 real(8) Fm
 real(8) x1
 x1 = x
 Fm = abs(m)/q(x1) - n 
 if (x1 >= 1.0d0) then
	Fm = abs(m)/q(1.0d0)/x1/x1 - n
 end if

end function Fm	

!=======  ne  10^19 m^-3 electron density  ============
function ne(x)
 implicit none
 real(8) ne
 real(8), intent(in) :: x 
 real(8) nsep, Delta, an0, xped, xmid, c1
 real(8) an1, aln1, aln2

 real(8) x1		 

 x1 = x - deltaX

 nsep = 1.2
 Delta = 0.086
 an0 = 3.4		 
 xmid = 1 - Delta/2.0d0		 
 c1 = 0 
 !parameters, which impact on the core profile
 xped = 1 - Delta
 an1 = 5
 aln1 = 0.8
 aln2 = 2.1
!black
 if (x <= xped) then
     ne=nsep+an0*(tanh(2*(1-xmid)/Delta)-tanh(2*((x1-c1)-xmid)/Delta)) + &
     HeavisideTheta(1 - (x)/xped)*an1*(1 - ((x)/xped)**aln1)**aln2 
     !HeavisideTheta(1 - (-c1 + x)/xped)*an1*(1 - ((-c1 + x)/xped)**aln1)**aln2 
 else
     ne=nsep+an0*(tanh(2*(1-xmid)/Delta)-tanh(2*((x1-c1)-xmid)/Delta))
 end if
end function ne
!==================  Te ev electron temp. ===========================================================
function Te(x)   ! eV
 implicit none
 real(8) Te
 real(8), intent(in) :: x
 real(8) tsep, Delta, at0, xped, xmid, c1, at1, alt1, alt2			
		 
 real(8) x1
		 		 
 x1 = x - deltaX - 0.008

 tsep = 100
 Delta = 0.07
 at0 = 370			
 xmid = 1 - Delta/2
 c1 = -0.029 
 !parameters, which impact on core profile
 xped = 1 - Delta
 at1 = 3500.0d0
 alt1 = 1.30d0
 alt2 = 1.50d0

 if (x <= xped) then
	Te=tsep+at0*(tanh(2*(1-xmid)/Delta)-tanh(2*((x1-c1)-xmid)/Delta)) + &
	HeavisideTheta(1 - (x)/xped)*at1*(1 - ((x)/xped)**alt1)**alt2
	!HeavisideTheta(1 - (-c1 + x)/xped)*at1*(1 - ((-c1 + x)/xped)**alt1)**alt2
 else 
	Te=tsep+at0*(tanh(2*(1-xmid)/Delta)-tanh(2*((x1-c1)-xmid)/Delta))
 end if

end function Te
!====== Ti eV, ion tempr. ===============================
function Ti(x)  !ev
 implicit none
 real(8) Ti
 real(8), intent(in) :: x
 real(8) tsep, Delta, at0, xped, xmid, c1, at1, alt1, alt2	
	
 real(8) x1	
		 
 x1 = x + 0.053  !- deltaX		 
		 

 tsep = 0.1
 Delta = 0.115
 at0 = 0.33			
 xmid = 1 - Delta/2
 c1 = 0.07
 !parameters, which impact on the core profile
 xped = 1 - Delta
 at1 = 5500.0
 alt1 = 0.80d0
 alt2 = 2.0
!black
 if (x <= xped) then
	Ti=1000.0d0*(tsep+at0*(tanh(2*(1-xmid)/Delta) - tanh(2*((x1-c1)-xmid)/Delta))) + &
	HeavisideTheta(1 - (x)/xped)*at1*(1 - ((x)/xped)**alt1)**alt2
	!HeavisideTheta(1 - (-c1 + x)/xped)*at1*(1 - ((-c1 + x)/xped)**alt1)**alt2
 else
	Ti=1000.0d0*(tsep+at0*(tanh(2*(1-xmid)/Delta) - tanh(2*((x1-c1)-xmid)/Delta)))
 end if 
	
end function Ti
!================= Pe, [Pa] electron pressure=========================
function Pe(x)   !Pa
 implicit none
 real(8) Pe
 real(8), intent(in) :: x

 Pe = 1.6*ne(x)*Te(x)		

end function Pe
!================= Pi, [Pa] ion pressure===============================
function Pi(x)  !Pa
 implicit none
 real(8) Pi
 real(8), intent(in) :: x

 Pi = 1.6*ne(x)*Ti(x)		

end function Pi
!================ Ef, [V/m] radial electric field ================== 
function Ef(x)
 implicit none
 real(8) Ef
 real(8), intent(in) :: x
 real(8) :: a1,r1,f1,be,ce,d1 !param. which define gradient 
 real(8) :: ga1, ga2, ga3,k1  !parameters, which define the Gausian 
 real(8) x1	
	
 x1 = x
	
 d1 = 7.32
 r1 = 0.59
 f1 = 1.042
 be = 7.71
 ce = 37
 a1 = 8
!parameters, which impact on the core profile
 ga1 = -673855.0d0
 ga2 = 57.3688
 ga3 = 1.4
 k1 = -0.6

if (x >=0.940d0) then
	Ef = 1000*(d1+r1*(x1-f1)*exp(be+ce*r1*(x1-f1))*exp(-(a1*r1*(x1-f1))**2))
else if (x < 0.940d0) then
	Ef = 1000*(ga1*exp(-ga2*(x1 - ga3)**2) + k1*x1)
end if
	!if (Ef >= 0 .and. x <= 1.0d0) Ef = 0.0d0
	
end function Ef 

!===================== speed of toroidal rotation cm/sec =======================

function V0(x)   !speed
 implicit none 
 real(8), intent(in) :: x
 real(8) V0
!V0 = -4500000	
 !V0 = -100000000    !cm/sec
 V0 = 0.0d0
end function V0	

!==============================conductivity===================================
function sigma(x)  !conductivity
 implicit none
 real(8), intent(in) :: x
 real(8) sigma	

 sigma = 1.2d+17 * (Te(x)/500.0d0)**1.5

end function sigma

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function HeavisideTheta(x)
 implicit none
 real(8), intent(in) :: x
 real(8) HeavisideTheta
 if (x < 0.0) then
	HeavisideTheta = 0.0
 else
	HeavisideTheta = 1.0
 end if

end function HeavisideTheta

end module Profiles
