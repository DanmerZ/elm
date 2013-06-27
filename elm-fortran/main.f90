program elm
 use Profiles
 use ArDef
 implicit none	  		  

 call CalculateFlux
 call ProfilesArrays 
 call CalculateQ 
 call Bootstrap		  
 !call SolveEq 
 call ProfilesToFiles 
 call Evans60

		 
end program elm
