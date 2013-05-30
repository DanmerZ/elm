program elm
 use Profiles
 use ArDef
 implicit none	  		  

 call CalculateFlux 
 call ProfilesArrays 
 call CalculateQ 		  
 !call SolveEq 
 call ProfilesToFiles 

		 
end program elm
