       subroutine numbernodeseqs(itype,ngen,ifg,nvar,np)
       
         implicit real *8 (a-h,o-z)
         integer itype,ngen,np

c
c--------Number of equations
c
	 if(itype.eq.1) then
	 
	   nvar = 3*ngen + 1
	   
	 elseif(itype.eq.2) then
	 
	   nvar = 3*ngen + 2
	   
	 elseif(itype.eq.3) then
	 
	   nvar = 3*ngen + 2
	   
	 elseif(itype.eq.4) then 
	 
	   nvar = 3*ngen + 3  
	      
 	 endif 
	 
c
c--------Number of points in quadrature
c
         if(itype.eq.1) then
	 
	   np = 12 + 60*(ifg+1)*ngen
	   
	 elseif(itype.eq.2) then
	 
	   np = 12 + 20 + 60*(ifg+1)*ngen
	   
	 elseif(itype.eq.3) then
	 
	   np = 12 + 30 + 60*(ifg+1)*ngen
	   
	 else
	 
	   np = 12 + 20 + 30 + 60*(ifg+1)*ngen
	   
	 endif	 
	 
  
         return
       
       end
