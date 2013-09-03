       double precision function dmacheps()
c===================================================================c
c        Computes machine epsilon, i.e., find the number such that  c
c         "1 + eps = 1"  and then * 2                               c
c===================================================================c
 
         implicit real *8 (a-h,o-z)
	 
	 real *8 eps
	 
	 eps = 0.5d0
	 
	 do while( (1.0d0 + eps) .gt. 1.0d0 ) 
	 
	   eps = 0.5d0 * eps
	 
	 enddo
		  
	 dmacheps = 2.0d0*eps 
	         
	 return
	 
       end 
