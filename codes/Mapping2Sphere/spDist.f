       double precision function spDist(p1,p2)
 
         implicit real *8 (a-h,o-z)
	 
	 real *8 p1(*), p2(*)
	 real *8 norm1, norm2, ip
c	 
c--------Calculate norms
c 
         norm1 = DNRM2(3,p1,1)
	 norm2 = DNRM2(3,p2,1)
c
c--------Inner product
c	 
         ip = p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3) 
c
c--------Spherical distance
c	  
	 spDist = dacos( ip / (norm1 * norm2) )
         
	 return
	 
       end 
