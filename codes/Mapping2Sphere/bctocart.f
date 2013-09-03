       subroutine bctocart(p1,p2,p3,x,idim,lambda)
c======================================================================c
c        Given barycentric coordinates lambda this routine returns the c
c        Cartesian coordinates of a point                              c
c                                                                      c
c        Input: lambda in  R^3 , idim = 2 or 3 dimensional             c
c                                                                      c
c        Output: x in R^2 or R^3                                       c
c                                                                      c
c======================================================================c
c
	 implicit real *8 (a-h,o-z)       
       	 real *8 p1(*),p2(*),p3(*),x(idim)
	 real *8 lambda(*)

         if(idim.eq.2) then
	   
	   x(1) = lambda(1)*p1(1) + lambda(2)*p2(1) + lambda(3)*p3(1)
	   
	   x(2) = lambda(1)*p1(2) + lambda(2)*p2(2) + lambda(3)*p3(2)
	 
	 else

	   x(1) = lambda(1)*p1(1) + lambda(2)*p2(1) + lambda(3)*p3(1)
	   
	   x(2) = lambda(1)*p1(2) + lambda(2)*p2(2) + lambda(3)*p3(2)
	 
	   x(3) = lambda(1)*p1(3) + lambda(2)*p2(3) + lambda(3)*p3(3)	 
	  
	 endif   
c
	 return
c	 
       end
