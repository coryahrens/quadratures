       subroutine findbcc3d(p1,p2,p3,x,lambda)
c======================================================================c
c        Given three points in the plane, p1, p2 and p3, this          c
c        routine returns the barycentric coordinates, lambda,          c
c        of the point x wrt p1-p2-p3 by using ratios of areas          c
c                                                                      c
c                                                                      c
c        Input: p1, p2, p3, pc, and x in R^3                           c
c                                                                      c
c        Output: lambda in R^3                                         c
c                                                                      c
c======================================================================c
c
	 implicit real *8 (a-h,o-z)       
       	 real *8 p1(*),p2(*),p3(*),x(*)
	 
	 real *8 v1a(3), v1b(3), n(3), na(3), nb(3), v(3), u(3)
	 real *8 v2a(3), v2b(3)
	 
	 real *8 lambda(3), meps, area 	 

         meps  = dmacheps()   
c
c--------Calculate area of large triangle p1-p2-p3
c
         v(1) = p2(1) - p1(1)
	 v(2) = p2(2) - p1(2)
	 v(3) = p2(3) - p1(3)
	 
	 u(1) = p3(1) - p1(1)
	 u(2) = p3(2) - p1(2)
	 u(3) = p3(3) - p1(3)
	 
	 call cross(v, u, n, 0)         
c
c--------Calculate area of "triangle a"
c
         v1a(1) = p3(1) - p2(1) 
	 v1a(2) = p3(2) - p2(2)
	 v1a(3) = p3(3) - p2(3)
	 
	 v2a(1) = x(1) - p2(1)
	 v2a(2) = x(2) - p2(2)
	 v2a(3) = x(3) - p2(3)
	 
	 call cross(v1a, v2a, na, 0)
	 	 
	 lambda(1) = DDOT(3,n,1,na,1) / DNRM2(3,n,1) ** 2
c
c--------Calculate area of "triangle b"
c
         v1b(1) = p1(1) - p3(1) 
	 v1b(2) = p1(2) - p3(2)
	 v1b(3) = p1(3) - p3(3)
	 
	 v2b(1) = x(1) - p3(1)
	 v2b(2) = x(2) - p3(2)
	 v2b(3) = x(3) - p3(3)
	 
	 call cross(v1b, v2b, nb, 0)
	 	 
	 lambda(2) = DDOT(3,n,1,nb,1) / DNRM2(3,n,1) ** 2	 
c
c--------Use normalization to get lambda(3)
c
         lambda(3) = 1.0d0 - lambda(1) - lambda(2)	 
c	 
	 return
c	 
       end
