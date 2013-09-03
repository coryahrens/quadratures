       subroutine findbcc2d(p1,p2,p3,x,lambda)
c======================================================================c
c        Given three points in the plane, p1, p2 and p3, this          c
c        routine returns the barycentric coordinates, lambda,          c
c        of the point x wrt p1-p2-p3 by solving the system             c
c                                                                      c
c           lambda(1)       + lambda(2)       + lambda(3)       = 1    c
c           lambda(1)*p1(1) + lambda(2)*p2(1) + lambda(3)*p3(1) = x(1) c
c           lambda(1)*p1(2) + lambda(2)*p2(2) + lambda(3)*p3(2) = x(2) c     
c                                                                      c
c        Input: p1, p2, p3 and x in R^2                                c
c                                                                      c
c        Output: lambda in R^3                                         c
c                                                                      c
c======================================================================c
c
	 implicit real *8 (a-h,o-z)       
       	 real *8 p1(2),p2(2),p3(2),x(2)
	 real *8 lambda(3), det, meps	 

         meps  = dmacheps()
c
c--------Determinant of linear system
c
         det = -p2(2)*p3(1) + p1(2)*(p3(1) - p2(1)) + 
     1  	p1(1)*(p2(2) - p3(2)) + p2(1)*p3(2)

         if(abs(det).lt.meps) then
	   write(*,*) 'Error: det = 0'
	   stop
	 endif
c
c--------Barycentric coordinates
c     
         lambda(1) = -(p2(2)*p3(1) - p2(1)*p3(2) - p2(2)*x(1) +
     1	               p3(2)*x(1) + p2(1)*x(2) - p3(1)*x(2))/det
	 
	 lambda(2) =  (p1(2)*p3(1) - p1(1)*p3(2) - p1(2)*x(1) + 
     1	               p3(2)*x(1) + p1(1)*x(2) - p3(1)*x(2))/det
	 
	 lambda(3) = -(p1(2)*p2(1) - p1(1)*p2(2) - p1(2)*x(1) + 
     1	               p2(2)*x(1) + p1(1)*x(2) - p2(1)*x(2))/det   
c	 
         if(abs(lambda(1)).lt.2.0d0*meps) then
	 
	   lambda(1) = 0.0d0
	   
	 elseif(abs(lambda(2)).lt.2.0d0*meps) then
	 
	   lambda(2) = 0.0d0
	   
	 elseif(abs(lambda(3)).lt.2.0d0*meps) then
	 
	   lambda(3) = 0.0d0
	   
	 endif
	 
	 return
c	 
       end
