       subroutine cross(p1,p2,p3,inorm)
c==============================================c
c        Vector cross product: p3 = p1 X p2    c
c        Input: p1, p2                         c
c               inorm = normalization flag     c
c                                              c
c        Output: p3                            c
c                                              c
c==============================================c
c
	 implicit real *8 (a-h,o-z)       
       	 real *8 p1(*),p2(*),p3(*), norm

         p3(1) = p1(2)*p2(3) - p2(2)*p1(3)
	 p3(2) = p2(1)*p1(3) - p1(1)*p2(3)
	 p3(3) = p1(1)*p2(2) - p2(1)*p1(2)

	 if(inorm.gt.0) then 
         
	   norm = DNRM2(3,p3,1)
	   
	   if(norm.gt.1.0d-12) then
	   
	     p3(1) = p3(1) / norm
	     p3(2) = p3(2) / norm
	     p3(3) = p3(3) / norm
	     
	   else
	     
	     call prinf('Error: norm close to zero*',0,0)
	     stop
	     
	   endif
           
	 endif
	 
	 return
c	 
       end
