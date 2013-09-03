       subroutine ico2plane(points,inp,ingens,plpoints)
c=====================================================================c
c                                                                     c
c        Description:  ...                                            c
c                                                                     c
c                                                                     c
c        Input:  points --    c
c                  inp  -- number of points to project                c
c                                                                     c
c                                                                     c
c        Output: points -- projected points                           c
c                                                                     c
c        Note: input array point is overwritten with new points       c
c                                                                     c
c=====================================================================c   
c
	 implicit real *8 (a-h,o-z)       
	 real *8 points(3,inp), plpoints(5,inp)
	 real *8 v1(3), v2(3), v3(3), v4(3), v5(3), v6(3)
	 real *8 v1i(3),v2i(3),v3i(3), v4i(3), v5i(3), v6i(3)
	 real *8 va1(3), va2(3), na(3), nb(3), n
	 real *8 lambda(3)
         real *8 sqrt3, m, mp, de
	 integer ingens
c
         sqrt3 = 3.0d0 ** 0.5d0	
c
	 m  = 0.5d0*(1.0d0        + 5.0d0**0.5d0)
	 mp = 0.5d0*(5.0d0**0.5d0 - 1.0d0       )
	 de = (mp*5.0d0**0.5)**0.5d0	
c
c--------Vertex one
c	 	 
	 v1i(1) =  0.0d0
	 v1i(2) = -mp / de
	 v1i(3) =  1.0d0 / de
c
c--------Vertex two
c	 
	 v2i(1) =  0.0d0
	 v2i(2) =  mp / de
	 v2i(3) =  1.0d0 / de  	 	 
c
c--------Vertex three
c	 
	 v3i(1) = 1.0d0 / de 
	 v3i(2) = 0.0d0
	 v3i(3) = mp / de 
c
c--------Vertex four
c
         v4i(1) =  mp / de 
	 v4i(2) = -1.0d0 / de
	 v4i(3) =  0.0d0
c
c--------Vertex five
c
         v5i(1) = mp / de 
	 v5i(2) = 1.0d0 / de
	 v5i(3) = 0.0d0
c
c--------Vertex six
c
         v6i(1) = -1.0d0 / de 
	 v6i(2) =  0.0d0
	 v6i(3) =  mp / de	 
c
c--------Verices in the plane

c
c--------Vertex one
c	 	 
	 v1(1) = 0.0d0
	 v1(2) = 0.0d0
	 v1(3) = 0.0d0
c
c--------Vertex two
c	 
	 v2(1) = 1.0d0
	 v2(2) = 0.0d0
	 v2(3) = 0.0d0 	 
c
c--------Vertex three
c	 
	 v3(1) = 0.5d0
	 v3(2) = 0.5d0 * sqrt3
	 v3(3) = 0.0d0 	 	 
c
c--------Vertex four
c
         v4(1) =  0.5d0 
	 v4(2) = -0.5d0 * sqrt3
	 v4(3) =  0.0d0
c
c--------Vertex five
c
         v5(1) = 1.5d0 
	 v5(2) = 0.5d0 * sqrt3
	 v5(3) = 0.0d0
c
c--------Vertex six
c
         v6(1) = -0.5d0
	 v6(2) =  0.5d0 * sqrt3
	 v6(3) =  0.0d0
c	 
c================================================================================
c
c
c--------Find barycentric coordinates w.r.t. v1i,v2i,v3i and project onto plane
c
         do i=1,3*ingens
	         
           call findbcc3d(v1i,v3i,v2i,points(1,i),lambda)
	   
	   plpoints(1,i) = lambda(1)*v1(1) + lambda(2)*v2(1) + 
     1      	           lambda(3)*v3(1)
	   plpoints(2,i) = lambda(1)*v1(2) + lambda(2)*v2(2) +
     1 	                   lambda(3)*v3(2)
	   plpoints(3,i) = lambda(1)
	   plpoints(4,i) = lambda(2)
	   plpoints(5,i) = lambda(3)
	        
         enddo
c
c--------Find barycentric coordinates w.r.t. v3i,v5i,v2i and project onto plane
c
         do i=1,3*ingens
	         
           call findbcc3d(v3i,v5i,v2i,points(1,3*ingens + i),lambda)
	   
	   plpoints(1,3*ingens + i) = lambda(1)*v2(1) + lambda(2)*v5(1) +
     1	                              lambda(3)*v3(1)
	   plpoints(2,3*ingens + i) = lambda(1)*v2(2) + lambda(2)*v5(2) + 
     1	                              lambda(3)*v3(2)
	   plpoints(3,3*ingens + i) = lambda(1)
	   plpoints(4,3*ingens + i) = lambda(2)
	   plpoints(5,3*ingens + i) = lambda(3)
	        
         enddo
c
c--------Find barycentric coordinates w.r.t. v2i,v6i,v1i and project onto plane
c
         do i=1,3*ingens
	         
           call findbcc3d(v2i,v6i,v1i,points(1,6*ingens + i),lambda)
	   
	   plpoints(1,6*ingens + i) = lambda(1)*v3(1) + lambda(2)*v6(1) + 
     1	                              lambda(3)*v1(1)
	   plpoints(2,6*ingens + i) = lambda(1)*v3(2) + lambda(2)*v6(2) +
     1	                              lambda(3)*v1(2)
	   plpoints(3,6*ingens + i) = lambda(1)
	   plpoints(4,6*ingens + i) = lambda(2)
	   plpoints(5,6*ingens + i) = lambda(3)
     
         enddo
c
c--------Find barycentric coordinates w.r.t. v1i,v4i,v3i and project onto plane
c
         do i=1,3*ingens
	         
           call findbcc3d(v1i,v4i,v3i,points(1,9*ingens + i),lambda)
	   
	   plpoints(1,9*ingens + i) = lambda(1)*v1(1) + lambda(2)*v4(1) +
     1	                              lambda(3)*v2(1)
	   plpoints(2,9*ingens + i) = lambda(1)*v1(2) + lambda(2)*v4(2) +
     1	                              lambda(3)*v2(2)
	   plpoints(3,9*ingens + i) = lambda(1)
	   plpoints(4,9*ingens + i) = lambda(2)
	   plpoints(5,9*ingens + i) = lambda(3)
	        
         enddo 
	 
	 return
c	 
       end
