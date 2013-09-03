       subroutine sphere2ico(points,inp,ingens)
c=====================================================================c
c                                                                     c
c        Description: Takes points on the sphere and projects them    c
c                     radial inward to the plane defined by ...       c
c                                                                     c
c                                                                     c
c        Input:  point -- points on the sphere to project onto plane  c
c                  inp -- number of points to project                 c
c                                                                     c
c                                                                     c
c        Output: point -- projected points                            c
c                                                                     c
c        Note: input array point is overwritten with new points       c
c                                                                     c
c=====================================================================c   
c
	 implicit real *8 (a-h,o-z)       
	 real *8 points(3,inp),normal(3),gamma, fc(3)
	 real *8 v1(3),v2(3),v3(3),v4(3),v5(3),v6(3)
	 real *8 vv1(3),vv2(3), nx, ny, nz
         real *8 m,mp,de,lambda(3)
	 integer ingens
c
	 m  = 0.5d0*(1.0d0        + 5.0d0**0.5d0)
	 mp = 0.5d0*(5.0d0**0.5d0 - 1.0d0       )
	 de = (mp*5.0d0**0.5)**0.5d0	
c
c--------Vertex one
c	 	 
	 v1(1) =  0.0d0
	 v1(2) = -mp / de
	 v1(3) =  1.0d0 / de
c
c--------Vertex two
c	 
	 v2(1) =  0.0d0
	 v2(2) =  mp / de
	 v2(3) =  1.0d0 / de  	 	 
c
c--------Vertex three
c	 
	 v3(1) = 1.0d0 / de 
	 v3(2) = 0.0d0
	 v3(3) = mp / de 
c
c--------Vertex four
c
         v4(1) =  mp / de 
	 v4(2) = -1.0d0 / de
	 v4(3) =  0.0d0
c
c--------Vertex five
c
         v5(1) = mp / de 
	 v5(2) = 1.0d0 / de
	 v5(3) = 0.0d0
c
c--------Vertex six
c
         v6(1) = -1.0d0 / de 
	 v6(2) =  0.0d0
	 v6(3) =  mp / de
	 
c=========================================================================
c
c--------Start projecting points from spherical triangle 1
c
         vv1(1) = v3(1) - v1(1)
         vv1(2) = v3(2) - v1(2)
	 vv1(3) = v3(3) - v1(3)
 
         vv2(1) = v2(1) - v1(1)
         vv2(2) = v2(2) - v1(2)
	 vv2(3) = v2(3) - v1(3)	 
c
c--------calculate normal vector to plane
c	 
	 call cross(vv1,vv2,normal,1)	 
c
c--------Project onto the plane
c
         do i=1,3*ingens
c
c----------Scaling factor
c	   
           gamma = DDOT(3,normal,1,v1,1)/DDOT(3,normal,1,points(1,i),1)          
c
c----------gamma must be < 1
c
           if(gamma.gt.1.0d0) then
	     call prinf('Problem with scaling: gamma > 1*', 0, 0)
	     stop
	   endif 
	   
	   points(1,i) = gamma*points(1,i)
	   points(2,i) = gamma*points(2,i)
	   points(3,i) = gamma*points(3,i)
	   
	 enddo   	 
c=========================================================================
c
c--------Start projecting points from spherical triangle 2
c
         vv1(1) = v5(1) - v3(1)
         vv1(2) = v5(2) - v3(2)
	 vv1(3) = v5(3) - v3(3)
 
         vv2(1) = v2(1) - v3(1)
         vv2(2) = v2(2) - v3(2)
	 vv2(3) = v2(3) - v3(3)
c
c--------calculate normal vector to plane
c	 
	 call cross(vv1,vv2,normal,1)
c
c--------Project onto the plane
c
         do i=1,3*ingens
c
c----------Scaling factor
c	   
           gamma = DDOT(3,normal,1,v3,1) / 
     1	           DDOT(3,normal,1,points(1,3*ingens+i),1)         
c
c----------gamma must be < 1
c	   
           if(gamma.gt.1.0d0) then
	     call prinf('Problem with scaling: gamma > 1*', 0, 0)
	     stop
	   endif 
	   
	   points(1,3*ingens+i) = gamma*points(1,3*ingens+i)
	   points(2,3*ingens+i) = gamma*points(2,3*ingens+i)
	   points(3,3*ingens+i) = gamma*points(3,3*ingens+i)
	   
	 enddo 
c
c=========================================================================
c
c--------Start projecting points from spherical triangle 3
c
         vv1(1) = v2(1) - v1(1)
         vv1(2) = v2(2) - v1(2)
	 vv1(3) = v2(3) - v1(3)
 
         vv2(1) = v6(1) - v1(1)
         vv2(2) = v6(2) - v1(2)
	 vv2(3) = v6(3) - v1(3)
c
c--------calculate normal vector to plane
c	 
	 call cross(vv1,vv2,normal,1)
c
c--------Project onto the plane
c
         do i=1,3*ingens
c
c----------Scaling factor
c	   
           gamma = DDOT(3,normal,1,v1,1) / 
     1	           DDOT(3,normal,1,points(1,6*ingens+i),1)         
c
c----------gamma must be < 1
c	   
           if(gamma.gt.1.0d0) then
	     call prinf('Problem with scaling: gamma > 1*', 0, 0)
	     stop
	   endif 
	   
	   points(1,6*ingens+i) = gamma*points(1,6*ingens+i)
	   points(2,6*ingens+i) = gamma*points(2,6*ingens+i)
	   points(3,6*ingens+i) = gamma*points(3,6*ingens+i)
	   
	 enddo 
c
c
c=========================================================================
c
c--------Start projecting points from spherical triangle 4
c
         vv1(1) = v4(1) - v1(1)
         vv1(2) = v4(2) - v1(2)
	 vv1(3) = v4(3) - v1(3)
 
         vv2(1) = v3(1) - v1(1)
         vv2(2) = v3(2) - v1(2)
	 vv2(3) = v3(3) - v1(3)
c
c--------calculate normal vector to plane
c	 
	 call cross(vv1,vv2,normal,1)
c
c--------Project onto the plane
c
         do i=1,3*ingens
c
c----------Scaling factor
c	   
           gamma = DDOT(3,normal,1,v1,1) / 
     1	           DDOT(3,normal,1,points(1,9*ingens+i),1)         
c
c----------gamma must be < 1
c	   
           if(gamma.gt.1.0d0) then
	     call prinf('Problem with scaling: gamma > 1*', 0, 0)
	     stop
	   endif 
	   
	   points(1,9*ingens+i) = gamma*points(1,9*ingens+i)
	   points(2,9*ingens+i) = gamma*points(2,9*ingens+i)
	   points(3,9*ingens+i) = gamma*points(3,9*ingens+i)
	   
	 enddo 

	 return
c	 
       end
