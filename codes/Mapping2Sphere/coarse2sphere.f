       subroutine coarse2sphere(allpnts,ing,gens)
c====================================================================c
c        Description:                                                c
c                                                                    c
c        Input: gens    -- (xyz) coordinates of generators           c
c               ing     -- number of generators                      c
c               allpnts -- all 12*ing points                         c
c                                                                    c
c        Output:                                                     c
c                                                                    c
c  BELOW ARE THE CORRECT VERTICES AND FACECENTER 5/1                 c
c                                                                    c
c====================================================================c   
c
	 implicit real *8 (a-h,o-z)       
	 real *8 gens(3,ing), rpnts(3,3*ing), allpnts(3,12*ing)
	 real *8 fc(3), v1(3), v2(3), v3(3)
	 real *8 r, rp, mp, m
	 integer ing
	 
	 pi = 4.0d0*atan(1.0d0)
	 m  = 0.5d0*(1.0d0        + 5.0d0**0.5d0)
	 mp = 0.5d0*(5.0d0**0.5d0 - 1.0d0       )
	 de = (mp*5.0d0**0.5)**0.5d0
c
c--------Face center
c
         fc(1) = mp / (3.0d0**0.5d0)
	 fc(2) =  0.0d0
	 fc(3) =  m / (3.0d0**0.5d0)	 
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
c===========Points on spherical triangle 1 	 
c
c--------Place generators in first "ing" spots
c	 
	 do i=1,ing
	   allpnts(1,i) = gens(1,i)
	   allpnts(2,i) = gens(2,i)
	   allpnts(3,i) = gens(3,i)
	 enddo
c
c--------Rotate generators about face center by 2Pi/3
c 
         call ScaleRotate(2.0d0*pi/3.0d0,gens,ing,rpnts,1.0d0,fc,3)
c	  
         do i=1,ing
	   allpnts(1,ing + i) = rpnts(1,i)
	   allpnts(2,ing + i) = rpnts(2,i)
	   allpnts(3,ing + i) = rpnts(3,i)
	 enddo 	 	 
c	         
         call ScaleRotate(4.0d0*pi/3.0d0,gens,ing,rpnts,1.0d0,fc,3)
c
         do i=1,ing
	   allpnts(1,2*ing + i) = rpnts(1,i)
	   allpnts(2,2*ing + i) = rpnts(2,i)
	   allpnts(3,2*ing + i) = rpnts(3,i)
	 enddo
c
c===========Points on spherical triangle 2	 
c
c--------Rotate around v3 by -2Pi/5 
c
         call ScaleRotate(-2.0d0*pi/5.0d0,allpnts,3*ing,rpnts,
     1	                   1.0d0,v3,3)	 
	 do i=1,3*ing
	   allpnts(1,3*ing + i) = rpnts(1,i)
	   allpnts(2,3*ing + i) = rpnts(2,i)
	   allpnts(3,3*ing + i) = rpnts(3,i)
	 enddo	 
c
c===========Points on spherical triangle 3
c
c--------Rotate around v1 by +2Pi/5 
c
         call ScaleRotate(2.0d0*pi/5.0d0,allpnts,3*ing,rpnts,
     1	                  1.0d0,v1,3)
c	 
	 do i=1,3*ing
	   allpnts(1,6*ing + i) = rpnts(1,i)
	   allpnts(2,6*ing + i) = rpnts(2,i)
	   allpnts(3,6*ing + i) = rpnts(3,i)
	 enddo	
c
c===========Points on spherical triangle 4 
c
c--------Rotate around v1 by -2Pi/5 
c
         call ScaleRotate(-2.0d0*pi/5.0d0,allpnts,3*ing,rpnts,
     1	                  1.0d0,v1,3)
c     	 
	 do i=1,3*ing
	   allpnts(1,9*ing + i) = rpnts(1,i)
	   allpnts(2,9*ing + i) = rpnts(2,i)
	   allpnts(3,9*ing + i) = rpnts(3,i)
	 enddo	 
c
	 return
c	 
       end
