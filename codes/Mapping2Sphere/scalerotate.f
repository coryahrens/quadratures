       subroutine ScaleRotate(theta,pnts,inpnts,rpnts,lambda,u,idim)
c=========================================================================c
c        Performs rotations and dilations in 2D and 3D. For 2D rotations, c
c        the rotation is in the "x,y" plane. For 3D rotations, the        c
c        counterclockwise rotation of is about the unit vector u.         c
c        Cartesian coordinates of a point                                 c
c                                                                         c
c        Input: theta  -- angle of rotation                               c
c               pnts   -- points to rotate                                c
c		inpnts -- number of points to rotate                      c
c		lambda -- dilation factor                                 c
c		u      -- unit vector about with 3D rotation is done      c 
c		idim   -- 2 or 3 D                                        c
c                                                                         c
c        Output: rpnts -- rotated and dilated points                      c
c                                                                         c
c=========================================================================c   

	 implicit real *8 (a-h,o-z)       
	 real *8 c, s, lambda, theta, u(*), r(3,3), y(3), ytmp(3), rt(3,3)
	 real *8 pnts(idim,inpnts), rpnts(idim,inpnts) 
	 integer inpnts
c  
c--------Rotation angle
c
         c = cos(theta)
         s = sin(theta)
c	 	 
	 if(idim .eq. 2) then 
c
c----------Apply rotation matrix
c	 
           do i=1,inpnts
	 
  	     rpnts(1,i) =  c*pnts(1,i) + s*pnts(2,i)
	     rpnts(2,i) = -s*pnts(1,i) + c*pnts(2,i)
         
	   end do
c
c----------Apply dilation matrix
c	 
           do i=1,inpnts
	 
	     rpnts(1,i) = lambda*rpnts(1,i)
	     rpnts(2,i) = lambda*rpnts(2,i)
         
	   end do
	   
	   
	 
	 else
c
c----------Rotation matrix
c	  	
	   r(1,1) = u(1)**2  + (1.0d0 - u(1)**2) * c
	   r(1,2) = u(1)*u(2)*(1.0d0 - c) - u(3) * s
	   r(1,3) = u(1)*u(3)*(1.0d0 - c) + u(2) * s
	   r(2,1) = u(1)*u(2)*(1.0d0 - c) + u(3) * s
	   r(2,2) = u(2)**2  + (1.0d0 - u(2)**2) * c
	   r(2,3) = u(2)*u(3)*(1.0d0 - c) - u(1) * s
	   r(3,1) = u(1)*u(3)*(1.0d0 - c) - u(2) * s
	   r(3,2) = u(2)*u(3)*(1.0d0 - c) + u(1) * s
	   r(3,3) = u(3)**2  + (1.0d0 - u(3)**2) * c		   
c
c----------Apply rotation matrix
c	
           do i=1,inpnts
	     
	     ytmp(1) = pnts(1,i)
	     ytmp(2) = pnts(2,i)
	     ytmp(3) = pnts(3,i)
	     
	     call DGEMV('N',3,3,1.0d0,r,3,pnts(1,i),1,0.0d0,y,1)

             rpnts(1,i) = y(1)
	     rpnts(2,i) = y(2)
	     rpnts(3,i) = y(3)
	          
           enddo	
c
c----------Apply dilation matrix
c	 	    
           do i=1,inpnts
	   
             rpnts(1,i) = lambda*rpnts(1,i)
	     rpnts(2,i) = lambda*rpnts(2,i) 
	     rpnts(3,i) = lambda*rpnts(3,i)
	     
           enddo
	 
	 endif	 
c
	 return
c	 
       end
