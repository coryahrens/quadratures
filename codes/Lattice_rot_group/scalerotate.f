       subroutine ScaleRotate(t,s,pnts,iint,rpnts,lambda)

	 implicit real *8 (a-h,o-z)       
	 real *8 ctheta,stheta,lambda,lsqrd,rise,run
	 real *8 pnts(*), rpnts(2*iint) 
	 integer t,s
c  
c--------Find rotation angle
c  
         lsqrd  = s*s + s*t + t*t
	 rise   = 3.0d0**0.5d0 * s / 2.0d0
	 run    = t + 1.0d0*s/2.0d0
         hypot  = lsqrd**0.5d0
         ctheta = run/hypot
         stheta = rise/hypot
	 
	 write(*,*) t,s,lsqrd,rise,run,hypot,ctheta,stheta
c
c--------Apply rotation matrix
c	 
         do i=1,iint
	 
	   rpnts(2*i-1) =  ctheta*pnts(2*i-1) + stheta*pnts(2*i)
	   rpnts(2*i)   = -stheta*pnts(2*i-1) + ctheta*pnts(2*i)
         
	 end do
c
c--------Apply dilation matrix
c	 
         do i=1,iint
	 
	   rpnts(2*i-1) = lambda*rpnts(2*i-1)
	   rpnts(2*i)   = lambda*rpnts(2*i)
         
	 end do
c
	 return
c	 
       end
