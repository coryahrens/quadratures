       subroutine cart2sphere(p,theta,phi,cphi,sphi)
c========================================================================c
c        Purpose: Convert from cartesian coordinates to spherical        c
c                                                                        c
c        Input :  p --- point in R^3, assumed to be on the unit sphere   c
c                                                                        c
c        Output:  theta -- polar angle (measured from north pole)        c
c                 phi   -- azimuthal angle                               c
c========================================================================c
c
         implicit real*8 (a-h,o-z)
         
	 real *8 p(*), theta, phi, phiTmp, length, stheta, pi
	 real *8 cphi, sphi
	 
	 pi  = 4.0d0*atan(1.0d0)
c
c--------Calculate lenght and check that it's not close to zero
c	 
	 length = (p(1)*p(1) + p(2)*p(2) + p(3)*p(3))**0.5d0
	 
	 if(length.le.1.0d-12) then
	 
	   write(*,*) 'Error, vector has length .lt. 1e-12'
	   stop
	   
	 endif
	 
	 p(1) = p(1) / length
	 p(2) = p(2) / length
	 p(3) = p(3) / length
	 	
         theta = acos( p(3) ) 
	 
	 stheta = (1.0d0 - p(3)*p(3)) ** 0.5d0

         if(stheta .le. 1.0d-12) then
c
c----------Too close to a pole where phi is indeterm. Set phi = 0	 
c
	   phi = 0.0d0
	   
	   cphi = 1.0d0
	   sphi = 0.0d0
	 
	 else
	 
	   cphi = p(1) / stheta
	   
	   sphi = p(2) / stheta
	 
C	   if( (p(2) / stheta) .gt. 1.0d0) then
C	   
C	     phiTmp = pi / 2.0d0
C	     
C	   elseif( (p(2) / stheta) .lt. -1.0d0) then
C	   
C	     phiTmp = -pi / 2.0d0
C	   
C	   else
C	 
C	     phiTmp = asin( p(2) / stheta )
C	   
C	   endif
C	   
C	   if( (p(1) .ge. 0.0d0)  .and. (p(2) .ge. 0.0d0) ) then
c
c------------1st quad
c	   
C             phi = phiTmp
C	     
C	   elseif( (p(1) .le. 0.0d0) .and. (p(2) .ge. 0.0d0) ) then
c
c------------2nd quad
c	   
C             phi = pi - abs(phiTmp)	
C	   
C	   elseif( (p(1) .le. 0.0d0) .and. (p(2) .le. 0.0d0) ) then
c
c------------3nd quad
c	   
C             phi = pi + abs(phiTmp)
C	     
C	   else 
c
c------------4th quad
c
C	     phi = 2.0d0 * pi - abs(phiTmp)
C   	   	   
C	   endif 
C	 
	 endif	 
	
         return
	
       end
