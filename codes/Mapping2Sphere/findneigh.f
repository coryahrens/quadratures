       subroutine findneigh(x,points,inpts,ipointer)
c=========================================================================c
c        Given a list of points in the plane         c
c                                                                         c
c        Input: points, idim, ipnts,  x                                   c
c                                                                         c
c        Output:                                                          c
c                                                                         c
c=========================================================================c
c
	 implicit real *8 (a-h,o-z) 
	 
	 integer inpts
       	 real *8 x(*), points(5,inpts),dist(inpts),diff(2)
	 integer ipointer(*), b
c
c--------Setup initial pointer array
c        
         do i=1,inpts
	 
	   ipointer(i) = i
	 
	 enddo	
c
c--------Chalcuate distances
c
         do i=1,inpts
	     
	   diff(1) = x(1) - points(1,i)
	   diff(2) = x(2) - points(2,i)
	     
	   dist(i) = DNRM2(2,diff,1)
	   
         enddo	  
c
c--------Sort distances and rearrange pointers
c
         do j=1,inpts
	 
	   a = dist(j)
	   b = ipointer(j)
	   
	   do i=j-1, 1, -1
	   
	     if(dist(i) .lt. a) then
	     
	       goto 10
	       
	     else
	       
	       dist(i+1)     = dist(i)
	       ipointer(i+1) = ipointer(i)
	       
	     endif
	   
	   enddo
	 
	   i = 0

10         dist(i + 1)   = a
           ipointer(i+1) = b
	 
	 enddo
	 
	 return
c	 
       end
