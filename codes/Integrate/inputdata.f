       subroutine inputdata(x,y,z,w,npoints)
       
         implicit real *8 (a-h,o-z)
	 
         real *8 x(*),y(*),z(*),w(*)
	 integer index(1000)
 
	 integer npoints
	 
	 character*50 qdata	 
c
c--------Read quadrature	 
c	 
c         call getfilename(qdata,'Quadrature data file: *')
	 
c	 call openfile(qdata,7,'old')
	 call openfile('ld-86.dat',7,'old')
	 
c	 read(6,*) npoints
         npoints = 86
	 do j=1,npoints
!	   read(6,*)  index(j), x(j), y(j), z(j), w(j)  
           read(7,*)  x(j), y(j), z(j), w(j)  
	 enddo
	 close(7)
c
         return
       
       end
       
