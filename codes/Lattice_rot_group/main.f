c
c       C.D.A. Feb. 2009 
c 
       program grid 
c       
         parameter(npnts = 5000)
	 implicit real *8 (a-h,o-z)
	 
	 real *8 bccoord(3*npnts),lobccoord(3*npnts)
	 real *8 rhbccoord(3*npnts),lhbccoord(3*npnts)
	 real *8 v1(3), v2(3), v3(3), fc(3)
	 real *8 t1(3), t2(3), t3(3), tfc(3)
	 real *8 x(npnts), y(npnts), z(npnts), mp, m
	 
	 
	 integer s, t, npnts, ios
	 
	 character *20 fname3D, fname2D
c
c--------Initialization for transcript file
c
         call prini(6,13)
         call getini

         do j=1,3*npnts
	   bccoord(j) = 0
	   rhbccoord(j)  = 0
	   lhbccoord(j)  = 0
	   lobccoord(j)  = 0
	 end do  	 
c
c--------Icosahedral vertices and face center
c
	 m  = 0.5d0*(1.0d0        + 5.0d0**0.5d0)
	 mp = 0.5d0*(5.0d0**0.5d0 - 1.0d0       )
	 de = (mp*5.0d0**0.5)**0.5d0
	 
c--------Face center
c
         fc(1) = mp / (3.0d0**0.5d0)
	 fc(2) =  0.0d0
	 fc(3) =  m / (3.0d0**0.5d0)	 
c
c--------Vertex three
c	 	 
	 v3(1) =  0.0d0
	 v3(2) = -mp / de
	 v3(3) =  1.0d0 / de
c
c--------Vertex two
c	 
	 v2(1) =  0.0d0
	 v2(2) =  mp / de
	 v2(3) =  1.0d0 / de  	 	 
c
c--------Vertex one
c	 
	 v1(1) = 1.0d0 / de 
	 v1(2) = 0.0d0
	 v1(3) = mp / de
	 	 
c--------Equalateral triangle in the plane	 
c
c--------Face center
c
         tfc(1) = 0.5d0
	 tfc(2) = (3.0d0 ** 0.5d0) / 6.0d0 
	 tfc(3) = 0.0d0	 
c
c--------Vertex one
c	 	 
	 t2(1) = 0.0d0
	 t2(2) = 0.0d0
	 t2(3) = 0.0d0
c
c--------Vertex two
c	 
	 t1(1) =  1.0d0
	 t1(2) =  0.0d0
	 t1(3) =  0.0d0  	 	 
c
c--------Vertex three
c	 
	 t3(1) =  0.5d0 
	 t3(2) = (3.0d0 ** 0.5d0) / 2.0d0
	 t3(3) =  0.0d0 	 	 
c
c======================================================c
c	 	
c
c--------Loop over lattice parameters
c
         do s=1,75
c
           do t=1,s  
c
c------------Get integer lattice positions
c
             call GetBaryCentricCoords(t,s,bccoord,lobccoord,
     1             rhbccoord,lhbccoord,iqtype,ipnts,irhs,ilhs,ilower)
               
	     ngen = ipnts + ilhs + ilower
c
c------------Change to get different quadrature types
c------------CHANGE FILENAME EXTENSION IN makefilename.f
c
             if((iqtype.eq.2)) then
c
c--------------Map to a fundamental domain on an icosahedral face and then to the sphere
c
	       call bctoxyz(bccoord,lobccoord,rhbccoord,
     1                      lhbccoord,iqtype,ipnts,irhs,
     2                       ilhs,ilower,v1,v3,fc,x,y,z)	       

               call makefilename(ngen,fname3D,fname2D,l)
 
	       open (unit = 3, file = fname3D,status = 'unknown', 
     1               iostat=ios, form = 'formatted')
               do i=1, ngen
c
c----------------Project radially outward to the unit sphere and write
c
	         rnorm = (x(i)*x(i) + y(i)*y(i) + z(i)*z(i))**0.5d0
	         write(3,1000) x(i)/rnorm, -y(i)/rnorm, z(i)/rnorm 
	   
	       enddo

	       close(3)
c
c--------------Map to the plane
c
	       call bctoxyz(bccoord,lobccoord,rhbccoord,
     1                      lhbccoord,iqtype,ipnts,irhs,
     2                      ilhs,ilower,t2,t1,tfc,x,y,z)
     	       open (unit = 3, file = fname2D,status = 'unknown',
     1               iostat=ios, form = 'formatted')
               do i=1, ngen

	         write(3,1000) x(i), y(i) 
	   
	       enddo

	       close(3)
	       
	       stop
	       
	     endif                              

CC           enddo
	   
CC	 enddo 

1000     format(3(2x,E20.10)) 
2000     format(2(2x,E20.10))
		 
       end
c
c======================================== 
c

