c
c       C.D.A. April 2009 
c 
       program grid 
c       
         parameter(npnts = 20000)
	 implicit real *8 (a-h,o-z)
	 
	 integer npnts, ngen0, ngen, ipointer(npnts)
	 
	 real *8 x(npnts), y(npnts), z(npnts)
	 	 
	 real *8 gen0(2,npnts), gen(2,npnts)
	 
	 real *8 quad(3,npnts), grid2D3D(5,npnts)
	 
	 real *8 lambda(3), weight(npnts)
	 
	 real *8 p1(2), p2(2), p3(2), pnt(2)
	 real *8 v1(3), v2(3), v3(3)
	 
	 real *8 m, mp, nrm
	 
	 character *25 pgenfile2D, pgenfile3D, ngenfile2D, ngenfile3D
	 
	 zero   = 0.0d0
	 one    = 1.0d0       
         pi     = 4.0d0*atan(1.0d0)
	 eps    = dmacheps()	 	 
c
c--------Initialization for transcript file
c
         call prini(6,13)
         call getini
c
c--------Read previous generators (in 2D)
c
         write(*,*) 'Input number of previous generators: '
         read(*,*) ngen0

         write(*,*) 'Input previous generator file name: '
         read(*,*) pgenfile2D	
	 
         open(unit = 7, file = pgenfile2D, status = 'unknown',
     1	      form = 'formatted')
c     
         do j=1,ngen0
	   read(7,*)  gen0(1,j), gen0(2,j) 
	 enddo
c
	 close(7)
c
c--------Read corresponding quadrature (in 3D)
c
         write(*,*) 'Input previous generator quadrature file name: '
         read(*,*) pgenfile3D
	 
         open(unit = 7, file = pgenfile3D,
     1	      status='unknown',form = 'formatted')
c     
         do j=1,ngen0
	   read(7,*)  quad(1,j), quad(2,j), quad(3,j), weight(j)  
	 enddo
c
	 close(7)

c===================================================	 
c--------Generate 2D grid based on lattice pnts and 
c        3D grid based on quadratures
	 	 
	 call gen2D3D(gen0,ngen0,quad,grid2D3D)
	 
c===================================================
         open(unit = 7, file = 'lattice.dat', status = 'unknown',
     1	      form = 'formatted')         
	 do j=1,(7*ngen0 + 3)
	   
	   write(7,*) grid2D3D(1,j),grid2D3D(2,j)
	   
	 enddo
	 
	 close(7)

         open(unit = 7, file = 'quadGrid.dat', status = 'unknown',
     1	      form = 'formatted')         
	 do j=1,(7*ngen0 + 3)
	   
	   write(7,*) grid2D3D(3,j),grid2D3D(4,j),grid2D3D(5,j)
	   
	 enddo
	 
	 close(7)



c
c--------Read new generators (in 2D)
c
         write(*,*) 'Input number of new generators: '
         read(*,*) ngen
	 
         write(*,*) 'Input new generator file name: '
         read(*,*) ngenfile2D	        

         open(unit = 7, file = ngenfile2D, status = 'unknown',
     1	      form = 'formatted')
c     
         do j=1,ngen
	   read(7,*)  gen(1,j), gen(2,j) 
	 enddo
c
	 close(7)

c===================================================
c        Start search for nearest neighbors and map 
c        to new points

         do j=1,ngen
c
c----------Find nearest neighbors
c
           call findneigh(gen(1,j),grid2D3D,7*ngen0+3,ipointer)
c
c----------Find barycentric coordinates of gen(1,j) wrt neighbors
c
           p1(1) = grid2D3D(1,ipointer(1))
	   p1(2) = grid2D3D(2,ipointer(1))
	   
	   p2(1) = grid2D3D(1,ipointer(2))
	   p2(2) = grid2D3D(2,ipointer(2))
	   
	   p3(1) = grid2D3D(1,ipointer(3))
	   p3(2) = grid2D3D(2,ipointer(3))
	   
	   pnt(1) = gen(1,j)
	   pnt(2) = gen(2,j)	
	     
	   
           call findbcc2d(p1,p2,p3,pnt,lambda)
	   
C	   call prinf('==========================*', 0, 0)
C	   call prinf('j = *', j, 1)
C	   call prin2x('point = *', pnt,2)   
C	     
C           call prin2x('lambda = *', lambda, 3)
	   slambda = lambda(1)+lambda(2)+lambda(3)
C	   call prin2x('sum l =*',slambda,1)
	   
C	   call prinf('pointer = *', ipointer,3) 
c
c----------Map to the sphere using neighbors and bc-coords
c
           v1(1) = grid2D3D(3,ipointer(1))
	   v1(2) = grid2D3D(4,ipointer(1))
	   v1(3) = grid2D3D(5,ipointer(1))

           v2(1) = grid2D3D(3,ipointer(2))
	   v2(2) = grid2D3D(4,ipointer(2))
	   v2(3) = grid2D3D(5,ipointer(2))

           v3(1) = grid2D3D(3,ipointer(3))
	   v3(2) = grid2D3D(4,ipointer(3))
	   v3(3) = grid2D3D(5,ipointer(3))
	   
C	   call prin2x('v1 = *', v1,3)
C	   call prin2x('v2 = *', v2,3)
C	   call prin2x('v3 = *', v3,3)
	   
C	   call prinf('==========================*', 0, 0)
	   	   	   

           x(j) = lambda(1)*v1(1) + lambda(2)*v2(1) + lambda(3)*v3(1)
           y(j) = lambda(1)*v1(2) + lambda(2)*v2(2) + lambda(3)*v3(2)
           z(j) = lambda(1)*v1(3) + lambda(2)*v2(3) + lambda(3)*v3(3)	  
	   
	   nrm  = (x(j)*x(j) + y(j)*y(j) + z(j)*z(j)) ** 0.5d0
	   
	   x(j) = x(j)/nrm
	   y(j) = y(j)/nrm
	   z(j) = z(j)/nrm
	  
	 enddo

         write(*,*) 'Input new mapped generator file name: '
         read(*,*) ngenfile3D

         open(unit = 7, file = ngenfile3D, status ='unknown',
     1	      form = 'formatted')
c     
         do j=1,ngen
	  write(7,2000)  x(j), y(j), z(j) 
	 enddo
c
	 close(7)
	 	 
	
1000     format(2(2x,E14.8))
2000     format(3(2x,E14.8))
	 
       end


