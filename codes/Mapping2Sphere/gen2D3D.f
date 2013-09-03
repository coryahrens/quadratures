       subroutine gen2D3D(gens,ngens,quad,grid2D3D)
c====================================================================c
c        Description:                                                c
c                                                                    c
c        Input:
c====================================================================c   
c
	 implicit real *8 (a-h,o-z)       
	 real *8 gens(2,ngens),sgens(2,ngens),rpnts2D(2,3*ngens)
	 real *8 rpnts3D(3,3*ngens), temp(3,3*ngens)
	         
	 real *8 grid2D3D(5,7*ngens + 3), quad(3,ngens)

	 real *8 fc(2), t1(2), t2(2), t3(2)
	 real *8 v1(3), v2(3), v3(3), ifc(3)
	 real *8 m, mp, de
	 
	 integer ngens
	 
	 pi = 4.0d0*atan(1.0d0)
c
c--------Icosahedral vertices and face center
c
	 m  =  0.5d0*(1.0d0        + 5.0d0**0.5d0)
	 mp =  0.5d0*(5.0d0**0.5d0 - 1.0d0       )
	 de = (mp*5.0d0**0.5)**0.5d0
	 
c--------Face center
c
         ifc(1) = mp / (3.0d0**0.5d0)
	 ifc(2) = 0.0d0
	 ifc(3) =  m / (3.0d0**0.5d0)	 
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
c--------Equalateral triangle in the plane	 
c
c--------Face center
c
         fc(1) =  0.5d0
	 fc(2) = (3.0d0 ** 0.5d0) / 6.0d0	 
c
c--------Vertex one
c	 	 
	 t1(1) = 0.0d0
	 t1(2) = 0.0d0
c
c--------Vertex two
c	 
	 t2(1) =  1.0d0
	 t2(2) =  0.0d0 	 	 
c
c--------Vertex three
c	 
	 t3(1) = 0.5d0 
	 t3(2) = (3.0d0 ** 0.5d0) / 2.0d0 
c===========================================================================
c--------Fundamental domain
         do j=1,ngens

	   grid2D3D(1,j) = gens(1,j)
	   grid2D3D(2,j) = gens(2,j)
	    
	   grid2D3D(3,j) = quad(1,j)
	   grid2D3D(4,j) = quad(2,j)
	   grid2D3D(5,j) = quad(3,j)

	 enddo	 
c	 
c--------shift face center to the origin
c
         do j=1,ngens

	   sgens(1,j) = gens(1,j) - fc(1)
	   sgens(2,j) = gens(2,j) - fc(2)

	 enddo
c===========================================================================
c--------Triangle 2
         call ScaleRotate(-2.0d0*pi/3.0d0,sgens,ngens,
     1 	                  rpnts2D,1.0d0,u,2)
         do j=1,ngens
	 
	   grid2D3D(1,ngens + j) = rpnts2D(1,j) + fc(1)
	   grid2D3D(2,ngens + j) = rpnts2D(2,j) + fc(2)
	 
	 enddo
c---------------------------------------------------------------------------	 
         call ScaleRotate(2.0d0*pi/3.0d0,quad,ngens,
     1 	                  rpnts3D,1.0d0,ifc,3)
         do j=1,ngens
	 
	   grid2D3D(3,ngens + j) = rpnts3D(1,j)
	   grid2D3D(4,ngens + j) = rpnts3D(2,j)
	   grid2D3D(5,ngens + j) = rpnts3D(3,j)
	 
	 enddo	        	    
c===========================================================================
c--------Triangle 3	 
         call ScaleRotate(-4.0d0*pi/3.0d0,sgens,ngens,
     1 	                  rpnts2D,1.0d0,u,2)
         do j=1,ngens
	 
	   grid2D3D(1,2*ngens + j) = rpnts2D(1,j) + fc(1)
	   grid2D3D(2,2*ngens + j) = rpnts2D(2,j) + fc(2)
	 
	 enddo
c---------------------------------------------------------------------------	 
         call ScaleRotate(4.0d0*pi/3.0d0,quad,ngens,
     1 	                  rpnts3D,1.0d0,ifc,3)
         do j=1,ngens
	 
	   grid2D3D(3,2*ngens + j) = rpnts3D(1,j)
	   grid2D3D(4,2*ngens + j) = rpnts3D(2,j)
	   grid2D3D(5,2*ngens + j) = rpnts3D(3,j)
	 
	 enddo	        	    
c===========================================================================
c--------Triangle 4
         call ScaleRotate(-pi/3.0d0,sgens,ngens,
     1 	                  rpnts2D,1.0d0,u,2)
     
         do j=1,ngens
	 
	   grid2D3D(1,3*ngens + j) = rpnts2D(1,j) 
	   grid2D3D(2,3*ngens + j) = rpnts2D(2,j) + 1.0d0/3.0d0**0.5d0
	 
	 enddo 
c---------------------------------------------------------------------------	
         do j=1,ngens
	 
	   temp(1,j) = grid2D3D(3,ngens + j)
	   temp(2,j) = grid2D3D(4,ngens + j)
	   temp(3,j) = grid2D3D(5,ngens + j)
	 
	 enddo

         call ScaleRotate(-2.0d0*pi/5.0d0,temp,ngens,
     1 	                  rpnts3D,1.0d0,v2,3)
         do j=1,ngens
	 
	   grid2D3D(3,3*ngens + j) = rpnts3D(1,j)
	   grid2D3D(4,3*ngens + j) = rpnts3D(2,j)
	   grid2D3D(5,3*ngens + j) = rpnts3D(3,j)
	 
	 enddo      	    
c===========================================================================
C--------Triangle 5
         call ScaleRotate(-pi,sgens,ngens,
     1 	                  rpnts2D,1.0d0,u,2)  
     
         do j=1,ngens
	 
	   grid2D3D(1,4*ngens + j) = rpnts2D(1,j) + fc(1)
	   grid2D3D(2,4*ngens + j) = rpnts2D(2,j) - fc(2) 
	 
	 enddo 
c---------------------------------------------------------------------------	
         do j=1,ngens
	 
	   temp(1,j) = grid2D3D(3,2*ngens + j)
	   temp(2,j) = grid2D3D(4,2*ngens + j)
	   temp(3,j) = grid2D3D(5,2*ngens + j)
	 
	 enddo

         call ScaleRotate(-2.0d0*pi/5.0d0,temp,ngens,
     1 	                  rpnts3D,1.0d0,v1,3)
         do j=1,ngens
	 
	   grid2D3D(3,4*ngens + j) = rpnts3D(1,j)
	   grid2D3D(4,4*ngens + j) = rpnts3D(2,j)
	   grid2D3D(5,4*ngens + j) = rpnts3D(3,j)
	 
	 enddo 	 	 
c===========================================================================
c--------Triangle 6
         call ScaleRotate(pi/3.0d0,sgens,ngens,
     1 	                  rpnts2D,1.0d0,u,2)
     	    	 	 
         do j=1,ngens
	 
	   grid2D3D(1,5*ngens + j) = rpnts2D(1,j) + 1.0d0
	   grid2D3D(2,5*ngens + j) = rpnts2D(2,j) + 1.0d0/3.0d0**0.5d0 
	 
	 enddo	
c---------------------------------------------------------------------------	
         do j=1,ngens
	 
	   temp(1,j) = grid2D3D(3,2*ngens + j)
	   temp(2,j) = grid2D3D(4,2*ngens + j)
	   temp(3,j) = grid2D3D(5,2*ngens + j)
	 
	 enddo

         call ScaleRotate(2.0d0*pi/5.0d0,temp,ngens,
     1 	                  rpnts3D,1.0d0,v2,3)
         do j=1,ngens
	 
	   grid2D3D(3,5*ngens + j) = rpnts3D(1,j)
	   grid2D3D(4,5*ngens + j) = rpnts3D(2,j)
	   grid2D3D(5,5*ngens + j) = rpnts3D(3,j)
	 
	 enddo 		  
c===========================================================================
c--------Triangle 7
         call ScaleRotate(pi/3.0d0,sgens,ngens,
     1 	                  rpnts2D,1.0d0,u,2)     	    	 	 
         do j=1,ngens
c	 
	   grid2D3D(1,6*ngens + j) = rpnts2D(1,j) + 0.5d0
	   grid2D3D(2,6*ngens + j) = rpnts2D(2,j) - 3.0d0**0.5d0 / 6.0d0 
c	 
	 enddo	
c---------------------------------------------------------------------------	
         do j=1,ngens
c	 
	   temp(1,j) = grid2D3D(3,j)
	   temp(2,j) = grid2D3D(4,j)
	   temp(3,j) = grid2D3D(5,j)
c	 
	 enddo
c
         call ScaleRotate(-2.0d0*pi/5.0d0,temp,ngens,
     1 	                  rpnts3D,1.0d0,v1,3)
         do j=1,ngens
c	 
	   grid2D3D(3,6*ngens + j) = rpnts3D(1,j)
	   grid2D3D(4,6*ngens + j) = rpnts3D(2,j)
	   grid2D3D(5,6*ngens + j) = rpnts3D(3,j)
c	 
	 enddo
c===========================================================================
c--------Vertices of large triangle and icosahedron
         grid2D3D(1,7*ngens + 2) = 0.0d0
	 grid2D3D(2,7*ngens + 2) = 0.0d0
	 grid2D3D(3,7*ngens + 2) = v1(1)
	 grid2D3D(4,7*ngens + 2) = v1(2)
	 grid2D3D(5,7*ngens + 2) = v1(3)
c
         grid2D3D(1,7*ngens + 1) = 1.0d0
	 grid2D3D(2,7*ngens + 1) = 0.0d0
	 grid2D3D(3,7*ngens + 1) = v3(1)
	 grid2D3D(4,7*ngens + 1) = v3(2)
	 grid2D3D(5,7*ngens + 1) = v3(3)
c
         grid2D3D(1,7*ngens + 3) = 0.5d0
	 grid2D3D(2,7*ngens + 3) = 3.0d0 ** 0.5d0 / 2.0d0
	 grid2D3D(3,7*ngens + 3) = v2(1)
	 grid2D3D(4,7*ngens + 3) = v2(2)
	 grid2D3D(5,7*ngens + 3) = v2(3)	 	 
c       
	 return
c	 
       end
