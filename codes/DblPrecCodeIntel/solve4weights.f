c--------------------------------------------------------------------------------------------c
c        C.A. Apr. 23 2009 
c        Last edit: 4/23/09
c
c        Purpose: Use initial node distribution to estimate weights
c
c        Input: ngen -- number of generators
c                gen -- quaternion representation of generators
c                ifg -- 0 = rotation subgroup; 1 = full icosahedral group
c              itype -- type of quadrature
c                       1 = vertex + general points, 2 = vertex + face center
c                       3 = vertex + side midpoints, 4 = all three
c                 np -- total number of points in quadrature
c               nmax -- maximum order and degree of subspace integrated
c             wvert0 -- vertex weights
c             wface0 -- face center weight
c             wside0 -- side center weight
c               wei0 -- generator weights
c           filename -- filename for output of weights     
c
c        Output: wvert0,wface0,wside,wei0 to file 'filename'
c--------------------------------------------------------------------------------------------c
c                   
       subroutine solve4weights(filename,ngen,gen,ifg,itype,np,nmax,
     1                   wvert0,wface0,wside0,wei0) 
c 
         implicit real *8 (a-h,o-z)

         integer np, nmax, jpvt(np)

	 real *8 rg(4,60),rgc(4,60),qtemp(4,60),qres(4,60),gen(4,ngen)
	 real *8 vert(4,12),fcent(4,20),scent(4,30)

         real *8 x(np),y(np),z(np)
	 real *8 A(np,np), b(np), work(20*np)
	 
	 real *8 p(0:nmax + 1), mu
	 real *8 wvert0,wface0,wside0,wei0(ngen)	 
	 
	 real *8 lnorm	 
	 
	 character*24 filename
	 
	 zero    = 0.0d0
	       
         pi      = 4.0d0*atan(1.0d0)
         fourpi  = 4.0d0*pi	
	 eps     = 1.0d-15
	 factor  = (nmax + 1.0d0) / fourpi	 
c
c--------load rotation group
c
         call rotg(rg)
c
c--------load vertices, face centers and side centers
c	 
	 call vertex(vert)
	 call fcenter(fcent)
	 call scenter(scent) 
c
c--------first 12 point are the verices
c	 	
	 do i=1,12
	   x(i) = vert(2,i)
	   y(i) = vert(3,i)
	   z(i) = vert(4,i)
	 enddo   	   
c
c--------others based on quadrature type
c 
         if(itype .eq. 1) then
	 
	   ioffset = 12
	 
         elseif(itype .eq. 2) then

	   do i=1,20
	     x(i+12) = fcent(2,i)
	     y(i+12) = fcent(3,i)
	     z(i+12) = fcent(4,i)
	   enddo 
	   
	   ioffset = 32

	 elseif(itype .eq. 3) then  
	   
	   do i=1,30
	     x(i+12) = scent(2,i)
	     y(i+12) = scent(3,i)
	     z(i+12) = scent(4,i)
	   enddo
	   
	   ioffset = 42
	 
	 elseif(itype .eq. 4) then
	 
	   do i=1,20
	     x(i+12) = fcent(2,i)
	     y(i+12) = fcent(3,i)
	     z(i+12) = fcent(4,i)
	   enddo	   

	   do i=1,30
	     x(i+32) = scent(2,i)
	     y(i+32) = scent(3,i)
	     z(i+32) = scent(4,i)
	   enddo
	   
	   ioffset = 62

	 endif	        
c
c--------compute conjugate quaternions
c
         do k=1,60
           rgc(1,k)   =  rg(1,k)
           do j=2,4
             rgc(j,k) = -rg(j,k)
           enddo 
         enddo	
	  
c---------------------------------------c
c--------loop over generators ----------c
c---------------------------------------c
         do j=1,ngen   
c
           if(ifg.eq.0) then
c	   	 	    
c------------use only the rotation subgroup 
c  
             do k=1,60
c	   
               call quatmult(gen(1,j),rgc(1,k),qtemp(1,k))
               call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
	       
	       indx = (j-1)*60 + k + ioffset
c	     
	       x(indx) = qres(2,k)
	       y(indx) = qres(3,k)
	       z(indx) = qres(4,k)
c
	       lnorm = x(indx)*x(indx) + y(indx)*y(indx) +
     1	               z(indx)*z(indx)
c     
               if(lnorm.le.0.99d0) then
	       
	         call prinf('Error in creating (x,y,z): stopping*',0,0)
                 stop
		 
	       else
	       
	         x(indx) = x(indx)/lnorm 
	         y(indx) = y(indx)/lnorm
	         z(indx) = z(indx)/lnorm
		 	       
	       endif	       
c	       
	     enddo
c	     
	   else
c	   	 	    
c------------Use the full group
c
c------------NEED TO PUT IN THE OFF SET LIKE ABOVE
c
	     do k=1,60
c	     
               call quatmult(gen(1,j),rgc(1,k),qtemp(1,k))
               call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
c	     
	       x((j-1)*120 + k) = qres(2,k)
	       y((j-1)*120 + k) = qres(3,k)
	       z((j-1)*120 + k) = qres(4,k)	     
c
	       x((j-1)*120 + 60 + k) = -qres(2,k)
	       y((j-1)*120 + 60 + k) = -qres(3,k)
	       z((j-1)*120 + 60 + k) = -qres(4,k)	

c--------------NEED A CHECK TO SEE IF UNIT VECTOR
	     
	     enddo
c	        
	   endif
c	   
         enddo
c
c--------Done creating the (x,y,z) quadrature points
c
c==========================================================
c==========================================================
c
c--------Form matrix and RHS
c
         t0 = second()
         do i=1,np
c
c----------RHS
c
           b(i) = 1.0d0
	   
	   do j=1,np
c
c------------Calculate mu = <omega_i, omega_j>
c
             mu = x(i)*x(j) + y(i)*y(j) + z(i)*z(j)
	     
	     if(mu.gt.1.0d0) then
               mu = 1.0d0
	     endif

	     call legendrep(nmax + 1,mu,p)
c
c------------Evaluate the kernel
c
             if((1.0d0 - mu).gt.eps) then 
	      
	       A(i,j) = factor*(p(nmax) - p(nmax+1))/(1.0d0 - mu)
	       
	     else
	     
	       A(i,j) = factor*(nmax + 1)
	     
	     endif	     	   
	   
	   
	   enddo
	   
	 enddo
	 
         t1 = second()
	 call prin2x('Time to create matrix: *', t1-t0,1)	 
c
c--------Solve linear system for initial weights
c	
         do i=1,np
	   jpvt(i) = 0
	 enddo
	 
	 rpara = 1.0d-12
   
         t0 = second()
        
         call DGELSY(np,np,1,A,np,b,np,jpvt,rpara,irank,
     $               work,20*np,info)

         t1 = second()
	 
	 if( (irank.lt.np) .or. (info.gt.0) ) then
	 
	   call prinf('Error in linear solve: rank of A=*',irank,1)
	   call prinf('Stopping...*',0,0)
	   
	 else
	 
 	   call prin2x('Time to solve linear system: *',t1-t0,1)
	   
         endif
c
c=====================================================================
c=====================================================================	 
c
c--------Parse solution b to get individual weights
c	 
c        itype -- type of quadrature
c                 1 = vertex + general points, 2 = vertex + face center
c                 3 = vertex + side midpoints, 4 = all three
         if(itype.eq.1) then
c
c----------Vertex weight
c
	   wvert0 = b(1)
c
c----------Generator weights
c	   
	   do i=1,ngen
             wei0(i) = b(ioffset + i*60)
	   enddo
	   
	   wface0 = zero
	   wside0 = zero
c	 
	 elseif(itype.eq.2) then	 
c
c----------Vertex weight
c
	   wvert0 = b(1)
c
c----------Face center weight
c	   
	   wface0 = b(13)	   
c
c----------Generator weights
c	   
	   do i=1,ngen
             wei0(i) = b(ioffset + i*60)
	   enddo	   

	   wside0 = zero	 
c	   	 
	 elseif(itype.eq.3) then
c
c----------Vertex weight
c
	   wvert0 = b(1)
c
c----------Side center weight
c	   
	   wside0 = b(13)	   
c
c----------Generator weights
c	   
	   do i=1,ngen
             wei0(i) = b(ioffset + i*60)
	   enddo	   

	   wface0 = zero
	   
	 else
c
c----------Vertex weight
c
	   wvert0 = b(1)
c
c----------Face center weight
c	   
	   wface0 = b(13)
c
c----------Side center weight
c	   
	   wside0 = b(43)	   	   	   
c
c----------Generator weights
c	   
	   do i=1,ngen
             wei0(i) = b(ioffset + i*60)
	   enddo	   
c	 
	 endif
c
c====================================================
c
c--------Write weights to a file
c
	 open(unit = 3, file = filename, status = 'unknown',
     1	      form = 'formatted')	 
	 
	 write(3,20) wvert0
	 write(3,20) wface0 
	 write(3,20) wside0
	 
	 do i=1,ngen
	   write(3,20) wei0(i)
	 enddo
	 
	 close(3)
	  
20       format(E19.12) 
 
         return
	 
       end
