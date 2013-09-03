c--------------------------------------------------------------------------------------------c
c        C.A. Dec. 17, 2008
c        Last edit: 12/17/08
c
c        Purpose: Solve the linear system Ax = b to get an initial estimate of the weights
c
c        Input: nfilename
c               wfilename
c		wvert
c		wface
c		wside
c		wei
c		ngen
c		itype      = 1 if vertex + general points; 2 if vertex + face center;
c                            3 if vertex + side midpoints; 4 if all three
c		ifg   
c
c        Output: weights
c--------------------------------------------------------------------------------------------c
c
       subroutine weights(filename,inp,inmax)
c
         implicit real *8 (a-h,o-z)
         real *8 x(inp), y(inp), z(inp)
	 real *8 mat(inp,inp)
	 real *8 b(inp)
	 real *8 work(20*inp)
	 real *8 tau(20*inp)
	 real *8 p(0:inmax+1), mu
	 real *8 nfact
	 integer jpvt(inp)
c	 real *8 wei(*), wface, wside, wvert
         character*24 filename
       
         zero    = 0.0d0
	 one     = 1.0d0      
         pi      = 4.0d0*atan(1.0d0)
	 twopi   = 2.0d0*pi
         fourpi  = 4.0d0*pi	
	 eps     = 1.0d-15
	 
	 factor  = (inmax + 1.0d0) / fourpi 
c	 	 
	 open (unit = 7, file = filename, status = 'unknown', iostat=ios,
     1         form = 'formatted')
c     
c--------Read in points
c     
         do i=1,inp
	   read(7,*) x(i),y(i),z(i)
	 enddo
	 close(7)      
c
c--------Form matrix and RHS
c
         t0 = second()
         do i=1,inp
c
c----------RHS
c
           b(i) = 1.0d0
	   
	   do j=1,inp
c
c------------Calculate mu = <omega_i, omega_j>
c
             mu = x(i)*x(j) + y(i)*y(j) + z(i)*z(j)
	     
	     if(mu.gt.1.0d0) then
	       call prinf('shit! *', 0,0)
	       stop
	     endif
c
c------------Evaluate the kernel
c
             if((1.0d0 - mu).gt.eps) then

	       call legendrep(inmax + 1,mu,p) 
	      
	       mat(i,j) = factor*(p(inmax) - p(inmax+1))/(1.0d0 - mu)
	       
	     else
	     
	       mat(i,j) = factor*(inmax + 1)
	     
	     endif	     	   
	   
	   
	   enddo
	   
	 enddo
         t1 = second()
	 call prin2x('time to create matrix: *', t1-t0,1)
c
c--------Solve linear system for initial weights
c	
         do i=1,inp
	   jpvt(i) = 0
	 enddo
	 
	 rpara = 1.0d-12
     
         t0 = second()
     
         call DGELSY(inp,inp,1,mat,inp,b,inp,jpvt,rpara,irank,
     $               work,20*inp,info)
     
     	 t1 = second()
	 
	 call prin2x('time to solve for weights: *', t1-t0,1)
c     
         if(irank .lt. inp) then
	   call prinf('Jacobian rank inp *',irank,1)
	   stop
	 endif 	 
	 
	 call prin2x('done with solve *',b,inp)
	 
	 return
	  	 
       end

