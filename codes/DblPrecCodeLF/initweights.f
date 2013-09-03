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
       subroutine initweights(filename,maxn,wvert,wface,wside,
     1                        wei,n,itype,ifg,ngen)
c
         implicit real *8 (a-h,o-z)
         real *8 x(n), y(n), z(n)
	 real *8 A(n,n)
	 real *8 b(n)
	 real *8 wrk(20*n)
	 real *8 tau(20*n)
	 real *8 p(0:maxn), mu
	 real *8 sigma(n), nfact
	 real *8 wei(ngen), wface, wside, wvert
         character*24 filename

         zero    = 0.0d0
	 one     = 1.0d0      
         pi      = 4.0d0*atan(1.0d0)
	 twopi   = 2.0d0*pi
         fourpi  = 4.0d0*pi
	 factor  = (maxn + 1)/fourpi	
	 eps     = 1.0d-15
c	 
         call prinf('inside *',0,0)
	 open (unit = 7, file = filename, status = 'unknown', iostat=ios,
     1         form = 'formatted')
c     
c--------Read in points
c     
         do i=1,n
	   read(7,*) x(i),y(i),z(i)
	 enddo
	 close(7)	 
c
c--------Form matrix and RHS
c
         do i=1,n
c
c----------RHS
c
           b(i) = 1.0d0
	   
	   do j=1,n
c
c------------Calculate mu = <omega_i, omega_j>
c
             mu = x(i)*x(j) + y(i)*y(j) + z(i)*z(j)
c
c------------Evaluate the kernel
c
             if((1.0d0 - mu).gt.eps) then
	     
	       call legendrep(maxn+1,mu,p) 
	      
	       A(i,j) = factor*(p(maxn) - p(maxn+1))/(1.0d0 - mu)
	      
	     else
	     
	       A(i,j) = factor*(maxn + 1)
	     
	     endif	     	   
	   
	   
	   enddo
	   
	 enddo	 	 
c
c--------Solve for weights
c
         call prinf('Start linear solve for weights *',0,0)
c	 
c--------Get QR factorization
c
         call dgeqrf(n,n,A,n,tau,wrk,20*n,info)
         if (info.ne.0) then
           print *, 'info dgeqrf',info
           stop
         endif
c
c--------Solve linear system using above QR factorization
c	
         call dgeqrs(n,n,1,A,n,tau,b,n,wrk,20*n,info)
c
         if (info.ne.0) then
           print *, 'info dgeqrs',info
           stop
	 else
	   call prinf('Done with solve for weights *',0,0)
         endif	
	 
	 
c
c--------Normalize solution	 
	 sum = zero
	 do i=1,n
	   sum = sum + b(i)
	 enddo  
	 
         nfact = fourpi/sum
	 
	 do i=1,n
	   b(i) = nfact*b(i)
	 enddo
c
c--------Parse weights vector to get generator weights
c	 
         if(itype.eq.1) then
c
c----------Vertex weight
c
	   wvert = b(1)
c
c----------Generator weights
c	   
	   do i=1,ngen
             wei(i) = b(12 + i*60)
	   enddo
	   
	   wface = zero
	   wside = zero
c	 
	 elseif(itype.eq.2) then
c	   wvert
c	   wface
	   	 
	 elseif(itype.eq.3) then
c	   wvert
c	   wside
	 else
c	   wvert
c          wface
c          wside	 
	 endif
         
	 	 
       end

