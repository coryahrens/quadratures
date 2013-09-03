c
c       C.A. July 2009
c 
        program driver 
	implicit real *8 (a-h,o-z)
	parameter(ngenmax = 2000)
	parameter(nnmax = 30)
c
	real *8 x(120*ngenmax), y(120*ngenmax), z(120*ngenmax)
	real *8 w(120*ngenmax), p1(3), p2(3), sR, sI, sumR, sumI
	real *8 theta, phi, cc, cphi, sphi, p(0:nnmax,0:nnmax)
	real *8 dp(0:nnmax,0:nnmax),phire(0:nnmax),phiim(0:nnmax),tol
	real *8 coefR(0:nnmax,0:nnmax),coefI(0:nnmax,0:nnmax)

        real *8 sum, pi, f_value

	integer npoints, nnmax, ngenmax, nmax

c
c     ----------------------------- 
        character*50 fname
c     ----------------------------- 
c
        tol = 1.0d-32
        pi  = dacos(-1.0d0)

        call prini(6,13)
        call getini	
c
c-------Input quadrature information 
c
c        call inputdata(x,y,z,w,npoints) 

!        open (unit=7,file='md008.00081.new.txt')
!        open(unit=7,file='qsph1-100-3432DP.dat')
        open(unit=7,file='qsph1-149-7512DP.dat')
!        open(unit=7,file='s6-tort-2.txt')
         
        sum = 0.0d0
        read(7,*) npoints
	do j=1,npoints  
          read(7,*)  x(j), y(j), z(j), w(j)
          sum = sum + w(j)  
	enddo
	close(7)

	do j=1,npoints  
          w(j) = w(j)*(4.0d0*pi/sum)  
	enddo



	sigma = 0.25d0
	p1(1) = 0.0d0 
	p1(2) = 0.0d0
	p1(3) = 1.0d0	
c
        nmax = 30

        do n=0,nmax
	 
	  do m=0,n
	  	
c
c-----------Integrate
c
            sumR = 0.0d0
	    sumI = 0.0d0
			
            do i=1,npoints
c
c-------------Integration node
c	
              p2(1) = x(i)
	      p2(2) = y(i)
	      p2(3) = z(i)	
c
c-------------value of function at integration node
c	
              f_value = f(p1,p2,sigma)
c
c-------------Spherical harmonic values
c
              call cart2sphere(p2,theta,phi,cphi,sphi)
	  
	      cc = cos(theta)
 
              call spharm(nnmax,nnmax,nnmax,cc,theta,sphi,cphi,p,dp,
     1                    phire,phiim,tol)
        
	      sR = p(m,n)*phire(m)
	      sI = p(m,n)*phiim(m)	 
	  
	      sumR = sumR + w(i) * sR * f_value
	      sumI = sumI + w(i) * sI * f_value     	      
c		
	    enddo

	    coefR(m,n) = sumR
            coefI(m,n) = sumI

C            if(sqrt(sumR**2 + sumI**2).gt.1e-11) then
C	    
C	      !write(*,*) n,m,sqrt(sumR**2 + sumI**2)	    	    
C
C	    endif
	  
	  enddo
	  
	  call prinf('Done with subspace of degee *', n, 1)
	  
	enddo
c
        open (unit=7,file='coefR_l.dat',status='NEW',iostat=ios,
     1         form='formatted')
        do n=0,nmax
	  write(7,*) coefR(0,n)
	enddo
	close(7)
c
        open (unit=7,file='coefR_lm.dat',status='NEW',iostat=ios,
     1         form='formatted')
        do n=1,nmax
	  do m=1,n
	    write(7,*) coefR(m,n)
	  enddo
	enddo
	close(7)
	
        open (unit=7,file='coefI.dat',status='NEW',iostat=ios,
     1         form='formatted')
        do n=0,nmax
	  do m=0,n
	    write(7,*) coefI(m,n)
	  enddo
	enddo
	close(7)	

c
c-------Done 
c        
	end




        function f(p1,p2,sigma)
c
c---------Gaussian on the sphere with decay away from the point p1
c	
	  implicit real *8 (a-h,o-z)

	  real *8 p1(*),p2(*)
	  
	  dist2 = (p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 +
     1	          (p1(3)-p2(3))**2 
	  
	  f = exp( -dist2 / (4.0d0 * sigma**2) )

	  return
	  
	end
	




























