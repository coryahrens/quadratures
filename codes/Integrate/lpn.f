       subroutine lpn(n,x,pn)
c        ===============================================
c        Purpose: Compute Legendre polynomials Pn(x)
c        Input :  x --- Argument of Pn(x)
c                 n --- Degree of Pn(x) ( n = 0,1,...)
c        Output:  Pn(x)
c        ===============================================
c
         implicit real*8 (a-h,o-z)
         dimension pn(0:n)
	
         pn(0) = 1.0d0
	 pn(1) = x
	
         p0 = 1.0d0
         p1 = x
	 
         do j = 2,n
	
           pf = (2.0d0*j - 1.0d0)/j*x*p1 - (j - 1.0d0)/j*p0
       	   pn(j) = pf
          
	   p0 = p1
           p1 = pf
	  
         enddo
	
         return
	
       end
