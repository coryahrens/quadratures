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
       subroutine iweight(filename,inp,inmax)
c
         implicit real *8 (a-h,o-z)
         real *8 x(inp), y(inp), z(inp)
	 real *8 A(inp,inp)
	 real *8 b(inp)
	 real *8 wrk(5*inp)
	 real *8 tau(5*inp)
	 real *8 p(0:inmax+1), mu
	 real *8 sigma(inp), nfact
c	 real *8 wei(*), wface, wside, wvert
         character*24 filename
       
         zero    = 0.0d0
	 one     = 1.0d0      
         pi      = 4.0d0*atan(1.0d0)
	 twopi   = 2.0d0*pi
         fourpi  = 4.0d0*pi	
	 eps     = 1.0d-15
c	 
         call prinf('inside *',0,0)
      
	 	 
       end

