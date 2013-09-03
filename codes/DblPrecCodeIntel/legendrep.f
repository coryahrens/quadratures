      subroutine legendrep(n, x, pnx)

c*********************************************************************72
c
cc LEGENDRE_POLY evaluates the Legendre polynomials P(N)(X) at X.
c
c  Discussion:
c
c    P(N)(1) = 1.
c    P(N)(-1) = (-1)**N.
c    | P(N)(X) | <= 1 in [-1,1].
c
c    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
c    function of the first kind and order N equals the Legendre polynomial
c    of the first kind and order N.
c
c    The N zeroes of P(N)(X) are the abscissas used for Gauss-Legendre
c    quadrature of the integral of a function F(X) with weight function 1
c    over the interval [-1,1].
c
c    The Legendre polynomials are orthonormal under the inner product defined
c    as integration from -1 to 1:
c
c      Integral ( -1 <= X <= 1 ) P(I)(X) * P(J)(X) dX 
c        = 0 if I =/= J
c        = 2 / ( 2*I+1 ) if I = J.
c
c    Except for P(0)(X), the integral of P(I)(X) from -1 to 1 is 0.
c
c    A function F(X) defined on [-1,1] may be approximated by the series
c      C0*P(0)(X) + C1*P(1)(X) + ... + CN*P(N)(X)
c    where
c      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I)(X) dx.
c
c  Differential equation:
c
c    (1-X*X) * P(N)(X)'' - 2 * X * P(N)(X)' + N * (N+1) = 0
c
c  First terms:
c
c    P( 0)(X) =       1
c    P( 1)(X) =       1 X
c    P( 2)(X) =  (    3 X**2 -       1)/2
c    P( 3)(X) =  (    5 X**3 -     3 X)/2
c    P( 4)(X) =  (   35 X**4 -    30 X**2 +     3)/8
c    P( 5)(X) =  (   63 X**5 -    70 X**3 +    15 X)/8
c    P( 6)(X) =  (  231 X**6 -   315 X**4 +   105 X**2 -     5)/16
c    P( 7)(X) =  (  429 X**7 -   693 X**5 +   315 X**3 -    35 X)/16
c    P( 8)(X) =  ( 6435 X**8 - 12012 X**6 +  6930 X**4 -  1260 X**2 +   35)/128
c    P( 9)(X) =  (12155 X**9 - 25740 X**7 + 18018 X**5 -  4620 X**3 +  315 X)/128
c    P(10)(X) =  (46189 X**10-109395 X**8 + 90090 X**6 - 30030 X**4 + 3465 X**2
c                 -63 ) /256
c
c  Recursion:
c
c    P(0)(X) = 1
c    P(1)(X) = X
c    P(N)(X) = ( (2*N-1)*X*P(N-1)(X)-(N-1)*P(N-2)(X) ) / N
c
c    P'(0)(X) = 0
c    P'(1)(X) = 1
c    P'(N)(X) = ( (2*N-1)*(P(N-1)(X)+X*P'(N-1)(X)-(N-1)*P'(N-2)(X) ) / N
c
c  Formula:
c
c    P(N)(X) = (1/2**N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X**(N-2*M)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to evaluate.
c    Note that polynomials 0 through N will be evaluated.
c
c    Input, double precision X, the point at which the polynomials
c    are to be evaluated.
c
c    Output, double precision CX(0:N), the values of the Legendre polynomials 
c    of order 0 through N at the point X.
c
c    Output, double precision CPX(0:N), the values of the derivatives of the
c    Legendre polynomials of order 0 through N at the point X.
c
        implicit real *8 (a-h,o-z) 
        real *8 pnx(0:*)
	integer i,n
	
        if ( n .lt. 0 ) then
          return
        end if

        pnx(0) = 1.0d0
  
        if (n.lt.1) then
          return
        end if

        pnx(1) = x
	
        do i=2, n
	
	  pnx(i) = ((2.0d0 * i - 1.0d0) * x * pnx(i-1) - 
     1	            (i - 1.0d0) * pnx(i-2))/(1.0d0 * i)

        end do

        return
	
      end
