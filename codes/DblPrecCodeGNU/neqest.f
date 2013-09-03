c
c         G.B. June 2008 
c 
      subroutine neqest(nnmax,neqmax)
      implicit real *8 (a-h,o-z)
      integer nnmax
      real *8 rr
c
      rr = ((nnmax+1.0d0)/30.0d0 + 1.0d0)*(nnmax+1.0d0)/2.0d0+2.0d0
      neqmax = rr
c
c
      return
      end
