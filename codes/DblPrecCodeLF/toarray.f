c
c     G.B. June 2008 
c 
c     This routine helps avoid errors in placing elements of a 2d array into 
c     a one-dimensional memory allocation 
c 
c     In the calling routine "array" is a one dimensional vector of size of at least 
c     nsize*nsize, but in this routine it is treated as a matrix
c     
c 
      subroutine toarray(m,n,array,i,j,val)
        implicit real *8 (a-h,o-z)
        real *8 array(m,n)
c
        array(i,j) = val
c
        return
      end
