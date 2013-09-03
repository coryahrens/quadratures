c--------------------------------------------------------------------------------------------c
c        C.A. June 27 2008 
c        Last edit: 6/27/08
c
c        Purpose: Calculate the condition number using SVD routine from Lapack
c
c        Input: Av -- vector of length n**2 -- represents square matrix of dimension n by n
c                n -- length of coloumn/row of matrix
c
c        Output: cnumber -- condition number
c--------------------------------------------------------------------------------------------c
c
         subroutine condnumber(Av,n,cnumber,sigma) 
c 
           implicit real *8 (a-h,o-z)
           real *8 Av(n*n), A(n,n), work(16*n), sigma(n)
c	 
	   lwork = 16*n
c
c----------Convert array into matrix for Lapack routine--------
c
           do i=1, n
             do j=1, n
 	       A(i,j) = Av( (i-1)*n + j )
 	     enddo
           enddo   
c
c----------Call SVD decomposition----------------
c              
           call dgesvd('N', 'N', n, n, A, n, sigma, u, n, v, n,
     1                 work, lwork, info)
           if (info.ne.0) then
             print *, 'Error in estimating cond #: info dgesvd',info           
           endif
c
c----------Calculate condition number----------
           cnumber = sigma(1)/sigma(n)
c   
           return
c	   
         end
