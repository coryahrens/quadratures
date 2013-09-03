       subroutine mvmult(intflag,m,n,A,x,y)
c-----------------------------------------------c
c        Matrix vector multiply
c        
c          A * x = y
c  
c        A is an m by n real matrix and x is an
c        n vector. The output is the m vector y
c
c        intflag = 0 <--> normal    --   A x = y
c        intflag = 1 <--> transpose -- A^T x = y
c 
c-----------------------------------------------c
         implicit real *8 (a-h,o-z)

         real *8 A(m,n), x(*), y(*)
c 
         if(intflag .eq.0 ) then
           do i=1,m
	     sum = 0.0d0
	     do j=1,n
	       sum = sum + A(i,j) * x(j)
	     enddo
	     y(i) = sum
	   enddo
	 else
           do i=1,n
	     sum = 0.0d0
	     do j=1,m
	       sum = sum + A(j,i) * x(j)
	     enddo
	     y(i) = sum
	   enddo	 
	 endif
         
	 return
      
       end
