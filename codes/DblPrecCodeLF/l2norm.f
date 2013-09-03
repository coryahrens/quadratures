       subroutine l2norm(n,x,norm)
c----------------------------c
c        l2 norm of a vector c
c----------------------------c
         implicit real *8 (a-h,o-z)

         real *8 x(n),norm
c 
	 sum = 0.0d0
	 do j=1,n
	   sum = sum + x(j)*x(j)
	 enddo
 
         norm = sqrt(sum)
 
	 return
      
       end
