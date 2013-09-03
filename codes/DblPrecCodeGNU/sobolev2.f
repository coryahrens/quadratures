c
c*********************************************************************************************       
c      Calculates the dimension of the subspace of degree n which is invariant under the     *
c      icosahedral group. The total number of equations needed for determining a quadrature  *
c      is then found by summing the number of equations for each subspace.                   *
c                                                                                            *
c	input:                                                                               *
c                                                                                            *
c	nmax	---	maximum degree of spherical harmonics to integrate                   *
c                                                                                            *
c	output:                                                                              *
c                                                                                            *
c	neqsub	---	vector with total number of equations for subspace of degree   n     *
c       ipn	---	pointer to contributing equation degree	                             *
c       ipm	---	pointer to contributing equation order	                             *
c	neq	---	total number of equations                                            *
c                                                                                            *
c                                                                                            *
c                                                                                            *
c      Reference: S. L. Sobolev, "Cubature formulas on the sphere invariant under finite     *
c                 groups of rotations," Dokl. Akad. Nauk SSR, 146 310-313 (1962)             *
c*********************************************************************************************      
         subroutine sobolev2(nmax,ipn,ipm,neq,ifg)
           integer neqsub(0:nmax), q(3), t(3)
           integer ipn(*), ipm(*), neq
           integer nmax, sum, M, ifg
c
c--------- the size of ipn and ipm, the total number of equations neq, should be estimated 
c--------- outside this routine.  Using nmax^2 is a good rough over-estimation
c--------- A tighter one is ((nnmax+1)/30+1)*(nnmax+1)/2+1 (Mathematica check up to nnmax=100000)  
c
c         
c          Data for the icosahedral group
c
	   M = 60
           data q(1), q(2), q(3) / 5,  3,  2/
	   data t(1), t(2), t(3) /12, 20, 30/
	
	 
	   if(ifg.eq.0) then
c
c------------Use only pure rotation subgroup of order 60
c	 
c------------Calculate s(j)
c
	     do j = 0, nmax
	     
	       neqsub(j) = int((1.0d0*j)/(1.0d0*q(1))) +
     1	                   int((1.0d0*j)/(1.0d0*q(2))) +
     2			   int((1.0d0*j)/(1.0d0*q(3))) - j + 1
c
c------------Calculate the sum over the set Q-star
c
c	       sum = 0
c	       do i = 1, 3
c	         if(mod(j,q(i)).ne.0) then
c	           sum = sum + t(i)
c                 endif
c	       enddo
c
c------------Check conditions of eq.(17)
c
c              if(2*j + 1 .le. sum) then
c	         neqsub(j) = int((2*j+1)/M)
c	       else
c	         neqsub(j) = int((2*j+1)/M) + 1
c	       endif
c	     	     
             enddo
c
c------------Construct pointers to indices n,m of spherical harmonics that 
c------------supply non-trivial equations
c
             icount = 0
	     
	     do j = 0, nmax
c
               if (neqsub(j).ne.0) then
c
c----------------Nontrivial subspace
c
                 do l=-j,j
c
                   icount = icount+1
		   
                   ipn(icount) = j
   
                   ipm(icount) = l
c
                 enddo
c
               endif
c
             enddo
c
             neq = icount
c
	   else
c
c------------Use full group of order 120--NEED TO BE REDONE !!!!! 5/22/09
c	 
c------------Calculate s(j)
c
	     do j = 0,nmax,2
c
c------------Calculate the sum over the set Q-star
c
	       sum = 0
	       do i = 1, 3
	         if(mod(j,q(i)).ne.0) then
	           sum = sum + t(i)
                 endif
	       enddo
c
c------------Check conditions of eq.(17)
c
               if(2*j + 1 .le. sum) then
	         neqsub(j) = int((2*j+1)/M)
	       else
	         neqsub(j) = int((2*j+1)/M) + 1
	       endif
c	     	     
             enddo
c
c------------Construct pointers to indices n,m of spherical harmonics that 
c------------supply non-trivial equations
c
             icount = 0
	     do j = 0, nmax
c
               if (neqsub(j).ne.0 .and. mod(j,2).eq.0) then
c
                 do l=1,neqsub(j)
c
                   icount = icount+1
                   ipn(icount) = j
                   if (mod(j,2).eq.0) then
                     ipm(icount) = 2*(l-1)
                   else
                     ipm(icount) = 2*l
                   endif
c
                 enddo
c
               endif
c
             enddo
c
             neq = icount
c	 
	   endif
c	 


           return
c
         end
