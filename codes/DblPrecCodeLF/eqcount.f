c
c*********************************************************************************************       
c      Calculates the dimension of the subspace of degree n which is invariant under the     *
c      icosahedral group. The total number of equations needed for determining a quadrature  *
c      is then found by summing the number of equations for each subspace.                   *
c                                                                                            *
c	input:                                                                               *
c                                                                                            *
c	nmax	---	maximum degree of spherical harmonics to integrate                   *
c	neq	---	total number of equations                                            *
c                                                                                            *
c	output:                                                                              *
c                                                                                            *
c	iflag	---	0 if a match and 1 otherwise                                         *
c                                                                                            *
c      Reference: S. L. Sobolev, "Cubature formulas on the sphere invariant under finite     *
c                 groups of rotations," Dokl. Akad. Nauk SSR, 146 310-313 (1962)             *
c*********************************************************************************************      
         subroutine eqcount(nmax,neq,icount,iflag)
           integer neqsub(0:nmax), q(3), t(3)
           integer neq
           integer nmax, sum, M
c
c----------the size and the total number of equations neq, should be estimated 
c----------outside this routine.  This is only a consitency check.  
c
c         
c          Data for the icosahedral group
c
	   M = 60
           data q(1), q(2), q(3) /5, 3, 2/
	   data t(1), t(2), t(3) /12, 20, 30/
c
c----------Calculate s(j)
c
	   do j = 0, nmax
c
c----------Calculate the sum over the set Q-star
c
	     sum = 0
	     do i = 1, 3
	       if(mod(j,q(i)).ne.0) then
	         sum = sum + t(i)
               end if
	     enddo
c
c----------Check conditions of eq.(17)
c
             if(2*j + 1 .le. sum) then
	         neqsub(j) = int((2*j+1)/M)
	     else
	         neqsub(j) = int((2*j+1)/M) + 1
	     end if	     
           enddo
c
           icount = 0
	   do j = 0, nmax
             if (neqsub(j).ne.0) then
               do l=1,neqsub(j)
                 icount = icount+1
               enddo
             endif
           enddo
c
           if (neq.eq.icount) then
             iflag = 0
           else
             iflag = 1
           endif
c	 
         return
       end
