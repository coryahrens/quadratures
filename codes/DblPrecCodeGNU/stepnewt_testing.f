        subroutine stepnewt_testing(neq,nrhs,nnmax,ngen,theta0,phi0,
     1             wei0,wvert0,wface0,wside0,theta,phi,wei,wvert,wface,
     2             wside,iter,ipn,ipm,itype,lambda,thres,
     3             ifg,isteps,gen,icflag,ipwflag)
         implicit real *8 (a-h,o-z)
         real *8 theta0(*), phi0(*), wei0(*), wvert0
         real *8 theta(*),  phi(*),  wei(*),  wvert

         real *8 gen(4,*)
         real *8 fval(neq),fvalNew(neq),jac(neq*nrhs),jacTmp(neq*nrhs)
         real *8 work(20*neq),rhs(nrhs), ynew(nrhs)
	 real *8 cnumber, sigma(neq), lambda
	 real *8 jacNew(neq*nrhs), ip

         integer jpvt(nrhs), itmax, maxdim
	 integer ipn(*), ipm(*), icflag, ipwflag
	 logical sdec, curve
	 
         itmax = 15
	 rpara = 2.0d0*dmacheps() 	 
         zero  = 0.0d0
c
c--------Calculate function values and jacobian
c
         if(ipwflag.eq.0) then
c
c----------weights can be both positive and negative
c	   
	   call eqquad(nrhs,nnmax,ngen,gen,wei0,wvert0,wface0,wside0,
     1                 fval,jac,neq,ipn,ipm,itype,ifg)
c     
         else
c
c----------weights are of the form Exp(wei0) > 0
c	 
	   call eqquadpw(nrhs,nnmax,ngen,gen,wei0,wvert0,wface0,wside0,
     1                   fval,jac,neq,ipn,ipm,itype,ifg)
c     
         endif
c	 
         do i=1,neq*nrhs
	   jacTmp(i) = jac(i)
	 enddo  	
c
c--------Norm at current step
c	
	 sum = zero
         do j=1,neq
           rhs(j) = fval(j)
           sum    = sum + fval(j)*fval(j)
         enddo
         sum = sqrt(sum)

         if(iter.eq.1) then	 
	   call prin2x('Initial rms = *', sum, 1)
	 endif

         itest = mod(iter,10)
	 
         if((itest.eq.0) .or. (iter.eq.1) .or. (icflag.eq.1) ) then
c
c----------Call SVD decomposition----------------
c              
           call dgesvd('N', 'N', neq, nrhs, jacTmp, neq, sigma, 
     1                  u, 1, v, 1, work, 20*neq, info)
           if (info.ne.0) then
             print *, 'Error in estimating cond #: info dgesvd',info           
           endif
c
c----------Calculate condition number----------
c
           cnumber = sigma(1)/sigma(min(neq,nrhs))
	   
	   if(icflag.eq.1) then
	     call prin2x('Singular values: *',sigma,min(neq,nrhs))
	     call prin2x('Condition number: *',cnumber,1)
	     call prinf('=====================================*',0,0)
             goto 10
	   else
	     call prinf('Iteration: *', iter, 1)
	     call prin2x('Condition number: *',cnumber,1)
	     call prinf('=====================================*',0,0)	  	     
	   endif

           do i=1,neq*nrhs
	     jacTmp(i) = jac(i)
	   enddo 
        	   	  
	 endif
c
c--------Parameters for Lapack QR routine
c
         do j=1,nrhs
	   jpvt(j) = 0
	 enddo
c	 
	 maxdim = max(neq,nrhs)
c
c--------Solve for Newton step using RR-QR
c 	 
         call DGELSY(neq,nrhs,1,jacTmp,neq,rhs,maxdim,jpvt,rpara,irank,
     1               work,20*neq,info)
c   
	 if(info.ne.0) then
	   call prinf('Error in linear solve: info = *',info,1)
	   stop
	 endif
c  
         if(irank .lt. min(neq,nrhs)) then
	   call prinf('Jacobian rank .lt. min(neq,nrhs) *',irank,1)
	   call prinf('Info from solve: *',info,1)
	 endif  
c	 	 	 
	 if(isteps.eq.1) then   
c
c----------Start iteration for line search using strong Wolfe conditions
c
           alpha = 1.0d0
c           
	   call prinf('Staring line search...*',0,0)
c	            	 	 
	   do i=1,itmax
c	
c------------Update angle and weight variables and generators
c	  
             do j=1,ngen
c	   
CC               theta(j) = theta0(j) - alpha*rhs(3*j-2)
CC               phi(j)   = phi0(j)   - alpha*rhs(3*j-1) 
CC               wei(j)   = wei0(j)   - alpha*rhs(3*j) 

               theta(j) = theta0(j) - alpha*rhs(         j)
               phi(j)   = phi0(j)   - alpha*rhs(  ngen + j) 
               wei(j)   = wei0(j)   - alpha*rhs(2*ngen + j)
c	     	     
	       gen(1,j) = zero
               gen(2,j) = sin(theta(j))*cos(phi(j))
               gen(3,j) = sin(theta(j))*sin(phi(j))
               gen(4,j) = cos(theta(j))
c	     
             enddo	   	  	  		  
c
c------------Update icosahedral and special points
c  
             wvert = wvert0 - alpha*rhs(nrhs) 
c	     
	     if(itype.eq.1) then
c	       
	       wface = zero
	       wside = zero
c
	     elseif(itype.eq.2) then
c	     
	       wface = wface0 - alpha*rhs(nrhs-1)	       
	       wside = zero
c	       
	     elseif(itype.eq.3) then
c	     
	       wside = wside0 - alpha*rhs(nrhs-1)
	       wface = zero
c	       
	     elseif((itype.eq.4).or.(itype.eq.5)) then
c	     
	       wside = wside0 - alpha*rhs(nrhs-1)
	       wface = wface0 - alpha*rhs(nrhs-2)
c	       
	     endif
c
c------------Calculate Function/Jacobian values at x_k - alpha*<(J_k)^-1, F_k>
c
             if(ipwflag.eq.0) then
c
c--------------weights can be both positive and negative
c	   
	       call eqquad(nrhs,nnmax,ngen,gen,wei,wvert,wface,wside,
     1                     fvalNew,jacNew,neq,ipn,ipm,itype,ifg)
c     
             else
c
c--------------weights are of the form Exp(wei0) > 0
c	 
               call eqquadpw(nrhs,nnmax,ngen,gen,wei,wvert,wface,wside,
     1                       fvalNew,jacNew,neq,ipn,ipm,itype,ifg)
c     
             endif 	   	   
c
c------------Test for sufficient decrease
c
             rhsnorm = DNRM2(neq,fvalNew,1)
c	   
	     factor = sqrt(1.0d0 -  1.0d-4 * alpha)
c 
             sdec = rhsnorm .lt. (factor * sum)
c
c------------Test for curvature 
c 	   
	     call DGEMV('T',neq,nrhs,1.0d0,jacNew,neq,fvalNew,
     1    	        1,0.0d0,ynew,1)
c            
	     ip = DDOT(nrhs,rhs,1,ynew,1)
c
             curve = (0.9d0*sum) .ge. abs(ip)
c	
	     if(sdec .and. curve) then
c	     
	       goto 10
c	     
             else
c	   
	       alpha = 0.5d0*alpha
c	   
	     endif
c	 	   
           enddo
c	 
	 else
c
c----------Use fixed step size of lambda for ||r|| > threshold
c     	 
           if (sum.le.thres) then
c
c------------Set damping factor to unity -- pure Newton
c
             alpha = 1.0d0
c	     
           else
c
c------------Use fixed size lambda
c   
	     alpha = lambda
c	   
	   endif
c	
c----------Update angle and weight variables and generators
c	   	  
           do j=1,ngen
c	   
CC             theta(j) = theta0(j) - alpha*rhs(3*j-2)
CC             phi(j)   = phi0(j)   - alpha*rhs(3*j-1) 
CC             wei(j)   = wei0(j)   - alpha*rhs(3*j) 

             theta(j) = theta0(j) - alpha*rhs(         j)
             phi(j)   = phi0(j)   - alpha*rhs(  ngen + j) 
             wei(j)   = wei0(j)   - alpha*rhs(2*ngen + j)
c	     
           enddo	   	  	  		  
c
c----------Update icosahedral and special points
c  
           wvert = wvert0 - alpha*rhs(nrhs) 
	   
	   if(itype.eq.1) then
	       
	     wface = zero
	     wside = zero
	   
	   elseif(itype.eq.2) then
	   
	     wface = wface0 - alpha*rhs(nrhs-1)
	     wside = zero
	     
	   elseif(itype.eq.3) then
	   
	     wside = wside0 - alpha*rhs(nrhs-1)
	     wface = zero
	     
	   elseif((itype.eq.4).or.(itype.eq.5)) then
	   
	     wside = wside0 - alpha*rhs(nrhs-1)
	     wface = wface0 - alpha*rhs(nrhs-2)
	     
	   endif
	   	 	 
	 endif
	 	 
10       continue

         if(isteps.eq.1) then
          call prin2x('Step size used *',alpha,1)
	 endif	
	 	 
	 return	 
	
       end


