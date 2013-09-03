c
c 
        program sphquad 
	 implicit real *8 (a-h,o-z)
        parameter(nneqmax=1750)
        parameter(ngenmax=1750)
c
        real *8 gen(4,ngenmax)
        real *8 theta0(ngenmax),phi0(ngenmax),wei0(ngenmax)

        real *8 theta(ngenmax),phi(ngenmax),wei(ngenmax)
        real *8 wvert0,wvert,wface,wside
	 real *8 wface0,wside0
	
	 real *8 f(nneqmax), jac(ngenmax*nneqmax)

 	 real *8 lambda, thres, zero, tRMS
	 real *8 x(120*ngenmax), y(120*ngenmax), z(120*ngenmax)
	
	 real *4 t1, t0, timeIter, ttime, avtime
	
	 integer ipn(nneqmax),ipm(nneqmax),itype,ifg,ngen	
	 integer inwflag, icflag, ipwflag
c
c     ----------------------------- 
        character*50 nfname, wfname, gfname
c     ----------------------------- 
c
        call prini(6,13)
        call getini
c      
        zero = 0.0d0
		
	 ctol = 1.0d-12
c
c------ Initial guess -------c
c
c       itype --- 1 = vertex + general points, 2 = vertex + face center
c                 3 = vertex + side midpoints, 4 = all three

	 isteps = 0
	
		
C        call inputfromfile(ngenmax,nmax,ngen,theta0,phi0,wei0,
C     1                  wvert0,wface0,wside0,itype,ifg,gen,nvar,nodes,
C     2                  io,nfname,wfname,gfname,lambda,thres,
C     3                  itest,isteps,ipwflag)  
     
        call inputdata(ngenmax,nmax,ngen,theta0,phi0,wei0,
     1                  wvert0,wface0,wside0,itype,ifg,gen,nvar,nodes,
     2                  io,nfname,wfname,gfname,lambda,thres,
     3                  itest,isteps,ipwflag)    
     
	 call sobolev(nmax,ipn,ipm,neq,ifg)
	
        call writesetup(nmax,ifg,nvar,nodes,neq,itype,lambda,
     1	                thres,itest,ctol,io,nfname,wfname,gfname,isteps)
         

	 call prinf('Starting Newton iterations...*', 0,0)  
	
	ttime = 0.0 	
c	
c-------Start Newton iteration--------c
c
	itermax = 50000
	icflag = 0
	
        do iter=1,itermax
c
	   call cpu_time(t0)

          call stepnewt_testing(neq,nvar,nmax,ngen,theta0,phi0,wei0,
     1                  wvert0,wface0,wside0,theta,phi,wei,wvert,wface,
     2                  wside,iter,ipn,ipm,itype,lambda,thres,
     3                  ifg,isteps,gen,icflag,ipwflag)
c
c---------Calculate RMS from previous iteration
c
          smth = 0.0d0
	   smph = 0.0d0
	   smwh = 0.0d0
	  
          do j=1,ngen
	  
	        smth = smth + (theta(j) - theta0(j))**2
	        smph = smph + (phi(j)   - phi0(j))**2
	        smwh = smwh + (wei(j)   - wei0(j))**2
	    
	    enddo
	    
	    smwh = smwh + (wvert - wvert0)**2 + (wface - wface0)**2 +
     1                  (wside - wside0)**2
c 
c---------Update angular variables and weights for generators
c
          do j=1,ngen
	  
            theta0(j) = theta(j)
            phi0(j)   = phi(j)
            wei0(j)   = wei(j)
	    
          enddo
c
c---------Update generators
c        
          do j=1,ngen
	  
            gen(1,j) = zero
            gen(2,j) = sin(theta0(j)) * cos(phi0(j))
            gen(3,j) = sin(theta0(j)) * sin(phi0(j))
            gen(4,j) = cos(theta0(j)) 
	      
          enddo
	  
	  wvert0 = wvert
	  wface0 = wface
	  wside0 = wside 
c
c---------Check for negative weights
c	  
	  if(ipwflag.eq.0) then
	    inwflag = 0
	    do j=1,ngen
	  
	      if(wei0(j).lt.0.0d0) then
	        inwflag = 1
	      endif
	    
	    enddo
	  endif	  	  
c
c---------Calculate function values and jacobian
c
          if(ipwflag.eq.0) then
c
c-----------weights can be both positive and negative
c
            call eqquad(nvar,nmax,ngen,gen,wei0,wvert0,wface0,wside0,
     1                  f,jac,neq,ipn,ipm,itype,ifg)
c     
          else
c
c-----------weights are of the form Exp(wei0) > 0
c	 
	    call eqquadpw(nvar,nmax,ngen,gen,wei0,wvert0,wface0,wside0,
     1                    f,jac,neq,ipn,ipm,itype,ifg)
c     
          endif	  
c
	   sum = zero
          do j=1,neq
            sum = sum + f(j)*f(j)
          enddo	  
	  
	  rms = sqrt(sum)	  
c
c---------Timing statistics
c
          call cpu_time(t1)	  
	  timeIter = t1 - t0	  
	  ttime  = ttime + timeIter
	  avtime = ttime / (1.0 * iter)	  
	    	  
c 	  if ((iter/ioutp)*ioutp.eq.iter) then
            call prinf('=====================================*',0,0)
	    call prinf('Iteration: *',iter,1)
	    call prinf('===RMS of variables in Newton step===*',0,0)
	    call prin2x('rms-theta = *', sqrt(smth),1)
	    call prin2x('rms-phi = *', sqrt(smph),1)
	    call prin2x('rms-weight = *', sqrt(smwh),1)
	    call prinf('=====================================*',0,0)	    
	    call prin2x('Function RMS = *',rms,1)
	    call prin('Average time/iteration = *', avtime,1)
	    if(ipwflag.eq.0) then
	      if((inwflag.eq.1).or.(wvert0.lt.0.0d0).or.
     1	         (wface0.lt.0.0d0).or.(wside0.lt.0.0d0)) then	    
	        call prin2x('g weights = *', wei0, ngen)
	        call prin2x('v weights = *', wvert0, 1)
	        call prin2x('f weight  = *', wface0, 1)
	        call prin2x('s weight  = *', wside0, 1)
	      endif
	    endif
c	  endif

c
c---------Check for convergence
c	  
          if(rms .lt. ctol) then
	    call prinf(' *',0,0)
	    call prinf('Convergence...*',0,0)
	    call prin2x('rms = *',rms,1)
	    call prinf('number of iterations = *',iter,1)
	    call prin('time to convergence = *',ttime,1)
	    call prinf('...exiting Newton iteration*',0,0)
c	
            icflag = 1
            call stepnewt_testing(neq,nvar,nmax,ngen,theta0,phi0,wei0,
     1                  wvert0,wface0,wside0,theta,phi,wei,wvert,wface,
     2                  wside,iter,ipn,ipm,itype,lambda,thres,
     3                  ifg,isteps,gen,icflag,ipwflag)
   
	    goto 10
c	    
	  endif
c
c---------Check for divergence
c
          tRMS = sqrt(smth) + sqrt(smph) + sqrt(smwh)
          if(tRMS.gt.1.0d2) then
	  
	    call prin2x('Divergence--RMS previous iterate:*',tRMS,1)
	    call prin2x('Stopping code*',0,0)
	    stop
	  
	  endif
c	  
        enddo
c
10      continue	

c
c-------End Newton iteration--------c
c
        if(ipwflag.eq.1) then
	
	  do j=1,ngen
	  
	    wei(j) = exp(wei(j))
	  
	  enddo
	  
	  wvert = exp(wvert)
	     
	  if(itype.eq.1) then
	 
	    wside = zero
	    wface = zero
	   
          elseif(itype.eq.2) then
	   
	    wface = exp(wface)
	    wside = zero
	   
	  elseif(itype.eq.3) then
	 
            wface = zero
	    wside = exp(wside)
	   
	  else
	 
            wface = exp(wface)
	    wside = exp(wside)	 
	 
	  endif
	
	endif
c
c-------Final results
c	
        call prinf('================Final result=================*',0,0)
        call prin2x('theta *', theta,ngen)
        call prin2x('phi *', phi,ngen)		
        call prin2x('generator weights *', wei,ngen)
	
	if(itype.eq.2) then
	
	  call prin2x('face weight *', wface,1)
	  
	elseif(itype.eq.3) then
	
	  call prin2x('side weight *', wside,1)
	  
	elseif(itype.eq.4) then
	
	  call prin2x('face weight *', wface,1)
	  call prin2x('side weight *', wside,1)
	  	
	endif	
        call prin2x('vertex weight *', wvert,1)
	call prinf('=============================================*',0,0)	
c
c-------Output of results
c

5000	if(io.eq.1) then

	  call prinf('============================== *',0,0)
	  call prinf('Preparing to output data files* ',0,0)

	  call toCart(ngen,gen,x,y,z,ifg)
	  
          call output(nfname,wfname,wvert,wface,wside,wei,
     1   	      ngen,itype,ifg,x,y,z,1)
     
     	  open(unit = 7, file = gfname, 
     1	       status = 'unknown', form = 'formatted')	  

	  do i=1,ngen
	  
	    write(7,2510) gen(2,i), gen(3,i), gen(4,i), wei(i)
	    
	  enddo

	  write(7,2520) wvert
	  
	  if(itype.eq.2) then
	
	    write(7,2520) wface
	  
	  elseif(itype.eq.3) then
	
	    write(7,2520) wside
	  
	  elseif(itype.eq.4) then
	
	    write(7,2520) wface
	    write(7,2520) wside
	  	
	  endif	
	    
	  close(7)
          
          call prinf('Done with data output* ',0,0)
	  call prinf('===================== *',0,0)
	  
	endif
c
c-------Check solution to see if it integrates correct subspace
c
5001    iflag = 0
	if(itest.eq.1) then
	
	  call prinf('Starting subspace test...*',0,0)
c  
          call quadtest(nmax,ngen,wei,wvert,wface,wside,gen,itype,
     1                  iflag,ifg)
     
	  if(iflag.eq.0) then
	  
	    call prinf('Test passed*',0,0)
	    
	  else
	  
	    call prinf('Test NOT passed*',0,0)
	    
	  endif
	
	else
	
	  call prinf('No testing of solution...*',0,0)
c	
	endif
c
		
1000    format(I5)
2000    format(E25.18) 
2500    format(3(2x,E25.18))
2510    format(4(2x,E25.18))
2520    format(1(2x,E25.18))

c
c
c-------Done for now
c        
	end

































