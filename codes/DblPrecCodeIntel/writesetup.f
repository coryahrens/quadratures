       subroutine writesetup(nmax,ifg,nvar,nodes,neq,itype,lambda,
     1                       thres,itest,ctol,io,nfname,wfname,gfname,
     2                       isteps)
         implicit real *8 (a-h,o-z)
	 real *8 lambda, thres, eff, ctol
	 integer itype, ifg, nmax, nvar, neq, io, isteps
	 
	 character*24 nfname, wfname, gfname
	 
	 eff = (nmax + 1)**2 / (3.0d0 * nodes)
	 
	 call prinf(' *',0,0)
         call prinf('======Summary of input data=============*', 0, 0)
	 	 	 
	 call prinf('Maximum order/degree: *',nmax,1)
	 call prinf('Number of unknowns: *',nvar,1)
	 call prinf('Number of equations: *',neq,1)
	 call prinf('Quadrature type: *', itype,1)
	 if(ifg.eq.1) then
	   call prinf('Using the full icosahedral group*',0,0)
	 else
           call prinf('Using the rotation sub-group*',0,0)	
	 endif
	 call prinf('Number of quadrature pnts: *',nodes,1)
	 if(isteps.eq.0) then
	   call prinf('Using fixed step size*',0,0)
	   call prin2x('  Damping: *',lambda,1)
	   call prin2x('  Threshold: *',thres,1)
	 else
	   call prinf('Using adaptive step size*',0,0)
	 endif
	 if(itest.eq.1) then
	   call prinf('Testing final results*',0,0)
	 else
           call prinf('No testing of final results*',0,0)	
	 endif
	 call prin2x('Convergence tolerance: *',ctol,1)	 
	 if(io.gt.0) then 
	   call prina('Quadrature node file name: *',nfname,24)
	   call prina('Quadrature weight file name: *',wfname,24)
	   call prina('Quadrature generator file name: *',gfname,24)
	 else
	   call prinf('No output of results*',0,0)
	 endif	 	 
	 call prin2x('Efficiency: *', eff,1)	
	  
         call prinf('======End summary of input data=========*', 0, 0)	
	 call prinf(' *',0,0)
	        
         return
       
       end
