       subroutine inputdata(ngenmax,nmax,ngen,theta0,phi0,wei0,
     1                  wvert0,wface0,wside0,itype,ifg,gen,nvar,nodes,
     2                  io,nfname,wfname,gfname,lambda,thres,itest,
     3                  isteps,ipwflag)
         implicit real*8 (a-h,o-z)
         real*8 theta0(ngenmax), phi0(ngenmax)
	 real*8 wei0(ngenmax), wvert0, wface0, wside0
	 real*8 gen(4,ngenmax), sr, sphi, lambda, thres
	 
	 integer itype, ifg, ngen, nmax, nvar, iweightsa
	 integer iweightsb, isteps, ipwflag
	 
	 character*50 nfname, wfname, gfname, gdata
	 
	 zero   = 0.0d0       
         pi     = 4.0d0*atan(1.0d0)
         fourpi = 4.0d0*pi	 
c
c--------Read type of quadrature,full or subgroup, # of generators, max degree 
c
         call prinf('Input quadrature type (1-4): *',0,0)
	 read(*,*) itype
	 
         call prinf('Input subgroup (0) or full group (1): *',0,0)
	 read(*,*) ifg
	 
	 call prinf('Input number of generators: *',0,0)
	 read(*,*) ngen
         
	 call prinf('Input maximum order/degree of subspace: *',0,0)
	 read(*,*) nmax
	 
	 call prinf('Output results: yes(1), no(0) *',0,0)
	 read(*,*) io
c
c--------Calculate number of nodes and equations
c	 
	 call numbernodeseqs(itype,ngen,ifg,nvar,nodes)	 
c	 	 
	 if(io.gt.0) then
	 
	   call makefilename(nmax,ngen,itype,nodes,nfname,l,0)	 
	   call prina('Quadrature node file name: *',nfname,l)
	   
           call makefilename(nmax,ngen,itype,nodes,wfname,l,2)	   
	   call prina('Quadrature node file name: *',wfname,l)

           call makefilename(nmax,ngen,itype,nodes,gfname,l,1)	   
	   call prina('Quadrature node file name: *',gfname,l)
	 
	 endif
	 	 
	 call prinf('Adaptive (1) or Fixed (0) step size: *',0,0)
	 read(*,*) isteps
	 
	 if(isteps.eq.0) then
	   call prinf('Damping parameter: *', 0,0)
	   read(*,*) lambda
	 
	   call prinf('Damping threshold: *', 0,0)
	   read(*,*) thres
	 endif
	 
	 call prinf('Force positive weights: yes(1), no(0)*',0,0)
	 read(*,*) ipwflag
	 
	 call prinf('Test final result: yes(1), no(0)*',0,0)
	 read(*,*) itest
	 
	 call prinf('Use uniform initial weights: yes(1), no(0)*',0,0)
	 read(*,*) iweightsa	 
c
c--------Read generator data	 
c	 

CCC         call getfilename(gdata,'Generator data file: *')

	 call prinf('Generator file name:*',0,0)
	 read(*,*) gdata
	 
	 call openfile(gdata,7,'old')
	 
	 if(iweightsa.eq.1) then	 
c
c----------Read only generator data
c
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j) 
	   enddo
	   close(7)	 

c
c----------Set uniform initial weights
c
           if(ipwflag.eq.0) then
	   
             do j=1,ngen
	       wei0(j) = fourpi/(1.0d0*nodes)
	     enddo
c 	 
             wvert0 = fourpi/(1.0d0*nodes)
	     if(itype.eq.1) then
	 
	       wside0 = zero
	       wface0 = zero
	   
             elseif(itype.eq.2) then
	   
	       wface0 = fourpi/(1.0d0*nodes)
	       wside0 = zero
	   
	     elseif(itype.eq.3) then
	 
               wface0 = zero
	       wside0 = fourpi/(1.0d0*nodes)
	   
	     else
	 
               wface0 = fourpi/(1.0d0*nodes)
	       wside0 = fourpi/(1.0d0*nodes)	 
	 
	     endif
	     
	   else
	   
             do j=1,ngen
	       wei0(j) = log(fourpi/(1.0d0*nodes))
	     enddo
c 	 
             wvert0 = log(fourpi/(1.0d0*nodes))
	     
	     if(itype.eq.1) then
	 
	       wside0 = zero
	       wface0 = zero
	   
             elseif(itype.eq.2) then
	   
	       wface0 = log(fourpi/(1.0d0*nodes))
	       wside0 = zero
	   
	     elseif(itype.eq.3) then
	 
               wface0 = zero
	       wside0 = log(fourpi/(1.0d0*nodes))
	   
	     else
	 
               wface0 = log(fourpi/(1.0d0*nodes))
	       wside0 = log(fourpi/(1.0d0*nodes))	 
	 
	     endif	   
	   
	   
	   endif
	 
	 else
	 
	   call prinf('Use new data format: yes(1), no(0)*',0,0)
	   read(*,*) iweightsb
	   
	   if(iweightsb.eq.1) then
	   
	     if(ipwflag.eq.0) then
c
c--------------Read generator and weight data in new format
c
	       do j=1,ngen
	         gen(1,j) = zero
	         read(7,*)  gen(2,j), gen(3,j), gen(4,j), wei0(j)  
	       enddo     	     
	   
	       read(7,*) wvert0

               if(itype.eq.1) then
	       
	         wface0 = 0.0d0
		 wside0 = 0.0d0

	       elseif(itype.eq.2) then
	
	         read(7,*) wface0
		 
		 wside0 = 0.0d0
	  
	       elseif(itype.eq.3) then
	
	         read(7,*) wside0
		 
		 wface0 = 0.0d0
	  
	       elseif(itype.eq.4) then
	
	         read(7,*) wface0
	         read(7,*) wside0
	  	
	       endif	
	     
	       close(7) 
	     
	     else
c
c--------------Read generator and weight data in new format
c
	       do j=1,ngen
	         gen(1,j) = zero
	         read(7,*)  gen(2,j), gen(3,j), gen(4,j), wei0(j)  
		 wei0(j) = log(wei0(j))
	       enddo     	     
	   
	       read(7,*) wvert0
	       wvert0 = log(wvert0)
	       
	       if(itype.eq.1) then
	       
	         wface0 = 0.0d0
		 wside0 = 0.0d0

	       elseif(itype.eq.2) then
	
	         read(7,*) wface0
		 wface0 = log(wface0)
		 
		 wside0 = 0.0d0
	  
	       elseif(itype.eq.3) then
	
	         read(7,*) wside0
		 wside0 = log(wside0)
		 
		 wface0 = 0.0d0
	  
	       elseif(itype.eq.4) then
	
	         read(7,*) wface0
	         read(7,*) wside0
		 
		 wface0 = log(wface0)
		 wside0 = log(wside0)
	  	
	       endif	
	     
	       close(7)	     
	     
	     endif     
   
	   else
c
c------------Read generator and weight data in old format
c
	     do j=1,ngen
	       gen(1,j) = zero
	       read(7,*)  gen(2,j), gen(3,j), gen(4,j) 
	     enddo
	     
	     close(7)

	     call prinf('Weight data file:*',0,0)
	     read(*,*) wfname

	     open (unit=3,file=wfname,status='old',iostat=ios,
     1             form='formatted')
             if(ios.ne.0) then
	     
	       write(*,*) 'Error opening file'
	       stop
	    
	     endif
	     
	     if(itype.eq.1) then
	   
	       do i=1,12
	         read(3,*) wvert0
	       enddo
	     	     
               do i=1, ngen
	         do j=1, 60
	           read(3,*)  wei0(i)
	         enddo
	       enddo
	       
	       wface0 = 0.0d0
	       wside0 = 0.0d0

	     elseif(itype.eq.2) then
	     
	       do i=1,12
	         read(3,*) wvert0
	       enddo
	       
	       do i=1,20
	         read(3,*) wface0
	       enddo
	
               do i=1, ngen
	         do j=1, 60
	           read(3,*)  wei0(i)
	         enddo
	       enddo
	       
	       wside0 = 0.0d0	       
	  
	     elseif(itype.eq.3) then
	
	       do i=1,12
	         read(3,*) wvert0
	       enddo
	       
	       do i=1,30
	         read(3,*) wside0
	       enddo
	
               do i=1, ngen
	         do j=1, 60
	           read(3,*)  wei0(i)
	         enddo
	       enddo
	       
	       wface0 = 0.0d0
	  
	     elseif(itype.eq.4) then
	
	       do i=1,12
	         read(3,*) wvert0
	       enddo

	       do i=1,20
	         read(3,*) wface0
	       enddo
	       
	       do i=1,30
	         read(3,*) wside0
	       enddo
	
               do i=1, ngen
	         do j=1, 60
	           read(3,*)  wei0(i)
	         enddo
	       enddo
	  	
	     endif	      
	   
	     close(3)	  	   
	   
	   endif	 

	 endif
c
c--------Calculate angles
c	   
	 do j=1,ngen
           theta0(j) = acos(gen(4,j))
           sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
           sphi      = gen(3,j)/sr
           phi0(j)   = dasin(sphi)
         enddo
       
         return
       
       end
       
