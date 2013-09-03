       subroutine quadtest(nmax,ngen,wei,wvert,wface,wside,
     1                     gen,itype,iflag,ifg) 
c---------------------------------------------------------------
c A check to see that the weights and nodes actually integrate
c the correct subspaces
c
c---------------------------------------------------------------
         implicit real *8 (a-h,o-z)
         real *8 wei(*),wvert,wface,wside
c	 real *8 theta(60,0:ngen),phi(60,0:ngen)

	 real *8 sumrv(0:nmax,0:nmax),sumiv(0:nmax,0:nmax)
	 real *8 sumrf(0:nmax,0:nmax),sumif(0:nmax,0:nmax)
	 real *8 sumrs(0:nmax,0:nmax),sumis(0:nmax,0:nmax)
	 real *8 sumrg(0:nmax,0:nmax),sumig(0:nmax,0:nmax)
	 
         real *8 p(0:nmax,0:nmax), dp(0:nmax,0:nmax)
         real *8 phire(0:nmax),phiim(0:nmax)
	 
	 real *8 vert(4,12),fcent(4,20),scent(4,30)
	 real *8 rg(4,60),rgc(4,60),qtemp(4,60),qres(4,60)
	 real *8 gen(4,ngen),rgen(4,ngen)
	 
	 real *8 tol, totalr, totali
	 
c	 real *8 theta0(ngen), phi0(ngen)
c
         zero  = 0.0d0
	 tol   = dlamch( 's' )
	 ctol  = 1.0d-13
c	 
	 if(ifg.eq.1) then
	   call prinf('Using full icosahedral group*',0,0)
	 else
           call prinf('Using rotation subgroup*',0,0)	
	 endif	 
c	 
c--------initialize sum
c        
         do l=0,nmax
	   do m=0,l
	     sumrv(l,m) = zero
	     sumiv(l,m) = zero
	     sumrf(l,m) = zero
	     sumif(l,m) = zero
	     sumrs(l,m) = zero
	     sumis(l,m) = zero
	     sumrg(l,m) = zero
	     sumig(l,m) = zero
	   enddo
	 enddo
c
c--------retrieve rotation group
c
         call rotg(rg)       
c
c--------compute conjugate quaternions
c
         do k=1,60
           rgc(1,k)   =  rg(1,k)
           do j=2,4
             rgc(j,k) = -rg(j,k)
           enddo 
         enddo
c
c--------retrieve vertices, face centers and side midpoints
c
         call vertex(vert)
	 call fcenter(fcent)
	 call scenter(scent)
c
c--------If using full group, reflect about the origin to get other generators
c
         if(ifg.eq.1) then
	 
           do j=1,ngen
	   
             rgen(1,j) =  zero
             rgen(2,j) = -gen(2,j)
             rgen(3,j) = -gen(3,j)
             rgen(4,j) = -gen(4,j)
	        
           enddo
	   
	 endif		 	 	 
c
c--------start loop over (l,m)
c
         iflag = 0
	 
         do l=0,nmax 
	   
	   do m=0,l
c
c------------compute contribution from vertices       
c
             do k=1,12	
	     
	       ctheta = vert(4,k)
               stheta = sqrt(vert(2,k)*vert(2,k) + 
     1	                       vert(3,k)*vert(3,k))
	       cphi   = vert(2,k)/stheta  
	       sphi   = vert(3,k)/stheta 
	              
               call spharm(nmax,nmax,nmax,ctheta,acos(ctheta),
     1	                   sphi,cphi,p,dp,phire,phiim,tol)
     
               sumrv(l,m) = sumrv(l,m) + p(m,l) * phire(m)
	       sumiv(l,m) = sumiv(l,m) + p(m,l) * phiim(m)
	       
             enddo
	      
             sumrv(l,m) = wvert * sumrv(l,m)
	     sumiv(l,m) = wvert * sumiv(l,m)
	     
c	     call prin2x('sumrv = *', sumrv,1)
c	     call prin2x('sumiv = *', sumiv,1)
c
             if(itype.eq.2) then
c--------------------------------------------------------------c
c--------------compute contribution from face centers----------c 
c--------------------------------------------------------------c
	       do k=1,20
	       	
	         ctheta = fcent(4,k)
                 stheta = sqrt(fcent(2,k)*fcent(2,k) + 
     1	                       fcent(3,k)*fcent(3,k))
	         cphi   = fcent(2,k)/stheta  
	         sphi   = fcent(3,k)/stheta
		    
c----------------compute values of spherical harmonics at face centers      
                 call spharm(nmax,nmax,nmax,ctheta,acos(ctheta),
     1	                     sphi,cphi,p,dp,phire,phiim,tol) 
     
                 sumrf(l,m) = sumrf(l,m) +  p(m,l) * phire(m)
	         sumif(l,m) = sumif(l,m) +  p(m,l) * phiim(m)
		 
               enddo 
	       
	       sumrf(l,m) = wface * sumrf(l,m)
	       sumif(l,m) = wface * sumif(l,m)
	       
	     elseif(itype.eq.3) then
c--------------------------------------------------------------c
c--------------compute contribution from side centers----------c      
c--------------------------------------------------------------c
               do k=1,30
	       
                 ctheta = scent(4,k)
		 
	         if((scent(2,k)*scent(2,k) + scent(3,k)*scent(3,k)).le.
     1		                                          1.0d-12) then
	           sphi = 0.0d0
		   
	           cphi = 1.0d0
		   
	         else
		 
                   sr = sqrt(scent(2,k)*scent(2,k) + 
     1		             scent(3,k)*scent(3,k))
                   sphi = scent(3,k)/sr
                   cphi = scent(2,k)/sr	
		        
	         endif	
		  
c----------------compute values of spherical harmonics at side centers   
                 call spharm(nmax,nmax,nmax,ctheta,acos(ctheta),
     1	                     sphi,cphi,p,dp,phire,phiim,tol)
     
                 sumrs(l,m) = sumrs(l,m) +  p(m,l) * phire(m)
	         sumis(l,m) = sumis(l,m) +  p(m,l) * phiim(m)
		 
               enddo
	       
	       sumrs(l,m) = wside * sumrs(l,m)
	       sumis(l,m) = wside * sumis(l,m)
	       	     
	     elseif(itype.eq.4) then
c--------------------------------------------------------------c
c--------------compute contribution from face centers----------c 
c--------------------------------------------------------------c
	       do k=1,20
	       	
	         ctheta = fcent(4,k)
                 stheta = sqrt(fcent(2,k)*fcent(2,k) + 
     1	                       fcent(3,k)*fcent(3,k))
	         cphi   = fcent(2,k)/stheta  
	         sphi   = fcent(3,k)/stheta  
		  
c----------------compute values of spherical harmonics at face centers      
                 call spharm(nmax,nmax,nmax,ctheta,acos(ctheta),
     1	                     sphi,cphi,p,dp,phire,phiim,tol) 
     
                 sumrf(l,m) = sumrf(l,m) +  p(m,l) * phire(m)
	         sumif(l,m) = sumif(l,m) +  p(m,l) * phiim(m)
		 
               enddo 
	       
	       sumrf(l,m) = wface * sumrf(l,m)
	       sumif(l,m) = wface * sumif(l,m)
	       	     
c--------------------------------------------------------------c
c--------------compute contribution from side centers----------c      
c--------------------------------------------------------------c
               do k=1,30
	       
                 ctheta = scent(4,k)
		 
	         if((scent(2,k)*scent(2,k) + scent(3,k)*scent(3,k)).le.
     1		                                          1.0d-10) then
	           sphi = 0.0d0
	           cphi = 1.0d0
		   
	         else
		 
                   sr = sqrt(scent(2,k)*scent(2,k) + 
     1		             scent(3,k)*scent(3,k))
                   sphi = scent(3,k)/sr
                   cphi = scent(2,k)/sr	
		        
	         endif	 
c----------------compute values of spherical harmonics at side centers   
                 call spharm(nmax,nmax,nmax,ctheta,acos(ctheta),
     1	                     sphi,cphi,p,dp,phire,phiim,tol)
     
                 sumrs(l,m) = sumrs(l,m) +  p(m,l) * phire(m)
	         sumis(l,m) = sumis(l,m) +  p(m,l) * phiim(m)
		 
               enddo
	       
	       sumrs(l,m) = wside * sumrs(l,m)
	       sumis(l,m) = wside * sumis(l,m)	
	            
	     endif
	     
c
             do j=1,60
	 
               qtemp(1,j) = zero
               qtemp(2,j) = zero
               qtemp(3,j) = zero
               qtemp(4,j) = zero   
	   
             enddo	     
c            
c------------compute contribution from generators       
c
             do j=1,ngen
c
               do k=1,60
                 call quatmult(gen(1,j),rgc(1,k),qtemp(1,k))	   
                 call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
               enddo
c
	       
               do k=1,60
		 
                 ct    = qres(4,k)
	         theta = acos(ct)
c	     	     
	         if((qres(2,k)*qres(2,k) + qres(3,k)*qres(3,k)) .le. 
     1	           1.0d-12) then
	           si = 0.0d0
	           co = 1.0d0
	         else
                   sr = sqrt(qres(2,k)*qres(2,k) + qres(3,k)*qres(3,k))
                   si = qres(3,k)/sr
                   co = qres(2,k)/sr	     
	         endif
c  
c--------------compute value of spherical harmonics and derivatives at image point for all orders and degrees
c	
	 
                 call spharm(nmax,nmax,nmax,ct,theta,si,co,
     1	                     p,dp,phire,phiim,tol)
         
                 sumrg(l,m) = sumrg(l,m) + wei(j)*p(m,l) * phire(m)
	         sumig(l,m) = sumig(l,m) + wei(j)*p(m,l) * phiim(m)
		 
	       enddo
c
c--------------If using the full icosahedral group, sum over other 60 points
c
               if(ifg.eq.1) then
c
                 do k=1,60
                   call quatmult(rgen(1,j),rgc(1,k),qtemp(1,k))
                   call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
                 enddo
c	       
                 do k=1,60
	           ctheta = qres(4,k)
                   stheta = sqrt(qres(2,k)*qres(2,k) +
     1  		         qres(3,k)*qres(3,k))
	           cphi   = qres(2,k)/stheta  
	           sphi   = qres(3,k)/stheta   
                   call spharm(nmax,nmax,nmax,ctheta,acos(ctheta),
     1	                       sphi,cphi,p,dp,phire,phiim,tol)
                   sumrg(l,m) = sumrg(l,m) + wei(j)*p(m,l) * phire(m)
	           sumig(l,m) = sumig(l,m) + wei(j)*p(m,l) * phiim(m)
	         enddo	
c		        
	       endif	       
c	       	       	       
	     enddo
c
c------------Check if quadrature is satisfied
c	 
             totalr = sumrv(l,m) + sumrf(l,m) + sumrs(l,m) + sumrg(l,m)
	     totali = sumiv(l,m) + sumif(l,m) + sumis(l,m) + sumig(l,m)  
	     if(((abs(totalr).gt.ctol).or. 
     1	        (abs(totali).gt.ctol)).and.
     2          (l.gt.0)) then
     
               call prinf('===Offending subspace===*',0,0)
               call prinf('========================*',0,0)
	       call prinf('l = *',l,1)
	       call prinf('m = *',m,1)
               call prin2x('TotalR = *',totalr,1)
	       call prin2x('TotalI = *',totali,1)
	       
	       call prin2x('sumrv  = *',sumrv(l,m),1)
	       call prin2x('sumiv  = *',sumiv(l,m),1)
	       call prin2x('sumrf  = *',sumrf(l,m),1)
	       call prin2x('sumif  = *',sumif(l,m),1)
	       call prin2x('sumrs  = *',sumrs(l,m),1)
	       call prin2x('sumis  = *',sumis(l,m),1)
	       call prin2x('sumrg  = *',sumrg(l,m),1)
	       call prin2x('sumig  = *',sumig(l,m),1)
	       iflag = 1
	       
	     else 
	     
	       iflag = 0
	       call prinf('===Subspace Okay===*',0,0)
	       call prinf('l = *',l,1)
	       call prinf('m = *',m,1)    	   	
c	       call prinf('===================*',0,0)
	       
	     endif
	     
	   enddo
	   
	 enddo
c

2500     format(E24.16)
	 
	 return
	 
       end
