       subroutine eqquad(nrhs,nnmax,ngen,gen,wei,wvert,wface,wside, 
     1               fval,jac,neq,ipn,ipm,itype,ifg)
c=======================================================================c 
c	 ngen   --- total number of generators
c        iqtype --- 1 = vertex + general points, 2 = vertex + face center
c                   3 = vertex + side midpoints, 4 = all three
c   
c 
c========================================================================c
         implicit real *8 (a-h,o-z)
         real *8 wei(ngen), wvert, gen(4,ngen), rgen(4,ngen)
         real *8 rg(4,60),rgc(4,60),qtemp(4,60),qres(4,60),rqres(4,60)
	 real *8 spnt(4,1)
         real *8 vert(4,12),fcent(4,20),scent(4,30)
         real *8 p(0:nnmax,0:nnmax),dp(0:nnmax,0:nnmax)
         real *8 phire(0:nnmax),phiim(0:nnmax)
	 real *8 cgenefg
c 
         real *8 fval(neq),vval(neq),cval(neq),sval(neq),qval(neq)

         real *8 pdth(neq),pdph(neq)
         real *8 ftheta1(60),fphi1(60)
         real *8 ftheta2(60),fphi2(60)
         real *8 jac(*), mlt

         integer ipn(*),ipm(*)
	 integer itype,ifg,ngen
c
c
         zero    =  0.0d0
         pi      =  4.0d0*atan(1.0d0)
c
c--------value of partial wrt icosahedral vertex weight
         cvert   =  6.0d0/sqrt(pi)
c
c--------value of partial wrt icosahedral face center weight	
         cface   = 10.0d0/sqrt(pi)
c
c--------value of partial wrt icosahedral side center weight	
         cside   = 15.0d0/sqrt(pi)
c
c--------value of partial wrt generator weights	for rotation subgroup
         cgene   = 30.0d0/sqrt(pi)
c
c--------value of partial wrt generator weights	for full group
         cgenefg = 60.0d0/sqrt(pi)	 
c
c--------normalization for (0,0)th equation
         crhs    =  2.0d0*sqrt(pi)
	 
         toll = dlamch( 's' )
	 
	 if(ifg.eq.0) then
	 
	   mlt = 1.0d0
	   
	 else
	 
	   mlt = 2.0d0
	   
	 endif
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
c--------initialization
c
         do i=1,neq
	 
           fval(i)    = zero
           vval(i)    = zero
	   cval(i)    = zero
	   sval(i)    = zero
	   
         enddo
	 
c----------------------------------------------------c
c--------compute contribution from vertices----------c      
c----------------------------------------------------c

         do k=1,12
c	 
c----------icosahedral vertex co-ordinates
c
           ct    = vert(4,k)
	   theta = acos(ct)
           sr    = sqrt(vert(2,k)*vert(2,k) + vert(3,k)*vert(3,k))
           si    = vert(3,k)/sr
           co    = vert(2,k)/sr	
c	    
c----------compute values of spherical harmonics at vertex  
c   
           call spharm(nnmax,nnmax,nnmax,ct,theta,si,co,p,
     1	                dp,phire,phiim,toll)         
c        
c----------sum over each equation     
c
           do i=1,neq
             vval(i) = vval(i) + p(ipm(i),ipn(i)) * phire(ipm(i))
           enddo
	   
         enddo		    
c
c--------multiply by the vertex weight 
c  
         do i=1,neq
	 
           fval(i) = vval(i)*wvert
	   
         enddo 
c
c--------Jacobian (derivative wrt "wvert" in the last column)
c	 
         call toarray(neq,nrhs,jac,1,nrhs,cvert)
c	 
         do i=2,neq
	 
           call toarray(neq,nrhs,jac,i,nrhs,vval(i))
	   
         enddo
c----------------------------------------------------------c
c--------check which type of quadrature is to be used------c
c----------------------------------------------------------c	  
         if(itype.eq.2) then
c----------------------------------------------------------c
c----------compute contribution from face centers----------c      
c----------------------------------------------------------c
           do k=1,20
c------------icosahedral face center coordinates
             ct    = fcent(4,k)
             theta = acos(ct)
             sr    = sqrt(fcent(2,k)*fcent(2,k) + fcent(3,k)*fcent(3,k))
             si    = fcent(3,k)/sr
             co    = fcent(2,k)/sr
	     	 
c------------compute values of spherical harmonics at face centers   
             call spharm(nnmax,nnmax,nnmax,ct,theta,si,co,p,
     1                     dp,phire,phiim,toll)
     
c------------sum over each equation      
             do i=1,neq
               cval(i) = cval(i) + p(ipm(i),ipn(i)) * phire(ipm(i))
             enddo
	     
           enddo
c
c----------multiply by the face center weight and add to running total
c
           do i=1,neq
	   
             fval(i) = fval(i) + cval(i) * wface
	     
           enddo
c
c--------Jacobian (derivative wrt "wface" in the second to last column)
c
           call toarray(neq,nrhs,jac,1,nrhs-1,cface)
c	 
           do i=2,neq
	   
             call toarray(neq,nrhs,jac,i,nrhs-1,cval(i))
	     
           enddo		      
	     	 
	 elseif(itype.eq.3) then
c----------------------------------------------------------c
c----------compute contribution from side centers----------c      
c----------------------------------------------------------c
           do k=1,30
	   
c------------icosahedral side center coordinates
             ct    = scent(4,k)
	     theta = acos(ct)
	     if((scent(2,k)*scent(2,k) + scent(3,k)*scent(3,k)) .le. 
     1 	         1.0d-10) then
	       si = 0.0d0
	       co = 1.0d0
	     else
               sr = sqrt(scent(2,k)*scent(2,k) + scent(3,k)*scent(3,k))
               si = scent(3,k)/sr
               co = scent(2,k)/sr	     
	     endif	
	      
c------------compute values of spherical harmonics at side centers   
             call   spharm(nnmax,nnmax,nnmax,ct,theta,si,co,p,
     1                     dp,phire,phiim,toll)
     
c------------sum over each equation    
             do i=1,neq	       
               sval(i) = sval(i) + p(ipm(i),ipn(i)) * phire(ipm(i))
             enddo
	     
           enddo
c
c----------multiply by the face center weight and add to running total
c
           do i=1,neq
             fval(i) = fval(i) + sval(i) * wside
           enddo
c
c----------Jacobian (derivative wrt "wside" in the second to last column)
c
           call toarray(neq,nrhs,jac,1,nrhs-1,cside)
c	 
           do i=2,neq
             call toarray(neq,nrhs,jac,i,nrhs-1,sval(i))
           enddo
	   
	 elseif(itype.eq.4) then
c
c-------------------------------------------------------------------c
c----------compute contribution from side and face centers----------c      
c-------------------------------------------------------------------c 
c
c----------Face centers-----------c
c
           do k=1,20
             ct    = fcent(4,k)
             theta = acos(ct)
             sr    = sqrt(fcent(2,k)*fcent(2,k) + fcent(3,k)*fcent(3,k))
             si    = fcent(3,k)/sr
             co    = fcent(2,k)/sr	 
c------------compute values of spherical harmonics at face centers   
             call   spharm(nnmax,nnmax,nnmax,ct,theta,si,co,p,
     1                     dp,phire,phiim,toll)
c------------sum over each equation      
             do i=1,neq
               cval(i) = cval(i) + p(ipm(i),ipn(i)) * phire(ipm(i))
             enddo
           enddo
c
c----------Side centers-----------c
c
           do k=1,30
             ct    = scent(4,k)
	     theta = acos(ct)
	     if((scent(2,k)*scent(2,k) + scent(3,k)*scent(3,k)) .le. 
     1	         1.0d-10) then
	       si = 0.0d0
	       co = 1.0d0
	     else
               sr = sqrt(scent(2,k)*scent(2,k) + scent(3,k)*scent(3,k))
               si = scent(3,k)/sr
               co = scent(2,k)/sr	     
	     endif	 
c------------compute values of spherical harmonics at side centers   
             call spharm(nnmax,nnmax,nnmax,ct,theta,si,co,p,
     1                     dp,phire,phiim,toll)
c------------sum over each equation    
             do i=1,neq	       
               sval(i) = sval(i) + p(ipm(i),ipn(i)) * phire(ipm(i))
             enddo
           enddo	   
c
c----------multiply by weights and add to running total
c
           do i=1,neq
             fval(i) = fval(i) + sval(i) * wside + cval(i) * wface
           enddo
c
c----------Jacobian (derivative wrt "wside" in the second to last column)
c
           call toarray(neq,nrhs,jac,1,nrhs-1,cside)
c	 
           do i=2,neq
             call toarray(neq,nrhs,jac,i,nrhs-1,sval(i))
           enddo
c
c----------Jacobian (derivative wrt "wface" in the third to last column)
c
           call toarray(neq,nrhs,jac,1,nrhs-2,cface)
c	 
           do i=2,neq
             call toarray(neq,nrhs,jac,i,nrhs-2,cval(i))
           enddo	   
	 elseif(itype.eq.5) then
c
c-------------------------------------------------------------------------------c
c----------compute contribution from 2 pnts on a side and face centers----------c      
c-------------------------------------------------------------------------------c 
c
c----------Face centers-----------c
c
           do k=1,20
             ct    = fcent(4,k)
             theta = acos(ct)
             sr    = sqrt(fcent(2,k)*fcent(2,k) + fcent(3,k)*fcent(3,k))
             si    = fcent(3,k)/sr
             co    = fcent(2,k)/sr	 
c
c------------compute values of spherical harmonics at face centers   
c
             call   spharm(nnmax,nnmax,nnmax,ct,theta,si,co,p,
     1                     dp,phire,phiim,toll)
c     
c------------sum over each equation      
c
             do i=1,neq
               cval(i) = cval(i) + p(ipm(i),ipn(i)) * phire(ipm(i))
             enddo
           enddo
c
c----------2 points on the side-------c
c
c----------Use quaternion multiplication to generate the 60 images of side point
c	 	       
           do k=1,60
             call quatmult(spnt(1,1),rgc(1,k),qtemp(1,k))
             call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
           enddo
c
c----------sum over the 60 images 
c
           do k=1,60
             ct    = qres(4,k)
	     theta = acos(ct)
             sr    = sqrt(qres(2,k)*qres(2,k)+qres(3,k)*qres(3,k))
             si    = qres(3,k)/sr
             co    = qres(2,k)/sr
c	 
c------------compute values of spherical harmonics at side centers 
c  
             call   spharm(nnmax,nnmax,nnmax,ct,theta,si,co,p,
     1                     dp,phire,phiim,toll)
c     
c------------sum over each equation 
   
             do i=1,neq	       
               sval(i) = sval(i) + p(ipm(i),ipn(i)) * phire(ipm(i))
             enddo
           enddo	   
c
c----------multiply by weights and add to running total
c
           do i=1,neq
             fval(i) = fval(i) + sval(i) * wside + cval(i) * wface
           enddo
c
c----------Jacobian (derivative wrt "wside" in the second to last column)
c
           call toarray(neq,nrhs,jac,1,nrhs-1,2.0d0*cside)
c	 
           do i=2,neq
             call toarray(neq,nrhs,jac,i,nrhs-1,sval(i))
           enddo
c
c----------Jacobian (derivative wrt "wface" in the third to last column)
c
           call toarray(neq,nrhs,jac,1,nrhs-2,cface)
c	 
           do i=2,neq
	   
             call toarray(neq,nrhs,jac,i,nrhs-2,cval(i))
	     
           enddo
	   	   	   
	 endif	 
	 
	 
c======================================================================c
c======================================================================c	  	 
c======================================================================c
	 

c---------------------------------------c
c--------loop over generators ----------c
c---------------------------------------c 
         do j=1,ngen		   
c
c----------Use quaternion multiplication to generate the 60 images of gen(j)
c	 	       
           do k=1,60
             call quatmult(gen(1,j),rgc(1,k),qtemp(1,k))	   
             call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
           enddo
c
c----------pre-compute coeffiecients needed for the jacobian
c
           call coefjac(gen(1,j),qres,ftheta1,fphi1,ftheta2,fphi2)
c
c----------initialize
c
           do i=1,neq
	   
             qval(i) = zero
             pdth(i) = zero
             pdph(i) = zero
	     
           enddo
c
c----------sum over the 60 images from pure rotation subgroup
c
           do k=1,60
	   
             ct    = qres(4,k)
	     theta = acos(ct)
	     	     
	     if((qres(2,k)*qres(2,k) + qres(3,k)*qres(3,k)) .le. 
     1	         1.0d-12) then
	       si = 0.0d0
	       co = 1.0d0
	     else
               sr = sqrt(qres(2,k)*qres(2,k) + qres(3,k)*qres(3,k))
               si = qres(3,k)/sr
               co = qres(2,k)/sr	     
	     endif	     
c
c------------compute value of spherical harmonics and derivatives at image point for all n,m
c               		 
             call spharm(nnmax,nnmax,nnmax,ct,theta,si,co,
     1	                 p,dp,phire,phiim,toll)        
c
c------------sum over equations
c
             do i=1,neq
c	     
               if (mod(ipn(i),2).eq.0) then
c
                 qval(i) = qval(i) + mlt*p(ipm(i),ipn(i))*phire(ipm(i))
c                    
                 pdth(i) = pdth(i) +
     1           mlt*(dp(ipm(i),ipn(i))*phire(ipm(i))*ftheta1(k) -
     2           ipm(i)*p(ipm(i),ipn(i))*phiim(ipm(i))*ftheta2(k))               
c                        
                 pdph(i) = pdph(i) +
     1           mlt*(dp(ipm(i),ipn(i))*phire(ipm(i))*fphi1(k) -
     2           ipm(i)*p(ipm(i),ipn(i))*phiim(ipm(i))*fphi2(k))
                
               else
c
                 qval(i) = qval(i) + p(ipm(i),ipn(i))*phiim(ipm(i))
c                 
                 pdth(i) = pdth(i) +
     1           dp(ipm(i),ipn(i))*phiim(ipm(i))*ftheta1(k) +
     2           ipm(i)*p(ipm(i),ipn(i))*phire(ipm(i))*ftheta2(k)
c                       
                 pdph(i) = pdph(i) +
     1           dp(ipm(i),ipn(i))*phiim(ipm(i))*fphi1(k) +
     2           ipm(i)*p(ipm(i),ipn(i))*phire(ipm(i))*fphi2(k)
c                     
               endif
c	       
             enddo
	     
           enddo	   
c
c----------Check to see if using the full icosahedral group
c----------If so, use the other 60 points generated by reflections
c     
	   if(ifg.eq.1) then
     


	   endif
c	   
c----------multiply by weight for generator
c
           do i=1,neq
	   
             fval(i) = fval(i) + qval(i)*wei(j)
	     
           enddo   
c
CCCCCCC----CHANGED VARIABLE ORDERING 2/3/2010
c
           if(ifg.eq.0) then

             call toarray(neq,nrhs,jac,1,2*ngen + j,cgene)
	     
	   else

             call toarray(neq,nrhs,jac,1,2*ngen + j,cgenefg)
	     
	   endif
c
c----------partial wrt weights of higher order equations
c
           do i=2,neq

	     call toarray(neq,nrhs,jac,i,2*ngen + j,qval(i))     
	     
           enddo
c	     
c----------partial wrt angles
c
           do i=1,neq
	   
             call toarray(neq,nrhs,jac,i,ngen + j,pdph(i)*wei(j))
             call toarray(neq,nrhs,jac,i,       j,pdth(i)*wei(j))
	     
           enddo	   
c
c----------end of loop over generators with distinct weights
c   
         enddo
c
c--------End loop over generators 
c 
	 
c
c--------For first equation to give zero, subtract a constant (integral of Y_0^0)
c
	 fval(1) = fval(1) - crhs

	 return
c	 
       end
