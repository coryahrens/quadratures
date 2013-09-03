c--------------------------------------------------------------------------------------------c
c        C.A. Dec. 8 2008 
c        Last edit: 12/8/08
c
c        Purpose: Convert from generators to cartesian coordinates of all nodes
c
c        Input: ngen -- number of generators
c               ifg  -- 0 = rotation subgroup; 1 = full icosahedral group
c                
c
c        Output: (x_i,y_i,z_i) -- coordinates of points i=1,2,...
c--------------------------------------------------------------------------------------------c
c
       subroutine toCart(ngen,gen,x,y,z,ifg) 
c 
         implicit real *8 (a-h,o-z)

	 real *8 rg(4,60),rgc(4,60),qtemp(4,60),qres(4,60),gen(4,ngen)
	 real *8 x(120*ngen), y(120*ngen), z(120*ngen)	 
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
c---------------------------------------c
c--------loop over generators ----------c
c---------------------------------------c
         do j=1,ngen   
c
           if(ifg.eq.0) then
c	   	 	    
c------------use only the rotation subgroup 
c  
             do k=1,60
c	   
               call quatmult(gen(1,j),rgc(1,k),qtemp(1,k))
               call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
c	     
	       x((j-1)*60 + k) = qres(2,k)
	       y((j-1)*60 + k) = qres(3,k)
	       z((j-1)*60 + k) = qres(4,k)
c
	     enddo
c	     
	   else
c	   	 	    
c------------use the full group
c
	     do k=1,60
c	     
               call quatmult(gen(1,j),rgc(1,k),qtemp(1,k))
               call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
c	     
	       x((j-1)*120 + k) = qres(2,k)
	       y((j-1)*120 + k) = qres(3,k)
	       z((j-1)*120 + k) = qres(4,k)	     
c
	       x((j-1)*120 + 60 + k) = -qres(2,k)
	       y((j-1)*120 + 60 + k) = -qres(3,k)
	       z((j-1)*120 + 60 + k) = -qres(4,k)	
	     
	     enddo
c	        
	   endif
c
         enddo	   
 
         return
	 
       end
