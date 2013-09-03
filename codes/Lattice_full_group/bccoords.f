       subroutine GetBaryCentricCoords(t,s,bccoord,lobccoord,rhbccoord,
     1                                 lhbccoord,iqtype,i,irhs,
     2                                 ilhs,ilower)
 
	 implicit real *8 (a-h,o-z)       
       	 real *8 bccoord(*),lobccoord(*),rhbccoord(*),lhbccoord(*)
	 real *8 t0, tc, tm	 
	 integer n,m,s,t


         iqtype = 1

         rdeter = s*s + t*s + t*t

         irhs = 0
	 ilhs = 0
	 ilower = 0
         i = 0  
	 
         do m=-s+1,t-1
	 
	   do n=1, t+s-1
c
c------------Lower left corner
c
             t0 = (rdeter - n*(2.0d0*s + t) - m*(s + 2.0d0*t)) / rdeter
c
c------------Midpoint between vertices
c     
             tm = 2.0d0*(n*(s - t) + m*(2.0d0*s + t)) / rdeter 
c
c------------Face-center
c
	     tc = 3.0d0*(n*t - m*s) / rdeter
     	     
c
c------------Check to see if point is strictly inside FD
c             
	     if( (t0.gt.0.0d0).and.(tc.gt.0.0d0).and.(tm.gt.0.0d0) ) 
     1	     then 

	       i = i + 1
	       bccoord(3*i-2) = t0
	       bccoord(3*i-1) = tc
	       bccoord(3*i)   = tm

             end if	
c
c------------Check to see if point is at the center of FD
c	          
             if ( (t0.eq.0.0d0).and.(tm.eq.0.0d0) ) then
	       
	       iqtype = 2
	     
	     end if
c
c------------Check to see if point is on lower boundary of larger triangle
c	     
	     if((tc.eq.0).and.(t0.gt.0).and.(tm.gt.0)) then
	     	     
	       ilower = ilower + 1
	       lobccoord(3*ilower - 2) = t0
	       lobccoord(3*ilower - 1) = tm
	       lobccoord(3*ilower)     = 0.0d0
	       	     
	     end if

c
c------------Check to see if point is on RHS boundary of larger triangle
c	     
	     if((t0.eq.0).and.(tc.gt.0).and.(tm.gt.0)) then
	     
	       irhs = irhs + 1
	       rhbccoord(3*irhs - 2) = 0.0d0
	       rhbccoord(3*irhs - 1) = tc
	       rhbccoord(3*irhs)     = tm
	     	     
	     end if
	     
c
c------------Check to see if point is on LHS boundary of larger triangle
c	     
	     if((tc.eq.0).and.(t0.gt.0).and.(tm.gt.0)) then

               ilhs = ilhs + 1
	       lhbccoord(3*ilhs - 2) = t0
	       lhbccoord(3*ilhs - 1) = 0.0d0
	       lhbccoord(3*ilhs)     = tm
	     
	     end if	     	     
	   	   
	   end do
	   
	 end do
	 
	 if( (iqtype.eq.2).and.(ilower.eq.1) ) then
	   
	   iqtype = 4 
	 
	 elseif( (iqtype.eq.1).and.(ilower.eq.1) ) then
	 
	   iqtype = 3
	   
	 elseif( (irhs.gt.0).or.(ilhs.gt.0).or.(ilower.gt.1) ) then
	 
	   iqtype = 5
	   
	 end if
c
	 return
c	 
       end
