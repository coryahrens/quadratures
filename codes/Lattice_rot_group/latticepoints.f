       subroutine GetLatticePoints(t,s,latind,loind,rhind,lhind,
     1                             iqtype,i,irhs,ilhs,ilower)
 
	 implicit real *8 (a-h,o-z)       
       	 integer latind(*),loind(*),rhind(*),lhind(*)	 
	 integer n,m,s,t,t0n,t1n,t2n


         iqtype = 1

         irhs = 0
	 ilhs = 0
	 ilower = 0
         i = 0  
	 
         do m=-s+1,t-1
	 
	   do n=1, t+s-1
	   
	     t0n = s*s + m*(s-t) + s*t + t*t - n*(s + 2*t)
	     t1n = (2*m + n)*s + (m - n)*t
	     t2n = 3*(n*t - m*s)
c
c------------Check to see if point is strictly inside FD
c             
	     if( (t0n.gt.0).and.(t1n.gt.0).and.(t2n.gt.0) ) then 

	       i = i + 1
	       latind(2*i-1) = m
	       latind(2*i)   = n

             end if	
c
c------------Check to see if point is at the center of FD
c	          
             if ( (t0n.eq.0).and.(t1n.eq.0) ) then
	       
	       iqtype = 2
	     
	     end if
c
c------------Check to see if point is on lower boundary of larger triangle
c	     
	     if((t2n.eq.0).and.(t0n.gt.0).and.(t1n.gt.0)) then
	     	     
	       ilower = ilower + 1
	       loind(2*ilower - 1) = m
	       loind(2*ilower)     = n
	       	     
	     end if

c
c------------Check to see if point is on RHS boundary of larger triangle
c	     
	     if((t0n.eq.0).and.(t1n.gt.0).and.(t2n.gt.0)) then
	     
	       irhs = irhs + 1
	       rhind(2*irhs - 1) = m
	       rhind(2*irhs)     = n
	     	     
	     end if
	     
c
c------------Check to see if point is on LHS boundary of larger triangle
c	     
	     if((t1n.eq.0).and.(t0n.gt.0).and.(t2n.gt.0)) then

               ilhs = ilhs + 1
	       lhind(2*ilhs - 1) = m
	       lhind(2*ilhs)     = n
	     
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
