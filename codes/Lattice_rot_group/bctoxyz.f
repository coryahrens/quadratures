       subroutine bctoxyz(bccoord,lobccoord,rhbccoord,
     1                             lhbccoord,iqtype,nn,irhs,
     2                             ilhs,ilower,p1,p2,p3,x,y,z,
     3                             ipointer,iflag,icount)
 
	 implicit real *8 (a-h,o-z)       
       	 real *8 bccoord(*),lobccoord(*),rhbccoord(*),lhbccoord(*)
	 real *8 p1(*),p2(*),p3(*),x(*),y(*),z(*)
         integer counter, ipointer(*)

         counter = 0	 


         if(iflag.eq.0) then
	   
	   write(*,*) 'nn=',nn

           do i=1,nn
	 
  	     x(i) = bccoord(3*i-2)*p1(1) + bccoord(3*i-1)*p2(1) +
     1	            bccoord(3*i)*p3(1)
	     y(i) = bccoord(3*i-2)*p1(2) + bccoord(3*i-1)*p2(2) +
     1	            bccoord(3*i)*p3(2) 
	     z(i) = bccoord(3*i-2)*p1(3) + bccoord(3*i-1)*p2(3) +
     1	            bccoord(3*i)*p3(3)
     
             write(*,*) 'x(i) = ',x(i)
     
               if(x(i).le.0.5d0) then
                 counter = counter + 1
                 ipointer(counter) = i
               end if
            
           end do

           icount = counter

         else

           do i=1,icount
	 
             jj   = ipointer(i)
            
  	     x(i) = bccoord(3*jj-2)*p1(1) + bccoord(3*jj-1)*p2(1) +
     1	            bccoord(3*jj)*p3(1)
	     y(i) = bccoord(3*jj-2)*p1(2) + bccoord(3*jj-1)*p2(2) +
     1	            bccoord(3*jj)*p3(2) 
	     z(i) = bccoord(3*jj-2)*p1(3) + bccoord(3*jj-1)*p2(3) +
     1	            bccoord(3*jj)*p3(3)
               
           end do


         end if

      
c           
         do i=1,ilower
	 
	   x(nn + i) = lobccoord(3*i-2)*p1(1) + lobccoord(3*i-1)*p2(1) +
     1	               lobccoord(3*i)*p3(1)
	   y(nn + i) = lobccoord(3*i-2)*p1(2) + lobccoord(3*i-1)*p2(2) +
     1	               lobccoord(3*i)*p3(2) 
	   z(nn + i) = lobccoord(3*i-2)*p1(3) + lobccoord(3*i-1)*p2(3) +
     1	               lobccoord(3*i)*p3(3)
         end do
c
c           
         do i=1,ilhs
	 
	   x(nn + ilower + i) = lhbccoord(3*i-2)*p1(1) + 
     1	               lhbccoord(3*i-1)*p2(1) + lhbccoord(3*i)*p3(1)
	   y(nn + ilower + i) = lhbccoord(3*i-2)*p1(2) + 
     1	               lhbccoord(3*i-1)*p2(2) + lhbccoord(3*i)*p3(2) 
	   z(nn + ilower+  i) = lhbccoord(3*i-2)*p1(3) + 
     1	               lhbccoord(3*i-1)*p2(3) + lhbccoord(3*i)*p3(3)
         end do

	 return
c	 
       end
