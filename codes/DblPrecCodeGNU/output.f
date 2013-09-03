c--------------------------------------------------------------------------------------------c
c        C.A. Dec. 8, 2008
c        Last edit: 12/8/08
c
c        Purpose: Write the results to a data file 
c
c        Input: nfilename
c               wfilename
c		wvert
c		wface
c		wside
c		wei
c		ngen
c		itype      = 1 if vertex + general points; 2 if vertex + face center;
c                            3 if vertex + side midpoints; 4 if all three
c		ifg   
c
c        Output: data files
c--------------------------------------------------------------------------------------------c
c
         subroutine output(nfilename,wfilename,wvert,wface,wside,wei,
     1   	           ngen,itype,ifg,x,y,z,io) 
c
           implicit real *8 (a-h,o-z)
           real *8 vert(4,12),fcent(4,20),scent(4,30)
	   real *8 xv(12), yv(12), zv(12), xf(20), yf(20), zf(20)
	   real *8 xs(30), ys(30), zs(30), x(*),y(*),z(*), wei(*)
	 
 	   character*50 nfilename, wfilename	   
c
c----------retrieve vertices, face centers and side midpoints
c
           call vertex(vert)
	   call fcenter(fcent)
	   call scenter(scent) 	
	   do i=1,12
	     xv(i) = vert(2,i)
	     yv(i) = vert(3,i)
	     zv(i) = vert(4,i)
	   enddo   	   
	   do i=1,20
	     xf(i) = fcent(2,i)
	     yf(i) = fcent(3,i)
	     zf(i) = fcent(4,i)
	   enddo 
	   do i=1,30
	     xs(i) = scent(2,i)
	     ys(i) = scent(3,i)
	     zs(i) = scent(4,i)
	   enddo 	   	      
c
c----------if using the full group, set multiplier to 2
c
	   mm = ifg + 1
c
c----------Write nodes in (x,y,z) format to nfilename
c
           if(itype.eq.1) then
c
c------------first write the 12 icosahedral vertices
c
             call xyzpair(nfilename,xv,yv,zv,12,ios,io)	     
c
c------------next write the general nodes
c
             call xyzpair(nfilename,x,y,z,mm*60*ngen,ios,io)
	     	   
	   elseif(itype.eq.2) then
c
c------------first write the 12 icosahedral vertices
c
             call xyzpair(nfilename,xv,yv,zv,12,ios,io)
c
c------------second write the 20 icosahedral face centers
c
             call xyzpair(nfilename,xf,yf,zf,20,ios,io)	     	     
c
c------------last write the general nodes
c
             call xyzpair(nfilename,x,y,z,mm*60*ngen,ios,io)	     
c	   
	   elseif(itype.eq.3) then
c
c------------first write the 12 icosahedral vertices
c
	     call xyzpair(nfilename,xv,yv,zv,12,ios,io)
c
c------------second write the 30 icosahedral side midpoints
c
	     call xyzpair(nfilename,xs,ys,zs,30,ios,io)	     	     
c
c------------last write the general nodes
c
	     call xyzpair(nfilename,x,y,z,mm*60*ngen,ios,io)	     
c	   
	   else
c
c------------first write the 12 icosahedral vertices
c
	     call xyzpair(nfilename,xv,yv,zv,12,ios,io)
c
c------------second write the 20 icosahedral face centers
c
	     call xyzpair(nfilename,xf,yf,zf,20,ios,io)	     
c
c------------third write the 30 icosahedral side midpoints
c
	     call xyzpair(nfilename,xs,ys,zs,30,ios,io)	     	     
c
c------------last write the general nodes
c
	     call xyzpair(nfilename,x,y,z,mm*60*ngen,ios,io)	     
c	   
	   endif
c	   	   
c=====================================c
c=====================================c
c=====================================c
	   open(unit = 12, file = wfilename, status = 'unknown',
     1	        form = 'formatted')
c
c----------Write weights to wfilename
c
           if(itype.eq.1) then
c
c------------first write the 12 icosahedral vertices
c 
	     do i=1, 12
	       write(12,20) wvert 
	     end do 
c
c------------next write the general nodes
c
             do i=1, ngen
	       do j=1, mm*60
	         write(12,20)  wei(i)
	       enddo
	     enddo  
c	     	     	     	   
	   elseif(itype.eq.2) then
c
c------------first write the 12 icosahedral vertices
c
	     do i=1, 12
	       write(12,20) wvert 
	     end do
c
c------------second write the 20 icosahedral face centers
c
	     do i=1, 20
	       write(12,20) wface  
	     end do	     	     
c
c------------last write the general nodes
c
             do i=1, ngen
	       do j=1, mm*60
	         write(12,20)  wei(i)
	       enddo
	     enddo	     
c	   
	   elseif(itype.eq.3) then
c
c------------first write the 12 icosahedral vertices
c
	     do i=1, 12
	       write(12,20) wvert  
	     end do
c
c------------second write the 30 icosahedral side midpoints
c
	     do i=1, 30
	       write(12,20) wside 
	     end do	     	     
c
c------------last write the general nodes
c
             do i=1, ngen
	       do j=1, mm*60
	         write(12,20)  wei(i)
	       enddo
	     enddo	     
c	   
	   else
c
c------------first write the 12 icosahedral vertices
c
	     do i=1, 12
	       write(12,20) wvert   
	     end do
c
c------------second write the 20 icosahedral face centers
c
	     do i=1, 20
	       write(12,20) wface 
	     end do	     
c
c------------third write the 30 icosahedral side midpoints
c
	     do i=1, 30
	       write(12,20) wside  
	     end do	     	     
c
c------------last write the general nodes
c
             do i=1, ngen
	       do j=1, mm*60
	         write(12,20)  wei(i)
	       enddo
	     enddo	     
c	   
	   endif

 	   close(12)
	   
 10        format(E23.16,2x,E23.16,2x,E23.16)  
 20        format(E23.16) 
	   	   
           return
         end
