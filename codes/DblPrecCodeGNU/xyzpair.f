c       $Log:   xypair.f,v $
c Revision 1.2  92/09/15  22:14:25  beylkin
c  (1) -> (*)
c
c Revision 1.1  92/08/11  10:46:14  beylkin
c Initial revision
c
c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c                                                                           c
c        This file is a part of                                             c
c        Double Precision Fast Wavelet Transform Library                    c
c        Contains proprietary information supplied by GB Consulting.        c
c        Copyright (C), 1992 GB Consulting. All rights reserved             c
c                                                                           c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c                                                                           c
c        This subr. outputs (x,y) pairs for ploting using Mathematica
c
c
c        input:
c
c        ax     --- one dimensional array of x-coordinates of length length 
c        ay     --- one dimensional array of y-coordinates of length length
c        az     --- one dimensional array of z-coordinates of length length  
c
c        output: 
c
c        ios       --- i/o status; 0 if o.k.
c
c
         subroutine xyzpair(fname,ax,ay,az,length,ios,io)
c    
           implicit real *8 (a-h,o-z)
	   integer eof
           real *8 ax(*),ay(*),az(*)
	   logical ex
           character*24 fname, stat
c
c check to see if file exists
c
           inquire(file = fname,exist = ex)

	   if(ex) then
	     stat = 'OLD'
	   else	   
	     stat = 'NEW'
	   endif
c
c open file
c
           open (unit = 3, file = fname, status = stat, iostat=ios,
     1           form = 'formatted')
           if (ios .ne. 0) then
             call prinf('----------------------*',0,0)
             call prina('unable to open file*',fname,8)
             call prinf('status=*',ios,1)
             return
           end if
c
c if old file, position at eof to write more data
c          
           if(ex) then
             eof = 0
             do while(eof.eq.0)
	       read(3,1000,iostat = eof)
	     end do 
	   endif
	   
	   backspace(3)	   
c
	   do 777 j=1,length
	     write(3,1000) ax(j),ay(j),az(j)
 777       continue
 
 1000      format(1x,3(2x,e40.32))
           close(3)
c
           if(io.eq.1) then
             call prinf('----------------------*',0,0)
             call prina('file *',fname,24)
	     if(ex) then
	       call prinf('has been appended*',0,0)
	     else
               call prinf('has been saved*',0,0)
	     endif
	   endif
	     
c
           return
	   
         end
