c	$Log:	hpsort1.f,v $
c Revision 1.1  92/08/04  11:12:20  beylkin
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
c
c
c
c	 modified  "hpsort" to handle tables with more than one entry
c        specifically, real "rentry".
c        keys to be sorted are in "item" 
c
c
         subroutine hpsort1(item,rentry,nn)
c
           implicit real *8 (a-h,o-z)
           integer item(1)
           real *8 rentry(1)
c
	   l  = nn/2+1
 200       continue  
  	   if (l.gt.1) then
	     l=l-1
             call sift1(item,rentry,l,nn)
             goto 200
	   endif
c
           ir=nn
 300       continue  
 
  	   if (ir.gt.1) then
	   
	     ix         = item(1)
             rxentry    = rentry(1)
             item(1)    = item(ir)
             rentry(1)  = rentry(ir)
             item(ir)   = ix
             rentry(ir) = rxentry
	     ir=ir-1
             call sift1(item,rentry,1,ir)
             goto 300
	   
	   endif
c
           return
	 
         end

c================================================================
c
c	 modified sift to handle tables with more than one entry
c        specifically, real "rentry". 
c        keys to be sorted are in "item" 
	 
         subroutine sift1(item,rentry,l,ir)
c
           implicit real *8 (a-h,o-z)
           integer item(1)
           real *8 rentry(1)
c
  	   i  = l
	   j  = 2*i
           ix   = item(i)
           rxentry = rentry(i)
	 
 200       continue  
 
  	   if (j.le.ir) then
	   
             if ((j.lt.ir).and.(item(j).lt.item(j+1))) j=j+1
	     if (ix.ge.item(j)) goto 100
	     item(i)   = item(j)
             rentry(i) = rentry(j)
	     i=j
             j=2*i
	     goto 200
	     
	   endif
 100       continue
 
	   item(i)   = ix
           rentry(i) = rxentry
c
           return
	   
         end	 
