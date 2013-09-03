       subroutine makefilename(nmax,ngen,iqtype,npoints,fname,l,io)

         implicit real*8 (a-h,o-z)	 
	 integer nmax, iqtype, npoints, l, ngen, io
	 
	 character*50  fname
	 character*20  fname1, fname2
	 character*10  chnmax
	 character*7   chnpoints, chngen
	 character*1   chiqtype	 
	 
	 if(io.eq.0) then
c
c----------Create standard filename for xyz output
c	 
	   write(chnmax,'(I4)') nmax     
           write(chnpoints, '(I7)') npoints
	   write(chiqtype, '(I1)') iqtype
	
  	   chnpoints = trim(adjustl(chnpoints))	 
	   chiqtype  = trim(adjustl(chiqtype))
	   chnmax    = trim(adjustl(chnmax))	 
  
           fname1 = 'qsph'//chiqtype//'-'//chnmax	 
	   fname2 = trim(adjustl(fname1))//'-'
	   fname2 = trim(adjustl(fname2))//chnpoints
	   fname  =  trim(adjustl(fname2))//'.dat'
	 
	   l = len_trim(fname)
	 
	 elseif(io.eq.1) then
c
c----------Create standard filename for generators
c
	   write(chnmax,'(I4)') nmax     
           write(chngen, '(I4)') ngen
	   write(chiqtype, '(I1)') iqtype
	   	 
  	   chngen   = trim(adjustl(chngen))	 
	   chiqtype = trim(adjustl(chiqtype))
	   chnmax   = trim(adjustl(chnmax))	 
  
           fname1 = 'gen'//chiqtype//'-'//chnmax	 
	   fname2 = trim(adjustl(fname1))//'-'
	   fname2 = trim(adjustl(fname2))//chngen
	   fname =  trim(adjustl(fname2))//'.dat'
	 
	   l = len_trim(fname)	 
	 
	 elseif(io.eq.2) then
c
c----------Create standard filename for weights
c
	   write(chnmax,'(I4)') nmax     
           write(chngen, '(I4)') ngen
	   write(chiqtype, '(I1)') iqtype
	   	 
  	   chngen   = trim(adjustl(chngen))	 
	   chiqtype = trim(adjustl(chiqtype))
	   chnmax   = trim(adjustl(chnmax))	 
  
           fname1 = 'weights'//chiqtype//'-'//chnmax	 
	   fname2 = trim(adjustl(fname1))//'-'
	   fname2 = trim(adjustl(fname2))//chngen
	   fname =  trim(adjustl(fname2))//'.dat'
	 
	   l = len_trim(fname)	
	 
	 else 

           write(chngen, '(I4)') ngen
	   fname = trim(adjustl(chngen))//'gen3D.dat'
	   l = len_trim(fname)
	 
	 endif
	 
         
         return
       
       end
