       subroutine makefilename(ngen,fname3D,fname2D,l)

         implicit real*8 (a-h,o-z)	 
	 integer  l, ngen
	 
	 character*20  fname1, fname2D, fname3D

	 character*4   chngen	
c
c--------Create standard filename for generators
c   
         write(chngen, '(I4)') ngen
 
  	 chngen   = trim(adjustr(chngen))	 	 
  
         fname1  = chngen//'gen'	 
	 fname3D = trim(adjustl(fname1))//'3D-1-fg.dat'
	 fname2D = trim(adjustl(fname1))//'2D-1-fg.dat'
	 
	 l = len_trim(fname2D)	 
	 	 
         return
       
       end
