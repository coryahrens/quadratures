       subroutine openfile(fname,iunit,stat)
       
         implicit real *8 (a-h,o-z)

	 character*24  fname
	 character*3   stat
	 integer iunit,ios
	
	 
         open (unit=iunit,file=fname,status=stat,iostat=ios,
     1         form='formatted')
         if (ios .ne. 0) then
	 
           write(*,*) '----------------------'
           write(*,*) 'unable to open file'
           write(*,*) 'status',ios
 
         end if
c
         return

       end
