       subroutine getfilename(fname,message)
       
         implicit real *8 (a-h,o-z)
         integer l
	 
	 character*50 filename
	 character*50 fname
	 character*50 message
	 character*4  ext
	 ext = '.dat'
        
	 call prinf(message,0,0)
	 
	 read(*,*) filename	
	 
	 l = len_trim(filename)
	 
	 fname = filename(1:l)//ext 
       
         return
       
       end
