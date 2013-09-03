       subroutine LatticeToXY(iint,latind,pnts)

	 implicit real *8 (a-h,o-z)       
       	 integer latind(*)
	 real *8 eonex,eoney,etwox,etwoy
	 real *8 pnts(2*iint) 

c
c--------Lattice vectors
c
	 eonex = 1.0d0 
	 eoney = 0.0d0
	 
	 etwox = 0.5d0
	 etwoy = 0.866025403784439d0

	 
         do i=1,iint

	   pnts(2*i-1) = latind(2*i-1)*eonex + latind(2*i)*etwox 
	   pnts(2*i)   = latind(2*i-1)*eoney + latind(2*i)*etwoy
         
	 end do
c
	 return
c	 
       end
