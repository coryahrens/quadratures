       subroutine setup(n,nmax,ngen,theta0,phi0,wei0,
     1                  wvert0,wface0,wside0,itype,ifg,gen)
         implicit real *8 (a-h,o-z)
         real *8 theta0(ngen), phi0(ngen)
	 real *8 wei0(ngen), wvert0, wface0, wside0
	 real *8 gen(4,ngen), sr, cphi
	 
	 zero   = 0.0d0
	 one    = 1.0d0
         two    = 2.0d0       
         pi     = 4.0d0*atan(1.0d0)
         fourpi = 4.0d0*pi

c----------------------------c	
         if(n.eq.206) then
           nmax  = 206
c          ngen  = 238
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0    
c
           open(unit=7,file='238gen3D.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	   
	 elseif(n.eq.155) then
c          c===========================================c
c          c===========================================c
           nmax  = 155
c          ngen  = 136
	   nodes = 60*ngen + 20 + 12
	   itype = 2
	   ifg   = 0 
c
           open(unit=7,file='136genSloan.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c


c          open(unit=7,file='N144theta.dat',status='OLD') 
c	   do j=1,ngen
c	     read(7,*)  theta0(j)  
c	   enddo
c	   close(7)
cc	   
c	   open(unit=7,file='N144phi.dat',status='OLD') 
c	   do j=1,ngen
c	     read(7,*)  phi0(j)  
c	   enddo
c	   close(7)
c	   
c	   do i=1,ngen
c	     gen(1,i) = zero
c            gen(2,i) = sin(theta0(i))*cos(phi0(i))
c            gen(3,i) = sin(theta0(i))*sin(phi0(i)) 
c            gen(4,i) = cos(theta0(i))
c	   enddo
cc 
c	   open(unit=7,file='N144weight.dat',status='OLD') 
c	   do j=1,ngen
c	     read(7,*)  wei0(j)  
c	   enddo
cc 	 
c          read(7,*) wvert0 
c          read(7,*) wface0 
c	   close(7)
c	   
c	   wside0 = zero
c
	 elseif(n.eq.154) then
c          c===========================================c
c          c===========================================c
           nmax  = 154
c          ngen  = 136
	   nodes = 60*ngen + 20 + 12
	   itype = 2
	   ifg   = 0 
c
           open(unit=7,file='N144theta.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  theta0(j)  
	   enddo
	   close(7)
	   
	   open(unit=7,file='N144phi.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  phi0(j)  
	   enddo
	   close(7)
	   
	   do i=1,ngen
	     gen(1,i) = zero
             gen(2,i) = sin(theta0(i))*cos(phi0(i))
             gen(3,i) = sin(theta0(i))*sin(phi0(i)) 
             gen(4,i) = cos(theta0(i))
	   enddo
c 
	   open(unit=7,file='N144weight.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  wei0(j)  
	   enddo
c 	 
           read(7,*) wvert0 
           read(7,*) wface0 
	   close(7)
	   
	   wside0 = zero
c

	 elseif(n.eq.151) then
c          c===========================================c
c          c===========================================c
           nmax  = 151
c          ngen  = 136
	   nodes = 60*ngen + 20 + 12
	   itype = 2
	   ifg   = 0 
c
           open(unit=7,file='N144theta.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  theta0(j)  
	   enddo
	   close(7)
	   
	   open(unit=7,file='N144phi.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  phi0(j)  
	   enddo
	   close(7)
	   
	   do i=1,ngen
	     gen(1,i) = zero
             gen(2,i) = sin(theta0(i))*cos(phi0(i))
             gen(3,i) = sin(theta0(i))*sin(phi0(i)) 
             gen(4,i) = cos(theta0(i))
	   enddo
c 
	   open(unit=7,file='N144weight.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  wei0(j)  
	   enddo
c 	 
           read(7,*) wvert0 
           read(7,*) wface0 
	   close(7)
	   
	   wside0 = zero
c
	 elseif(n.eq.147) then
           nmax  = 147
c           ngen  = 136
	   nodes = 60*ngen + 12 + 20
	   itype = 2
	   ifg   = 0 
c
           open(unit=7,file='N144theta.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  theta0(j)  
	   enddo
	   close(7)
	   
	   open(unit=7,file='N144phi.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  phi0(j)  
	   enddo
	   close(7)
	   
	   do i=1,ngen
	     gen(1,i) = zero
             gen(2,i) = sin(theta0(i))*cos(phi0(i))
             gen(3,i) = sin(theta0(i))*sin(phi0(i)) 
             gen(4,i) = cos(theta0(i))
	   enddo
c 
	   open(unit=7,file='N144weight.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  wei0(j)  
	   enddo
c 	 
           read(7,*) wvert0 
           read(7,*) wface0 
	   close(7)
	   
	   wside0 = zero
c
	 elseif(n.eq.146) then
           nmax  = 146
c           ngen  = 120
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	   
           open(unit=7,file='120genQuad145.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
           enddo
	   close(7)
c	   
	    do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           open(unit=7,file='weights.dat',status='OLD')
           do j=1,ngen
	     read(7,*) wei0(j) 
	   enddo
 	 
           read(7,*) wvert0 
	   
	   close(7)
	   
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.145) then
           nmax  = 145
c          ngen  = 120
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='120gen3DmapB.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c   
	 elseif(n.eq.148) then
           nmax  = 148
c           ngen  = 143
	   nodes = 60*ngen + 12
	   mm    = 0
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='N148genxyz.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c
c          open(unit=7,file='InAnglesN148.dat',status='OLD') 
c	   do j=1,ngen
c	     
c	     read(7,*) theta0(j), phi0(j) 
c	     gen(1,j) = 0.0d0 
c	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
c	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
c	     gen(4,j) = cos(theta0(j)) 
c	   enddo
c	   close(7)
c 
c          call prinf('Reading initial weights...*',0,0)
c          open(unit=7,file='initweight.dat',status='OLD')
c	   read(7,*) wvert0
c          ead(7,*) wei0(j)
c	   close(7)

           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	 
         elseif(n.eq.144) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax  = 144
c           ngen  = 136
	   nodes = 60*ngen + 20 + 12
	   mm    = 0
	   itype = 2
	   ifg   = 0
c
c           open(unit=7,file='InAnglesN144.dat',status='OLD') 
c	   do j=1,ngen
c	     read(7,*) theta0(j), phi0(j) 
c	   enddo
c	   close(7)
	   
	   open(unit=7,file='N144theta.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*) theta0(j) 
	   enddo
	   close(7)
	   
	   
	   open(unit=7,file='N144phi.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*) phi0(j) 
	   enddo
	   close(7)
	   
	   open(unit=7,file='N144weight.dat',status='OLD')    
           do j=1,ngen
	     read(7,*) wei0(j) 
           enddo
	   
	   read(7,*) wvert0
	   read(7,*) wface0  
	   
	   close(7)	
	   
	   wside0 = zero
c	  
         elseif(n.eq.143) then
c          c======================c
c          c= does not converges =c
c          c======================c
           nmax  = 143
c          ngen  = 115
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0
           open(unit=7,file='115gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	   	
	 elseif(n.eq.130) then
           nmax  = 130
c           ngen  = 113
	   nodes = 60*ngen + 12
	   mm    = 0
	   itype = 1
	   ifg   = 0
c	
           open(unit=7,file='sphN130b.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*) theta0(j), phi0(j)  
	   enddo
	   close(7)
c	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wside0 = zero
           wface0 = zero
           wvert0    = fourpi/(1.0d0*nodes)
c
         elseif(n.eq.127) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax  = 127
c           ngen  = 108
	   nodes = 60*ngen + 30 + 12
 	   mm    = 2
  	   itype = 3
	   ifg   = 0
c	
           open(unit=7,file='N127.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*) theta0(j), phi0(j)  
	   enddo
	   close(7)
c	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wside0  = fourpi/(1.0d0*nodes)
	   wvert0  = fourpi/(1.0d0*nodes)
	   wface0  = zero
c
	 elseif(n.eq.124) then
c          c===========================================c
c          c= 
c          c===========================================c
           nmax  = 124
c          ngen  = 103
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
c          t=5, s=22
c
           open(unit=7,file='103gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.103) then
c          c===========================================c
c          c= does not converge
c          c===========================================c
           nmax  = 103
c           ngen  = 60
	   nodes = 60*ngen + 12
	   mm    = 0
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='60gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.108) then
           nmax  = 108
c           ngen  = 80
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
c----------Sloan's initial data-----------------------
c           open(unit=7,file='N108.dat',status='OLD') 
c	   do j=1,ngen
c	     
c	     read(7,*) theta0(j), phi0(j) 
c	     
c	     gen(1,j) = 0.0d0 
c	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
c	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
c	     gen(4,j) = cos(theta0(j)) 
c	   enddo
c	   close(7)
c-----------------------------------------------------
	
           open(unit=7,file='xyz108.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
   
	   do j=1,ngen
            theta0(j) = acos(gen(4,j))
            sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
            cphi      = gen(2,j)/sr
            phi0(j)   = acos(cphi)
           enddo

c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.100) then
           nmax  = 100
c           ngen  = 57
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
	   open(unit=7,file='57gen3Dmap.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j)   
	   enddo	   
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = log(fourpi/(1.0d0*nodes))
	   enddo
c 	 
           wvert0 = log(fourpi/(1.0d0*nodes))

	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.99) then
           nmax  = 99
c           ngen  = 57
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
	   open(unit=7,file='57gen3Dmap.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j)  
	   enddo	   
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = log(fourpi/(1.0d0*nodes))
	   enddo
c 	 
           wvert0 = log(fourpi/(1.0d0*nodes))
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.98) then
           nmax  = 98
c           ngen  = 57
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
	   open(unit=7,file='57genQuad97_new.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
	        enddo
	   
	   read(7,*) wvert0	   
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
c           do j=1,ngen
c	     wei0(j) = fourpi/(1.0d0*nodes)
c	   enddo
c 	 
c           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.97) then
c          c===========================================c
c          c===========================================c
           nmax  = 97
c          ngen  = 57
	   nodes = 60*ngen  + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='57gen3D.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
             read(7,*)  gen(2,j),gen(3,j),gen(4,j)  
	   enddo	   	   
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo	   
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.95) then
c          c===========================================c
c          c===========================================c
           nmax  = 95
c          ngen  = 51
	   nodes = 60*ngen  + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='51gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
             read(7,*)  gen(2,j),gen(3,j),gen(4,j)  
	   enddo	   	   
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo	   
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.92) then
c	 
           nmax  = 92
c          ngen  = 57
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='57gen91QuadSol.dat',status='OLD')
c           open(unit=7,file='57gen3Dmap.dat',status='OLD')
c           open(unit=7,file='57genQuadSol.dat',status='OLD') 	
c           open(unit=7,file='57genQuad96_new2.dat',status='OLD')
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
c	   read(7,*) wvert0
	   close(7)
   
	   do j=1,ngen
            theta0(j) = acos(gen(4,j))
            sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
            cphi      = gen(2,j)/sr
            phi0(j)   = acos(cphi)
           enddo
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	   
	 elseif(n.eq.91) then
           nmax  = 91
c          ngen  = 57
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c 	
           open(unit=7,file='57gen3Dmap.dat',status='OLD')
c           open(unit=7,file='57genQuadSol.dat',status='OLD') 	
c           open(unit=7,file='57genQuad96_new2.dat',status='OLD')
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
c	   read(7,*) wvert0
	   close(7)
   
	   do j=1,ngen
            theta0(j) = acos(gen(4,j))
            sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
            cphi      = gen(2,j)/sr
            phi0(j)   = acos(cphi)
           enddo
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	   
	 elseif(n.eq.90) then
c          c===========================================c
c          c= Sloan's initial data converges          =c
c          c= =c
c          c===========================================c
           nmax  = 90
c          ngen  = 57
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c----------Sloan's initial data------------------------
c           open(unit=7,file='N90_b.dat',status='OLD') 
c	   do j=1,ngen
c	     
c	     read(7,*) theta0(j), phi0(j) 
c	     
c	     gen(1,j) = 0.0d0 
c	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
c	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
c	     gen(4,j) = cos(theta0(j))
c	      
c          enddo
c	  close(7)
c------------------------------------------------------	   	
           open(unit=7,file='57gen3Dmap.dat',status='OLD') 
c           open(unit=7,file='N90_c.dat',status='OLD')
	   do j=1,ngen
             gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
c   
	   do j=1,ngen
            theta0(j) = acos(gen(4,j))
            sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
            cphi      = gen(2,j)/sr
            phi0(j)   = acos(cphi)
           enddo
	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.85) then
           nmax  = 85
c          ngen  = 51
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
c----------Sloan's initial data-----------------------
           open(unit=7,file='N85_b.dat',status='OLD') 
	   do j=1,ngen
	     
	     read(7,*) theta0(j), phi0(j) 
	     
	     gen(1,j) = 0.0d0 
	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
	     gen(4,j) = cos(theta0(j))
	      
	   enddo
	   close(7)
c-----------------------------------------------------	   
	   do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
	   
c	   do j=1,ngen
c             theta0(j) = acos(gen(4,j))
c             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
c             cphi      = gen(2,j)/sr
c             phi0(j)   = acos(cphi)
c           enddo 
c 
c 	
c          open(unit=7,file='N90.dat',status='OLD') 
c	   do j=1,ngen
c	     read(7,*)  theta0(j), phi0(j)  
c	   enddo
c	   close(7)

           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
c	 elseif(n.eq.84) then
c           nmax  = 84
cc           ngen  = 49
c	   nodes = 60*ngen + 20 + 30 + 12
c	   mm    = 0
c	   itype = 4
c	   ifg   = 0 
cc	
c           open(unit=7,file='N84.dat',status='OLD') 
c	   do j=1,ngen
c	     read(7,*)  theta0(j), phi0(j)  
c	   enddo
c	   close(7)
cc	  	  
c           do j=1,ngen
c	     wei0(j) = fourpi/(1.0d0*nodes)
c	   enddo
cc 
c           wvert0 = fourpi/(1.0d0*nodes)
c	   wside0 = fourpi/(1.0d0*nodes)
c           wface0 = fourpi/(1.0d0*nodes)
c
c	 elseif(n.eq.79) then
c	   nmax  = 79
c          ngen  = 43
c	   nodes = 60*ngen + 12
c	   itype = 1
c	   ifg   = 0
c 
c           open(unit=7,file='theta.dat',status='OLD') 
c	   do j=1,ngen	     
c	     read(7,*) theta0(j) 
c	   enddo
c	   close(7)
c
c	   open(unit=7,file='phi.dat',status='OLD') 
c	   do j=1,ngen	     
c	     read(7,*) phi0(j) 
c	   enddo
c	   close(7)
c	   
c	   do j=1,ngen	     
c	     gen(1,j) = 0.0d0 
c	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
c	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
c	     gen(4,j) = cos(theta0(j))  
c	   enddo
c	   
c           do j=1,ngen
c	     wei0(j) = fourpi/(1.0d0*nodes)
c	   enddo
c	   
c           wvert0 = fourpi/(1.0d0*nodes)
c	   wside0 = zero
c           wface0 = zero
c
	 elseif(n.eq.87) then
	   nmax  = 87
c          ngen  = 43
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0
c 
           open(unit=7,file='theta.dat',status='OLD') 
	   do j=1,ngen	     
	     read(7,*) theta0(j) 
	   enddo
	   close(7)
c
	   open(unit=7,file='phi.dat',status='OLD') 
	   do j=1,ngen	     
	     read(7,*) phi0(j) 
	   enddo
	   close(7)
	   
	   do j=1,ngen	     
	     gen(1,j) = 0.0d0 
	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
	     gen(4,j) = cos(theta0(j))  
	   enddo
	   
	   open(unit=7, file='weights.dat',status='OLD')
           do j=1,ngen
	     read(7,*) wei0(j)
	     wei0(j) = log(wei0(j))
	   enddo
	   
	   close(7)
	   
           wvert0 = log(0.3625061829535985d-02)
	   wside0 = zero
           wface0 = zero
c
c	 elseif(n.eq.77) then
c          c=============c
c          c= Converges =c
c          c=============c	 
c	   nmax  = 77
c          ngen  = 43
c	   nodes = 60*ngen + 12
c	   itype = 1
c	   ifg   = 0
c 
c           open(unit=7,file='xyz77.dat',status='OLD') 
c	   do j=1,ngen
c	     gen(1,j) = zero
c	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
c	   enddo
c	   close(7)
c
c           open(unit=7,file='tTest.dat',status='OLD') 
c	   do j=1,ngen	     
c	     read(7,*) theta0(j) 
c	   enddo
c	   close(7)
c
c	   open(unit=7,file='pTest.dat',status='OLD') 
c	   do j=1,ngen	     
c	     read(7,*) phi0(j) 
c	   enddo
c	   close(7)
c	   do j=1,ngen
c	     
c	     gen(1,j) = 0.0d0 
c	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
c	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
c	     gen(4,j) = cos(theta0(j))  
c	   enddo
c	   
c
c           do j=1,ngen
c	     wei0(j) = fourpi/(1.0d0*nodes)
c	   enddo
c	   
c---------Sloan's initial distribution---------------
c           open(unit=7,file='N77.dat',status='OLD') 
c	   do j=1,ngen
c	     
c	     read(7,*) theta0(j), phi0(j) 
c	     
c	     gen(1,j) = 0.0d0 
c	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
c	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
c	     gen(4,j) = cos(theta0(j)) 
c	   enddo
c
c	   close(7)
c
c-----------------------------------------------------	   
c	   do j=1,ngen
c             theta0(j) = acos(gen(4,j))
c             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
c             cphi      = gen(2,j)/sr
c             phi0(j)   = acos(cphi)
c           enddo
c
c          do j=1,ngen
c	     wei0(j) = fourpi/(1.0d0*nodes)
c	   enddo
c
c 	   wside0 = zero
c           wface0 = zero
c	   wvert0 = fourpi/(1.0d0*nodes)
c
	 elseif(n.eq.83) then
           nmax  = 83
c          ngen  = 40
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c 	
           open(unit=7,file='40gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.69) then
c          c===========================================c
c          c= does not converge
c          c===========================================c
           nmax  = 69
c           ngen  = 27
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='27gen3Dmap.dat',status='OLD') 
c	   open(unit=7,file='27genSloane.dat',status='OLD')
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
c           do j=1,ngen
c	    wei0(j) = fourpi/(1.0d0*nodes)
c	   enddo
c 	 
c           wvert0 = fourpi/(1.0d0*nodes)
c	   wside0 = zero
c           wface0 = zero
c
	 elseif(n.eq.66) then
           nmax  = 66
c          ngen  = 32
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='xyz66_1.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
         elseif(n.eq.65) then
           nmax  = 65
           ngen  = 24
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='24genQuad64.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
	   enddo
	   
	   read(7,*) wvert0
	   
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c          
	   wside0 = zero
           wface0 = zero
c
         elseif(n.eq.64) then
           nmax  = 64
           ngen  = 24
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='24genQuad58.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
	   enddo
	   
	   read(7,*) wvert0
	   
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c          
	   wside0 = zero
           wface0 = zero
c
         elseif(n.eq.63) then
           nmax  = 63
           ngen  = 24
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='24genQuad58.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
	   enddo
	   
	   read(7,*) wvert0
	   
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c          
	   wside0 = zero
           wface0 = zero
c

	 elseif(n.eq.61) then
           nmax  = 61
c          ngen  = 21
	   nodes = 60*ngen + 12 + 20
	   itype = 2
	   ifg   = 0 
c
           open(unit=7,file='21gen3D_2.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = fourpi/(1.0d0*nodes)
c
	 elseif(n.eq.60) then
           nmax  = 60
c           ngen  = 25
	   nodes = 60*ngen + 20 + 30 + 12
	   mm    = 0
	   itype = 4
	   ifg   = 0 
c	
           open(unit=7,file='N60.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  theta0(j), phi0(j)  
	   enddo
	   close(7)
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = fourpi/(1.0d0*nodes)
           wface0 = fourpi/(1.0d0*nodes)
c
	 elseif(n.eq.59) then
           nmax  = 59
c          ngen  = 25
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='xyz59.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.58) then
           nmax  = 58
c          ngen  = 24
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
c           open(unit=7,file='xyzgenN58.dat',status='OLD') 
c	   do j=1,ngen
c	     gen(1,j) = zero
c	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
c	   enddo
c	   close(7)
	    
  	theta0(1) = 0.8438912552709054E+00
  	theta0(2) = 0.1500887469614057E+01
  	theta0(3) = 0.1150758191068289E+01
  	theta0(4) = 0.8370351232828056E+00
  	theta0(5) = 0.6632391633482384E+00
  	theta0(6) = 0.1012243470404262E+01
  	theta0(7) = 0.1407710889903340E+01
  	theta0(8) = 0.3483185855343895E+00
  	theta0(9) = 0.1252304209049151E+01
  	theta0(10) = 0.1363830118576745E+01
  	theta0(11) = 0.9662352085074769E+00
  	theta0(12) = 0.8886661283812348E+00
  	theta0(13) = 0.1329055807324737E+01
  	theta0(14) = 0.8223893089716992E+00
  	theta0(15) = 0.1346082865428686E+01
  	theta0(16) = 0.7092248177096380E+00
  	theta0(17) = 0.1048421615254036E+01
  	theta0(18) = 0.4426065926297103E+00
  	theta0(19) = 0.1511261601974008E+01
  	theta0(20) = 0.1526425195802320E+01
  	theta0(21) = 0.8580989162155506E+00
  	theta0(22) = 0.9796456665060419E+00
  	theta0(23) = 0.9765683856605565E+00
  	theta0(24) = 0.1608563917597192E+01
 
  	phi0(1) = 0.3598058288519538E+00
  	phi0(2) = 0.8798651455542741E+00
  	phi0(3) = 0.1130999001035378E+01
  	phi0(4) = 0.6726578496472480E+00
 	phi0(5) = 0.8670445621202446E+00
 	phi0(6) = 0.6495418288660540E+00
  	phi0(7) = 0.1317915547625973E+01
  	phi0(8) = 0.4738378656591005E+00
  	phi0(9) = 0.1126663270124580E+01
  	phi0(10) = 0.1415186780488475E+01
  	phi0(11) = 0.1054869975402722E+01
  	phi0(12) = 0.1202301773149220E+01
  	phi0(13) = 0.1003368186899938E+01
  	phi0(14) = 0.5116726620986390E+00
  	phi0(15) = 0.1113522026136739E+01
  	phi0(16) = 0.1186493038177151E+01
  	phi0(17) = 0.1384701037070480E+01
  	phi0(18) = 0.1270053184801110E+01
  	phi0(19) = 0.3566285372169135E+00
  	phi0(20) = 0.1537046291563691E+01
  	phi0(21) = 0.1403678097436229E+01
  	phi0(22) = 0.9454097544454161E-01
  	phi0(23) = 0.1162070434268147E+01
  	phi0(24) = 0.6302457102821801E+00  
	
	  do i=1,ngen
	    gen(1,i) = zero
            gen(2,i) = sin(theta0(i))*cos(phi0(i))
            gen(3,i) = sin(theta0(i))*sin(phi0(i)) 
            gen(4,i) = cos(theta0(i))
	  enddo 
	   
c	   do j=1,ngen
c             theta0(j) = acos(gen(4,j))
c             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
c             cphi      = gen(2,j)/sr
c             phi0(j)   = acos(cphi)
c           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.56) then
           nmax  = 56
c          ngen  = 23
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='xyz56_1.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.55) then
c          c===========================================c
c          c= Converges                               =c
c          c===========================================c
           nmax  = 55
c          ngen  = 22
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='22gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.53) then
c          c===========================================c
c          c= Converges                               =c
c          c===========================================c
           nmax  = 53
c          ngen  = 16
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='16genSloane.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
           wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.52) then
           nmax  = 52
c          ngen  = 16
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c
           open(unit=7,file='16gen3D.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.50) then
           nmax  = 50
c           ngen  = 18
	   nodes = 60*(ngen + 1) + 20 + 12
	   mm    = 0
	   itype = 5
	   ifg   = 0 
c
           open(unit=7,file='xyz50.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 	 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = fourpi/(1.0d0*nodes)
           wface0 = fourpi/(1.0d0*nodes)
c
	 elseif(n.eq.49) then
           nmax  = 49
c           ngen  = 10
	   nodes = 120*ngen + 12
	   mm    = 0
	   itype = 1
	   ifg   = 1 
c	
           open(unit=7,file='N49.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  theta0(j), phi0(j)  
	   enddo
	   close(7)
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c 
	 elseif(n.eq.47) then
           nmax  = 47
c          ngen  = 12
	   nodes = 60*ngen + 20 + 30 + 12
	   itype = 4
	   ifg   = 0 
c	
           open(unit=7,file='12gen3D_4.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = fourpi/(1.0d0*nodes)
           wface0 = fourpi/(1.0d0*nodes)
c 
	 elseif(n.eq.46) then
           nmax  = 46
c           ngen  = 16
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 
           open(unit=7,file='16gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c 
c	 elseif(n.eq.44) then
c           nmax  = 44
cc          ngen  = 11
c	   nodes = 60*ngen + 12
c	   itype = 1
c	   ifg   = 0 
c	                   
c           open(unit=7,file='11QuadSol_old.dat',status='OLD') 
c	   
c	   do j=1,ngen
c	     gen(1,j) = zero
c	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
c	   enddo
c
c	   do j=1,ngen
c	   
c	     read(7,*)  theta0(j), phi0(j)
c	     
c	    gen(1,j) = zero
c            gen(2,j) = sin(theta0(j))*cos(phi0(j))
c            gen(3,j) = sin(theta0(j))*sin(phi0(j)) 
c            gen(4,j) = cos(theta0(j))	     
c	       
c	   enddo
c
c	   close(7)
c	   
c	   do j=1,ngen
c             theta0(j) = acos(gen(4,j))
c             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
c             cphi      = gen(2,j)/sr
c             phi0(j)   = acos(cphi)
c          enddo 
c          do j=1,ngen
c	     wei0(j) = fourpi/(1.0d0*nodes)
c	   enddo
c 
c           wei0(1)  = 0.1753971860566286E-01
c           wei0(2)  = 0.1807289030819800E-01
c           wei0(3)  = 0.1609086142081690E-01
c           wei0(4)  = 0.1895598860991113E-01
c           wei0(5)  = 0.1477572208213813E-01
c           wei0(6)  = 0.2286003658763882E-01
c           wei0(7)  = 0.2170264253743665E-01
c           wei0(8)  = 0.1236451111086792E-01
c           wei0(9)  = 0.1689999656201042E-01
c           wei0(10) = 0.2039056406403785E-01
c           wei0(11) = 0.2519742066734817E-01 
c           wvert0   = 0.2294578841626432E-01	   
c	   
c	   wside0 = zero
c           wface0 = zero
c 
         elseif(n.eq.43) then
           nmax  = 43
           ngen  = 11
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='11gen3Dmap.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j)  
	   enddo
	   
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c    
           do j=1,ngen
	     wei0(j) = log(fourpi/(1.0d0*nodes))
	   enddo
	   
           wvert0 = log(fourpi/(1.0d0*nodes))        
	   wside0 = zero
           wface0 = zero
c
c         elseif(n.eq.43) then
c           nmax  = 43
c           ngen  = 11
c	   nodes = 60*ngen + 12
c	   itype = 1
c	   ifg   = 0 
c	 	   
c	   open(unit=7,file='11QuadSol_old.dat',status='OLD')
c   
c	   do j=1,ngen
c	     gen(1,j) = zero
c	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
c	   enddo
c	   close(7)
c	   
c	   do j=1,ngen
c             theta0(j) = acos(gen(4,j))
c             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
c             cphi      = gen(2,j)/sr
c             phi0(j)   = acos(cphi)
c           enddo 
c    
c           wvert0 = 0.2294578841626432D-01
c	   wside0 = zero
c           wface0 = zero
c
         elseif(n.eq.42) then
           nmax  = 42
           ngen  = 11
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='11QuadSol_old.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c    
           wvert0 = 0.2294578841626432D-01
	   wside0 = zero
           wface0 = zero
c
         elseif(n.eq.41) then
           nmax  = 41
           ngen  = 11
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='11QuadSol_old.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c    
           wvert0 = 0.2294578841626432D-01
	   wside0 = zero
           wface0 = zero
c
         elseif(n.eq.40) then
           nmax  = 40
           ngen  = 11
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='11gen3Dmap.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j)  
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c    
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.39) then
           nmax  = 39
           ngen  = 11
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 	   
	   open(unit=7,file='11QuadSol_old.dat',status='OLD')
c   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j),gen(3,j),gen(4,j),wei0(j)  
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c    
           wvert0 = 0.2294578841626432D-01
	   wside0 = zero
           wface0 = zero
c 
	 elseif(n.eq.38) then
           nmax  = 38
           ngen  = 11
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0 
c	 
           open(unit=7,file='11genSloane.dat',status='OLD') 
	   
c	   open(unit=7,file='11QuadSol_old.dat',status='OLD')
   
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j) 
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
           wvert0 = fourpi/(1.0d0*nodes)
	   
c           wvert0 = 0.2294578841626432D-01
	   wside0 = zero
           wface0 = zero
c  
	 elseif(n.eq.37) then
	   nmax  = 37
c	   ngen  = 8
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0
c 
c           open(unit=7,file='xyz34.dat',status='OLD') 
c           open(unit=7,file='iSldata.dat',status='OLD')
	   open(unit=7,file='8gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
	 elseif(n.eq.36) then
c          c=============c
c          c= Converges =c
c          c=============c
           nmax  = 36
c           ngen  = 9
	   nodes = 60*ngen + 20 + 12
	   mm    = 0
	   itype = 2
	   ifg   = 0 
c	 
           open(unit=7,file='xyz36.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = fourpi/(1.0d0*nodes)
c
	 elseif(n.eq.35) then
           nmax  = 35
c           ngen  = 7
	   nodes = 60*ngen +  12
	   itype = 1
	   ifg   = 0 
c	 
           open(unit=7,file='7gen3D.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	   	  
	 elseif(n.eq.34) then
	   nmax  = 34
c	   ngen  = 8
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0
c 
c           open(unit=7,file='xyz34.dat',status='OLD') 
c           open(unit=7,file='iSldata.dat',status='OLD')
	   open(unit=7,file='8gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
c	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo 


c           theta0(1) =  0.3880292627287484E+00
c           theta0(2) =  0.5542690200310353E+00
c           theta0(3) =  0.2411298482781721E+00
c           theta0(4) =  0.9429646217008456E-01
c           theta0(5) =  0.2974067127481551E+00
c           theta0(6) =  0.5505926093627633E+00
c           theta0(7) =  0.7838832271728686E+00
c           theta0(8) =  0.4244876821614893E+00
c 
c           phi0(1) =  -.2684367584822231E+00
c           phi0(2) =   -.1613853047050257E+00
c           phi0(3) =   0.8125771443445494E+00
c           phi0(4) =   0.1737659051412132E+00
c           phi0(5) =   -.1131425707649321E+01
c           phi0(6) =   0.7630191473666869E+00
c           phi0(7) =   0.1546249602895481E+00
c           phi0(8) =   -.1414508379499374E+01
c
c---------Sloan's initial distribution---------------
c           open(unit=7,file='N34.dat',status='OLD') 
c	   do j=1,ngen
c	     
c	     read(7,*) theta0(j), phi0(j) 
cc	     
c	     gen(1,j) = 0.0d0 
c	     gen(2,j) = sin(theta0(j))*cos(phi0(j))
c	     gen(3,j) = sin(theta0(j))*sin(phi0(j))
c	     gen(4,j) = cos(theta0(j)) 
c	   enddo
c	   close(7)
c----------------------------------------------------- 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	  
	 elseif(n.eq.32) then
c          c=============c
c          c= Converges =c
c          c=============c	 
	   nmax  = 32
c          ngen  = 7
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0
c 
           open(unit=7,file='xyz32.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c
 	   wside0 = zero
           wface0 = zero
	   wvert0 = fourpi/(1.0d0*nodes)
c	         	
	 elseif(n.eq.31) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax  = 31
c           ngen  = 6
	   nodes = 60*ngen + 12
	   itype = 1
	   ifg   = 0

c
           open(unit=7,file='6gen3Dmap.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero	   
	   
c
c           theta0(1) = 0.177683577730264d+01 
c           theta0(2) = 0.1300307975063925d+01
c           theta0(3) = 0.1059272963452778d+01
c           theta0(4) = 0.8026312137796407d+00
c           theta0(5) = 0.2941596161927831d+01
c           theta0(6) = 0.1983712243622630d+01
cc
c           phi0(1) = 0.3460277225335893d+01
c           phi0(2) = 0.6097932457843251d+01
c           phi0(3) = 0.1777437150516778d+01
c           phi0(4) = 0.6136775160787267d+01
c           phi0(5) = 0.8867243967924299d+00 
c           phi0(6) = 0.2787690864874015d+01
cc
c           wei0(1) = 0.3991277824630497d-01
c           wei0(2) = 0.3563392770902429d-01
c           wei0(3) = 0.3485019633098942d-01
c           wei0(4) = 0.3379777203184547d-01
c           wei0(5) = 0.2760879566031245d-01 
c           wei0(6) = 0.3223176368718324d-01
cc
c           wvert0 = 0.2702138286829945d-01
c  	    wside0 = zero
c           wface0 = zero
c	   	 
 	 elseif(n.eq.30) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax  = 30
c           ngen  = 5
	   nodes = 60*ngen + 30 + 20 + 12
	   itype = 4
	   ifg   = 0
c	  
	   theta0(1) = 0.41924280511374873585279406345494d0
           theta0(2) = 0.54627676914414576418832324222353d0
           theta0(3) = 0.34468410722531180393407851295636d0
	   theta0(4) = 0.54879596719275900100840858561553d0
	   theta0(5) = 0.52396973423025603548037080873643d0
c	 
           phi0(1)   = 0.51341322858052919743144263097697d0
           phi0(2)   = 0.20268651662470488394928802517922d0
           phi0(3)   = 0.99292343007973107985722221923123d0
	   phi0(4)   = -0.87038953969988002954664024592077d0
	   phi0(5)   = -1.2473152265312774525516437189896d0
c	   	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	    
           wei0(1) = fourpi/(1.0d0*nodes)
           wei0(2) = fourpi/(1.0d0*nodes)
	   wei0(3) = fourpi/(1.0d0*nodes)
           wei0(4) = fourpi/(1.0d0*nodes)
	   wei0(5) = fourpi/(1.0d0*nodes)
           wside0  = fourpi/(1.0d0*nodes)
	   wface0  = fourpi/(1.0d0*nodes) 
	   wvert0  = fourpi/(1.0d0*nodes)
c	 
	 elseif(n.eq.29) then	 
	   nmax  = 29
c	   ngen  = 5
	   nodes = 60*ngen + 12
	   itype = 1
 	   ifg   = 0
	   
           open(unit=7,file='5gen3D.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	
         elseif(n.eq.271) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax  = 27
c           ngen  = 4
	   nodes = 60*ngen + 30 + 12
	   mm    = 0
	   itype = 3
	   ifg   = 0
c	  
	   theta0(1) = 0.476074175265880923002205391549631305624787541d0
           theta0(2) = 0.422733902178372561408733496790331374713357794d0
           theta0(3) = 0.226421122343199283177695078504384246534740006d0
	   theta0(4) = 0.455313411663907116830681694491432493104495849d0
c	 
           phi0(1) = 0.18021421790926451981239097993472059435697140430d0
           phi0(2) = 0.68825027307369203697251883279276704348208128121d0
           phi0(3) = 1.14807540122775448526666118837266799247537387989d0
	   phi0(4) = 1.20229322030577344858231364022293527708627722917d0
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	     
           wei0(1) = fourpi/(1.0d0*nodes)
           wei0(2) = fourpi/(1.0d0*nodes)
	   wei0(3) = fourpi/(1.0d0*nodes)
           wei0(4) = fourpi/(1.0d0*nodes)
           wside0  = fourpi/(1.0d0*nodes)
	   wvert0  = fourpi/(1.0d0*nodes)
	   wface0  = zero 
c
	 elseif(n.eq.272) then
           nmax  = 27
c           ngen  = 3
	   nodes = 120*ngen + 12
	   mm    = 0
	   itype = 1
	   ifg   = 1 
c	
           open(unit=7,file='N27b.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*)  theta0(j), phi0(j)  
	   enddo
	   close(7)
c	  	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c	   
        elseif(n.eq.26) then
c          c=============c
c          c= Converges =c
c          c=============c	
           nmax     = 26
c           ngen      = 4
	   nodes     = 60*ngen + 12
	   itype     = 1
	   ifg       = 0
c
           open(unit=7,file='genxyz.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero	
c	  	  
         elseif(n.eq.24) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax = 24
c           ngen  = 3
	   nodes = 60*ngen + 20 + 12
  	   itype = 2
	   ifg   = 0
c
           gen(1,1) = zero
           gen(2,1) = 0.656509d0
           gen(3,1) = 0.476982d0 
           gen(4,1) = sqrt(1-gen(2,1)*gen(2,1) - gen(3,1)*gen(3,1))
c	  
           gen(1,2) = zero
           gen(2,2) = 0.166604d0
           gen(3,2) = 0.0690097d0
           gen(4,2) = sqrt(1-gen(2,2)*gen(2,2) - gen(3,2)*gen(3,2))
c	  
           gen(1,3) = zero
           gen(2,3) = 0.370986d0
           gen(3,3) = 0.510618d0
           gen(4,3) = sqrt(1-gen(2,3)*gen(2,3) - gen(3,3)*gen(3,3))
c
           do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo  	
c	  	  
           wei0(1) = fourpi/(1.0d0*nodes)
           wei0(2) = fourpi/(1.0d0*nodes)
           wei0(3) = fourpi/(1.0d0*nodes)
           wface0  = fourpi/(1.0d0*nodes)
	   wvert0  = fourpi/(1.0d0*nodes)
	   wside0  = zero	 	 
c
         elseif(n.eq.23) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax = 23
c           ngen  = 3
	   nodes = 60*ngen + 12
  	   itype = 1
	   ifg   = 0
c	   
           open(unit=7,file='3gen3D.dat',status='OLD') 
	   do j=1,ngen
	     gen(1,j) = zero
	     read(7,*)  gen(2,j), gen(3,j), gen(4,j)  
	   enddo
	   close(7)
	   
	   do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
 
c 
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c	   
           wvert0 = fourpi/(1.0d0*nodes)
	   wside0 = zero
           wface0 = zero
c
CCCCC
c           gen(1,1) = zero
c           gen(2,1) = 0.6205666579941735d0  
c 	    gen(3,1) = 0.6286715278136170d0
c           gen(4,1) = 0.4686887379726923d0
c	  
c           gen(1,2) = zero
c           gen(2,2) =  -.1037954894679029d0
c           gen(3,2) =  -.9934404366822338d0
c           gen(4,2) =  0.04798536371364673d0
c	  
c           gen(1,3) = zero
c           gen(2,3) = 0.2578794980492158d0
c           gen(3,3) = 0.6356620242079164d0
c           gen(4,3) = 0.7276207497493312d0
c
c           do j=1,ngen
c             theta0(j) = acos(gen(4,j))
c             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
c             cphi      = gen(2,j)/sr
c             phi0(j)   = acos(cphi)
c           enddo  	
c	  	  
c           wei0(1) = fourpi/(1.0d0*nodes)
c           wei0(2) = fourpi/(1.0d0*nodes)
c           wei0(3) = fourpi/(1.0d0*nodes)
c           wface0  = zero
c	    wvert0  = fourpi/(1.0d0*nodes)
c	    wside0  = zero	 	 
c
         elseif(n.eq.21) then
c          c=============c
c          c= Converges =c
c          c= w_f < 0   =c
c          c=============c	 
           nmax = 21
c           ngen  = 2
	   nodes = 60*ngen + 30 + 20+ 12
 	   mm    = 0
  	   itype = 4
	   ifg   = 0
c  
 	   gen(1,1) = zero
           gen(2,1) = 0.83920167645432497405266758505604d0
           gen(3,1) = 0.44709478033866528745576829351194d0 
           gen(4,1) = sqrt(1 - gen(2,1)*gen(2,1) - gen(3,1)*gen(3,1))
c	  
     	   gen(1,2) = zero
           gen(2,2) = 2.1718654858637769332069607149340d-9
           gen(3,2) = 0.95642852487143150239465591021752d0
           gen(4,2) = sqrt(1 - gen(2,2)*gen(2,2) - gen(3,2)*gen(3,2))	
   
           do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo	   		   
c	  	  
           wei0(1) = fourpi/(1.0d0*nodes)
           wei0(2) = 0.8d0*fourpi/(1.0d0*nodes)
           wside0  = fourpi/(1.0d0*nodes)
	   wvert0  = fourpi/(1.0d0*nodes)
	   wface0  = 0.8d0*fourpi/(1.0d0*nodes)
c	  	
         elseif(n.eq.20) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax  = 20
c           ngen  = 2
	   nodes = 60*ngen + 30 + 12
 	   mm    = 0
  	   itype = 3
	   ifg   = 0
c
	   gen(1,1) = zero
           gen(2,1) =  0.83920167645432497405266758505604d0
           gen(3,1) = -0.44709478033866528745576829351194d0 
           gen(4,1) = sqrt(1 - gen(2,1)*gen(2,1) - gen(3,1)*gen(3,1))
c	  
           gen(1,2) = zero
           gen(2,2) = 2.1718654858637769332069607149340d-9
           gen(3,2) = 0.95642852487143150239465591021752d0
           gen(4,2) = sqrt(1 - gen(2,2)*gen(2,2) - gen(3,2)*gen(3,2))	
   
           do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo	   		   
c	  	  
           wei0(1) = fourpi/(1.0d0*nodes)
           wei0(2) = fourpi/(1.0d0*nodes)
           wside0  = fourpi/(1.0d0*nodes)
	   wvert0  = fourpi/(1.0d0*nodes)
	   wface0  = zero
c	  
	 elseif(n.eq.19) then
c          c=============c
c          c= Converges =c
c          c=============c
           nmax  = 19
c           ngen  = 2
	   nodes = 60*ngen + 12
	   mm    = 0
	   itype = 1
	   ifg   = 0
c
	   gen(1,1) = zero
           gen(2,1) = 0.5882973389992774696d0
           gen(3,1) = 0.4181238515553936275d0 
           gen(4,1) = 0.6923330007135656164d0
           gen(1,2) = zero
           gen(2,2) = 0.2988140811238982425d0
           gen(3,2) = 0.4844155031786993688d0
           gen(4,2) = 0.8222245632603329231d0
c
           do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo
c	 
	   wei0(1) = 0.0986995175669551111d0
           wei0(2) = 0.0947571240643577593d0
           wvert0  = 0.0799143430400334012d0
	   wface0  = 0.0d0
	   wside0  = 0.0d0
c	   
         elseif(n.eq.17) then
c          c=============c
c          c= Converges =c
c          c=============c	 
           nmax  = 17
c           ngen  = 1
	   nodes = 60*ngen + 30 + 20 + 12
 	   mm    = 0
  	   itype = 4
	   ifg   = 0
c  
 	   gen(1,1) = zero
           gen(2,1) =  0.83920167645432497405266758505604d0
           gen(3,1) = -0.44709478033866528745576829351194d0 
           gen(4,1) = sqrt(1 - gen(2,1)*gen(2,1) - gen(3,1)*gen(3,1))
c   
           do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo	   		   
c	  	  
           wei0(1) = fourpi/(1.0d0*nodes)
           wside0  = fourpi/(1.0d0*nodes)
	   wvert0  = fourpi/(1.0d0*nodes)
	   wface0  = fourpi/(1.0d0*nodes)
c
         elseif(n.eq.15) then	 
           nmax  = 15
c           ngen  = 1
	   nodes = 120*ngen + 12
 	   mm    = 0
  	   itype = 1
	   ifg   = 1
c  
 	   gen(1,1) = zero
           gen(2,1) =  0.83920167645432497405266758505604d0
           gen(3,1) = -0.44709478033866528745576829351194d0 
           gen(4,1) = sqrt(1 - gen(2,1)*gen(2,1) - gen(3,1)*gen(3,1))
c   
           do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo	   		   
c	  	  
           wei0(1) = fourpi/(1.0d0*nodes)
           wside0  = zero
	   wvert0  = fourpi/(1.0d0*nodes)
	   wface0  = zero	   
c	
         elseif(n.eq.14) then	 
           nmax  = 14
c           ngen  = 1
	   nodes = 60*ngen + 12
  	   itype = 1
	   ifg   = 0
 
  	   gen(1,1) = zero
           gen(2,1) = 0.25d0
           gen(3,1) = 0.75d0
           gen(4,1) = sqrt(1.0d0 - gen(2,1)**2 - gen(3,1)**2)   		   

           do j=1,ngen
             theta0(j) = acos(gen(4,j))
             sr        = sqrt(gen(2,j)*gen(2,j)+gen(3,j)*gen(3,j))
             cphi      = gen(2,j)/sr
             phi0(j)   = acos(cphi)
           enddo

c          theta0(1) = 0.218372511747457d0
c	   phi0(1)   = 8.625891663086355d0
	  	  
c          wei0(1) = 0.1782729164537028d0   
c	   wvert0  = 0.1558329689280721d0

           wei0(1) = 5.0d0*pi/72.0d0 
	   wvert0  = 5.0d0*pi/72.0d0
	   
	   wface0  = zero
	   wside0  = zero
	   
c    	    gen(1,1) = zero
c           gen(2,1) = sin(theta0(1))*cos(phi0(1))
c           gen(3,1) = sin(theta0(1))*sin(phi0(1))
c           gen(4,1) = cos(theta0(1))
	   
         endif
	 	 
	 return
	 
       end
