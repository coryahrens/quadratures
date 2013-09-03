       subroutine input(n,nmax,ngen,theta0,phi0,wei0,wvert0,
     1                  wface0,wside0,mm,itype,ifg)
         implicit real *8 (a-h,o-z)
         real *8 theta0(ngen), phi0(ngen)
	 real *8 wei0(ngen), wvert0, wface0, wside0

c----------------------------c	
         if(n.eq.223) then
           nmax  = 223
           ngen  = 312
	   nodes = 60*ngen + 20 + 12
	   mm    = 2
	   itype = 2
	   ifg   = 0    
c
           open(unit=7,file='InAnglesN223.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*) theta0(j), phi0(j) 
	   enddo
	   close(7)

           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
           enddo	
c	  	  
           wface0 = fourpi/(1.0d0*nodes)
	   wvert0 = 0.8d0*fourpi/(1.0d0*nodes)
	   wside0 = zero
	   
c
         elseif(n.eq.177) then
           nmax = 177
           ngen  = 200
	   nodes = 60*ngen + 30 + 20 + 12
 	   mm    = 0
  	   itype = 4
	   ifg   = 0
c	 
	 
         elseif(n.eq.144) then
           nmax  = 144
           ngen  = 136
	   nodes = 60*ngen + 20 + 12
	   mm    = 0
	   itype = 2
	   ifg   = 0
c
           open(unit=7,file='InAnglesN144.dat',status='OLD') 
	   do j=1,ngen
	     read(7,*) theta0(j), phi0(j) 
	   enddo
	   close(7)

           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
           enddo	
c	  	  
           wface0 = fourpi/(1.0d0*nodes)
	   wvert0 = 0.8d0*fourpi/(1.0d0*nodes)
	   wside0 = zero
c	  
         elseif(n.eq.131) then
           nmax  = 131
           ngen  = 115
	   nodes = 60*ngen + 12
	   mm    = 2
	   itype = 1
	   ifg   = 0
c	
           smt = 0.0d0
	   smp = 0.0d0
           open(unit=7,file='n130theta.dat',status='OLD') 
	   do j=1,ngen-2
	     read(7,*) theta0(j)
	     smt = smt + theta0(j)  
	   enddo
	   close(7)
	   open(unit=7,file='n130phi.dat',status='OLD') 
	   do j=1,ngen-2
	     read(7,*) phi0(j)
	     smp = smp + phi0(j)  
	   enddo
	   close(7)
c	  
           theta0(ngen - 1) = smt/(1.0d0*(ngen-2)) 
	   theta0(ngen)     = theta0(ngen - 1) - 1.0d-2
	   phi0(ngen - 1)   = smp/(1.0d0*(ngen-2))
	   phi0(ngen)       = phi0(ngen - 1) - 1.0d-3
	  
           do j=1,ngen
	     wei0(j) = fourpi/(1.0d0*nodes)
	   enddo
c 
           wside0 = zero
           wface0 = zero
           wvert0 = fourpi/(1.0d0*nodes)
c	   	
	 elseif(n.eq.130) then
           nmax  = 130
           ngen  = 113
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
	elseif(n.eq.90) then
          nmax  = 90
          ngen  = 56
	  nodes = 60*ngen + 12
	  mm    = 0
	  itype = 1
	  ifg   = 0 
c	
          open(unit=7,file='N90.dat',status='OLD') 
	  do j=1,ngen
	    read(7,*)  theta0(j), phi0(j)  
	  enddo
	  close(7)
c
c	  do j=1,ngen
c	    phi0(j) = -1.0d0*phi0(j)
c	  enddo
	  	  
          do j=1,ngen
	    wei0(j) = fourpi/(1.0d0*nodes)
	  enddo
c 
          wvert0 = fourpi/(1.0d0*nodes)
	  wside0 = zero
          wface0 = zero
c  
        elseif(n.eq.40) then
          nmax  = 40
          ngen  = 12
	  nodes = 60*ngen + 12
	  mm    = 0 
	  itype = 1
	  ifg   = 0
c     
	  theta0(1)  = 1.13715650743740398506901862314542d0 
	  theta0(2)  = 0.9840339496337197407738334715349d0
	  theta0(3)  = 1.13557875539228816691596929245292d0
	  theta0(4)  = 0.7031563530579348031264780934891d0
	  theta0(5)  = 1.28221684837322503473817035092264d0
	  theta0(6)  = 1.2101173103341434013645793932959d0
	  theta0(7)  = 0.9928229983247869448696959425478d0
	  theta0(8)  = 1.1397451875172817038471692357282d0
	  theta0(9)  = 1.45433999663722357702969057353621d0
	  theta0(10) = 1.2899464624470464601401988567817d0
	  theta0(11) = 1.0805175172100070900712180951641d0
	  theta0(12) = 1.4970260301003631920995983250363d0 
c	                   	 
	  phi0(1)  =  -0.00923485040026692419367708473503d0
	  phi0(2)  = 0.26186385276474526637215185819509d0 
	  phi0(3)  = 0.25636807601276525715043893739223d0 
	  phi0(4)  = 0.1698780867663118744928476959756d0
	  phi0(5)  = 0.23965249682185714382995882392140d0 
	  phi0(6)  = 0.37530487879516096374536850477115d0
	  phi0(7)  = 0.5347062572812928469583609598306d0
	  phi0(8)  = 0.5091386905870656382763035314890d0 
	  phi0(9)  = -0.18363595549448831849727344858203d0 
	  phi0(10) = 0.487989550990427550821553504693d0
	  phi0(11) = 0.6551778920174871199640828676548d0
	  phi0(12) = 0.3325853043207359747824584140511d0
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
	  nmax = 34
	  ngen  = 8
	  nodes = 60*ngen + 12
	  mm    = 0
	  itype = 1
	  ifg   = 0
c 
          open(unit=7,file='N34b.dat',status='OLD') 
	  do j=1,ngen
	    read(7,*) theta0(j), phi0(j)  
	  enddo
	  close(7)
c	  
          do j=1,ngen
	    wei0(j) = fourpi/(1.0d0*nodes)
	  enddo
c 
          wvert0    = fourpi/(1.0d0*nodes)
	  wside0 = zero
          wface0 = zero
c
c          theta0(1) = 0.2609917649433072d+01
c          theta0(2) = 0.1315658192744376d+01
c          theta0(3) = 0.1092209431815446d+01
c          theta0(4) = 0.1320606274272225d+01
c          theta0(5) = 0.4113612870800398d+00
c          theta0(6) = 0.2943634105538610d+01
c          theta0(7) = 0.2548162781691452d+01
c          theta0(8) = 0.1818066848674111d+01
c 
c          phi0(1)   = 0.2386687751345448d+01
c          phi0(2)   = 0.2976210101251069d+01
c          phi0(3)   = 0.2228010916977576d+01
c          phi0(4)   = 0.4956319532124677d+00
c          phi0(5)   = 0.6056323961348414d+01
c          phi0(6)   = 0.5562917824983825d+01
c          phi0(7)   = 0.4416437611987314d+01
c          phi0(8)   = 0.1148761998328532d+01
c 
c          wei0(1)   = 0.2409920694582632d-01
c          wei0(2)   = 0.2758636828427778d-01
c          wei0(3)   = 0.2332880455499617d-01
c          wei0(4)   = 0.2407030324181014d-01
c          wei0(5)   = 0.2239007407446654d-01
c          wei0(6)   = 0.2806388234292413d-01
c          wei0(7)   = 0.2754223856217822d-01
c          wei0(8)   = 0.2755197145550283d-01
c
c          wvert0    = 0.2403330388668789d-01	  
c
	elseif(n.eq.32) then
	  nmax  = 32
          ngen  = 7
	  nodes = 60*ngen + 12
	  neq   = 21
	  mm    = 2
	  itype = 1
	  ifg   = 0
c
          theta0(1) = 0.177683577730264d+01 
          theta0(2) = 0.1300307975063925d+01
          theta0(3) = 0.1059272963452778d+01
          theta0(4) = 0.8026312137796407d+00
          theta0(5) = 0.2941596161927831d+01
          theta0(6) = 0.1983712243622630d+01
	  sm = 0.0d0
	  do j=1,6
	    sm = sm + theta0(j)
	  enddo
	  theta0(7) = sm/6.0d0
c
          phi0(1) = 0.3460277225335893d+01
          phi0(2) = 0.6097932457843251d+01
          phi0(3) = 0.1777437150516778d+01
          phi0(4) = 0.6136775160787267d+01
          phi0(5) = 0.8867243967924299d+00 
          phi0(6) = 0.2787690864874015d+01
	  sm = 0.0d0
	  do j=1,6
	    sm = sm + phi0(j)
	  enddo
	  phi0(7)   =  sm/6.0d0
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
          nmax  = 31
          ngen  = 6
	  nodes = 60*ngen + 12
	  mm    = 0
	  itype = 1
	  ifg   = 0
c
          theta0(1) = 0.177683577730264d+01 
          theta0(2) = 0.1300307975063925d+01
          theta0(3) = 0.1059272963452778d+01
          theta0(4) = 0.8026312137796407d+00
          theta0(5) = 0.2941596161927831d+01
          theta0(6) = 0.1983712243622630d+01
c
          phi0(1) = 0.3460277225335893d+01
          phi0(2) = 0.6097932457843251d+01
          phi0(3) = 0.1777437150516778d+01
          phi0(4) = 0.6136775160787267d+01
          phi0(5) = 0.8867243967924299d+00 
          phi0(6) = 0.2787690864874015d+01
c
          wei0(1) = 0.3991277824630497d-01
          wei0(2) = 0.3563392770902429d-01
          wei0(3) = 0.3485019633098942d-01
          wei0(4) = 0.3379777203184547d-01
          wei0(5) = 0.2760879566031245d-01 
          wei0(6) = 0.3223176368718324d-01
c
          wvert0 = 0.2702138286829945d-01
	  wside0 = zero
          wface0 = zero
c	   	 
	elseif(n.eq.30) then
          nmax  = 30
          ngen  = 5
	  nodes = 60*ngen + 30 + 20 + 12
	  mm    = 2
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
	  ngen  = 5
	  nodes = 60*ngen + 12
	  mm    = 2
	  itype = 1
	  ifg   = 0
c
          theta0(1) = 0.8238091031919634205354069791068
          theta0(2) = 0.48539096450974506456390720001268
          theta0(3) = 0.41896134224315124671948841045490
          theta0(4) = 1.2352006762033047081854300819211
          theta0(5) = 1.1911489463508967267792422716752
c
	  phi0(1) = 0.7739051194765860442474598029292
          phi0(2) = 1.16400760393468220062242292289504
          phi0(3) = 0.6683019620108908215579623510807
          phi0(4) = 0.9881003719304176748396822098814
          phi0(5) = 0.4891397532128760033166122600101
c	
          do j=1,ngen-1
	    wei0(j) = fourpi/(1.0d0*nodes)
	  enddo
c 
          wvert0 = fourpi/(1.0d0*nodes)
	  wside0 = zero
          wface0 = zero
c	
        elseif(n.eq.27) then
          nmax  = 27
          ngen  = 4
	  nodes = 60*ngen + 30 + 12
	  mm    = 0
	  itype = 3
	  ifg   = 0
c	  
	  theta0(1) = 0.476074175265880923002205391549631305624787541664d0
          theta0(2) = 0.422733902178372561408733496790331374713357794138d0
          theta0(3) = 0.226421122343199283177695078504384246534740006404d0
	  theta0(4) = 0.455313411663907116830681694491432493104495849093d0
c	 
          phi0(1) = 0.180214217909264519812390979934720594356971404301d0
          phi0(2) = 0.688250273073692036972518832792767043482081281251d0
          phi0(3) = 1.148075401227754485266661188372667992475373879893d0
	  phi0(4) = 1.202293220305773448582313640222935277086277229171d0
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
        elseif(n.eq.26) then
           nmax     = 26
           ngen      = 4
	   nodes     = 60*ngen + 12
	   mm        = 0
	   itype     = 1
	   ifg       = 0
c
           theta0(1) = 0.28802E+01
           theta0(2) = 0.72996E+00
           theta0(3) = 0.13605E+01
           theta0(4) = 0.26032E+01
c	 
           phi0(1)   = 0.14599E+01
           phi0(2)   = 0.27211E+01
           phi0(3)   = 0.52064E+01
           phi0(4)   = 0.33729E+01	
c	  	  
           wei0(1)   = 0.49867E-01
           wei0(2)   = 0.49867E-01
           wei0(3)   = 0.49867E-01
           wei0(4)   = 0.49867E-01
	   wvert0    = 0.49867E-01
	   wside0    = zero
	   wface0    = zero
c
         elseif(n.eq.24) then
           nmax = 24
           ngen  = 3
	   nodes = 60*ngen + 20 + 12
 	   mm    = 0
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
         elseif(n.eq.21) then
           nmax = 21
           ngen  = 2
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
           nmax  = 20
           ngen  = 2
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
           nmax  = 19
           ngen  = 2
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
           nmax = 17
           ngen  = 1
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
         endif
	 
	 return
	 
       end
