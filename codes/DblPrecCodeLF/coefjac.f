c=========================================================================c
c                                                                         c
c Coefficients for calculating Jacobian --- values derived in Mathematica c
c file XXX.nb. See CA & GB Proc.Royal Soc. paper appendix for details     c
c                                                                         c
c                                                                         c
c         G.B. June 2008                                                  c
c=========================================================================c
c 
          subroutine coefjac(q,qres,ftheta1,fphi1,ftheta2,fphi2)
c 
          implicit real *8 (a-h,o-z)
          real *8 coeftheta1(3,60)
          real *8 coefphi1(2,60)
          real *8 coeftheta2(2,60)
          real *8 coefphi2(3,60)
          real *8 q(4),qres(4,60)
          real *8 ftheta1(60),fphi1(60)
          real *8 ftheta2(60),fphi2(60)
c
c
	  coeftheta1(1,1)=1.000000000000000000000000000000000d0
	  coeftheta1(2,1)=0.0d0
	  coeftheta1(3,1)=0.0d0
	  coeftheta1(1,2)=0.0d0
	  coeftheta1(2,2)=1.000000000000000000000000000000000d0
	  coeftheta1(3,2)=0.0d0
	  coeftheta1(1,3)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,3)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,3)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,4)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,4)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,4)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,5)=-1.000000000000000000000000000000000d0
	  coeftheta1(2,5)=0.0d0
	  coeftheta1(3,5)=0.0d0
	  coeftheta1(1,6)=0.0d0
	  coeftheta1(2,6)=-1.000000000000000000000000000000000d0
	  coeftheta1(3,6)=0.0d0
	  coeftheta1(1,7)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,7)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,7)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,8)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,8)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,8)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,9)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,9)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,9)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,10)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,10)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,10)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,11)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,11)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,11)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,12)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,12)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,12)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,13)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,13)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,13)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,14)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,14)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,14)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,15)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,15)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,15)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,16)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,16)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,16)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,17)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,17)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,17)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,18)=0.0d0
	  coeftheta1(2,18)=0.0d0
	  coeftheta1(3,18)=1.000000000000000000000000000000000d0
	  coeftheta1(1,19)=0.0d0
	  coeftheta1(2,19)=0.0d0
	  coeftheta1(3,19)=-1.000000000000000000000000000000000d0
	  coeftheta1(1,20)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,20)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,20)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,21)=0.0d0
	  coeftheta1(2,21)=0.0d0
	  coeftheta1(3,21)=1.000000000000000000000000000000000d0
	  coeftheta1(1,22)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,22)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,22)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,23)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,23)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,23)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,24)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,24)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,24)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,25)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,25)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,25)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,26)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,26)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,26)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,27)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,27)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,27)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,28)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,28)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,28)=0.5000000000000000000000000000000000d0
	  coeftheta1(1,29)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,29)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,29)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,30)=1.000000000000000000000000000000000d0
	  coeftheta1(2,30)=0.0d0
	  coeftheta1(3,30)=0.0d0
	  coeftheta1(1,31)=0.0d0
	  coeftheta1(2,31)=-1.000000000000000000000000000000000d0
	  coeftheta1(3,31)=0.0d0
	  coeftheta1(1,32)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,32)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,32)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,33)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,33)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,33)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,34)=-1.000000000000000000000000000000000d0
	  coeftheta1(2,34)=0.0d0
	  coeftheta1(3,34)=0.0d0
	  coeftheta1(1,35)=0.0d0
	  coeftheta1(2,35)=1.000000000000000000000000000000000d0
	  coeftheta1(3,35)=0.0d0
	  coeftheta1(1,36)=0.0d0
	  coeftheta1(2,36)=0.0d0
	  coeftheta1(3,36)=-1.000000000000000000000000000000000d0
	  coeftheta1(1,37)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,37)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,37)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,38)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,38)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,38)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,39)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,39)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,39)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,40)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,40)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,40)=-0.8090169943749474241022934171828191d0
	  coeftheta1(1,41)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,41)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,41)=-0.5000000000000000000000000000000000d0
	  coeftheta1(1,42)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,42)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,42)=-0.5000000000000000000000000000000000d0
	  coeftheta1(1,43)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,43)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,43)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,44)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,44)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,44)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,45)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,45)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,45)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,46)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,46)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,46)=-0.3090169943749474241022934171828191d0
	  coeftheta1(1,47)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,47)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,47)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,48)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,48)=0.3090169943749474241022934171828191d0
	  coeftheta1(3,48)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,49)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,49)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,49)=-0.5000000000000000000000000000000000d0
	  coeftheta1(1,50)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,50)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,50)=-0.5000000000000000000000000000000000d0
	  coeftheta1(1,51)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,51)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,51)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,52)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,52)=-0.5000000000000000000000000000000000d0
	  coeftheta1(3,52)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,53)=0.5000000000000000000000000000000000d0
	  coeftheta1(2,53)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,53)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,54)=-0.5000000000000000000000000000000000d0
	  coeftheta1(2,54)=-0.3090169943749474241022934171828191d0
	  coeftheta1(3,54)=0.8090169943749474241022934171828191d0
	  coeftheta1(1,55)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,55)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,55)=-0.5000000000000000000000000000000000d0
	  coeftheta1(1,56)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,56)=-0.8090169943749474241022934171828191d0
	  coeftheta1(3,56)=-0.5000000000000000000000000000000000d0
	  coeftheta1(1,57)=-0.8090169943749474241022934171828191d0
	  coeftheta1(2,57)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,57)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,58)=0.8090169943749474241022934171828191d0
	  coeftheta1(2,58)=0.5000000000000000000000000000000000d0
	  coeftheta1(3,58)=0.3090169943749474241022934171828191d0
	  coeftheta1(1,59)=-0.3090169943749474241022934171828191d0
	  coeftheta1(2,59)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,59)=-0.5000000000000000000000000000000000d0
	  coeftheta1(1,60)=0.3090169943749474241022934171828191d0
	  coeftheta1(2,60)=0.8090169943749474241022934171828191d0
	  coeftheta1(3,60)=-0.5000000000000000000000000000000000d0
c
c
	  coefphi1(1,1)=0.0d0
	  coefphi1(2,1)=-1.000000000000000000000000000000000d0
	  coefphi1(1,2)=1.000000000000000000000000000000000d0
	  coefphi1(2,2)=0.0d0
	  coefphi1(1,3)=0.5000000000000000000000000000000000d0
	  coefphi1(2,3)=0.8090169943749474241022934171828191d0
	  coefphi1(1,4)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,4)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,5)=0.0d0
	  coefphi1(2,5)=1.000000000000000000000000000000000d0
	  coefphi1(1,6)=-1.000000000000000000000000000000000d0
	  coefphi1(2,6)=0.0d0
	  coefphi1(1,7)=0.8090169943749474241022934171828191d0
	  coefphi1(2,7)=0.3090169943749474241022934171828191d0
	  coefphi1(1,8)=0.8090169943749474241022934171828191d0
	  coefphi1(2,8)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,9)=0.3090169943749474241022934171828191d0
	  coefphi1(2,9)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,10)=0.3090169943749474241022934171828191d0
	  coefphi1(2,10)=0.5000000000000000000000000000000000d0
	  coefphi1(1,11)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,11)=0.5000000000000000000000000000000000d0
	  coefphi1(1,12)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,12)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,13)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,13)=0.3090169943749474241022934171828191d0
	  coefphi1(1,14)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,14)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,15)=0.5000000000000000000000000000000000d0
	  coefphi1(2,15)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,16)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,16)=0.8090169943749474241022934171828191d0
	  coefphi1(1,17)=0.5000000000000000000000000000000000d0
	  coefphi1(2,17)=0.8090169943749474241022934171828191d0
	  coefphi1(1,18)=0.0d0
	  coefphi1(2,18)=0.0d0
	  coefphi1(1,19)=0.0d0
	  coefphi1(2,19)=0.0d0
	  coefphi1(1,20)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,20)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,21)=0.0d0
	  coefphi1(2,21)=0.0d0
	  coefphi1(1,22)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,22)=0.5000000000000000000000000000000000d0
	  coefphi1(1,23)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,23)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,24)=0.3090169943749474241022934171828191d0
	  coefphi1(2,24)=0.5000000000000000000000000000000000d0
	  coefphi1(1,25)=0.8090169943749474241022934171828191d0
	  coefphi1(2,25)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,26)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,26)=0.3090169943749474241022934171828191d0
	  coefphi1(1,27)=0.8090169943749474241022934171828191d0
	  coefphi1(2,27)=0.3090169943749474241022934171828191d0
	  coefphi1(1,28)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,28)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,29)=0.3090169943749474241022934171828191d0
	  coefphi1(2,29)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,30)=0.0d0
	  coefphi1(2,30)=-1.000000000000000000000000000000000d0
	  coefphi1(1,31)=-1.000000000000000000000000000000000d0
	  coefphi1(2,31)=0d0
	  coefphi1(1,32)=0.5000000000000000000000000000000000d0
	  coefphi1(2,32)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,33)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,33)=0.8090169943749474241022934171828191d0
	  coefphi1(1,34)=0.0d0
	  coefphi1(2,34)=1.000000000000000000000000000000000d0
	  coefphi1(1,35)=1.000000000000000000000000000000000d0
	  coefphi1(2,35)=0.0d0
	  coefphi1(1,36)=0.0d0
	  coefphi1(2,36)=0.0d0
	  coefphi1(1,37)=0.3090169943749474241022934171828191d0
	  coefphi1(2,37)=0.5000000000000000000000000000000000d0
	  coefphi1(1,38)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,38)=0.5000000000000000000000000000000000d0
	  coefphi1(1,39)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,39)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,40)=0.3090169943749474241022934171828191d0
	  coefphi1(2,40)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,41)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,41)=0.3090169943749474241022934171828191d0
	  coefphi1(1,42)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,42)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,43)=0.5000000000000000000000000000000000d0
	  coefphi1(2,43)=0.8090169943749474241022934171828191d0
	  coefphi1(1,44)=0.5000000000000000000000000000000000d0
	  coefphi1(2,44)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,45)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,45)=0.8090169943749474241022934171828191d0
	  coefphi1(1,46)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,46)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,47)=0.3090169943749474241022934171828191d0
	  coefphi1(2,47)=0.5000000000000000000000000000000000d0
	  coefphi1(1,48)=0.3090169943749474241022934171828191d0
	  coefphi1(2,48)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,49)=0.8090169943749474241022934171828191d0
	  coefphi1(2,49)=0.3090169943749474241022934171828191d0
	  coefphi1(1,50)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,50)=0.3090169943749474241022934171828191d0
	  coefphi1(1,51)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,51)=0.8090169943749474241022934171828191d0
	  coefphi1(1,52)=-0.5000000000000000000000000000000000d0
	  coefphi1(2,52)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,53)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,53)=-0.5000000000000000000000000000000000d0
	  coefphi1(1,54)=-0.3090169943749474241022934171828191d0
	  coefphi1(2,54)=0.5000000000000000000000000000000000d0
	  coefphi1(1,55)=0.8090169943749474241022934171828191d0
	  coefphi1(2,55)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,56)=-0.8090169943749474241022934171828191d0
	  coefphi1(2,56)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,57)=0.5000000000000000000000000000000000d0
	  coefphi1(2,57)=0.8090169943749474241022934171828191d0
	  coefphi1(1,58)=0.5000000000000000000000000000000000d0
	  coefphi1(2,58)=-0.8090169943749474241022934171828191d0
	  coefphi1(1,59)=0.8090169943749474241022934171828191d0
	  coefphi1(2,59)=0.3090169943749474241022934171828191d0
c
c
	  coeftheta2(1,1)=0.0d0
	  coeftheta2(2,1)=-1.000000000000000000000000000000000d0
	  coeftheta2(1,2)=1.000000000000000000000000000000000d0
	  coeftheta2(2,2)=0.0d0
	  coeftheta2(1,3)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,3)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,4)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,4)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,5)=0.0d0
	  coeftheta2(2,5)=1.000000000000000000000000000000000d0
	  coeftheta2(1,6)=-1.000000000000000000000000000000000d0
	  coeftheta2(2,6)=0.0d0
	  coeftheta2(1,7)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,7)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,8)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,8)=-0.3090169943749474241022934171828191d0
	  coeftheta2(1,9)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,9)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,10)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,10)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,11)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,11)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,12)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,12)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,13)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,13)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,14)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,14)=-0.3090169943749474241022934171828191d0
	  coeftheta2(1,15)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,15)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,16)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,16)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,17)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,17)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,18)=0.0d0
	  coeftheta2(2,18)=0.0d0
	  coeftheta2(1,19)=0.0d0
	  coeftheta2(2,19)=0.0d0
	  coeftheta2(1,20)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,20)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,21)=0.0d0
	  coeftheta2(2,21)=0.0d0
	  coeftheta2(1,22)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,22)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,23)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,23)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,24)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,24)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,25)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,25)=-0.3090169943749474241022934171828191d0
	  coeftheta2(1,26)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,26)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,27)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,27)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,28)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,28)=-0.3090169943749474241022934171828191d0
	  coeftheta2(1,29)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,29)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,30)=0.0d0
	  coeftheta2(2,30)=-1.000000000000000000000000000000000d0
	  coeftheta2(1,31)=-1.000000000000000000000000000000000d0
	  coeftheta2(2,31)=0.0d0
	  coeftheta2(1,32)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,32)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,33)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,33)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,34)=0.0d0
	  coeftheta2(2,34)=1.000000000000000000000000000000000d0
	  coeftheta2(1,35)=1.000000000000000000000000000000000d0
	  coeftheta2(2,35)=0.0d0
	  coeftheta2(1,36)=0.0d0
	  coeftheta2(2,36)=0.0d0
	  coeftheta2(1,37)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,37)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,38)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,38)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,39)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,39)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,40)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,40)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,41)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,41)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,42)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,42)=-0.3090169943749474241022934171828191d0
	  coeftheta2(1,43)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,43)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,44)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,44)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,45)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,45)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,46)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,46)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,47)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,47)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,48)=0.3090169943749474241022934171828191d0
	  coeftheta2(2,48)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,49)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,49)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,50)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,50)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,51)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,51)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,52)=-0.5000000000000000000000000000000000d0
	  coeftheta2(2,52)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,53)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,53)=-0.5000000000000000000000000000000000d0
	  coeftheta2(1,54)=-0.3090169943749474241022934171828191d0
	  coeftheta2(2,54)=0.5000000000000000000000000000000000d0
	  coeftheta2(1,55)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,55)=-0.3090169943749474241022934171828191d0
	  coeftheta2(1,56)=-0.8090169943749474241022934171828191d0
	  coeftheta2(2,56)=-0.3090169943749474241022934171828191d0
	  coeftheta2(1,57)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,57)=0.8090169943749474241022934171828191d0
	  coeftheta2(1,58)=0.5000000000000000000000000000000000d0
	  coeftheta2(2,58)=-0.8090169943749474241022934171828191d0
	  coeftheta2(1,59)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,59)=0.3090169943749474241022934171828191d0
	  coeftheta2(1,60)=0.8090169943749474241022934171828191d0
	  coeftheta2(2,60)=-0.3090169943749474241022934171828191d0
	  coefphi1(1,60)=0.8090169943749474241022934171828191d0
	  coefphi1(2,60)=-0.3090169943749474241022934171828191d0
c
c
	  coefphi2(1,1)=-1.0000000000000000000000000000000000d0
	  coefphi2(2,1)=0.0d0
	  coefphi2(3,1)=0.0d0
	  coefphi2(1,2)=0.0d0
	  coefphi2(2,2)=-1.000000000000000000000000000000000d0
	  coefphi2(3,2)=0.0d0
	  coefphi2(1,3)=0.8090169943749474241022934171828191d0
	  coefphi2(2,3)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,3)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,4)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,4)=0.5000000000000000000000000000000000d0
	  coefphi2(3,4)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,5)=1.000000000000000000000000000000000d0
	  coefphi2(2,5)=0.0d0
	  coefphi2(3,5)=0.0d0
	  coefphi2(1,6)=0.0d0
	  coefphi2(2,6)=1.000000000000000000000000000000000d0
	  coefphi2(3,6)=0.0d0
	  coefphi2(1,7)=0.3090169943749474241022934171828191d0
	  coefphi2(2,7)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,7)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,8)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,8)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,8)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,9)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,9)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,9)=0.8090169943749474241022934171828191d0
	  coefphi2(1,10)=0.5000000000000000000000000000000000d0
	  coefphi2(2,10)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,10)=0.8090169943749474241022934171828191d0
	  coefphi2(1,11)=0.5000000000000000000000000000000000d0
	  coefphi2(2,11)=0.3090169943749474241022934171828191d0
	  coefphi2(3,11)=0.8090169943749474241022934171828191d0
	  coefphi2(1,12)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,12)=0.3090169943749474241022934171828191d0
	  coefphi2(3,12)=0.8090169943749474241022934171828191d0
	  coefphi2(1,13)=0.3090169943749474241022934171828191d0
	  coefphi2(2,13)=0.8090169943749474241022934171828191d0
	  coefphi2(3,13)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,14)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,14)=0.8090169943749474241022934171828191d0
	  coefphi2(3,14)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,15)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,15)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,15)=0.3090169943749474241022934171828191d0
	  coefphi2(1,16)=0.8090169943749474241022934171828191d0
	  coefphi2(2,16)=0.5000000000000000000000000000000000d0
	  coefphi2(3,16)=0.3090169943749474241022934171828191d0
	  coefphi2(1,17)=0.8090169943749474241022934171828191d0
	  coefphi2(2,17)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,17)=0.3090169943749474241022934171828191d0
	  coefphi2(1,18)=0.0d0
	  coefphi2(2,18)=0.0d0
	  coefphi2(3,18)=-1.000000000000000000000000000000000d0
	  coefphi2(1,19)=0.0d0
	  coefphi2(2,19)=0.0d0
	  coefphi2(3,19)=1.000000000000000000000000000000000d0
	  coefphi2(1,20)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,20)=0.5000000000000000000000000000000000d0
	  coefphi2(3,20)=0.3090169943749474241022934171828191d0
	  coefphi2(1,21)=0.0d0
	  coefphi2(2,21)=0.0d0
	  coefphi2(3,21)=-1.000000000000000000000000000000000d0
	  coefphi2(1,22)=0.5000000000000000000000000000000000d0
	  coefphi2(2,22)=0.3090169943749474241022934171828191d0
	  coefphi2(3,22)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,23)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,23)=0.3090169943749474241022934171828191d0
	  coefphi2(3,23)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,24)=0.5000000000000000000000000000000000d0
	  coefphi2(2,24)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,24)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,25)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,25)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,25)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,26)=0.3090169943749474241022934171828191d0
	  coefphi2(2,26)=0.8090169943749474241022934171828191d0
	  coefphi2(3,26)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,27)=0.3090169943749474241022934171828191d0
	  coefphi2(2,27)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,27)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,28)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,28)=0.8090169943749474241022934171828191d0
	  coefphi2(3,28)=-0.5000000000000000000000000000000000d0
	  coefphi2(1,29)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,29)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,29)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,30)=-1.000000000000000000000000000000000d0
	  coefphi2(2,30)=0.0d0
	  coefphi2(3,30)=0.0d0
	  coefphi2(1,31)=0.0d0
	  coefphi2(2,31)=1.000000000000000000000000000000000d0
	  coefphi2(3,31)=0d0
	  coefphi2(1,32)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,32)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,32)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,33)=0.8090169943749474241022934171828191d0
	  coefphi2(2,33)=0.5000000000000000000000000000000000d0
	  coefphi2(3,33)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,34)=1.000000000000000000000000000000000d0
	  coefphi2(2,34)=0.0d0
	  coefphi2(3,34)=0.0d0
	  coefphi2(1,35)=0.0d0
	  coefphi2(2,35)=-1.000000000000000000000000000000000d0
	  coefphi2(3,35)=0.0d0
	  coefphi2(1,36)=0.0d0
	  coefphi2(2,36)=0.0d0
	  coefphi2(3,36)=1.000000000000000000000000000000000d0
	  coefphi2(1,37)=0.5000000000000000000000000000000000d0
	  coefphi2(2,37)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,37)=0.8090169943749474241022934171828191d0
	  coefphi2(1,38)=0.5000000000000000000000000000000000d0
	  coefphi2(2,38)=0.3090169943749474241022934171828191d0
	  coefphi2(3,38)=0.8090169943749474241022934171828191d0
	  coefphi2(1,39)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,39)=0.3090169943749474241022934171828191d0
	  coefphi2(3,39)=0.8090169943749474241022934171828191d0
	  coefphi2(1,40)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,40)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,40)=0.8090169943749474241022934171828191d0
	  coefphi2(1,41)=0.3090169943749474241022934171828191d0
	  coefphi2(2,41)=0.8090169943749474241022934171828191d0
	  coefphi2(3,41)=0.5000000000000000000000000000000000d0
	  coefphi2(1,42)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,42)=0.8090169943749474241022934171828191d0
	  coefphi2(3,42)=0.5000000000000000000000000000000000d0
	  coefphi2(1,43)=0.8090169943749474241022934171828191d0
	  coefphi2(2,43)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,43)=0.3090169943749474241022934171828191d0
	  coefphi2(1,44)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,44)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,44)=0.3090169943749474241022934171828191d0
	  coefphi2(1,45)=0.8090169943749474241022934171828191d0
	  coefphi2(2,45)=0.5000000000000000000000000000000000d0
	  coefphi2(3,45)=0.3090169943749474241022934171828191d0
	  coefphi2(1,46)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,46)=0.5000000000000000000000000000000000d0
	  coefphi2(3,46)=0.3090169943749474241022934171828191d0
	  coefphi2(1,47)=0.5000000000000000000000000000000000d0
	  coefphi2(2,47)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,47)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,48)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,48)=-0.3090169943749474241022934171828191d0
	  coefphi2(3,48)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,49)=0.3090169943749474241022934171828191d0
	  coefphi2(2,49)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,49)=0.5000000000000000000000000000000000d0
	  coefphi2(1,50)=0.3090169943749474241022934171828191d0
	  coefphi2(2,50)=0.8090169943749474241022934171828191d0
	  coefphi2(3,50)=0.5000000000000000000000000000000000d0
	  coefphi2(1,51)=0.8090169943749474241022934171828191d0
	  coefphi2(2,51)=0.5000000000000000000000000000000000d0
	  coefphi2(3,51)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,52)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,52)=0.5000000000000000000000000000000000d0
	  coefphi2(3,52)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,53)=-0.5000000000000000000000000000000000d0
	  coefphi2(2,53)=0.3090169943749474241022934171828191d0
	  coefphi2(3,53)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,54)=0.5000000000000000000000000000000000d0
	  coefphi2(2,54)=0.3090169943749474241022934171828191d0
	  coefphi2(3,54)=-0.8090169943749474241022934171828191d0
	  coefphi2(1,55)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,55)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,55)=0.5000000000000000000000000000000000d0
	  coefphi2(1,56)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,56)=0.8090169943749474241022934171828191d0
	  coefphi2(3,56)=0.5000000000000000000000000000000000d0
	  coefphi2(1,57)=0.8090169943749474241022934171828191d0
	  coefphi2(2,57)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,57)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,58)=-0.8090169943749474241022934171828191d0
	  coefphi2(2,58)=-0.5000000000000000000000000000000000d0
	  coefphi2(3,58)=-0.3090169943749474241022934171828191d0
	  coefphi2(1,59)=0.3090169943749474241022934171828191d0
	  coefphi2(2,59)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,59)=0.5000000000000000000000000000000000d0
	  coefphi2(1,60)=-0.3090169943749474241022934171828191d0
	  coefphi2(2,60)=-0.8090169943749474241022934171828191d0
	  coefphi2(3,60)=0.5000000000000000000000000000000000d0
c
c
          do k=1,60
	  
            sr  = dsqrt(q(2)*q(2) + q(3)*q(3))
            sth = dsqrt(1.0d0 - q(4)*q(4))
	      
            qresdenom = qres(2,k)*qres(2,k) + qres(3,k)*qres(3,k)
c
            ftheta1(k) = (coeftheta1(1,k)*q(2)          +
     1                    coeftheta1(2,k)*q(3))*q(4)/sr +
     2                    coeftheta1(3,k)*sth
c
            fphi1(k)   = coefphi1(1,k)*q(2) + coefphi1(2,k)*q(3) 

c
            ftheta2(k) = (coeftheta2(1,k)*q(2)+
     1                    coeftheta2(2,k)*q(3))/(sr*qresdenom)
c

            fphi2(k)   = ((coefphi2(1,k)*q(2) + 
     1                     coefphi2(2,k)*q(3))*q(4) + 
     2                     coefphi2(3,k)*sth**2 )/qresdenom

c
          enddo
c
          return
        
	end
