c
c         G.B. June 2008 
c 
      real *8 function qnorm(qq)
      implicit real *8 (a-h,o-z)
      real *8 qq(4)
c     
      qnorm=sqrt(qq(1)*qq(1)+qq(2)*qq(2)+qq(3)*qq(3)+qq(4)*qq(4))
c
c
      return
      end
