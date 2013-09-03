c
c         G.B. June 2008 
c 
          subroutine quatmult(q1,q2,qres)
c 
          implicit real *8 (a-h,o-z)
          real *8 q1(4),q2(4),qres(4)
c
c         multiplication of quaternions with no attempt to optimize the code 
c
          qres(1) = q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4)
          qres(2) = q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3)
          qres(3) = q1(1)*q2(3)+q1(3)*q2(1)+q1(4)*q2(2)-q1(2)*q2(4)
          qres(4) = q1(1)*q2(4)+q1(4)*q2(1)+q1(2)*q2(3)-q1(3)*q2(2)
c
          return
          end
