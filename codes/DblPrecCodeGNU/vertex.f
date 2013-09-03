c
c         G.B. June 2008 
c 
c         quaternions representing vertices of an icosohedron
c
          subroutine vertex(vert)
c 
            implicit real *8 (a-h,o-z)
            real *8 vert(4,12)
c
c           constant below is (1-Sqrt[5])/2 
c
            data grm/-0.6180339887498948482045868343656381d0/
            data five/5.0d0/
            data one/1.0d0/
            data zero/0.0d0/
c 
            fact=sqrt(-grm*sqrt(five))
c
            vert(1,1) = zero
            vert(2,1) = zero
            vert(3,1) = grm/fact
            vert(4,1) = one/fact
c
            vert(1,2) = zero
            vert(2,2) = zero
            vert(3,2) = grm/fact
            vert(4,2) = -one/fact
c
            vert(1,3) = zero
            vert(2,3) = one/fact
            vert(3,3) = zero
            vert(4,3) = grm/fact
c
            vert(1,4) = zero
            vert(2,4) = one/fact
            vert(3,4) = zero
            vert(4,4) = -grm/fact
c
            vert(1,5) = zero
            vert(2,5) = grm/fact
            vert(3,5) = one/fact
            vert(4,5) = zero 
c
            vert(1,6) = zero
            vert(2,6) = grm/fact
            vert(3,6) = -one/fact
            vert(4,6) = zero 
c
            do k=1,6
              do j=1,4
                vert(j,k+6) = -vert(j,k)
              enddo
            enddo
c
c
            return
          end
