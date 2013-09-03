          subroutine fcenter(cntr)
c 
            implicit real *8 (a-h,o-z)
            real *8 cntr(4,20)
	    real *8 m, mp

c
            data  m/1.6180339887498948482045868343656381d0/
            data mp/-0.6180339887498948482045868343656381d0/
            data one/1.0d0/
            data zero/0.0d0/
c 
            fact = sqrt(3.0d0)	    	    
c
            cntr(1,1) = zero
            cntr(2,1) = one/fact
            cntr(3,1) = one/fact
            cntr(4,1) = one/fact
c
            cntr(1,2) = zero
            cntr(2,2) = -one/fact
            cntr(3,2) = -one/fact
            cntr(4,2) = one/fact
c
            cntr(1,3) = zero
            cntr(2,3) = one/fact
            cntr(3,3) = -one/fact
            cntr(4,3) = -one/fact
c
            cntr(1,4) = zero
            cntr(2,4) = -one/fact
            cntr(3,4) = one/fact
            cntr(4,4) = -one/fact
c
            cntr(1,5) = zero
            cntr(2,5) = zero
            cntr(3,5) = m/fact
            cntr(4,5) = mp/fact 
c
            cntr(1,6) = zero
            cntr(2,6) = zero
            cntr(3,6) = m/fact
            cntr(4,6) = -mp/fact 
c
            cntr(1,7) = zero
            cntr(2,7) = mp/fact
            cntr(3,7) = zero
            cntr(4,7) = m/fact
c
            cntr(1,8) = zero
            cntr(2,8) = -mp/fact
            cntr(3,8) = zero
            cntr(4,8) = m/fact
c
            cntr(1,9) = zero
            cntr(2,9) = m/fact
            cntr(3,9) = mp/fact
            cntr(4,9) = zero
c
            cntr(1,10) = zero
            cntr(2,10) = m/fact
            cntr(3,10) = -mp/fact
            cntr(4,10) = zero
c
            do k=1,10
              do j=1,4
                cntr(j,k+10) = -cntr(j,k)
              enddo
            enddo
c
c
            return
          end
