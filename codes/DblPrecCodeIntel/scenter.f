          subroutine scenter(cntr)
c 
            implicit real *8 (a-h,o-z)
            real *8 cntr(4,30)
	    real *8 m, mp

c
            data    m/1.6180339887498948482045868343656381d0/
            data  mp/-0.6180339887498948482045868343656381d0/
            data  one/1.0d0/
	    data  two/2.0d0/
            data zero/0.0d0/
    	    
c
            cntr(1,1) = zero
            cntr(2,1) = one
            cntr(3,1) = zero
            cntr(4,1) = zero
c
            cntr(1,2) = zero
            cntr(2,2) = zero
            cntr(3,2) = one
            cntr(4,2) = zero
c
            cntr(1,3) = zero
            cntr(2,3) = zero
            cntr(3,3) = zero
            cntr(4,3) = one
c
            cntr(1,4) = zero
            cntr(2,4) = one/two
            cntr(3,4) = mp/two
            cntr(4,4) = m/two
c
            cntr(1,5) = zero
            cntr(2,5) = one/two
            cntr(3,5) = mp/two
            cntr(4,5) = -m/two 
c
            cntr(1,6) = zero
            cntr(2,6) = one/two
            cntr(3,6) = -mp/two
            cntr(4,6) = m/two 
c
            cntr(1,7) = zero
            cntr(2,7) = one/two
            cntr(3,7) = -mp/two
            cntr(4,7) = -m/two
c
            cntr(1,8) = zero
            cntr(2,8) = m/two
            cntr(3,8) = one/two
            cntr(4,8) = mp/two
c
            cntr(1,9) = zero
            cntr(2,9) = m/two
            cntr(3,9) = one/two
            cntr(4,9) = -mp/two
c
            cntr(1,10) = zero
            cntr(2,10) = -m/two
            cntr(3,10) = one/two
            cntr(4,10) = mp/two
c
            cntr(1,11) = zero
            cntr(2,11) = -m/two
            cntr(3,11) = one/two
            cntr(4,11) = -mp/two
c
            cntr(1,12) = zero
            cntr(2,12) = mp/two
            cntr(3,12) = m/two
            cntr(4,12) = one/two	    
c
            cntr(1,13) = zero
            cntr(2,13) = mp/two
            cntr(3,13) = -m/two
            cntr(4,13) = one/two	   
c
            cntr(1,14) = zero
            cntr(2,14) = -mp/two
            cntr(3,14) = m/two
            cntr(4,14) = one/two	   
c
            cntr(1,15) = zero
            cntr(2,15) = -mp/two
            cntr(3,15) = -m/two
            cntr(4,15) = one/two
c
            do k=1,15
              do j=1,4
                cntr(j,k+15) = -cntr(j,k)
              enddo
            enddo
c
c
            return
          end
