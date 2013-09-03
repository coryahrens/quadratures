c
c         G.B. June 2008 
c 
        program sphquad 
	implicit real *8 (a-h,o-z)
        parameter(maxterms=500)
        parameter(neqmax=maxterms*maxterms)
c
        real *8 q1(4),q2(4),qq(4)
        real *8 rg(4,60),rgc(4,60),qtemp(4,60),qres(4,60)
        real *8 pols(0:maxterms,60),der(0:maxterms,60)
        real *8 vert(4,12)
        real *8 vval(0:maxterms,12),dummy(0:maxterms,12)
        real *8 vsum(0:maxterms),fsum(0:maxterms)
        real *8 dsum(0:maxterms)
        real *8 p(0:maxterms,0:maxterms),dp(0:maxterms,0:maxterms)
        real *8 phire(0:maxterms),phiim(0:maxterms)
        real *8 ptemp(0:maxterms,0:maxterms)
        integer neqsub(0:maxterms),ipn(neqmax),ipm(neqmax)
c
c     ------------------------- 
         character*12 fname
c     ------------------------- 
c
        call prini(6,13)
        call getini
c
c
        two    = 2.0d0
        one    = 1.0d0
        pi     = 4.0d0*atan(1.0d0)
        fourpi = 4.0d0*pi
c
        nnn =100
        call sobolev(nnn,neqsub,ipn,ipm,neq)
cc        call prinf('neqsub *',neqsub,nnn)
cc        call prinf('neq *',neq,1)
cc        call prinf('ipn *',ipn,neq)
cc        call prinf('ipm *',ipm,neq)
cc        stop
c
c  set nn to test
c
c
        nn = 5
        mm = 2
c
        qq(1)=0
        qq(2)=pi
        qq(3)=pi*pi
        qq(4)=pi*pi*pi
        qn=qnorm(qq)
        qq(1)=qq(1)/qn
        qq(2)=qq(2)/qn
        qq(3)=qq(3)/qn
        qq(4)=qq(4)/qn
c
        call prin2('norm qq *',qnorm(qq),1)
c
c
         sr=sqrt(qq(2)*qq(2)+qq(3)*qq(3))
         si=qq(3)/sr
         co=qq(2)/sr
         ct=qq(4)
         call prin2x('ct *',ct,1)
         call prin2x('si *',si,1)
         call prin2x('co *',co,1)
c
c
         call spharm(maxterms,nn,nn,ct,si,co,p,dp,phire,phiim)
c
         call prin2('phire *',phire,nn)
         call prin2('phiim *',phiim,nn)
c
ccc        call assleg(maxterms,nn,mm,c,p,dp)
c
        do n=0,nn
          call prinf('==========n===========-*',n,1)
          do m=0,nn
          ptemp(m,n)=p(m,n)*phire(m)
          enddo
          call prinf('------------real-------------*',0,0)
          call prin2('p *',ptemp(0,n),nn+1)
          call prinf('------------imag-------------*',0,0)
          do m=0,nn
            ptemp(m,n)=p(m,n)*phiim(m)
          enddo
          call prin2('p *',ptemp(0,n),nn+1)
        enddo
c
c
        call prinf('------------------------*',0,0)
        do n=0,nn
          call prinf('----*',n,1)
          do m=0,nn
          ptemp(m,n)=dp(m,n)*phire(m)
          enddo
          call prinf('------------real-------------*',0,0)
          call prin2('dp *',ptemp(0,n),nn+1)
          call prinf('------------imag-------------*',0,0)
          do m=0,nn
            ptemp(m,n)=dp(m,n)*phiim(m)
          enddo
          call prin2('dp *',ptemp(0,n),nn+1)
        enddo
c
        stop
c
        call vertex(vert)
c
cc        do k=1,12
cc         call prinf('k= *',k,1)
cc         call prin2('vert *',vert(1,k),4)
cc        enddo
c
c       test multiplication of quaternions 
c
c        ee=exp(1.0d0)
c        call prin2('pi : *',pi, 1)
c        call prin2('ee : *',ee, 1)
c        q1(1)=pi
c        q1(2)=pi*pi
c        q1(3)=pi*pi*pi
c        q1(4)=pi*pi*pi*pi
c
c        q2(1)=ee
c        q2(2)=ee*ee
c        q2(3)=ee*ee*ee
c        q2(4)=ee*ee*ee*ee
c
c        call quatmult(q1,q2,qres(1,1))
c        call prin2('qres(1,1): *',qres(1,1),4)
c
c
c        test of rotation group
c
         call rotg(rg)       
c
         call prin2('rg *',rg(1,10),4)
c
        qq(1)=0
        qq(2)=pi
        qq(3)=pi*pi
        qq(4)=pi*pi*pi
        qn=qnorm(qq)
        qq(1)=qq(1)/qn
        qq(2)=qq(2)/qn
        qq(3)=qq(3)/qn
        qq(4)=qq(4)/qn
c
        call prin2('norm qq *',qnorm(qq),1)
c
         do k=1,60
             rgc(1,k) = rg(1,k)
             do j=2,4
             rgc(j,k)=-rg(j,k)
             enddo 
         enddo 
c
         do k=1,60
            call quatmult(qq,rgc(1,k),qtemp(1,k))
            call quatmult(rg(1,k),qtemp(1,k),qres(1,k))
         enddo 
         call prin2('qres *',qres(1,10),4)
c
c --------- at this point we have 60 images of the generator qq
c
c zero degree sp. harm. Sqrt[(2 n +1)/(4 Pi)] LegendreP[n,Cos[t]] 
c
c  set nn to test
c
c
         nn=16
         do k=1,12
         call lege1d(nn,vert(4,k),vval(0,k),dummy(0,k))
             do j=0,nn
             vval(j,k)= sqrt((two*j +one)/fourpi)*vval(j,k)
             enddo
         enddo
         do j=0,nn
            vsum(j)=zero
            do k=1,12
            vsum(j)=vsum(j)+ vval(j,k)
            enddo
         enddo
         call prin2('vsum *',vsum,nn)
c
c                 
c
         do k=1,60
c
cxxxxxxxxx          call lege1d(nn,qres(4,k),pols(0,k),der(0,k))
c
         sr=sqrt(qres(2,k)*qres(2,k)+qres(3,k)*qres(3,k))
         si=qres(3,k)/sr
         co=qres(2,k)/sr
         ct=qres(4,k)
         call spharm(maxterms,nn,nn,ct,si,co,p,dp,phire,phiim)
         do j=0,nn/2
         pols(2*j,k)=p(0,2*j)
         pols(2*j+1,k)=p(2,2*j+1)*phiim(2)
         enddo
cxxxxxxxxx             do j=0,nn
cxxxxxxxxx             pols(j,k)= sqrt((two*j +one)/fourpi)*pols(j,k)
cxxxxxxxxx             der(j,k) = sqrt((two*j +one)/fourpi)*der(j,k)
cxxxxxxxxx             enddo
         enddo 
         do j=0,nn
            fsum(j) = zero
            dsum(j) = zero
            do k=1,60
            fsum(j)=fsum(j)+ pols(j,k)
            dsum(j)=dsum(j)+ der(j,k)
            enddo
         enddo
c
         call prin2('fsum *',fsum,nn)
         call prin2('dsum *',dsum,nn)
ccc         do k=1,60
ccc         call prin2('pols *',pols(15,k),1)
ccc         enddo
c
c
         stop
         end

































