      program Dcr2d
c****************************************************************************
c
c                                D C R 2 D
c
c     This program solves the system of equations
c
c              -  (K_ij U,j),i + A_i U,i + S U = F
c
c     where U = (u1,u2,un), F = (f1,f2,fn) and K_ij, A_i and S are  n x n
c     given matrices (n = 2 or 3).
c
c     Ramon Codina
c     Last revision: October, 1997
c
c****************************************************************************
      implicit none
      integer
     .  maxpoin,maxelem,maxnode,maxgaus,maxdofn,maxevab
      parameter (maxpoin=  5000, maxelem=  5000, maxnode= 9, 
     .           maxgaus=     9, maxdofn=  3,    maxevab=27)
      real*8   time1,time2
      integer
     .  lnods(maxnode,maxelem),   ifpre(maxdofn*maxpoin)
      real*8
     .  coord(2,maxpoin),         elcod(2,maxnode),
     .  unkno(maxdofn,maxpoin),   rhsgl(maxdofn,maxpoin),
     .  rhslo(maxdofn,maxnode)
      real*8
     .  posgx(maxgaus),           posgy(maxgaus), 
     .  weigp(maxgaus),           derst(2,maxnode,maxgaus),
     .  shape(maxnode,maxgaus),   hesst(3,maxnode,maxgaus),   
     .  derxy(2,maxnode,maxgaus), hesxy(3,maxnode,maxgaus),
     .  pertu(maxevab,maxdofn),   resid(maxdofn,maxevab),
     .  amate(maxevab,maxevab)

      call cputim(
     .  time1)
      call opfile
      call inputd(
     .  coord,ifpre,lnods,posgx,posgy,weigp,unkno)
      call solver(
     .  coord,lnods,ifpre,rhsgl,rhslo,shape,derxy,
     .  hesxy,unkno,amate,posgx,posgy,elcod,derst,
     .  hesst,weigp,resid,pertu)
      call cputim(
     .  time2)
      call output(
     .  coord,lnods,unkno,time1,time2)

      stop
      end
      subroutine cputim(rtime)
c***********************************************************************
c*
c*   This routine finds out the CPU time in seconds
c*
c***********************************************************************
      implicit none
      real*8 rtime
      real*4 stime(2)
      real*4 ETIME
      external ETIME

      rtime = ETIME(stime)

      end
      subroutine opfile
c*****************************************************************************
c
c     Open files
c
c*****************************************************************************
      implicit none
      integer nr,nw,np,iout
      character*20 arx_lect,arx_escr,arx_post
      common/inpout/nr,nw,np,iout
      data    nr,nw,np /8,9,10/
      
      open(7,file='Dcr2d.que',status='old')
      read(7,*) arx_lect
      read(7,*) arx_escr
      read(7,*) arx_post
      open(nr,file=arx_lect,status='old')
      open(nw,file=arx_escr,status='unknown')
      open(np,file=arx_post,form='unformatted', status='unknown')

      end
      subroutine inputd(coord,ifpre,lnods,posgx,posgy,weigp,unkno)
c*************************************************************************
c
c     Read most of the input data, renumber the equations and allocate
c     memory for the algebraic system      
c
c*************************************************************************
      implicit none
      integer
     .  npoin,nelem,nnode,ngaut,ndofn,nevab,ksoty,
     .  kprec,kstab,ktaum,ntotv,nr,nw,np,iout,idumy
      real*8
     .  hnatu,patau,difma(3,3,2,2),conma(3,3,2),reama(3,3),
     .  force(3),dummy,cputi
      real*8
     .  young,poiss,thick,zer
      integer
     .  ifpre(*), lnods(*),plate
      real*8
     .  coord(*), posgx(*), posgy(*), unkno(*), weigp(*)
      common/contro/npoin,nelem,nnode,ngaut,ndofn,nevab
      common/proper/difma,conma,reama,force
      common/numert/hnatu,patau,ksoty,kprec,kstab,ktaum
      common/inpout/nr,nw,np,iout

      zer = 0.0
      call initia(difma,36,zer)
      call initia(conma,18,zer)
      call initia(reama,9 ,zer)
      call initia(force,3 ,zer)

      read(nr,1) npoin,nelem,nnode,ngaut,ndofn
      if(ndofn.eq.2) then
        read(nr,2) difma(1,1,1,1),difma(1,2,1,1),
     .             difma(2,1,1,1),difma(2,2,1,1)
        read(nr,2) difma(1,1,1,2),difma(1,2,1,2),
     .             difma(2,1,1,2),difma(2,2,1,2)
        read(nr,2) difma(1,1,2,2),difma(1,2,2,2),
     .             difma(2,1,2,2),difma(2,2,2,2)
        read(nr,2) conma(1,1,1)  ,conma(1,2,1),
     .             conma(2,1,1)  ,conma(2,2,1)
        read(nr,2) conma(1,1,2)  ,conma(1,2,2),
     .             conma(2,1,2)  ,conma(2,2,2)
        read(nr,2) reama(1,1)    ,reama(1,2),
     .             reama(2,1)    ,reama(2,2)
        read(nr,3) force(1)      ,force(2)
      else if(ndofn.eq.3) then                              
        read(nr,5)
     .    difma(1,1,1,1),difma(1,2,1,1),difma(1,3,1,1),
     .    difma(2,1,1,1),difma(2,2,1,1),difma(2,3,1,1),
     .    difma(3,1,1,1),difma(3,2,1,1),difma(3,3,1,1)
        read(nr,5)
     .    difma(1,1,1,2),difma(1,2,1,2),difma(1,3,1,2),
     .    difma(2,1,1,2),difma(2,2,1,2),difma(2,3,1,2),
     .    difma(3,1,1,2),difma(3,2,1,2),difma(3,3,1,2)
        read(nr,5)
     .    difma(1,1,2,2),difma(1,2,2,2),difma(1,3,2,2),
     .    difma(2,1,2,2),difma(2,2,2,2),difma(2,3,2,2),
     .    difma(3,1,2,2),difma(3,2,2,2),difma(3,3,2,2)
        read(nr,5)
     .    conma(1,1,1)  ,conma(1,2,1)  ,conma(1,3,1), 
     .    conma(2,1,1)  ,conma(2,2,1)  ,conma(2,3,1),
     .    conma(3,1,1)  ,conma(3,2,1)  ,conma(3,3,1)
        read(nr,5)
     .    conma(1,1,2)  ,conma(1,2,2)  ,conma(1,3,2), 
     .    conma(2,1,2)  ,conma(2,2,2)  ,conma(2,3,2),
     .    conma(3,1,2)  ,conma(3,2,2)  ,conma(3,3,2)
        read(nr,5)
     .    reama(1,1)    ,reama(1,2)    ,reama(1,3),
     .    reama(2,1)    ,reama(2,2)    ,reama(2,3),
     .    reama(3,1)    ,reama(3,2)    ,reama(3,3)
        read(nr,3)
     .    force(1)      ,force(2)      ,force(3)
      end if

      plate = 0
      if(ndofn.eq.-3) then
        read(nr,3) young,poiss,thick,force(3)
        ndofn=3
        plate=1
      end if                                                

      read(nr,4) ksoty,kprec,hnatu,kstab,ktaum,patau,iout
      call geodat(coord,ifpre,lnods,posgx,posgy,weigp,unkno)

      if(plate.eq.1)
     .  call plamat(young,poiss,thick,difma,conma,reama)

      nevab = ndofn*nnode

c The slash / descriptor begins a new line (record) on output and skips to the next line on input, ignoring any unread information on the current record format(6/) o 6/
    1 format(6(/),5(39x,i10,/))
    2 format(39x,2(e15.5),/,39x,2(e15.5))
    3 format(39x,4(e15.5))
    4 format(/,2(39x,i10,/),(39x,e15.5,/),2(39x,i10,/),
     .  (39x,e15.5,/),(39x,i10,/)/)
    5 format(39x,3(e15.5),/,39x,3(e15.5),/,39x,3(e15.5))

      difma(1,1,2,1)=difma(1,1,1,2)
      difma(1,2,2,1)=difma(2,1,1,2)
      difma(2,1,2,1)=difma(1,2,1,2)
      difma(2,2,2,1)=difma(2,2,1,2)
      if(ndofn.eq.3) then
        difma(1,3,2,1)=difma(3,1,1,2)
        difma(2,3,2,1)=difma(3,2,1,2)
        difma(3,1,2,1)=difma(1,3,1,2)
        difma(3,2,2,1)=difma(2,3,1,2)
        difma(3,3,2,1)=difma(3,3,1,2)
      end if

      ntotv=ndofn*npoin
      call inisol
c task 1 is renumber of the equation      
      call MatMan(
     .  lnods,ifpre,dummy,dummy,ndofn,nnode,nnode,
     .  idumy,nelem,npoin,ntotv,idumy,   nw,idumy,
     .  idumy,    1,idumy,idumy,cputi,    1)
c task 2 is memory allocation of iprob
      call MatMan(
     .  lnods,idumy,dummy,dummy,ndofn,nnode,nnode,
     .  nnode,nelem,npoin,ntotv,    1,   nw,idumy,
     .  1,    1,idumy,idumy,cputi,    2)
      end
      subroutine geodat(
     .  coord,ifpre,lnods,posgx,posgy,weigp,unkno)
c***************************************************************************
c
c     Read geometrical data
c
c***************************************************************************
      implicit none
      integer
     .  npoin,nelem,nnode,ngaut,ndofn,nevab,nr,nw,np,iout,
     .  ipoin,numel,nofix,iofix,ielem,inode,jpoin,idime
      common/contro/npoin,nelem,nnode,ngaut,ndofn,nevab
      common/inpout/nr,nw,np,iout
      integer  ifpre(ndofn,npoin), lnods(nnode,nelem)  
      real*8
     .  coord(2,npoin), unkno(ndofn,npoin),  posgx(ngaut),
     .  posgy(ngaut),   weigp(ngaut) 
c
c***  Initializations
c
      do ipoin=1,npoin
        ifpre(    1,ipoin)=0
        ifpre(    2,ipoin)=0
        ifpre(ndofn,ipoin)=0
        unkno(    1,ipoin)=0.0
        unkno(    2,ipoin)=0.0
        unkno(ndofn,ipoin)=0.0
      end do
c
c***  Geometry
c      
      do ielem=1,nelem                                      
        read(nr,*)  numel,(lnods(inode,numel),inode=1,nnode)
      end do
      do jpoin=1,npoin
        read(nr,*)  ipoin,(coord(idime,ipoin),idime=1,2)
      end do
c
c***  Fixity conditions
c
      if(ndofn.eq.2) then
        nofix=1
        do while(nofix.lt.npoin) 
          read(nr,*)  iofix,ifpre(1,iofix),ifpre(2,iofix),
     .      unkno(1,iofix),unkno(2,iofix)
          nofix=iofix
        end do
      else if (ndofn.eq.3) then
        nofix=1
        do while(nofix.lt.npoin) 
          read(nr,*)
     .      iofix,ifpre(1,iofix),ifpre(2,iofix),
     .      ifpre(3,iofix),unkno(1,iofix),unkno(2,iofix),
     .      unkno(3,iofix)
          nofix=iofix
        end do
      end if
        
      call gaussq(posgx,posgy,weigp)

      end
      subroutine gaussq(posgx,posgy,weigp)
c**************************************************************************
c
c     Local coordinates and weights for integration points
c
c**************************************************************************
      implicit none
      integer  npoin,nelem,nnode,ngaut,ndofn,nevab
      real*8   ex1,et1,ez1,ex2,et2,ez2,a,b
      real*8   posgx(ngaut), posgy(ngaut), weigp(ngaut)
      common/contro/npoin,nelem,nnode,ngaut,ndofn,nevab

      if(nnode.eq.4. or .nnode.eq.9) then
        if(ngaut.eq.1) then
          posgx(1)=0.0
          posgy(1)=0.0
          weigp(1)=2.0
        else if(ngaut.eq.4) then
          posgx(1)=-sqrt(1.0/3.0)
          posgx(2)= sqrt(1.0/3.0)
          posgx(3)=-sqrt(1.0/3.0)
          posgx(4)= sqrt(1.0/3.0)
          posgy(1)=-sqrt(1.0/3.0)
          posgy(2)=-sqrt(1.0/3.0)
          posgy(3)= sqrt(1.0/3.0)
          posgy(4)= sqrt(1.0/3.0)
          weigp(1)=1.0
          weigp(2)=1.0
          weigp(3)=1.0
          weigp(4)=1.0
        else if(ngaut.eq.9) then
          posgx(1)=-sqrt(3.0/5.0)
          posgx(2)= 0.0
          posgx(3)= sqrt(3.0/5.0)
          posgx(4)=-sqrt(3.0/5.0)
          posgx(5)= 0.0
          posgx(6)= sqrt(3.0/5.0)
          posgx(7)=-sqrt(3.0/5.0)
          posgx(8)= 0.0
          posgx(9)= sqrt(3.0/5.0)
          posgy(1)=-sqrt(3.0/5.0)
          posgy(2)=-sqrt(3.0/5.0)
          posgy(3)=-sqrt(3.0/5.0)
          posgy(4)= 0.0
          posgy(5)= 0.0
          posgy(6)= 0.0
          posgy(7)= sqrt(3.0/5.0)
          posgy(8)= sqrt(3.0/5.0)
          posgy(9)= sqrt(3.0/5.0)
          weigp(1)= 25.0/81.0
          weigp(2)= 40.0/81.0
          weigp(3)= 25.0/81.0
          weigp(4)= 40.0/81.0
          weigp(5)= 64.0/81.0
          weigp(6)= 40.0/81.0
          weigp(7)= 25.0/81.0
          weigp(8)= 40.0/81.0
          weigp(9)= 25.0/81.0
        end if
      else if(nnode.eq.3. or .nnode.eq.6) then
        if(ngaut.eq.1) then
          posgx(1)=1.0/3.0
          posgy(1)=1.0/3.0
          weigp(1)=0.5
        else if(ngaut.eq.3) then
          posgx(1)=0.0
          posgx(2)=0.5
          posgx(3)=0.5
          posgy(1)=0.5
          posgy(2)=0.0
          posgy(3)=0.5
          weigp(1)=1.0/6.0
          weigp(2)=1.0/6.0
          weigp(3)=1.0/6.0
        else if(ngaut.eq.4) then
          posgx(1)=  1.0/3.0
          posgx(2)=  1.0/5.0
          posgx(3)=  1.0/5.0
          posgx(4)=  3.0/5.0
          posgy(1)=  1.0/3.0
          posgy(2)=  3.0/5.0
          posgy(3)=  1.0/5.0
          posgy(4)=  1.0/5.0
          weigp(1)=-27.0/96.0
          weigp(2)= 25.0/96.0
          weigp(3)= 25.0/96.0
          weigp(4)= 25.0/96.0
        else if(ngaut.eq.6) then
          ex1 = 0.81684 75729 80459
          et1 = 0.09157 62135 09771
          ez1 = 0.09157 62135 09771
          ex2 = 0.10810 30181 68070
          et2 = 0.44594 84909 15965
          ez2 = 0.44594 84909 15965
          posgx(1)= ex1
          posgx(2)= et1
          posgx(3)= ez1
          posgx(4)= ex2
          posgx(5)= et2
          posgx(6)= ez2
          posgy(1)= et1
          posgy(2)= ez1
          posgy(3)= ex1
          posgy(4)= et2
          posgy(5)= ez2
          posgy(6)= ex2
          a = 0.054975870996713638
          b = 0.1116907969117165    
          weigp(1)  = a
          weigp(2)  = a
          weigp(3)  = a
          weigp(4)  = b
          weigp(5)  = b
          weigp(6)  = b
        end if
      end if

      end
      subroutine initia(a,n,c)
c**********************************************************************
c
c     Initialize the vector a(n) to the value c (all components)
c
c**********************************************************************
      implicit none
      integer  n,i
      real*8   c
      real*8 a(n)

      do i=1,n
        a(i)=c
      end do

      end
      subroutine jacob2(
     .  djacb,elcod,shape,derst,hesst,derxy,
     .  hesxy,xjacm,xjaci,nnode)
c**************************************************************************
c
c     Compute the Jacobian matrix, its determinant and its inverse,
c     the Cartesian derivatives and the Hessian matrix
c
c**************************************************************************
      implicit none
      integer  nr,nw,np,iout,idime,jdime,inode,nnode
      real*8   djacb
      real*8
     .  elcod(2,nnode), shape(nnode),   derst(2,nnode),
     .  derxy(2,nnode), xjacm(2,2),     xjaci(2,2),
     .  hesst(3,nnode), hesxy(3,nnode)
      common/inpout/nr,nw,np,iout
      
      do idime=1,2
        do jdime=1,2
          xjacm(idime,jdime)=0.0
          do inode=1,nnode
            xjacm(idime,jdime)=xjacm(idime,jdime)+
     .       derst(idime,inode)*elcod(jdime,inode)
          end do
        end do
      end do

      djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
      if(djacb.lt.1.e-8)
     .  call runend('Element with non-positive Jacobian')
      xjaci(1,1)= xjacm(2,2)/djacb
      xjaci(2,2)= xjacm(1,1)/djacb
      xjaci(1,2)=-xjacm(1,2)/djacb
      xjaci(2,1)=-xjacm(2,1)/djacb

      do idime=1,2
        do inode=1,nnode
          derxy(idime,inode)=0.0
          do jdime=1,2
            derxy(idime,inode)=derxy(idime,inode)+
     .       xjaci(idime,jdime)*derst(jdime,inode)
          end do
        end do
      end do

      do inode=1,nnode                                      
        hesxy(1,inode)
     .    =xjaci(1,1)*xjaci(1,1)*hesst(1,inode)+
     .   2.0*xjaci(1,1)*xjaci(2,1)*hesst(2,inode)+
     .   xjaci(2,1)*xjaci(2,1)*hesst(3,inode)
        hesxy(2,inode)
     .    =xjaci(1,1)*xjaci(1,2)*hesst(1,inode)+
     .    (xjaci(1,1)*xjaci(2,2)+xjaci(2,1)*xjaci(1,2))
     .    *hesst(2,inode)+
     .    xjaci(2,1)*xjaci(2,2)*hesst(3,inode)
        hesxy(3,inode)=xjaci(2,2)*xjaci(2,2)*hesst(3,inode)+
     .   2.0*xjaci(2,2)*xjaci(1,2)*hesst(2,inode)+
     .   xjaci(1,2)*xjaci(1,2)*hesst(1,inode)
      end do

      end
      subroutine output(coord,lnods,unkno,time1,time2)
c**********************************************************************
c
c     Output of results
c
c**********************************************************************
      implicit none
      integer  npoin,nelem,nnode,ngaut,ndofn,nevab,nr,
     .         nw,np,iout,ipoin
      real*8   time1,time2
      integer  lnods(nnode,nelem)
      real*8   coord(2,npoin),  unkno(ndofn*npoin)
      common/contro/npoin,nelem,nnode,ngaut,ndofn,nevab
      common/inpout/nr,nw,np,iout
c
c*** Output to formatted file
c
      time2=time2-time1
      write(nw,800) time2
      if(iout.eq.1) then
        write(nw,900)
        if(ndofn.eq.2) then
          do ipoin=1,npoin
            write(nw,910)
     .        ipoin,coord(1,ipoin),coord(2,ipoin),
     .        unkno(2*ipoin-1),unkno(2*ipoin)
          end do
        else if(ndofn.eq.3) then
          do ipoin=1,npoin
            write(nw,910)                                   
     .        ipoin,coord(1,ipoin),coord(2,ipoin),
     .        unkno(3*ipoin-2),unkno(3*ipoin-1),
     .        unkno(3*ipoin)
          end do
        end if
      end if                                                
  800 format(/,5x,67('#'),//,5x,'Total CPU time :  ', e12.5)
  900 format(//,5x,'Point   x-coord.   y-coord.      ',
     .  'x-unknown   y-unknown ')
  910 format(5x,i5,2(2x,f9.4),3x,3f12.5)

c.......................... modificar      
      do ipoin = 1,21
        write(45,*) coord(1,ipoin), unkno(3*ipoin-2)
      end do
c
c*** Output to unformatted file (for post-process)
c
      write (np) nelem, npoin, nnode, ndofn
      write (np) coord
      write (np) lnods
      write (np) unkno

      end
      subroutine shape2(s,t,nnode,shape,derst,hesst)
c***********************************************************************
c
c     Evaluate the shape functions and its natural first derivatives 
c     and Hessian matrix
c
c***********************************************************************
      implicit none
      integer  nnode,n_ini
      real*8
     .  s,t,st,s2,t2,a1,a2,a3,ss,tt,s1,t1,s9,t9,v_ini
      real*8
     .  shape(nnode), derst(2,nnode), hesst(3,nnode)

      st = s*t
      s2 = 2.0*s
      t2 = 2.0*t
      n_ini = 3*nnode
      v_ini = 0.0
      call initia(hesst,n_ini,v_ini)

      if(nnode.eq.3) then
        shape(1)=1.0-s-t
        shape(2)=s
        shape(3)=t

        derst(1,1)=-1.0
        derst(1,2)= 1.0
        derst(1,3)= 0.0
        derst(2,1)=-1.0
        derst(2,2)= 0.0
        derst(2,3)= 1.0

      else if(nnode.eq.4) then
        shape(1)=(1.-t-s+st)*0.25
        shape(2)=(1.-t+s-st)*0.25
        shape(3)=(1.+t+s+st)*0.25
        shape(4)=(1.+t-s-st)*0.25

        derst(1,1)=(-1.+t)*0.25
        derst(1,2)=(+1.-t)*0.25
        derst(1,3)=(+1.+t)*0.25
        derst(1,4)=(-1.-t)*0.25
        derst(2,1)=(-1.+s)*0.25
        derst(2,2)=(-1.-s)*0.25
        derst(2,3)=(+1.+s)*0.25
        derst(2,4)=(+1.-s)*0.25

        hesst(2,1)= 0.25
        hesst(2,2)=-0.25
        hesst(2,3)= 0.25
        hesst(2,4)=-0.25

      else if(nnode.eq.6) then
        a1=1.0-s-t
        a2=s
        a3=t
        shape(1)=(2.0*a1-1.0)*a1
        shape(2)=(2.0*a2-1.0)*a2
        shape(3)=(2.0*a3-1.0)*a3
        shape(4)= 4.0*a1*a2
        shape(5)= 4.0*a2*a3
        shape(6)= 4.0*a1*a3

        derst(1,1)=1.0-4.0*a1
        derst(1,2)=4.0*a2-1.0
        derst(1,3)=0.0
        derst(1,4)=4.0*(a1-a2)
        derst(1,5)=4.0*a3
        derst(1,6)=-4.0*a3
        derst(2,1)=1.0-4.0*a1
        derst(2,2)=0.0
        derst(2,3)=4.0*a3-1.0
        derst(2,4)=-4.0*a2
        derst(2,5)=4.0*a2
        derst(2,6)=4.0*(a1-a3)

        hesst(1,1)= 4.0
        hesst(1,2)= 4.0
        hesst(1,4)=-8.0
        hesst(2,1)= 4.0
        hesst(2,4)=-4.0
        hesst(2,5)= 4.0
        hesst(2,6)=-4.0
        hesst(3,1)= 4.0
        hesst(3,3)= 4.0
        hesst(3,6)=-8.0

      else if(nnode.eq.9) then
        ss=s*s
        st=s*t
        tt=t*t
        s1=s+1.0
        t1=t+1.0
        s9=s-1.0
        t9=t-1.0
        shape(1)=0.25*s9*st*t9
        shape(2)=0.25*s1*st*t9
        shape(3)=0.25*s1*st*t1  
        shape(4)=0.25*s9*st*t1
        shape(5)=0.5*(1.0-ss)*t*t9
        shape(6)=0.5*s*s1*(1.0-tt)
        shape(7)=0.5*(1.0-ss)*t*t1
        shape(8)=0.5*s*s9*(1.0-tt)
        shape(9)=(1.0-ss)*(1.0-tt)

        derst(1,1)=0.25*t*t9*(-1.0+s2)
        derst(1,2)=0.25*(1.0+s2)*t*t9
        derst(1,3)=0.25*(1.0+s2)*t*t1
        derst(1,4)=0.25*(-1.0+s2)*t*t1
        derst(1,5)=-st*t9
        derst(1,6)=0.5*(1.0+s2)*(1.0-tt)
        derst(1,7)=-st*t1
        derst(1,8)=0.5*(-1.0+s2)*(1.0-tt)
        derst(1,9)=-s2*(1.0-tt)
        derst(2,1)=0.25*(-1.0+t2)*s*s9
        derst(2,2)=0.25*s*s1*(-1.0+t2)
        derst(2,3)=0.25*s*s1*(1.0+t2)
        derst(2,4)=0.25*s*s9*(1.0+t2)
        derst(2,5)=0.5*(1.0-ss)*(-1.0+t2)
        derst(2,6)=-st*s1
        derst(2,7)=0.5*(1.0-ss)*(1.0+t2)
        derst(2,8)=-st*s9
        derst(2,9)=-t2*(1.0-ss)

        hesst(1,1)= 0.5*t*t9
        hesst(1,2)= 0.5*t*t9
        hesst(1,3)= 0.5*t*t1
        hesst(1,4)= 0.5*t*t1
        hesst(1,5)=-t*t9
        hesst(1,6)= 1.0-tt
        hesst(1,7)=-t*t1
        hesst(1,8)= 1.0-tt
        hesst(1,9)=-2.0*(1.0-tt)
        hesst(2,1)= 0.25*(-1.0+t2)*(s9+s)
        hesst(2,2)= 0.25*(-1.0+t2)*(s1+s)
        hesst(2,3)= 0.25*( 1.0+t2)*(s1+s)
        hesst(2,4)= 0.25*( 1.0+t2)*(s9+s)
        hesst(2,5)=-s*(-1.0+t2)
        hesst(2,6)=-t*s1-st
        hesst(2,7)=-s*(1.0+t2)
        hesst(2,8)=-t*s9-st
        hesst(2,9)= s2*t2
        hesst(3,1)= 0.5*s*s9
        hesst(3,2)= 0.5*s*s1
        hesst(3,3)= 0.5*s*s1
        hesst(3,4)= 0.5*s*s9
        hesst(3,5)= 1.0-ss
        hesst(3,6)=-s*s1
        hesst(3,7)= 1.0-ss
        hesst(3,8)=-s*s9
        hesst(3,9)=-2.0*(1.0-ss)
      end if

      end
      subroutine solver(
     .  coord,lnods,ifpre,rhsgl,rhslo,shape,
     .  derxy,hesxy,unkno,amate,posgx,posgy,
     .  elcod,derst,hesst,weigp,resid,pertu)
c************************************************************************
c
c     Main routine
c
c************************************************************************
      implicit none
      integer                                               
     .  nr,nw,np,iout,npoin,nelem,nnode,ngaut,ndofn,nevab,
     .  kstab,ktaum,kprec,igaus,ielem,idumy,ntotv,ksoty,
     .  n_ini
      real*8
     .  hnatu,patau,difma(3,3,2,2),conma(3,3,2),reama(3,3),
     .  force(3),s,t,djacb,workm(2,2),hx,hy,hmaxi,cputi,
     .  dvolu,v_ini
      common/numert/hnatu,patau,ksoty,kprec,kstab,ktaum
      common/inpout/nr,nw,np,iout
      common/contro/npoin,nelem,nnode,ngaut,ndofn,nevab
      common/proper/difma,conma,reama,force 
      integer
     .  lnods(nnode,nelem),       ifpre(ndofn,npoin)
      real*8
     .  coord(2,npoin),           unkno(ndofn*npoin),
     .  shape(  nnode,ngaut),     derst(2,nnode,ngaut), 
     .  derxy(2,nnode,ngaut),     hesst(3,nnode,ngaut), 
     .  hesxy(3,nnode,ngaut),     posgx(ngaut),            
     .  posgy(ngaut),             weigp(ngaut),            
     .  xjacm(2,2),               xjaci(2,2),              
     .  elcod(2,nnode),           amate(nevab,nevab),  
     .  rhslo(nevab),             rhsgl(ndofn*npoin),           
     .  resid(ndofn,nevab),       pertu(nevab,ndofn),
     .  tauma(3,3)
c
c***  Construct the shape functions
c
      do igaus=1,ngaut
        s=posgx(igaus)
        t=posgy(igaus)
        call shape2(s,t,nnode,shape(1,igaus),
     .    derst(1,1,igaus),hesst(1,1,igaus))
      end do
      ntotv=ndofn*npoin
c
c***  Construct the matrix of the linear system
c
      v_ini = 0.0
      call initia(rhsgl,ntotv,v_ini)
      do ielem=1,nelem
        call gather(coord,elcod,lnods(1,ielem),2,nnode)
        n_ini = nevab*nevab
        call initia(amate,n_ini,v_ini)
        call initia(rhslo,nevab,v_ini)
        do igaus=1,ngaut
          call jacob2(djacb,elcod,shape(1,igaus),
     .      derst(1,1,igaus),hesst(1,1,igaus),
     .      derxy(1,1,igaus),hesxy(1,1,igaus),
     .      xjacm,xjaci,nnode)
          dvolu=djacb*weigp(igaus)
          hx=sqrt(xjaci(1,1)**2+xjaci(2,1)**2)
          hy=sqrt(xjaci(1,2)**2+xjaci(2,2)**2)
          hmaxi=hnatu/(min(hx,hy))
          call taumat(tauma,hmaxi,ndofn)
c
c***  Galerkin contribution to the system matrix and RHS          
c
          call galcon(
     .      nnode,ndofn,nevab,dvolu,shape(1,igaus),
     .      derxy(1,1,igaus),amate,rhslo)
c
c***  Contribution from the stabilization term          
c
          call stacon(
     .      nnode,ndofn,nevab,dvolu,pertu,workm,
     .      shape(1,igaus),derxy(1,1,igaus),
     .      hesxy(1,1,igaus),resid,tauma,amate,
     .      rhslo,kstab)
        end do                                              !igaus=1,ngaut
c
c***  Assembly
c
        call scater(
     .    rhsgl,rhslo,lnods(1,ielem),ndofn,nnode)      
        call assemb(
     .    lnods(1,ielem),amate,nnode,ndofn,nevab,
     .    ntotv,ksoty,kprec,nelem,ielem)
      end do                                                !ielem=1,nelem
c
c***  Solution of the linear system
c
c task 4 is the solution of linear system 
      call MatMan(
     .  lnods,ifpre,rhsgl,unkno,ndofn,nnode,nnode,
     .  nnode,nelem,npoin,ntotv,    1,   nw,    1,
     .  idumy,    1,    0,    1,cputi,    4)
  
      end 
      subroutine assemb(lnods,elmat,nnode,ndofn,nevab,ntotv,
     .                  ksoty,kprec,nelem,ielem)
c*****************************************************************************
c
c     Assembly of element matrices
c
c*****************************************************************************
      implicit  none
      integer   nnode,ndofn,nevab,ntotv,ksoty,kprec,nelem,ielem,
     .          inode,ipoin,idofn,ievab,itotv,jnode,jpoin,jdofn,
     .          jevab,jtotv,idumy,idest(2)
      real*8    value,dummy
      integer   lnods(nnode)
      real*8    elmat(nevab,nevab)
c
c***  Direct & SK iterative solver
c
      if((ksoty.eq.0).or.(kprec.ge.0)) then
        do inode=1,nnode
          ipoin=lnods(inode)
          do idofn=1,ndofn
            ievab=(inode-1)*ndofn+idofn
            itotv=(ipoin-1)*ndofn+idofn
            do jnode=1,nnode
              jpoin=lnods(jnode)
              do jdofn=1,ndofn
                jevab=(jnode-1)*ndofn+jdofn
                jtotv=(jpoin-1)*ndofn+jdofn
                value=elmat(ievab,jevab)
                idest(1)=itotv
                idest(2)=jtotv
c task 3 is the construction of linear system                                              
                call MatMan(                                
     .            idumy,idumy,dummy,dummy,idumy,idumy,idumy,
     .            idumy,idumy,idumy,idumy,idumy,idumy,idumy,
     .            idumy,    1,    1,idest,value,    3)
              end do
            end do
          end do
        end do
c
c***  Iterative solver
c task 3 is the construction of linear system     
      else if(kprec.lt.0) then
        call MatMan(
     .    idumy,idumy,dummy,dummy,ndofn,nnode,nnode,
     .    nnode,nelem,idumy,idumy,idumy,idumy,idumy,
     .    idumy,    1,nevab,ielem,elmat,    3)
      end if
      
      end
      subroutine scater(vecgl,veclo,lnods,ndofn,nnode)
c*****************************************************************************
c
c     Scater of vectors
c
c*****************************************************************************
      implicit none
      integer  inode,nnode,idofn,ndofn,ipoin,ievab,itotv
      integer  lnods(nnode)
      real*8   veclo(*), vecgl(*)

      do inode=1,nnode
        ipoin=lnods(inode)
        do idofn=1,ndofn
          ievab=(inode-1)*ndofn+idofn
          itotv=(ipoin-1)*ndofn+idofn
          vecgl(itotv)=vecgl(itotv)+veclo(ievab)
        end do
      end do

      end
      subroutine gather(vecgl,veclo,lnods,ndofn,nnode)
c*****************************************************************************
c
c     Gather operations
c
c*****************************************************************************
      implicit none
      integer   inode,nnode,idofn,ndofn,ipoin,ievab,itotv
      integer   lnods(nnode)
      real*8    veclo(ndofn*nnode), vecgl(*)

      do inode=1,nnode
        ipoin=lnods(inode)
        do idofn=1,ndofn
          ievab=(inode-1)*ndofn+idofn
          itotv=(ipoin-1)*ndofn+idofn
          veclo(ievab)=vecgl(itotv)
        end do
      end do

      end
      subroutine pertur(
     .  kstab,idofn,jdofn,workm,derxy,shape,pertu)
c*****************************************************************************
c
c     Peturbation of the test function according to the type of method:
c
c       SUPG    :                 A_i   V,i           (kstab=1)
c       GLS     : -(K_ij V,j),i + A_i   V,i + S   V   (kstab=2)
c       SGS, TG :  (K_ij V,j),i + A^t_i V,i - S^t V   (kstab=3,5)
c       CG      :            diag(A_i)  V,i           (kstab=4)
c
c*****************************************************************************
      implicit none
      integer
     .  kstab,idofn,jdofn,k,l
      real*8
     .  difma(3,3,2,2),conma(3,3,2),reama(3,3),force(3),
     .  workm(2,2),derxy(2),shape,pertu,prod1,prod2,prod3
      common/proper/difma,conma,reama,force 
c
c*** SUPG
c      
      if(kstab.eq.1) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        pertu=prod1
c
c*** Galerkin least squares
c      
      else if(kstab.eq.2) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            prod2=prod2+difma(jdofn,idofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(jdofn,idofn)*shape
        pertu=-prod2+prod1+prod3
c
c*** Subgrid scale & Taylor Galerkin
c      
      else if((kstab.eq.3).or.(kstab.eq.5)) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(idofn,jdofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            prod2=prod2+difma(idofn,jdofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(idofn,jdofn)*shape
        pertu=prod2+prod1-prod3
c
c*** Characteristic Galerkin
c      
      else if(kstab.eq.4) then
        prod1=0.0
        if(idofn.eq.jdofn) then
          do k=1,2
            prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
          end do
        end if
        pertu=prod1
      end if

      end 
      subroutine taumat(tauma,hmaxi,ndofn)
c*****************************************************************************
c
c     Matrix of intrinsic time scales, computed as
c
c     TAU = PATAU * [ 4 K / h^2 + 2 A / h + S ]^{-1}
c
c*****************************************************************************
      implicit none
      integer
     .  ndofn,kstab,ktaum,kprec,ksoty,i,j,k
      real*8
     .  hnatu,patau,difma(3,3,2,2),conma(3,3,2),reama(3,3),
     .  force(3),tauma(3,3),hmaxi,chadi(3,3),chaco(3,3),
     .  chare(3,3),tauin(3,3),a,b,c,tau,det,v_ini
      common/numert/hnatu,patau,ksoty,kprec,kstab,ktaum
      common/proper/difma,conma,reama,force

      v_ini = 0.0
      call initia(tauma,9,v_ini)
      if(kstab.eq.0) return
      call initia(tauin,9,v_ini)
      call initia(chaco,9,v_ini)
      call initia(chadi,9,v_ini)
c
c***  Characteristic convection matrix: A = sqrt | A_i A_i |
c
      do i=1,ndofn
        do j=1,ndofn
          chaco(i,j)=0.0
          do k=1,ndofn
            chaco(i,j) = chaco(i,j)
     .        + conma(i,k,1)*conma(k,j,1)
     .        + conma(i,k,2)*conma(k,j,2)
          end do
        end do
      end do
      call sqrtma(chaco,chaco,ndofn)
c
c***  Characteristic diffusion matrix: K = sqrt( K_ij K_ij )
c
      do i=1,ndofn
        do j=1,ndofn
          chadi(i,j)=0.0
          do k=1,ndofn
            chadi(i,j)=chadi(i,j)
     .        + difma(i,k,1,1)*difma(k,j,1,1)
     .        + difma(i,k,1,2)*difma(k,j,1,2)*2.0
     .        + difma(i,k,2,2)*difma(k,j,2,2)
          end do
        end do
      end do
      call sqrtma(chadi,chadi,ndofn)
c      
c***  Characteristic reaction matrix: S = | S |
c
      do i=1,ndofn
        do j=1,ndofn
          chare(i,j)=0.0
          do k=1,ndofn
            chare(i,j)=chare(i,j) + reama(i,k)*reama(k,j)
          end do
        end do
      end do
      call sqrtma(chare,chare,ndofn)
c
c***  Invers of the matrix of characteristic times
c
      do i=1,ndofn
        do j=1,ndofn
          tauin(i,j) = 4.0*chadi(i,j)/(hmaxi*hmaxi)
     .      + 2.0*chaco(i,j)/hmaxi + chare(i,j)
        end do
      end do
c
c***  Matrix tau, corresponding to:
c     KTAUM = 0: T = t I, where t is the minimum of all the admissible tau's
c           = 1: T = diag(t1,t2,t3), where ti is the minimum of the admissible
c                    tau's for the i-th row (equation)
c           = 2: T = [ 4 K / h^2 + 2 A / h + S ]^{-1}      
c
      if(ktaum.eq.0) then       
        tau = 0.0
        do i=1,ndofn
          do j=1,ndofn
            tau = max(tau,abs(tauin(i,j)))
          end do
        end do
        tau = patau/tau
        do i=1,ndofn
          tauma(i,i) = tau
        end do
      else if(ktaum.eq.1) then 
        a = 0.0
        b = 0.0
        c = 0.0
        do j=1,ndofn
          a = max(a,abs(tauin(    1,j)))
          b = max(b,abs(tauin(    2,j)))
          c = max(c,abs(tauin(ndofn,j)))
        end do
        a = patau/a
        b = patau/b
        c = patau/c
        tauma(    1,    1) = a
        tauma(    2,    2) = b
        tauma(ndofn,ndofn) = c
      else if(ktaum.eq.2) then
        call invmtx(tauin,tauma,det,ndofn)
        do i = 1,ndofn
          do j = 1,ndofn
            tauma(i,j) = tauma(i,j)*patau
          end do
        end do
        tauma(ndofn,ndofn) = 0.0
      else if(ktaum.eq.3) then                              
        a = 1.0/(patau*difma(1,1,1,1)
     .    /(hmaxi*hmaxi) + reama(1,1))
        tauma(1,1) = a
        tauma(2,2) = a
        a = (hmaxi*hmaxi*hmaxi*hmaxi)/(patau*patau)
        a = a*(patau/(hmaxi*hmaxi*reama(1,1))
     .    + 1.0d0/(difma(1,1,1,1)))
        tauma(3,3) = a
      end if

      end
      subroutine inisol
c*****************************************************************************
c
c     Solver parameters
c
c*****************************************************************************
      implicit none
      include 'MatMan.h'
      ! include '/Users/molivag/Dropbox/1.Doctorado/1.Research/Computing/Fortran/DCR/matman.free/Sources/include'
      integer ksoty,kprec,kstab,ktaum
      real*8  hnatu,patau
      common/numert/hnatu,patau,ksoty,kprec,kstab,ktaum

      lengp = 4
      ksymm(1)=1
c
c***  Direct solver (LDU decomposition with skyline storage)
c
      if(ksoty.eq.0) then
        krenu(1)=1
        kmsip(1)=1   
        kites(1)=11  
c
c***  Iterative solver
c
      else if(ksoty.ge.1) then
        krenu(1)=0
        if(kprec.ge.0) then
          kmsip(1)=kprec   
          kites(1)=30 + ksoty
          thres(1)= 1.0e-5
          lfill(1)= 5
        else if(kprec.lt.0) then
          kmsip(1)=kprec
          kites(1)=20 + ksoty
          if(ksoty.eq.1) ksymm(1) = 0
        end if
        itmax(1)=300
        toler(1)=1.0e-6
        kkryl(1)=50
      end if
      
      wosol(1)='DCR2D SOLVER   '
      conde(1)= .false.
      
      end
      subroutine runend(message)
c*****************************************************************************
c
c     Stop the run
c
c*****************************************************************************
      implicit none
      integer nr,nw,np,iout
      character message*(*)
      common/inpout/nr,nw,np,iout

      write(nw,'(x,a)') message
      stop
      
      end
      subroutine sqrtma(mainp,maout,ndofn)
c*****************************************************************************
c
c     Square root of matrix. In the case NDOFN = 3, it is assumed that this
c     matrix has the form diag(A,a), where A is a 2 x 2 matrix.
c
c*****************************************************************************
      implicit none
      integer  ndofn,i
      real*8   mainp(3,3), maout(3,3)

      do i = 1,2
        mainp(i,i)=abs(mainp(i,i))
      end do
      call sqrtm2(mainp,maout)
      if(ndofn.eq.3)
     .  maout(3,3) = sqrt(abs(mainp(3,3)))

      end
      subroutine sqrtm2(mainp,maout)
c*****************************************************************************
c
c     Square root of a 2 x 2 matrix (from a 3 x 3 matrix)
c
c*****************************************************************************
      implicit none
      real*8   mainp(3,3), maout(3,3)
      real*8   a,b,c,d,aux1,aux2,vap1,vap2,det,sq1,sq2

      a = mainp(1,1)
      b = mainp(1,2)
      c = mainp(2,1)
      d = mainp(2,2)
      aux1 =  0.5*(a+d)
      aux2 = 0.25*(a-d)*(a-d) + b*c
      if(aux2.lt.0.0) then
        call runend('SQRTMA: Non real eigenvalue in A')
      else if(aux2.lt.1.0e-10) then                         ! b or c = 0, a = d
        maout(1,1) = sqrt(a) 
        maout(1,2) = 0.0
        maout(2,1) = 0.0
        maout(2,2) = sqrt(d)
        if(abs(b).gt.1.0e-10) then
          maout(1,2) = b/(sqrt(a) + sqrt(d))
        else if(abs(c).gt.1.0e-10) then
          maout(2,1) = c/(sqrt(a) + sqrt(d))
        end if
      else
        vap1 = aux1 + sqrt(aux2)                            ! vep1 = ( b ,vap1 -a )
        vap2 = aux1 - sqrt(aux2)                            ! vep2 = ( vap2 -d, c )
        if(abs(b)+abs(vap1-a).lt.1.0e-15) then
          sq1  = vap1
          vap1 = vap2
          vap2 = sq1
        end if
        sq1 = sqrt(vap1)
        sq2 = sqrt(vap2)
        vap1 = vap1 - a  
        vap2 = vap2 - d  
        det = b*c - vap1*vap2
        maout(1,1) = (b*c*sq1-vap1*vap2*sq2)/det
        maout(1,2) =  b*vap2*(-sq1+sq2)/det
        maout(2,1) =  c*vap1*( sq1-sq2)/det
        maout(2,2) =(-vap1*vap2*sq1+b*c*sq2)/det
      end if
      
      end
      subroutine galcon(
     .  nnode,ndofn,nevab,dvolu,shape,derxy,
     .  amate,rhslo)
c*****************************************************************************
c
c     Galerkin contribution to the system matrix and RHS
c
c*****************************************************************************
      implicit none
      integer
     .  nnode,ndofn,nevab,ievab,inode,idofn,jevab,jnode,
     .  jdofn,i,j
      real*8
     .  dvolu,prod1,prod2,prod3
      real*8
     .  shape(nnode),derxy(2,nnode),amate(nevab,nevab),
     .  rhslo(nevab),difma(3,3,2,2),conma(3,3,2),reama(3,3),
     .  force(3)
      common/proper/difma,conma,reama,force
      
      ievab=0
      do inode=1,nnode
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nnode
            do jdofn=1,ndofn
              jevab=jevab+1
              prod1=0.0
              do i=1,2
                do j=1,2
                  prod1=prod1+derxy(i,inode)
     .              *difma(idofn,jdofn,i,j)*derxy(j,jnode)
                end do
              end do
              prod2=0.0
              do i=1,2
                prod2=prod2+shape(inode)
     .            *conma(idofn,jdofn,i)*derxy(i,jnode)
              end do
              prod3=shape(inode)*reama(idofn,jdofn)*shape(jnode)
              amate(ievab,jevab)=amate(ievab,jevab)
     .          +(prod1+prod2+prod3)*dvolu
            end do
          end do
          rhslo(ievab)
     .      = rhslo(ievab)+shape(inode)*force(idofn)*dvolu
        end do
      end do

      end
      subroutine stacon(
     .  nnode,ndofn,nevab,dvolu,pertu,workm,shape,
     .  derxy,hesxy,resid,tauma,amate,rhslo,kstab)
c*****************************************************************************
c
c     Contribution to the system matrix and RHS from the stabilization term
c
c*****************************************************************************
      implicit none
      integer
     .  ievab,inode,idofn,jdofn,jevab,jnode,k,l
      integer
     .  nnode,ndofn,nevab,kstab,n_ini
      real*8
     .  dvolu,prod1,prod2,prod3,v_ini
      real*8
     .  pertu(nevab,ndofn), workm(2,2), shape(nnode), 
     .  derxy(2,nnode), hesxy(3,nnode), resid(ndofn,nevab),
     .  tauma(3,3), amate(nevab,nevab), rhslo(nevab),
     .  difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3)
      common/proper/difma,conma,reama,force
      
      ievab = 0
      n_ini = ndofn*nevab
      v_ini = 0.0
      call initia(pertu,n_ini,v_ini)
      do inode=1,nnode
        workm(1,1)=hesxy(1,inode)
        workm(2,2)=hesxy(3,inode)
        workm(1,2)=hesxy(2,inode)
        workm(2,1)=workm(1,2)
        do idofn=1,ndofn
          ievab=ievab+1
          do jdofn=1,ndofn
            prod1=reama(jdofn,idofn)*shape(inode)
            prod2=0.0
            do k=1,2                                        
              prod2=prod2
     .          +conma(jdofn,idofn,k)*derxy(k,inode)
            end do
            prod3=0.0
            do k=1,2
              do l=1,2
                prod3=prod3
     .            +difma(jdofn,idofn,k,l)*workm(k,l)
              end do
            end do
            resid(jdofn,ievab)=prod1+prod2-prod3
            call pertur(
     .        kstab,idofn,jdofn,workm,derxy(1,inode),
     .        shape(inode),pertu(ievab,jdofn))
          end do
        end do
      end do
      ievab=0
      do inode=1,nnode
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nnode
            do jdofn=1,ndofn
              jevab=jevab+1
              prod1=0.0
              do k=1,ndofn
                do l=1,ndofn
                  prod1=prod1
     .              +pertu(ievab,k)*tauma(k,l)*resid(l,jevab)
                end do
              end do
              amate(ievab,jevab)=amate(ievab,jevab)
     .          +prod1*dvolu
            end do
          end do
          prod1=0.0
          do k=1,ndofn
            do l=1,ndofn
              prod1=prod1+pertu(ievab,k)*tauma(k,l)*force(l)
            end do
          end do
          rhslo(ievab)=rhslo(ievab)+prod1*dvolu
        end do
      end do

      end
      subroutine invmtx(a,b,deter,nsize)
c***********************************************************************
c
c**** This routine inverts a square matrix A -> Mat(nsize,nsize). The
c**** inverse is stored in B. Its determinant is DETER
c
c***********************************************************************
      implicit none
      integer  nsize
      real*8   deter,t1,t2,t3,denom
      real*8   a(3,3), b(3,3)
c
c*** Invers of a 1*1 matrix
c
      if(nsize.eq.1) then
        deter=a(1,1)
        if(deter.eq.0.0) return
        b(1,1) = 1.0/a(1,1)
        return
      endif
c
c*** Invers of a 2*2 matrix
c
      if(nsize.eq.2) then
        deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
        if(deter.eq.0.) return
        denom=1.0/deter
        b(1,1) = a(2,2)*denom
        b(2,2) = a(1,1)*denom
        b(2,1) =-a(2,1)*denom
        b(1,2) =-a(1,2)*denom
        return
      endif
c
c*** Inverse of a 3*3 matrix
c
      if(nsize.eq.3) then
        t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
        t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
        t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
        deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
        if(deter.eq.0.0) return
        denom = 1./deter
        b(1,1) = t1*denom
        b(2,1) = t2*denom
        b(3,1) = t3*denom                                   
        b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
        b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
        b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
        b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
        b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
        b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom
        return
      endif

      end
      subroutine plamat(young,poiss,thick,difma,conma,reama)
c***********************************************************************
c
c**** This routine computes the coefficient matrices for plates
c
c***********************************************************************
      implicit none
      real*8
     .  young,poiss,thick,difma(3,3,2,2),
     .  conma(3,3,2),reama(3,3)
      real*8  k1,k2,epsi,v_ini
      
      k1 = young/(24.0*(1.0+poiss))
      k2 = young/(24.0*(1.0-poiss))
      epsi = (thick*thick*2.0*(1.0 + poiss)/young)*(6.0/5.0)
      epsi = 1.0/epsi
      
      call initia(difma,36,v_ini)
      call initia(conma,18,v_ini)
      call initia(reama, 9,v_ini)
        
      difma(1,1,1,1) = k1 + k2 
      difma(2,2,1,1) = k1
      difma(3,3,1,1) = epsi
      difma(1,2,1,2) = 0.5*k2
      difma(2,1,1,2) = 0.5*k2
      difma(1,1,2,2) = k1
      difma(2,2,2,2) = k1 + k2
      difma(3,3,2,2) = epsi
      conma(1,3,1) = -epsi
      conma(3,1,1) =  epsi
      conma(2,3,2) = -epsi
      conma(3,2,2) =  epsi
      reama(1,1) = epsi
      reama(2,2) = epsi
      
      difma(1,2,2,1)=difma(2,1,1,2)
      difma(2,1,2,1)=difma(1,2,1,2)

      end
      
