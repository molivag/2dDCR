      subroutine buildi(lpntn,iprob,value,nvalu,idest)
c**********************************************************************
c                 
c**** This routine constructs the matrix for problem IPROB when direct
c**** methods are used. This routine assumes that only one value has
c**** to be assembled.
c                 
c**********************************************************************
      implicit  none
      include  'MatMan.h'
      integer   iprob,itotv,jtotv,nvalu
      integer   lpntn(*),idest(2)
      integer
     .  p_maf(mprob),p_map(mprob),p_lpo(mprob),
     .  p_apc(mprob),p_rof(mprob)                           ! pointer
      real*8    value
      save      p_maf,p_map,p_lpo,p_apc,p_rof 
c
c*** Case of prescribed or condensed degree of freedom
c
      itotv=idest(1)
      jtotv=idest(2)
      if((lpntn(itotv).le.0).or.(lpntn(jtotv).eq.0)) return
c
c*** Initializations
c
      if(kmadi(iprob).gt.0) then
        call spoias(%val(ialpo+(iprob-1)*lengp),p_lpo(iprob))
        call spoias(%val(iamat+(iprob-1)*lengp),p_maf(iprob))
        call svecze(kmadi(iprob),%val(p_maf(iprob)))
        call spoias(%val(iapma+(iprob-1)*lengp),p_map(iprob))
        call spoias(%val(iapii+(iprob-1)*lengp),p_rof(iprob))
        call spoias(%val(iapij+(iprob-1)*lengp),p_apc(iprob))
        if(p_map(iprob).ne.0)
     .    call svecze(nzecp(iprob),%val(p_map(iprob)))
        kmadi(iprob)=-kmadi(iprob)
      end if
c
c*** Coefficient of matrix A_p
c        
      if(lpntn(jtotv).lt.0) then
        if(p_map(iprob).ne.0) then
          call macoam(value,itotv,jtotv,nzecp(iprob),
     .      %val(p_rof(iprob)),%val(p_apc(iprob)),
     .      %val(p_map(iprob)))
        end if
c
c*** Coefficient of the system matrix
c
      else if(lpntn(jtotv).gt.0) then
c
c*** Band matrix storage
c
        if(kmsip(iprob).eq.0) then
          call macoba(value,lpntn(itotv),lpntn(jtotv),nband(iprob),      
     .                %val(p_maf(iprob)),ksymm(iprob))
c
c*** Profile matrix storage
c
        else if(kmsip(iprob).eq.1) then
          call macopr(value,lpntn(itotv),lpntn(jtotv),
     .      %val(p_maf(iprob)),%val(p_lpo(iprob)),
     .      nequa(iprob),ksymm(iprob))
c
c*** Sparse matrix storage
c
        else if(kmsip(iprob).eq.2) then
          call macoba(value,lpntn(itotv),lpntn(jtotv),
     .      nband(iprob),%val(p_maf(iprob)),ksymm(iprob))
c
c*** If a real sparse matrix storage is to be used (thus complicating the
c*** process of the assembly) the routine that must be called to construct
c*** the matrix is MACOSP. The call is as follows
c     
c         call macosp(value,lpntn(itotv),lpntn(jtotv),%val(p_maf),
c    .                %val(incol(iprob)),%val(inrow(iprob)),nonze(iprob))
        end if
      end if
      
      end
      subroutine builma(iprob,value,nvalu,idest,nelem,ndofn,
     .                  nnode,nnodw)
c**********************************************************************
c                 
c**** This routine constructs the matrix for problem IPROB 
c                 
c**********************************************************************
      implicit  none
      include  'MatMan.h'
      integer   iprob,nvalu,nevaw,nelem,ndofn,nnode,nnodw
      integer   idest(*)
      integer   p_mat,p_lpn                                 ! pointer
      real*8    value(*)

      if(sofam(iprob).eq.1) then
        call spoias(%val(ialpn+(iprob-1)*lengp),p_lpn)
        call buildi(%val(p_lpn),iprob,value,nvalu,idest)    ! Direct methods
      else if(sofam(iprob).eq.2) then
        call spoias(%val(iamat+(iprob-1)*lengp),p_mat)
        nevaw=ndofn*nnode+nnodw-nnode
        call stomat(value,nvalu,nevaw,idest,nelem,          ! EE Iterative methods
     .              ksymm(iprob),p_mat,lengp,kmadi(iprob))
      else if(sofam(iprob).eq.3) then
        call spoias(%val(ialpn+(iprob-1)*lengp),p_lpn)
        call spamat(%val(p_lpn),iprob,value,nvalu,idest)    ! SK Iterative methods
      end if

      end
      subroutine macoam(value,itotv,jtotv,nzeco,roffs,amcol,amatr)
c************************************************************************
c                 
c**** This routine assembles a scalar VALUE to the proper position of 
c**** a sparse matrix A.
c                 
c************************************************************************
      implicit  none
      integer   itotv,jtotv,izeco,jcomp,nzeco
      integer   amcol(nzeco), roffs(*)
      real*8    value
      real*8    amatr(nzeco)

      izeco=roffs(itotv)
      jcomp=amcol(izeco)
      do while(jcomp.ne.jtotv)
        izeco=izeco+1
        jcomp=amcol(izeco)
      end do

      amatr(izeco)=amatr(izeco)+value
      
      end
      subroutine macoba(value,itotv,jtotv,nband,xmatr,ksymm)
c************************************************************************
c                 
c**** This routine assembles a scalar VALUE to the proper position of 
c**** a banded matrix XMATR
c                 
c************************************************************************
      implicit   real*8(a-h,o-z)
      real*8     xmatr(nband*(ksymm+1)+1,*)

      if(itotv.ge.jtotv) then                       ! Lower matrix
        xmatr(ksymm*nband+1+itotv-jtotv,jtotv)=
     .    xmatr(ksymm*nband+1+itotv-jtotv,jtotv)
     .    +value
      else if(ksymm.eq.1) then                      ! Upper matrix
        xmatr(nband+1+itotv-jtotv,itotv)=
     .    xmatr(nband+1+itotv-jtotv,itotv)+value
      end if

      end
      subroutine macopr(value,idest,jdest,xmatr,lpont,nequa,ksymm)
c************************************************************************
c                 
c**** This routine assembles a scalar VALUE to the proper position of 
c**** a skyline matrix XMATR
c                 
c************************************************************************
      implicit   real*8(a-h,o-z)
      integer    lpont(*)
      real*8     xmatr(*)              ! xmatr(nequa+lpont(nequa)*(1+ksymm))

      if (idest.eq.jdest) then                              ! Diagonal terms
        xmatr(idest)=xmatr(idest)+value
      else if (idest.lt.jdest) then                         ! Upper matrix
        jloca=idest+lpont(jdest)-jdest+1+nequa
        xmatr(jloca)=xmatr(jloca)+value
      else if (ksymm.eq.1) then                             ! Lower matrix
        jloca=jdest+lpont(idest)-idest+1+nequa+lpont(nequa)
        xmatr(jloca)=xmatr(jloca)+value
      end if

      end
      subroutine macosp(value,idest,jdest,xmatr,lncol,lnrow,nonze)
c************************************************************************
c                 
c**** This routine assembles a scalar VALUE to the proper position of 
c**** a skyline matrix XMATR
c                 
c************************************************************************
      implicit   real*8(a-h,o-z)
      integer    lncol(*), lnrow(*)
      real*8     xmatr(*)
      logical*2 newva

      newva=.true.
      do iasse=1,nonze
        if (idest.eq.lnrow(iasse)) then
          if (jdest.eq.lncol(iasse)) then
            xmatr(iasse)=xmatr(iasse)+value
            newva=.false.
          end if
        end if
      end do
      if (newva) then
        nonze=nonze+1
        xmatr(nonze)=value
        lnrow(nonze)=idest
        lncol(nonze)=jdest
      end if
      
      end
      subroutine putmat(amatr,elmat,nvalo,nevaw,ksymm)
c*****************************************************************************
c
c**** This routine adds the elemental matrices for iterative solvers
c
c*****************************************************************************
      implicit none
      integer   nvalo,nevaw,ksymm,ievaw,jevaw,icolu
      real*8    amatr(nevaw*nevaw), elmat(nvalo,nvalo)

      icolu=0
      if (ksymm.eq.1) then                                  ! Non symmetric
        do jevaw=1,nevaw
          do ievaw=1,nevaw
            amatr(ievaw+icolu)=amatr(ievaw+icolu)+elmat(ievaw,jevaw)
          end do
          icolu=icolu+nevaw
        end do
      else                                                  ! Symmetric
        do jevaw=1,nevaw
          do ievaw=jevaw,nevaw
            amatr(ievaw+icolu)=amatr(ievaw+icolu)+elmat(ievaw,jevaw)
          end do
          icolu=icolu+nevaw-jevaw
        end do
      end if

      end
      subroutine spamat(lpntn,iprob,value,nvalu,idest)
c**********************************************************************
c                 
c**** This routine constructs the matrix for problem IPROB when it is
c**** a sparse matrix and an iterative solver is used (SOFAM = 3). This
c**** routine assumes that only one value has to be assembled.
c                 
c**********************************************************************
      implicit  none
      include  'MatMan.h'
      integer   iprob,itotv,jtotv,nvalu,itotn,jtotn
      integer   lpntn(*),idest(2)
      integer   p_mat,p_apc,p_rof                           ! pointer
      real*8    value
c
c*** Initializations
c
      itotv=idest(1)
      jtotv=idest(2)
      if(kmadi(iprob).gt.0) then
        call spoias(%val(iamat+(iprob-1)*lengp),p_mat)
        call svecze(nzecf(iprob),%val(p_mat))
        call spoias(%val(iapma+(iprob-1)*lengp),p_mat)
        if(p_mat.ne.0) call svecze(nzecp(iprob),%val(p_mat))
        kmadi(iprob)=-kmadi(iprob)
      end if
c
c*** Case of prescribed or condensed degree of freedom
c
      if((lpntn(itotv).le.0).or.(lpntn(jtotv).eq.0)) then
        return
c
c*** Coefficient of matrix A_p
c        
      else if(lpntn(jtotv).lt.0) then
        call spoias(%val(iapma+(iprob-1)*lengp),p_mat)
        if(p_mat.ne.0) then
          call spoias(%val(iapii+(iprob-1)*lengp),p_rof)
          call spoias(%val(iapij+(iprob-1)*lengp),p_apc)
          call macoam(value,itotv,jtotv,nzecp(iprob),
     .                %val(p_rof),%val(p_apc),%val(p_mat))
        end if
c
c*** Coefficient of the system matrix
c
      else if(lpntn(jtotv).gt.0) then
        itotn=lpntn(itotv)
        jtotn=lpntn(jtotv)
        call spoias(%val(iamat+(iprob-1)*lengp),p_mat)
        call spoias(%val(iafii+(iprob-1)*lengp),p_rof)
        call spoias(%val(iafij+(iprob-1)*lengp),p_apc)
        call macoam(value,itotn,jtotn,nzecf(iprob),
     .              %val(p_rof),%val(p_apc),%val(p_mat))
      end if
      
      end
      subroutine stomat(elmat,nvalo,nevaw,ielem,nelem,ksymm,
     .                  iamat,lengp,kmadi)
c*****************************************************************************
c
c**** This routine stores the elemental matrices for iterative solvers
c
c*****************************************************************************
      implicit none
      integer   nvalo,nevaw,ksymm,ielem,nelem,lengp,
     .          kmadi,ndimt,jelem
      integer   iamat,p_mat                                 ! pointer
      real*8    elmat(nvalo,nvalo)
c
c*** If needed, clear the space for the element matrices
c
      if(kmadi.eq.1) then
        ndimt=(1+ksymm)*(nevaw*nevaw-nevaw)/2+nevaw           
        do jelem=1,nelem                            
          call spoias(%val(iamat+(jelem-1)*lengp),p_mat)
          call svecze(ndimt,%val(p_mat))
        end do
        kmadi=-1
      end if
c
c*** Store the matrix of element IELEM
c
      call spoias(%val(iamat+(ielem-1)*lengp),p_mat)
      call putmat(%val(p_mat),elmat,nvalo,nevaw,ksymm)

      end



      subroutine asbddd(nequa,nsban,amatr,cmatr,nsist,bvect,
     .                  xvect,iflag,lunit)
c****************************************************************************
c
c*** This routine solves a linear system  A x = b using the Cholesky
c*** decomposition and a forward and a backward substitution. 
c
c    nequa: number of equations
c    nsban: half bandwidth
c    amatr: matrix of the system, assumed to be symmetric and  positive-
c           definite, stored by rows (amatr(nsban+1,nequa))
c    cmatr: working space (decomposition of amatr, cmatr(nsban+1,nequa)
c    nsist: number of different right-hand-sides
c    bvect: right-hand-side
c    xvect: unknown
c    iflag: decomposition is needed (=1) or not (=0)
c    lunit: logical unit to print errors
c
c***************************************************************************
      implicit real*8(a-h,o-z)
      real*8  amatr(nsban+1,nequa),bvect(nequa,nsist),
     .        cmatr(nsban+1,nequa),xvect(nequa,nsist)

      if(iflag.eq.1)
     .  call chobde(nequa,nsban,amatr,cmatr,lunit)           !Cholesky decom.
      call fosbsu(nequa,nsban,cmatr,nsist,bvect,xvect,lunit) !forward subst.
      call basbsu(nequa,nsban,cmatr,nsist,xvect,xvect,lunit) !backward subst.

      end
      subroutine aubidd(nequa,nsban,amatr,cmatr,nsist,bvect,
     .                  xvect,iflag,lunit)
c****************************************************************************
c
c*** This routine solves a linear system  A x = b using the Crout
c*** decomposition and a forward and a backward substitution. 
c
c    nequa: number of equations
c    nsban: half bandwidth
c    amatr: matrix of the system, stored by rows (amatr(2*nsban+1,nequa))
c    cmatr: working space (decomposition of amatr, cmatr(2*nsban+1,nequa)
c    nsist: number of different right-hand-sides
c    bvect: right-hand-side
c    xvect: unknown
c    iflag: decomposition is needed (=1) or not (=0)
c    lunit: logical unit to print errors
c
c***************************************************************************
      implicit real*8(a-h,o-z)
      real*8  amatr(2*nsban+1,nequa),bvect(nequa,nsist),
     .        cmatr(2*nsban+1,nequa),xvect(nequa,nsist)

      if(iflag.eq.1)
     .  call crobde(nequa,nsban,amatr,cmatr,lunit)           !Crout decom.
      call foubsu(nequa,nsban,cmatr,nsist,bvect,xvect,lunit) !forward subst.
      call baubsu(nequa,nsban,cmatr,nsist,xvect,xvect,lunit) !backward subst.

      end

      subroutine basbsu(n,l,c,m,z,x,iu)
c****************************************************************************
c
c*** This routine performs a backward substitution to solve the system `Cx=z',
c*** `n' being the number of equations, `l' half the bandwidth,`m' the number
c*** of substitutions to be performed and `iu' the logical unit where errors
c*** have to be printed
c
c****************************************************************************
      implicit real*8(a-h,o-z)
      real*8  c(l+1,n),z(n,m),x(n,m)

      do is = 1,m
        x(n,is) = z(n,is)/c(1,n)
        do i = n-1,1,-1
          x(i,is) = z(i,is)
          k2 = i+l
          if (k2.gt.n) k2 = n
          if (i+1.le.k2) then 
            do k = i+1,k2
              x(i,is) = x(i,is)-c(k-i+1,i)*x(k,is)
            end do
          end if
          if (c(1,i).le.0.) 
     .    call srunen('BASBSU: NON-POSITIVE DEFINITE MATRIX',iu)
          x(i,is) = x(i,is)/c(1,i)
        end do
      end do
 
      return
      end
      subroutine batosp(nequa,nband,banma,spama,lncol,lnrow,
     .                  nonze,multi,ksymm,bvect)
c****************************************************************************
c
c*** This routine stores a banded matrix BANMA in sparse form SPAMA
c
c***************************************************************************
      implicit real*8(a-h,o-z)
      real*8    banma(nband*(1+ksymm)+1,nequa), spama(multi*nonze),
     .          bvect(nequa)
      integer   lncol(multi*nonze), lnrow(multi*nonze)
      
      iasse=0
      do iband=ksymm*nband+1,ksymm*nband+nband+1
        do iequa=1,nequa
          if (abs(banma(iband,iequa)).gt.1.0d-12) then
            jdest=iequa
            idest=iband-ksymm*nband-1+jdest
            iasse=iasse+1
            spama(iasse)=banma(iband,iequa)
            lncol(iasse)=jdest
            lnrow(iasse)=idest
            if (ksymm.eq.0) then
              if (idest.ne.jdest) then
                iasse=iasse+1
                spama(iasse)=banma(iband,iequa)
                lncol(iasse)=idest
                lnrow(iasse)=jdest
              end if
            end if
           end if 
         end do
       end do
       if (ksymm.eq.1) then
         do iband=1,nband
           do iequa=1,nequa
             if (abs(banma(iband,iequa)).gt.1.0d-12) then
               idest=iequa
               jdest=-iband+nband+1+idest
               iasse=iasse+1
               spama(iasse)=banma(iband,iequa)
               lncol(iasse)=jdest
               lnrow(iasse)=idest
             end if
           end do
         end do
       end if

       do iasse=nonze+1,multi*nonze
         spama(iasse)=0.0d0
         lncol(iasse)=0
         lnrow(iasse)=0
       end do

      end
      subroutine baubsu(n,l,c,m,z,x,iu)
c****************************************************************************
c
c*** This routine performs a backward substitution to solve the system `Cx=z',
c*** `n' being the number of equations, `l' half the bandwidth,`m' the number
c*** of substitutions to be performed and `iu' the logical unit where errors
c*** have to be printed. 'C' is an upper triangular matrix with C(i,i)=1.
c
c****************************************************************************
      implicit real*8(a-h,o-z)
      real*8  c(2*l+1,n),z(n,m),x(n,m)
      
      id=l+1
      do is = 1,m
        x(n,is) = z(n,is)
        do i = n-1,1,-1
          x(i,is) = z(i,is)
          k2 = i+l
          if (k2.gt.n) k2 = n
          do k = i+1,k2
            x(i,is) = x(i,is)-c(id+i-k,i)*x(k,is)
          end do
        end do
      end do

      end
      subroutine canonz(nequa,nband,banma,nonze,ksymm)
c****************************************************************************
c
c*** This routine calculates the number of non-zero values in BANMA
c
c***************************************************************************
      implicit real*8(a-h,o-z)
      real*8    banma(nband*(1+ksymm)+1,nequa)

      nonze=0
      if (ksymm.eq.0) then
        do iequa=1,nequa
          do iband=2,nband+1
            if (abs(banma(iband,iequa)).gt.1.0d-12) nonze=nonze+2
          end do
          if (abs(banma(1,iequa)).gt.1.0d-12) nonze=nonze+1
        end do
      else
        do iband=1,2*nband+1
          do iequa=1,nequa
            if (abs(banma(iband,iequa)).gt.1.0d-12) nonze=nonze+1
          end do
        end do
      end if

      end
      subroutine chobde(n,l,a,c,iu)
c****************************************************************************
c
c*** This routine performs the cholesky decomposition of a matrix `A'
c*** of bandwidth `l' and `n' rows. The resulting matrix is  `C'. Errors
c*** are printed in unit `iu'
c
c****************************************************************************
      implicit real*8(a-h,o-z)
      real*8  a(l+1,n),c(l+1,n)

      do i = 1,n
        c(1,i) = a(1,i)
        k1 = i-l
        if (k1.lt.1) k1 = 1
        if (k1.le.i-1) then 
          do k = k1,i-1
            c(1,i) = c(1,i) - c(i-k+1,k)*c(i-k+1,k)
          end do
        end if
        if (c(1,i).le.0.) 
     .  call srunen('CHOBDE: NON-POSITIVE DEFINITE MATRIX',iu)
        c(1,i) = sqrt(c(1,i))
        j2 = i + l
        if (j2.gt.n) j2 = n
        if (i+1.le.j2) then 
          do j = i+1,j2
            c(j-i+1,i) = a(j-i+1,i)
            k1 = j - l
            if (k1.lt.1) k1 = 1
            if (k1.le.i-1) then 
              do k = k1,i-1
                c(j-i+1,i) = c(j-i+1,i) - c(i-k+1,k)*c(j-k+1,k)
              end do
            end if
            c(j-i+1,i) = c(j-i+1,i)/c(1,i)
          end do
        end if
      end do

      end
      subroutine crobde(n,l,a,c,iu)
c****************************************************************************
c
c*** This routine performs the Crout descomposition of a matrix `A'
c*** of bandwidth `l' and `n' rows. The resulting matrix is  `C'. Errors
c*** are printed in unit `iu'
c
c****************************************************************************
      implicit real*8(a-h,o-z)
      real*8  a(2*l+1,n),c(2*l+1,n)

      nb=l+1
      id=l+1
      do m = 1,n
        rleup=0.0
        l1=1
        if (m-l.gt.1) l1=m-l
        do lj=l1,m-1
          rleup=rleup+c(id+m-lj,lj)*c(id+lj-m,lj)
        end do
        c(id,m)=a(id,m)-rleup
        if (c(id,m).eq.0.0)
     .    call srunen('CROBDE: PIVOTING IS NEEDED',iu)
        do i=m+1,nb+m-1
          rleup=0.0
          l1=1
          if (i-l.gt.1) l1=i-l
          do lj=l1,m-1
            rleup=rleup+c(id+i-lj,lj)*c(id+lj-m,lj)
          end do
          c(id+i-m,m)=a(id+i-m,m)-rleup
        end do
        do j=m+1,nb+m-1
          rleup=0.0
          l1=1
          if(j-l.gt.1) l1=j-l
          do lj=l1,m-1
            rleup=rleup+c(id+m-lj,lj)*c(id+lj-j,lj)
          end do
          c(id+m-j,m)=1/c(id,m)*(a(id+m-j,m)-rleup)
        end do
      end do

      end





      subroutine dirsol(
     .  rhsid,unkno,ntotv,nsist,kfact,iprob,outso)
c****************************************************************************
c
c**** This routine drives the library to solve a linear system of equations
c**** by direct methods using the LDU decomposition. This decomposition is
c**** done only when KFACT = 1. Otherwise, only a backward and a forward
c**** substituion is performed.       
c
c****************************************************************************
      implicit none
      include   'MatMan.h'
      integer    ntotv,nsist,iprob,multi,ifail,kfact,
     .           iflag
      real*8     rhsid(ntotv,nsist), unkno(ntotv,nsist)
      integer    p_mat,iwoso,imasp,p_col,p_row,p_lpn,       ! pointer
     .           p_lpo,p_rof,p_apc,p_rhs
      logical*1  outso
      integer    lbyta,lbyts,svomt
      data       lbyta,svomt/0,0/
      save       lbyta,iwoso,svomt
c
c*** Allocate a work space for the redefinition of RHSID and UNKNO
c
      lbyts=8*ntotv*nsist
      if (lbyta.eq.0) then
        call sadmem(0,iwoso,lbyts,svomt)
        lbyta=lbyts
      else
        if (lbyts.gt.lbyta) then
          svome(iprob)=max(svome(iprob),svomt)        
          call sadmem(2,iwoso,lbyta,svomt)
          call sadmem(0,iwoso,lbyts,svomt)
          lbyta=lbyts
        end if
      end if
      call spoias(%val(ialpn+(iprob-1)*lengp),p_lpn)
c
c*** Compute A_p x_p, where x_p are the prescribed values of the unknown
c
      iflag=1
      call spoias(%val(iapma+(iprob-1)*lengp),p_mat)
      if(p_mat.eq.0) iflag=0
      call spoias(%val(irhsx+(iprob-1)*lengp),p_rhs)
      if((p_mat.ne.0).and.(kfact.le.1)) then
        call spoias(%val(iapii+(iprob-1)*lengp),p_rof)
        call spoias(%val(iapij+(iprob-1)*lengp),p_apc)
        call modrhs(
     .    %val(p_lpn),%val(p_rof),%val(p_apc),
     .    unkno,%val(p_rhs),%val(p_mat),ntotv,
     .    nzecp(iprob),nsist,nequa(iprob),iflag,
     .    zeroc)
      end if
c
c*** Redefine RHSID according to the renumbering strategy
c
      if(iflag.eq.0) p_rhs=iwoso                            ! not used 
      call redefr(
     .  %val(iwoso),%val(p_lpn),%val(p_rhs),rhsid,ntotv,
     .  nequa(iprob),nsist,iflag)
c
c*** Solve the final algebraic system using a band storage
c
      call spoias(%val(iamat+(iprob-1)*lengp),p_mat)
      if(kmsip(iprob).eq.0) then
        if(ksymm(iprob).eq.0) then
          if(kites(iprob).eq.11) then
            call asbddd(nequa(iprob),nband(iprob),%val(p_mat),
     .        %val(p_mat),nsist,%val(iwoso),%val(iwoso),
     .        kfact,lusol)
          else if(kites(iprob).eq.12) then
            call srunen(
     .        'DIRSOL: PIVOTING NOT YET IMPLEMENTED',lusol)
          end if            
        else if(ksymm(iprob).eq.1) then
          if(kites(iprob).eq.11) then
            call aubidd(nequa(iprob),nband(iprob),%val(p_mat),
     .        %val(p_mat),nsist,%val(iwoso),%val(iwoso),
     .        kfact,lusol)
          else if(kites(iprob).eq.12) then
            call srunen(
     .        'DIRSOL: PIVOTING NOT YET IMPLEMENTED',lusol)
          end if            
        end if
c
c*** Solve the final algebraic system using a profile storage
c
      else if (kmsip(iprob).eq.1) then
        call spoias(%val(ialpo+(iprob-1)*lengp),p_lpo)
        if(kites(iprob).eq.11) then
          call skylin(nequa(iprob),%val(p_mat),%val(iwoso),
     .      %val(p_lpo),nsist,ksymm(iprob),kfact,
     .      lusol)
        else if(kites(iprob).eq.12) then
          call srunen('DIRSOL: PIVOTING NOT YET IMPLEMENTED',lusol)
        end if          
c
c*** Solve the final algebraic system using a sparse storage
c
      else if (kmsip(iprob).eq.2) then
        if(kites(iprob).eq.11) then
          call srunen(
     .      'DIRSOL: ONLY THE PIVOTING VERSION AVAILABLE',lusol)
        else if(kites(iprob).eq.12) then
          call canonz(nequa(iprob),nband(iprob),
     .      %val(p_mat),
     .      nonze(iprob),ksymm(iprob))
          multi=nmult(iprob)
          ifail=5
          do while (ifail.eq.5)
            call sadmem(0,imasp,multi*nonze(iprob)*8    ,svomt)         
            call sadmem(0,p_col,multi*nonze(iprob)*lengp,svomt)
            call spoias(p_col,%val(incol+(iprob-1)*lengp))
            call sadmem(0,p_row,multi*nonze(iprob)*lengp,svomt)
            call spoias(p_row,%val(inrow+(iprob-1)*lengp))
            call batosp(nequa(iprob),nband(iprob),
     .        %val(p_mat),%val(imasp),
     .        %val(p_col),%val(p_row),
     .        nonze(iprob),multi,ksymm(iprob),%val(iwoso))
            call sparse(
     .        nequa(iprob),nonze(iprob),multi,%val(imasp),
     .        %val(p_col),%val(p_row),%val(iwoso),
     .        nsist,ifail,kfact,lusol,svomt,svome(iprob),
     .        lengp)
            svome(iprob)=max(svome(iprob),svomt)        
            call sadmem(2,imasp,multi*nonze(iprob)*8    ,svomt) ! Deal. temp.
            call sadmem(2,p_col,multi*nonze(iprob)*lengp,svomt) ! memory
            call sadmem(2,p_row,multi*nonze(iprob)*lengp,svomt)
            if (ifail.ne.0)
     .        call redefr(%val(iwoso),%val(p_lpn),%val(p_rhs),
     .        rhsid,ntotv,nequa(iprob),nsist,iflag)
            multi=multi+1
          end do
          nmult(iprob)=multi-1
        end if
      end if
c
c*** Redefine UNKNO according to the renumbering strategy
c
      call redunk(%val(iwoso),%val(p_lpn),unkno,
     .            ntotv,nequa(iprob),nsist)
c
c*** Write temporary memory requirements 
c
      if(outso) write(lusol,110) svome(iprob)
c
c*** Formats
c
  110 format(11x,'*** VOLATILE MEMORY (BYTES): ',i12)

      end
      subroutine fosbsu(n,l,c,m,b,z,iu)
c****************************************************************************
c
c*** This routine performs a forward substitution to solve the system `Cz=b',
c*** `n' being the number of equations, `l' half the bandwidth, `m' the number
c*** of substitutions to be performed and `iu' the logical unit where errors
c*** have to be printed.
c
c****************************************************************************
      implicit real*8(a-h,o-z)
      real*8 c(l+1,n),b(n,m),z(n,m)

      do is = 1,m
        z(1,is) = b(1,is)/c(1,1)
        do i = 2,n
          z(i,is) = b(i,is)
          k1 = i-l
          if (k1.lt.1) k1 = 1
          if (k1.le.i-1) then 
            do k = k1,i-1
              z(i,is) = z(i,is)-c(i-k+1,k)*z(k,is)
            end do
          end if
          if (c(1,i).le.0.) 
     .    call srunen('FOSBSU: NON-POSITIVE DEFINITE MATRIX',iu)
          z(i,is) = z(i,is)/c(1,i)
        end do
      end do

      end
      subroutine foubsu(n,l,c,m,b,z,iu)
c****************************************************************************
c
c*** This routine performs a forward substitution to solve the system `Cz=b',
c*** `n' being the number of equations, `l' half the bandwidth, `m' the number
c*** of substitutions to be performed and `iu' the logical unit where errors
c*** have to be printed
c
c****************************************************************************
      implicit real*8(a-h,o-z)
      real*8 c(2*l+1,n),b(n,m),z(n,m)

      id=l+1
      do is = 1,m
        z(1,is) = b(1,is)/c(id,1)
        do i = 2,n
          z(i,is) = b(i,is)
          k1 = i-l
          if (k1.lt.1) k1 = 1
          do k = k1,i-1
            z(i,is) = z(i,is)-c(id+i-k,k)*z(k,is)
          end do
          if (c(id,i).eq.0.) 
     .    call srunen('FOUBSU: NON DEFINITE MATRIX',iu)
          z(i,is) = z(i,is)/c(id,i)
        end do
      end do

      end
      subroutine modrhs(lpntn,roffs,apcol,unkno,rhsax,apmat,
     .                  ntotv,nzeco,nsist,nequa,iflag,zeroc)
c****************************************************************************
c
c**** This routine computes the product A_p x_p when a direct solver
c**** is used.       
c
c****************************************************************************
      implicit none
      integer   ntotv,iflag,itotv,itotn,jtotv,izec0,izec1,izeco,
     .          nzeco,isist,nsist,nequa,izecn
      integer   lpntn(ntotv), roffs(ntotv+1), apcol(nzeco)
      logical*1 first
      real*8    xcomp,zeroc
      real*8    apmat(nzeco),       unkno(ntotv,nsist),
     .          rhsax(ntotv,nsist)
c
c*** Initialization
c
      iflag=0
      first=.false.
      call svecze(nsist*ntotv,rhsax)
      do isist=1,nsist
c
c*** Check if x_p is zero or not 
c
        do itotv=1,ntotv
          if (lpntn(itotv).lt.0) then
            itotn=-lpntn(itotv)
            xcomp=abs(unkno(itotv,isist))
            if(xcomp.gt.zeroc) iflag=1
          end if
        end do
c
c*** Compute the product A_p x_p
c
        if(iflag.eq.1) then
          do itotv=1,ntotv
            izec0=roffs(itotv)
            if(izec0.gt.0) then
              if(.not.first) then
                jtotv=itotv+1
                izecn=roffs(jtotv)
                do while(izecn.eq.0)
                  jtotv=jtotv+1
                  izecn=roffs(jtotv)
                end do
                if(izecn.gt.izec0) first=.true.
              end if
              if(first) then
                jtotv=itotv+1
                izec1=roffs(jtotv)-1
                do while(izec1.le.0)
                  jtotv=jtotv+1
                  izec1=roffs(jtotv)-1
                end do
                do izeco=izec0,izec1
                  jtotv=apcol(izeco)
                  rhsax(itotv,isist)=rhsax(itotv,isist)
     .              +apmat(izeco)*unkno(jtotv,isist)
                end do
              end if
            end if
          end do
        end if
      end do

      end
      subroutine redefr(rwork,lpntn,rhsax,rhsid,ntotv,nequa,
     .                  nsist,iflag)
c****************************************************************************
c
c**** This routine redefines RHSID according to the renumbering strategy
c
c****************************************************************************
      implicit none
      integer  ntotv,nsist,iflag,itotv,itotn,isist,nequa
      integer  lpntn(ntotv)
      real*8   rhsax(ntotv,*),     rhsid(ntotv,nsist),
     .         rwork(nequa,nsist)

      do isist=1,nsist
c
c*** External RHS
c
        do itotv=1,ntotv
          if (lpntn(itotv).gt.0)
     .      rwork(lpntn(itotv),isist)=rhsid(itotv,isist)
        end do
c
c*** Contribution from A_p x_p
c
        if(iflag.eq.1) then
          do itotv=1,ntotv
            itotn=lpntn(itotv)
            if(itotn.gt.0)
     .        rwork(itotn,isist)=rwork(itotn,isist)
     .        -rhsax(itotv,isist)
          end do
        end if
      end do
      
      end
      subroutine redunk(rwork,lpntn,unkno,ntotv,nequa,nsist)
c****************************************************************************
c
c**** This routine redefines the vector of unknowns according to the 
c**** renumbering strategy
c
c****************************************************************************
      implicit real*8(a-h,o-z)
      real*8      unkno(ntotv,nsist), rwork(nequa,nsist)
      integer     lpntn(ntotv)

      do isist=1,nsist
        do itotv=1,ntotv
          if (lpntn(itotv).gt.0)
     .      unkno(itotv,isist)=rwork(lpntn(itotv),isist)
        end do
      end do
      
      end
      subroutine skybak(gstdi,gstlo,gstup,lpont,neqns,refor,lures)
c***********************************************************************
c
c**** This routine performs the back-solution of symmetric equations
c     stored in profile form .
c     coefficient matrix must be decomposed into its triangular
c     factors using skylin before using skybak.
c
c.... input parameters
c
c         gstdi(neqns)       -reciprocal of diagonal of triangular 
c                               factor
c         gstlo(lpont(neqns) -lower triangular factor of matrix
c         gstup(lpont(neqns) -upper triangular factor of matrix
c                               (gstup and gstlo have same calling 
c                               address for symmetric matrices)
c         refor(neqns)       -right hand side vector in equations
c         lpont(neqns)       -pointer array to bottom of columns 
c                               of gstlo and gstup
c         neqns              -number of equations to be solved
c
c.... output parameter
c
c         refor(neqns)       -solution of equations
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      integer  lures
      real*8   gstlo(*), gstup(*), gstdi(*), refor(*)
      integer  lpont(*)
      zeros=0.0d0
c
c*** Find the first nonzero entry in the right hand side
c
      do 100 keqns=1,neqns
        ieqns=keqns
        if(refor(ieqns).ne.zeros) go to 200
  100 continue
      if(neqns.gt.0) write(lures,2000)
      return
  200 if(ieqns.lt.neqns) then
c
c*** Reduce the right hand side
c
        do 300 jeqns=ieqns+1,neqns
          jrows=lpont(jeqns-1)
          jheig=lpont(jeqns)-jrows
          if(jheig.gt.0) then
            refor(jeqns)=refor(jeqns)-
     .        svecdo(gstlo(jrows+1),refor(jeqns-jheig),jheig)
          endif
  300   continue
      endif
c
c*** Multiply by inverse of diagonal elements
c
      do 400 jeqns=ieqns,neqns
        refor(jeqns)=refor(jeqns)*gstdi(jeqns)
  400 continue
c
c*** Backsubstitution
c
      if(neqns.gt.1) then
        do 500 jeqns=neqns,2,-1
          jrows=lpont(jeqns-1)
          jheig=lpont(jeqns)-jrows
          if(jheig.gt.0) then
            kterm=jeqns-jheig
            do 600 iterm=1,jheig
              refor(kterm)=refor(kterm)-gstup(jrows+iterm)*refor(jeqns)
              kterm=kterm+1
  600       continue
          endif
  500   continue
      endif
      return
 2000 format(11x,'*** WARNING: ZERO RIGHT-HAND-SIDE VECTOR')
      end
      
      subroutine skycek(gstup,nterm,saval)
c***********************************************************************
c
c*** This routine tests for rank
c
c.... inputs
c
c       gstup(nterm) -  column to of unreduced elements in array
c       nterm        -  number of elements in column
c
c.... outputs
c
c       saval        -  sum of absolute values
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      dimension gstup(*)

      saval=0.
      do 100 iterm=1,nterm
  100 saval=saval+abs(gstup(iterm))

      return
      end
      subroutine skydia(gstlo,gstup,gstdi,nterm,isymm,pivot)
c***********************************************************************
c
c*** This routine diagonalizes element in triangular decomposition
c
c.... input parameters
c         gstup(nterm) - column of upper triangular part of matrix
c         gstlo(nterm) - row of lower triangular part of matrix
c         gstdi(nterm) - reciprocal of diagonals in triangular factors
c         nterm        - number of terms in vectors
c         isymm        - if 0 -> equations are unsymmetric
c                        if 1 -> equations are symmetric
c         pivot        - diagonal in matrix to be factored
c
c.... output parameter
c         pivot        - reduced diagonal of factor
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      dimension gstlo(*),gstup(*),gstdi(*)
c
c      <<<<<<<<<<<<<<<<<<<<< cute trick to avoid convex error
c
      if(isymm.eq.0) then
        do 100 iterm= 1,nterm
        pivot=pivot-gstlo(iterm)*gstup(iterm)*gstdi(iterm)
        gstup(iterm)=gstup(iterm)*gstdi(iterm)
  100   continue
      else if(isymm.eq.1) then
        do 105 iterm= 1,nterm
        pivot=pivot-gstup(iterm)*gstup(iterm)*gstdi(iterm)
        gstup(iterm)=gstup(iterm)*gstdi(iterm)
  105   continue
      endif
c
c***finish computation of column of gstlo for unsymmetric matrices
c
      if(isymm.eq.0) then
        do 200 iterm=1,nterm
  200   gstlo(iterm)=gstlo(iterm)*gstdi(iterm)
      endif

      return
      end
      subroutine skyli0(neqns,gstdi,gstup,gstlo,eload,
     .                  lpont,nsist,ksymm,kresl,lures)
c***********************************************************************
c
c**** This routine solves the set of linear equations 
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      integer   lures
      integer   lpont(*)
      real*8    gstdi(neqns),
     .          gstlo(*),                                   !gstlo(lpont(neqns)),
     .          gstup(*),                                   !gstup(lpont(neqns)),
     .          eload(neqns,nsist)
c
c*** Factorize matrix if neccesary
c
      isymm=0                           ! Non symmetric (ksymm=1)
      if (ksymm.eq.0) isymm=1           ! Symmetric
      if(kresl.eq.1)
     .  call skytri(gstdi,gstlo,gstup,isymm,lpont,
     .  neqns,lures)
c
c*** Solve the equations 
c
      do isist=1,nsist
        call skybak(gstdi,gstlo,gstup,lpont,neqns,
     .    eload(1,isist),lures)
      end do
      
      end
      subroutine skylin(neqns,amatr,bvect,lpont,nsist,ksymm,
     .                  kfact,lures)
c****************************************************************************
c
c*** This routine solves a linear system  A x = b using the Skyline
c*** solver
c
c    neqns: number of equations
c    amatr: matrix of the system (GSTDI, GSTUP, GSTLO) see skygra.f
c    nsist: number of different right-hand-sides
c    bvect: right-hand-side
c
c***************************************************************************
      implicit real*8(a-h,o-z)
      integer   lures,kfact
      integer   lpont(neqns)
      real*8    amatr(*),    bvect(neqns,nsist)

      igstd=1                            ! GSTDI(NEQNS)
      igstu=igstd+neqns                  ! GSTUP(LPONT(NEQNS))
      igstl=igstu+ksymm*lpont(neqns)     ! GSTLO(LPONT(NEQNS))
      call skyli0(neqns,amatr(igstd),amatr(igstu),
     .            amatr(igstl),bvect,lpont,nsist,ksymm,
     .            kfact,lures)

      end
      subroutine skytri(gstdi,gstlo,gstup,isymm,lpont,neqns,lures)
c***********************************************************************
c
c*** This routine performs the triangular decomposition 
c    of a matrix stored in profile form
c
c....input parameters
c
c       gstlo(lpont(neqns)) - lower triangular part of matrix
c       gstup(lpont(neqns)) - upper part of triangular matrix
c       gstdi(neqns)        - diagonals of triangular matrix
c       lpont(neqns)        - pointers to bottom of colums of 
c                             gstlo and gstup arrays
c       neqns               - number of equations to be solved
c       isymm               - =0 equations are unsymmetric
c                             =1 equations are symmetric and calling
c                                 address of gstlo may be same 
c                                 as that for gstup (i.e., gstup 
c                                 and gstlo share same memory)
c
c.... output parameters
c
c       gstlo(lpont(neqns)) - lower triangular factor of matrix
c       gstup(lpont(neqns)) - upper triangular factor of matrix
c       gstdi(neqns)        - inverse of diagonal matrix in triangular 
c                              factor
c
c.... local  parameters
c
c       tolle               - tollerance used to check for null pivot
c                              it should be set to approximate half-word
c                              precision 
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      integer  lpont(neqns)
      real*8   gstlo(*),                                    !gstlo(lpont(neqns))
     .         gstup(*),                                    !gstup(lpont(neqns))
     .         gstdi(neqns)

      nillc=0                                               ! Number of ILL-C
      zeros=0.0d0
      unity=1.0d0
      tolle=0.5d-07
c
c*** If isymm=0 or isymm=1 perform triangular decomposition
c
c*** Loop through the columns to perform the triangular decomposition
c
      kpivo=0
      idiag=1
      do 200 ieqns=1,neqns
        irows=idiag+1
        idiag=lpont(ieqns)
        iheig=idiag-irows
        if(iheig.gt.0) then
          istar=ieqns-iheig
          iends=ieqns-1
c
c*** If diagonal is zero compute a norm for singularity test
c
          if(gstdi(ieqns).eq.zeros) 
     .      call skycek(gstup(irows),iheig,saval)
          do 100 jeqns=istar,iends
            irows=irows+1
            jdiag=lpont(jeqns)
            jheig=min0(jdiag-lpont(jeqns-1),jeqns-istar+1)
            if(jheig.gt.0) then
              irhei=irows-jheig
              jdhei=jdiag-jheig+1
              gstup(irows)=gstup(irows)-
     .          svecdo(gstup(irhei),gstlo(jdhei),jheig) 
              if(isymm.eq.0) 
     .          gstlo(irows)=gstlo(irows)-
     .          svecdo(gstlo(irhei),gstup(jdhei),jheig) 
            end if
  100     continue
        end if
c
c*** Reduce the diagonal
c
        if(iheig.ge.0) then
          pivot=gstdi(ieqns)
          irows=idiag-iheig
          irhei=ieqns-iheig-1
          call skydia(gstlo(irows),gstup(irows),gstdi(irhei),iheig+1,
     .      isymm,gstdi(ieqns))
c
c*** Check for possible errors and print warnings
c
          if(pivot.lt.zeros) kpivo=kpivo+1
          if(abs(gstdi(ieqns)).lt.tolle*abs(pivot))
     .      nillc=nillc+1
          if(abs(gstdi(ieqns)).lt.1.0d-15)
     .      write(lures,2001) ieqns
c
c*** Complete rank test for a zero diagonal case
c
          if(pivot.eq.zeros.and.iheig.gt.0) then
            if(abs(gstdi(ieqns)).lt.tolle*saval) then
              write(lures,2003) ieqns
              stop
            end if
          end if
        end if
c
c*** Store reciprocal of diagonal
c
        if(gstdi(ieqns).ne.zeros) gstdi(ieqns)=unity/gstdi(ieqns)

  200 continue
      if (nillc.ne.0) write(lures,2000) nillc
c
c*** Formats
c
 2000 format(11x,'*** ILL-CONDITIONING.',
     .  ' LOSS OF AT LEAST 7 DIGITS FOR ',I8,' EQS.')
 2001 format(11x,'*** SINGULAR MATRIX.',
     .  ' REDUCED DIAGONAL IS ZERO FOR EQ.',I5)
 2003 format(11x,'*** SINGULAR MATRIX.',
     .  ' RANK FAILURE FOR ZERO UNREDUCED DIAGONAL IN EQ. ',I5)

      end
      subroutine solsys(
     .  iffix,lnods,rhsid,unkno,ntotv,nelem,nnode,
     .  nnodt,nnodw,ndofn,npoin,nsist,kfact,iprob,
     .  istep,iiter,timei) 
c****************************************************************************
c
c**** This routine drives the library to solve a linear system
c**** of equations. 
c
c****************************************************************************
      implicit none
      include   'MatMan.h'
      integer    ntotv,nelem,nnode,nnodt,nnodw,ndofn,npoin,
     .           nsist,iprob,kfact,iiter,istep
      integer    lnods(nnodt,nelem), iffix(ntotv)
      real*8     rhsid(ntotv,nsist), unkno(ntotv,nsist)
      real*8     elcpu,timei,timef
      logical*1  outso
c
c*** Write titles, if required
c
      outso=.false.
      if(iiter+istep.gt.0) outso=.true.
      if(outso) write(lusol,100) wosol(iprob),iiter,istep
c
c*** Direct to appropriate solver scheme
c      
      if (sofam(iprob).eq.1) then                           ! Direct solvers
        call dirsol(rhsid,unkno,ntotv,nsist,kfact,iprob,
     .              outso)
      else if(sofam(iprob).eq.2) then                       ! EE Iterative solvers
        call ite_ee(iffix,lnods,unkno,rhsid,nnode,nnodt,
     .              nnodw,nelem,ndofn,ntotv,npoin,nsist,
     .              iprob,kfact,outso)
      else if (sofam(iprob).eq.3) then                      ! SK Iterative solvers
        call ite_sk(rhsid,unkno,ntotv,nsist,kfact,iprob,
     .              outso)
      end if
      kmadi(iprob) = abs(kmadi(iprob))                      ! Solution done
      svome(iprob) = 0
      if(outso) then
        call scputi(timef)
        elcpu=timef-timei
        write(lusol,110) elcpu
      end if
c
c*** Formats
c
  100 format(//,5x,67('#'),//,5x,
     .  '>>>> SOLVER INFORMATION FOR:  ',a15,/,5x,
     .  '     ----------------------   ',//,
     .  11x,'*** ITERATION NUM.         : ',i12/,
     .  11x,'*** TIME STEP NUM.         : ',i12)
  110 format(11x,'*** ELAPSED CPU TIME       : ',f12.2)

      end
      subroutine sparse(nequa,nonze,multi,amatr,lncol,lnrow,bvect,
     .                  nsist,ifail,kfact,lunit,svomt,svome,lengp)
c****************************************************************************
c
c*** This routine solves a linear system  A x = b using
c*** Gaussian elimination. A is stored in a sparse way.
c
c    nequa: number of equations
c    nonze: number of non-zero elements in A
c    multi: multi*nonze is the dimension of amatr
c    amatr: matrix of the system (multi*nonze)
c    lncol: lncol(j) contains number of the column of amatr(j)
c    lnrow: lnrow(j) contains number of the row of amatr(j)
c    bvect: rigth hand side (nequa,nsist)
c    nsist: number of systems to be solved
c    lunit: unity to write errors
c    kfact: (=1) solve, (=0) not solve (only substitution is done)
c
c***************************************************************************
      implicit real*8(a-h,o-z)
      integer  ip                                           ! pointer
      real*8   amatr(multi*nonze), bvect(nequa,nsist)
      integer  lncol(multi*nonze), lnrow(multi*nonze)
      integer  svomt,svome,lengp
      
      iha=nequa+1
      ipoi1=0                                               ! pivot(nequa)
      ipoi2=ipoi1+nequa*8                                   ! ha(iha,11)
      ipoi3=ipoi2+iha*11*lengp                              ! aflag(8)
      ipoi4=ipoi3+8*8                                       ! iflag(10)
      ipoi5=ipoi4+10*lengp
      call sadmem(0,ip,ipoi5,svomt)
      call svecze(ipoi5/8,%val(ip))
      ifail=0
      call y12maf(nequa,nonze,amatr,lncol,multi*nonze,lnrow,
     .            multi*nonze,%val(ip+ipoi1),%val(ip+ipoi2),
     .            iha,%val(ip+ipoi3),%val(ip+ipoi4),bvect,
     .            ifail,kfact,nsist)
      if (ifail.ne.0.and.ifail.ne.5)
     .  call srunen('SPARSE: THE PARAMETER ifail IS NOT 0 OR 5',
     .  lunit)

      svome=max(svome,svomt)
      call sadmem(2,ip,ipoi5,svomt)
      
      end
      subroutine y12maf(n, z, a, snr, nn, rnr, nn1, pivot, ha,
     .  iha,aflag,iflag,b,ifail,kfact,nsist)
      implicit real*8 (a-b,g,p,t-y), integer   (c,f,h-n,r-s,z)
      real*8 a(nn), pivot(n), aflag(8),b(n,nsist)
      integer   snr(nn), rnr(nn1), ha(iha,11), iflag(10)

      aflag(1)=16.0d0
      aflag(2)=1.d-12
      aflag(3)=1.d+16
      aflag(4)=1.d-12
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1
      if(kfact.eq.1) then
        call y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
        if(ifail.eq.0)
     .    call y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     .    ifail)
      end if
      if(ifail.eq.0) then
        do isist=1,nsist
          call y12mdf(n,a,nn,b(1,isist),pivot,snr,ha,iha,iflag,ifail)
        end do
      end if

      end
      subroutine y12mbf(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag,
     .  iflag,ifail)
c
c
c  the non-zero elements of a sparse matrix a are prepared  in order to
c  solve the system ax=b by use of sparse matrix technique/
c
c
      implicit real*8 (a-b,g,p,t-y),integer   (c,f,h-n,r-s,z)
      real*8 a(nn), aflag(8)
      integer   snr(nn), rnr(nn1), ha(iha,11), iflag(10)

      mode=iflag(4)
      ifail=0
      if(n.lt.2)ifail=12
      if(z.le.0)ifail=13
      if(nn.lt.2*z)ifail=5
      if(nn1.lt.z)ifail=6
      if(ifail.eq.0.and.n.gt.z)ifail=14
      if(iha.lt.n)ifail=15
      if(mode.lt.0)ifail=16
      if(mode.gt.2)ifail=16
      if(ifail.ne.0) go to 22
      gt1=0.0d0
      do 10 i=1,n
      ha(i,2)=0
      ha(i,3)=0
   10 ha(i,6)=0
c
c  find the number of the non-zero elements in each row and column;move
c  the non-zero elements in the end of the arrays a and snr;find the
c  largest non-zero element in a(in absolute value).
c
      do 20 i=1,z
      t=abs(a(i))
      l3=rnr(i)
      l4=snr(i)
      if(l4.gt.n.or.l4.lt.1)ifail=24
      if(l3.gt.n.or.l3.lt.1)ifail=25
      ha(l3,3)=ha(l3,3)+1
      ha(l4,6)=ha(l4,6)+1
      if(t.gt.gt1)gt1=t
      a(z+i)=a(i)
   20 snr(z+i)=snr(i)
      if(ifail.gt.0)go to 22
c
c  store the information of the row starts(in ha(i,1))and of the column
c  starts(in ha(i,4)).
c
      l1=1
      l2=1
      do 40 i=1,n
      l3=ha(i,3)
      l4=ha(i,6)
      if(l3.gt.0)go to 21
      ifail=17
      go to 22
   21 if(l4.gt.0)go to 23
      ifail=18
      go to 22
   23 if(mode.eq.2)go to 30
      ha(i,9)=l3
      ha(i,10)=l4
      ha(i,11)=0
      ha(l3,2)=ha(l3,2)+1
      ha(i,5)=l3
   30 ha(i,1)=l1
      ha(i,4)=l2
      l1=l1+l3
      l2=l2+l4
      ha(i,3)=0
   40 ha(i,6)=0
c
c  store the non-zero elements of matrix a(ordered in rows) in the
c  first z locations of the array a.do the same for their column numbers
c
      do 50 i=1,z
      l1=z+i
      l3=rnr(i)
      l2=ha(l3,1)+ha(l3,3)
      a(l2)=a(l1)
      snr(l2)=snr(l1)
   50 ha(l3,3)=ha(l3,3)+1
c
c  store the row numbers of the non-zero elements ordered by columns in
c  the first z locations of the array rnr. store information about row
c  ends(in ha(i,3)).
c
      l4=1
      do 70 i=1,n
      if(mode.eq.2)go to 60
      if(ha(i,2).eq.0)go to 60
      ha(i,11)=l4
      l4=l4+ha(i,2)
      ha(i,2)=ha(i,11)
   60 ha(i,3)=ha(i,1)+ha(i,3)-1
      l1=ha(i,1)
      l2=ha(i,3)
      do 70 j=l1,l2
      l3=snr(j)
      r=ha(l3,6)
      index=ha(l3,4)+r
      rnr(index)=i
      if(r.eq.0)go to 70
      if(j.eq.l1)go to 70
      if(rnr(index-1).ne.i)go to 70
      ifail=11
      go to 22
   70 ha(l3,6)=r+1
      do 90 i=1,n
      if(mode.eq.2)go to 80
      l3=ha(i,5)
      l5=ha(l3,2)
      ha(l5,8)=i
      ha(i,7)=l5
      ha(l3,2)=ha(l3,2)+1
   80 continue
   90 ha(i,6)=ha(i,4)+ha(i,6)-1
      aflag(6)=gt1
      iflag(6)=0
      iflag(7)=0
      iflag(8)=z
      iflag(1)=-1
   22 return
     
      end
      subroutine y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha, aflag,
     .  iflag,ifail)
c
c  systens of linear equations are solved by use of sparse matrix tech-
c  nique and by gaussian elimination.
c
      implicit real*8(a-b,g,p,t-y),integer  (c,f,h-n,r-s,z)
      real*8 a(nn),b(n),pivot(n),aflag(8)
c
c  information which is necessary to begin the elimination is stored.
c
      integer   snr(nn),rnr(nn1),ha(iha,11), iflag(10)

      ifail=0
      if(iflag(1).ne.-1)ifail=2
      if(aflag(1).lt.1.0d0)aflag(1)=1.0005 d0
      if(aflag(3).lt.1.0d+5)aflag(3)=1.0d+5
      if(aflag(4).lt.0.0d0)aflag(4)=-aflag(4)
      if(iflag(2).lt.1)ifail=19
      if(iflag(3).lt.0.or.iflag(3).gt.2)ifail=20
      if(iflag(5).lt.1.or.iflag(5).gt.3)ifail=21
      if(iflag(5).eq.3)ifail=22
      if(ifail.gt.0)go to 1110
      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
c
c  use the information about fill-ins if it is possible.
c
      zz=z
      nr=n*n
      if(iflag(4).ne.2)go to 100
      if(iflag(10).gt.nn)go to 50
      l1=iflag(10)
      l5=l1+1
      if(l5.le.nn)snr(l5)=0
      do 40 i=1,n
      l=n8-i
      l2=ha(l,3)+1
      l3=l2-ha(l,1)
      do 10 j=1,l3
      snr(l5-j)=snr(l2-j)
   10 a(l5-j)=a(l2-j)
      ha(l,3)=l1
      ha(l,1)=l5-l3
      l6=l1-l3
      l5=l5-ha(l,9)
      if(l5.gt.l6)go to 30
      do 20 j=l5,l6
   20 snr(j)=0
   30 continue
   40 l1=l5-1
   50 if(iflag(9).gt.nn1)go to 100
      l2=iflag(9)
      l5=l2+1
      if(l5.le.nn1)rnr(l5)=0
      do 90 i=1,n
      l=n8-i
      l1=ha(l,6)+1
      l4=l1-ha(l,4)
      do 60 j=1,l4
   60 rnr(l5-j)=rnr(l1-j)
      ha(l,4)=l5-l4
      ha(l,6)=l2
      l6=l2-l4
      l5=l5-ha(l,10)
      if(l5.gt.l6)go to 80
      do 70 j=l5,l6
   70 rnr(j)=0
   80 continue
   90 l2=l5-1
  100 r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)
      do 110 i=1,n
      pivot(i)=0.0 d0
      ha(i,2)=ha(i,1)
  110 ha(i,5)=ha(i,4)
      index=ha(n,8)
c
c  start of gaussian elimination.
c
      slut=ha(index,3)-ha(index,2)+1
      do 950 i=1,n7
      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0)go to 350
      if(iflag(4).ne.2)go to 120
      rrow=ha(i,7)
      rcoll=ha(i,8)
      go to 220
  120 l4=ha(i,8)
      if(iflag(3).eq.1)go to 130
      rrow=l4
      rcoll=rrow
      rpivot=i
      go to 170
  130 r=nr
      v=0.0 d0
      index=iflag(2)
      do 160 kk=1,index
      l1=i-1+kk
      if(l1.gt.n)go to 170
      j=ha(l1,8)
      r7=ha(j,2)
      r8=ha(j,3)
      r9=r8-r7
      t=0.0 d0
      do 140 k=r7,r8
      td=abs(a(k))
  140 if(t.lt.td)t=td
      t=t/u
      do 160 k=r7,r8
      td=abs(a(k))
      if(td.lt.t)go to 150
      r6=snr(k)
      r3=r9*(ha(r6,6)-ha(r6,5))
      if(r3.gt.r)go to 150
      if(r3.lt.r)go to 151
      if(v.ge.td)go to 150
  151 v=td
      rrow=j
      rcoll=r6
      r=r3
      rpivot=l1
  150 continue
  160 continue
  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
c
c  remove the pivot row of the list where the rows are ordered by
c  increasing numbers of non-zero elements.
c
      ha(i,9)=r3
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1)ha(l2,11)=l
      if(l.gt.n)go to 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1
      if(rpivot.lt.l)go to 190
      ha(l4,7)=l
      ha(l,8)=l4
      l4=l5
      l1=l2
      l2=l3
      l3=n8
      go to 180
  190 if(l2.eq.l1)go to 200
      if(l3.eq.l2)go to 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i)go to 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
c
c  row interchanges.
c
      ha(i,8)=rcoll
  220 if(rrow.eq.i)go to 290
      t=b(rrow)
      b(rrow)=b(i)
      b(i)=t
      do 250 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
      r10=ha(l1,6)
  240 r=r+1
      if(rnr(r).ne.i)go to 240
      rnr(r)=rnr(r10)
  250 rnr(r10)=rrow
      rr3=ha(rrow,2)
      rr4=ha(rrow,3)
      do 270 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
  260 r=r+1
      if(rnr(r).ne.rrow)go to 260
  270 rnr(r)=i
      do 280 j=1,3
      r3=ha(rrow,j)
      ha(rrow,j)=ha(i,j)
c
c  column interchanges.
c
  280 ha(i,j)=r3
  290 if(rcoll.eq.i)go to 350
      do 310 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
      r10=ha(l1,3)
  300 r=r+1
      if(snr(r).ne.i)go to 300
      t=a(r10)
      a(r10)=a(r)
      a(r)=t
      snr(r)=snr(r10)
  310 snr(r10)=rcoll
      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)
      do 330 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
  320 r=r+1
      if(snr(r).ne.rcoll)go to 320
  330 snr(r)=i
      do 340 j=4,6
      r3=ha(rcoll,j)
      ha(rcoll,j)=ha(i,j)
c
c end of the interchanges.
c the row ordered list and the column ordered list are prepared to
c begin step i of the elimination.
c
  340 ha(i,j)=r3
  350 r9=rr4-rr3
      do 360 rr=rr3,rr4
      if(snr(rr).eq.i)go to 370
  360 continue
      ifail=9
      go to 1110
  370 v=a(rr)
      pivot(i)=v
      td=abs(v)
      if(td.lt.aflag(8))aflag(8)=td
      if(td.ge.grmin)go to 380
      ifail=3
      go to 1110
  380 r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0)go to 431
      do 430 j=rr3,rr4
      index=snr(j)
  430 pivot(index)=a(j)
  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i)go to 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i)go to 390
      t=a(l)/v
      if(iflag(5).eq.2)go to 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      go to 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      b(r1)=b(r1)-b(i)*t
      if(r9.le.0)go to 669
      r=rr1
      if(r.gt.rr2)go to 470
      do 460 l=r,rr2
      l1=snr(l)
      td=pivot(l1)
      if(td.eq.0.0d0)go to 450
      pivot(l1)=0.0 d0
      td=a(l)-td*t
      a(l)=td
      td1=abs(td)
      if(td1.gt.aflag(7))aflag(7)=td1
c
c  too small element is created.remove it from the lists.
c
      if(td1.gt.aflag(2))go to 450
      z=z-1
      a(l)=a(rr1)
      snr(l)=snr(rr1)
      a(rr1)=a(i1)
      snr(rr1)=snr(i1)
      snr(i1)=0
      rr1=rr1+1
      i1=i1+1
      ha(r1,2)=rr1
      ha(r1,1)=i1
      r3=ha(l1,5)
      r2=r3-1
      l4=ha(l1,4)
      l5=rnr(l4)
      l6=rnr(r3)
  440 r2=r2+1
      if(rnr(r2).ne.r1)go to 440
      rnr(r2)=l6
      rnr(r3)=l5
      rnr(l4)=0
      ha(l1,5)=r3+1
      ha(l1,4)=l4+1
  450 continue
  460 continue
  470 continue
      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0d0)go to 740
      tol3=-tol2*t
      tol1=abs(tol3)
      if(tol1.lt.aflag(2))go to 740
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2
      if(iflag(4).ne.1)go to 480
      if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
      if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
  480 if(i1.eq.1)go to 490
      if(snr(i1-1).eq.0)go to 600
  490 if(rr2.eq.nn)go to 500
      if(snr(rr2+1).eq.0)go to 580
c
c  collection in row ordered list.
c
  500 r10=nn-lfr
      if(r10.ge.r4)go to 560
      iflag(6)=iflag(6)+1
      do 520 jj=1,n
      l1=ha(jj,3)
      if(l1.lt.ha(jj,1))go to 510
      ha(jj,3)=snr(l1)
      snr(l1)=-jj
  510 continue
  520 continue
      l3=0
      l4=1
      do 550 jj=1,r4
      if(snr(jj).eq.0)go to 540
      l3=l3+1
      if(snr(jj).gt.0)go to 530
      l5=-snr(jj)
      snr(jj)=ha(l5,3)
      ha(l5,3)=l3
      l6=l4+ha(l5,2)-ha(l5,1)
      ha(l5,2)=l6
      ha(l5,1)=l4
      l4=l3+1
  530 a(l3)=a(jj)
      snr(l3)=snr(jj)
  540 continue
  550 continue
      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j
      if(r10.ge.r4)go to 560
      ifail=5
c
c fill-in takes place in the row ordered list.
c
      go to 1110
  560 r8=lfr-1
      rr2=r4+lfr
      if(r8.le.0)go to 579
      l3=i1-1
      do 570 ll=1,r8
      l4=r4+ll
      l5=l3+ll
      a(l4)=a(l5)
      snr(l4)=snr(l5)
  570 snr(l5)=0
  579 rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      go to 590
  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4)go to 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      go to 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=abs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z) iflag(8)=z
      if(c2.eq.1)go to 620
      if(rnr(c2-1).eq.0)go to 720
  620 if(cr2.eq.nn1)go to 630
      if(rnr(cr2+1).eq.0)go to 700
c
c  collection in column ordered list.
c
  630 r10=nn1-lfc
      if(r10.ge.r5)go to 680
      iflag(7)=iflag(7)+1
      do 640 jj=i,n
      l1=ha(jj,6)
      ha(jj,6)=rnr(l1)
  640 rnr(l1)=-jj
      l3=0
      l4=1
      do 670 jj=1,r5
      if(rnr(jj).eq.0)go to 660
      l3=l3+1
      if(rnr(jj).gt.0)go to 650
      l5=-rnr(jj)
      rnr(jj)=ha(l5,6)
      ha(l5,6)=l3
      l6=l4+ha(l5,5)-ha(l5,4)
      ha(l5,5)=l6
      ha(l5,4)=l4
      l4=l3+1
  650 rnr(l3)=rnr(jj)
  660 continue
  670 continue
      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)
      if(r10.ge.r5)go to 680
      ifail=6
c
c fill-in takes place in the column ordered list.
c
      go to 1110
  680 r8=lfc-1
      cr2=r5+lfc
      if(r8.le.0)go to 699
      l3=c2-1
      do 690 l=1,r8
      l4=r5+l
      l5=l3+l
      rnr(l4)=rnr(l5)
  690 rnr(l5)=0
  699 cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      go to 710
  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5)go to 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      go to 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2)go to 760
      ifail=7
c
c  update the information in the list where the rows are ordered by
c  increasing numbers of the non-zero elements.
c
      go to 1110
  760 if(iflag(4).eq.2)go to 870
      if(iflag(3).eq.0)go to 870
      l1=rr2-rr1+1
      if(l1.eq.l2)go to 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2)go to 820
      if(l6.gt.l4)go to 780
      if(l4.eq.n)go to 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 790
  770 ha(l2,11)=0
      go to 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1)go to 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2)go to 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      go to 870
  820 if(l6.gt.l4)go to 840
      if(l4.eq.n)go to 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut)go to 850
      l3=n
      slut=l1
      l2=l1
      go to 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1)go to 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1)go to 840
  870 continue
  880 continue
      if(r9.le.0)go to 882
      do 881 j=rr3,rr4
      index=snr(j)
  881 pivot(index)=0.0 d0
  882 continue
      cr3=ha(i,4)
      do 890 j=cr3,cr4
  890 rnr(j)=0
      if(r9.le.0)go to 930
      l2=ha(i,2)-1
      do 920 ll=1,r9
      r=snr(l2+ll)
      r1=ha(r,5)
      r2=ha(r,6)
      if(r2.gt.r1)go to 900
      ifail=8
      go to 1110
  900 ha(r,5)=r1+1
      r3=r1-1
  910 r3=r3+1
      if(rnr(r3).ne.i)go to 910
      rnr(r3)=rnr(r1)
  920 rnr(r1)=i
  930 aflag(5)=aflag(7)/aflag(6)
      if(aflag(5).lt.aflag(3))go to 940
      ifail=4
      go to 1110
  940 continue
c
c  preparation to begin the back substitution.
c
  950 continue
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0 d0
      td=abs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td
      if(td.gt.grmin)go to 960
      ifail=3
      go to 1110
  960 if(iflag(4).ne.1)go to 1060
      iflag(10)=ha(n,9)
      iflag(9)=ha(n,10)
      do 990 i=1,n7
      r1=n-i
      iflag(10)=iflag(10)+ha(r1,9)
      iflag(9)=iflag(9)+ha(r1,10)
      if(iflag(3).eq.0)go to 980
      do 970 j=9,10
      r2=ha(r1,j-2)
      r6=ha(r2,j)
      ha(r2,j)=ha(r1,j)
  970 ha(r1,j)=r6
  980 continue
  990 continue
1060  continue
      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end
      subroutine y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
      implicit real*8 (a-b,g,p,t-y),integer   (c,f,h-n,r-s,z)
      real*8 a(nn), pivot(n), b(n)
      integer   snr(nn), ha(iha,11), iflag(10)

      ifail=0
      if(iflag(1).eq.-2)go to 1000
      ifail=1
      go to 1110
 1000 mode=iflag(4)
      ipiv=iflag(3)
      n8=n+1
      n7=n-1
      state=iflag(5)
c
c  solve the system with lower triangular matrix  l  (if the
c  lu-factorization is available).
c
      if(state.ne.3)go to 1051
      if(ipiv.eq.0)go to 1020
      do 1010 i=1,n7
      l1=ha(i,7)
      t=b(l1)
      b(l1)=b(i)
      b(i)=t
 1010 continue
 1020 continue
      do 1050 i=1,n
      rr1=ha(i,1)
      rr2=ha(i,2)-1
      if(rr1.gt.rr2)go to 1040
      do 1030 j=rr1,rr2
      l1=snr(j)
 1030 b(i)=b(i)-a(j)*b(l1)
 1040 continue
 1050 continue
c
c  solve the system with upper triagular matrix.
c
 1051 continue
      do 1090 i=1,n
      r1=n8-i
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      if(rr2.lt.rr1)   go to 1080
      do 1070 j=rr1,rr2
      r2=snr(j)
 1070 b(r1)=b(r1)-a(j)*b(r2)
 1080 continue
 1090 b(r1)=b(r1)/pivot(r1)
c
c if interchanges were used during the  elimination then a reordering in
c lution vector is made.
c
      if(ipiv.eq.0)go to 1110
      do 1100 i=1,n7
      r1=n-i
      r2=ha(r1,8)
      t=b(r2)
      b(r2)=b(r1)
 1100 b(r1)=t
 1110 return
      end
      subroutine assdia(lnods,elmat,adiag,nnode,nnodt,nnodw,
     .                  nevaw,ndofn,npoin,npois,nequa,ksymm)
c****************************************************************************
c
c**** This routine assembles the diagonal for pre-conditioning
c
c****************************************************************************
      implicit none
      integer   nnode,nnodt,nnodw,nevaw,ndofn,npoin,npois,
     .          nequa,ksymm,icolu,ntotm,ievat,inode,itott,
     .          idofn,ipoin
      integer   lnods(nnodt)
      real*8    adiag(nequa), elmat(nevaw*nevaw)

      icolu=0
      ntotm=nequa-npois                                     ! = ndofn*npoin
      if (ksymm.eq.1) then                                  ! Non symmetric
        ievat=0                                               
        do inode=1,nnode                                    ! Main assembly 
          ipoin=lnods(inode)
          if(ipoin.ne.0) then
            itott=(ipoin-1)*ndofn   
            do idofn=1,ndofn               
              ievat=ievat+1
              itott=itott+1
              adiag(itott)=adiag(itott)+elmat(ievat+icolu)
              icolu=icolu+nevaw
            end do
          end if
        end do
        do inode=nnode+1,nnodw                              ! Secondary assemb.
          ipoin=lnods(inode)
          if(ipoin.ne.0) then
            ievat=ievat+1 
            itott=ntotm+ipoin-npoin    
            adiag(itott)=adiag(itott)+elmat(ievat+icolu)
            icolu=icolu+nevaw
          end if
        end do
      else                                                  ! Symmetric case
        ievat=0
        do inode=1,nnode                                    ! Main assembly 
          ipoin=lnods(inode)
          if(ipoin.ne.0) then
            itott=(ipoin-1)*ndofn 
            do idofn=1,ndofn             
              ievat=ievat+1
              itott=itott+1
              adiag(itott)=adiag(itott)+elmat(ievat+icolu)
              icolu=icolu+nevaw-ievat
            end do
          end if
        end do
        do inode=nnode+1,nnodw                              ! Secondary asse. 
          ipoin=lnods(inode)
          if(ipoin.ne.0) then
            ievat=ievat+1   
            itott=ntotm+ipoin-npoin  
            adiag(itott)=adiag(itott)+elmat(ievat+icolu)
            icolu=icolu+nevaw-ievat
          end if
        end do
      end if

      end
      
      subroutine cgprec(lnods,eamat,lengp,adiag,x_k,bvect,toler,kmaxi,
     .                  r_k_1,p_k,z_k_1,w_k,elpve,elwve,nnode,nnodt,
     .                  nnodw,nelem,nequa,nevaw,ndofn,npoin,npois,
     .                  iffix,lures,outso)
c****************************************************************************
c
c**** This routine solves a linear system of equations by the CG method.
c**** A diagonal preconditioning is used.
c
c****************************************************************************
      implicit none
      integer   eamat                                       ! pointer
      integer   nnode,nnodt,nnodw,nelem,nequa,nevaw,ndofn,npoin,
     .          npois,iequa,kiter,kmaxi,lengp,lures
      real*8    toler,bnor2,rnor2,tole1,beta1,beta2,betak,alfak,
     .          reini,refin,rerel
      integer   lnods(nnodt,nelem), iffix(nequa)
      real*8    adiag(nequa), bvect(nequa), x_k(nequa),
     .          r_k_1(nequa), p_k(nequa),   z_k_1(nequa), 
     .          w_k(nequa),   elpve(nevaw), elwve(nevaw)
      logical*1 outso

      call mabvec(lnods,eamat,lengp,x_k,w_k,elpve,elwve,    ! w_k=A*x_k
     .            nequa,nevaw,nnode,nnodt,nnodw,ndofn,
     .            npoin,npois,nelem, 0)
      do iequa=1,nequa
        if(abs(adiag(iequa)).lt.1.d-15) adiag(iequa)=1.0
        r_k_1(iequa)=bvect(iequa)-w_k(iequa)                ! r_0=b-A*x
      end do
      bnor2=0.0d0
      rnor2=0.0d0
      do iequa=1,nequa
        if(iffix(iequa).ne.1) then
          bnor2=bnor2+bvect(iequa)*bvect(iequa)
          rnor2=rnor2+r_k_1(iequa)*r_k_1(iequa)
        end if
      end do 
      reini=sqrt(rnor2/bnor2)
      kiter=0                                               ! k=0
      tole1=toler*toler*bnor2

      do while ((rnor2.gt.tole1.and.kiter.lt.kmaxi).or.kiter.eq.0)
        kiter=kiter+1                                       ! k=k+1
        do iequa=1,nequa
          z_k_1(iequa)=r_k_1(iequa)/adiag(iequa)            ! Solve M*z_k=r_k
        end do
        if (kiter.eq.1) then
          beta1=0.0d0
          do iequa=1,nequa
            p_k(iequa)=z_k_1(iequa)                         ! p_1=z_0
            beta1=beta1+r_k_1(iequa)*z_k_1(iequa)           ! r_0*z_0
          end do
        else
          beta2=beta1
          beta1=0.0d0
          do iequa=1,nequa
            beta1=beta1+r_k_1(iequa)*z_k_1(iequa)           ! r_k-1*z_k-1
          end do
          betak=beta1/beta2                                 ! B_k
          do iequa=1,nequa
            p_k(iequa)=z_k_1(iequa)+betak*p_k(iequa)        ! p_k=z_k-1+B_k*p_k-1
          end do
        end if
        call mabvec(lnods,eamat,lengp,p_k,w_k,elpve,elwve,  ! w_k=A*p_k
     .              nequa,nevaw,nnode,nnodt,nnodw,ndofn,
     .              npoin,npois,nelem, 0)
        alfak=0.0d0
        do iequa=1,nequa
          alfak=alfak+p_k(iequa)*w_k(iequa)                 ! p_k*w_k
        end do
        alfak=beta1/alfak                                   ! Alfa_k
        rerel=rnor2
        rnor2=0.0
        do iequa=1,nequa
          x_k(iequa)=x_k(iequa)+alfak*p_k(iequa)      ! x_k=x_k-1 + Alfa_k*p_k
          r_k_1(iequa)=r_k_1(iequa)-alfak*w_k(iequa)  ! r_k=r_k-1 + Alfa_k*w_k
          rnor2=rnor2+r_k_1(iequa)*r_k_1(iequa)             ! rnor2=r_k*r_k
        end do
      end do

      refin=sqrt(rnor2/bnor2)
      rerel=sqrt(rnor2/rerel)
      if(outso) write(lures,100) reini,refin,rerel,kiter

  100 format(11x,'*** INITIAL RESIDUAL NORM  : ',e12.5,/,
     .       11x,'*** FINAL RESIDUAL NORM    : ',e12.5,/,
     .       11x,'*** LAST RESIDUAL FACTOR   : ',e12.5,/,
     .       11x,'*** NUMBER OF ITERATIONS   : ',i12)

      end
      subroutine conjgr(
     .  lnods,eamat,xvect,bvect,toler,kmaxi,
     .  nnode,nnodt,nnodw,nelem,nequa,ndofn,
     .  npoin,lengp,iffix,lures,svome,svomt,
     .  kprec,outso)
c****************************************************************************
c
c**** This routine solves a linear system of equations by the conjugate 
c**** gradient. A diagonal preconditioning is used.
c
c****************************************************************************
      implicit none
      integer   a,eamat,p_mat                               ! pointer
      integer
     .  nnode,nnodt,nnodw,nelem,nequa,nevaw,ndofn,npoin,
     .  npois,ielem,adiag,r_k_1,p_k  ,z_k_1,lengp,w_k  ,
     .  elpve,elwve,lbyts,ndimt,kmaxi,lures,svome,
     .  svomt,kprec
      integer   lnods(nnodt,nelem),iffix(nequa)
      real*8    toler
      real*8    xvect(nequa),       bvect(nequa)
      logical*1 outso
      integer   lbyta
      data      lbyta/0/
      save      lbyta,a
c
c*** Auxiliar memory pointers
c      
      npois=nequa-npoin*ndofn
      nevaw=nnode*ndofn+nnodw-nnode
      adiag=0   
      r_k_1=adiag+nequa*8
      p_k  =r_k_1+nequa*8
      z_k_1=p_k  +nequa*8
      w_k  =z_k_1+nequa*8
      elpve=w_k  +nequa*8
      elwve=elpve+nevaw*8
      lbyts=elwve+nevaw*8
c
c*** Allocate memory and initialize auxiliar memory
c
      if (lbyta.eq.0) then
        call sadmem(0,a,lbyts,svomt)
        lbyta=lbyts
      else
        if (lbyts.gt.lbyta) then
          svome=max(svome,svomt)
          call sadmem(2,a,lbyta,svomt)
          call sadmem(0,a,lbyts,svomt)
          lbyta=lbyts
        end if
      end if
      ndimt=lbyts/8 
      call svecze(ndimt,%val(a)) 
c
c*** Assembly of the diagonal, if needed
c
      if(kprec.eq.-1) then                                  ! No precond.
        call svecpu(1.0,nequa,%val(a+adiag))
      else if(kprec.eq.-2) then                             ! Diagonal precond.
        do ielem=1,nelem
          call spoias(%val(eamat+(ielem-1)*lengp),p_mat)
          call assdia(lnods(1,ielem),%val(p_mat),       
     .      %val(a+adiag),nnode,nnodt,nnodw,nevaw,
     .      ndofn,npoin,npois,nequa,0)
        end do
      end if
c
c*** Solve the algebraic system using the preconditionned conjugate
c*** gradient method      
c
      call cgprec(lnods,eamat,lengp,%val(a+adiag),xvect,   
     .            bvect,toler,kmaxi,                       
     .            %val(a+r_k_1),%val(a+p_k),%val(a+z_k_1), 
     .            %val(a+w_k),%val(a+elpve),%val(a+elwve),
     .            nnode,nnodt,nnodw,nelem,nequa,nevaw,ndofn,
     .            npoin,npois,iffix,lures,outso)

      end
      subroutine elmbve(elmat,elvec,elrhs,nevaw,ksymm)
c****************************************************************************
c
c**** This routine multiplies a elemental matrix by a vector
c
c****************************************************************************
      implicit none
      integer    nevaw,ksymm,icolu,ievaw,jevaw
      real*8     elmat(nevaw*nevaw), elvec(nevaw), elrhs(nevaw)
        
      if (ksymm.eq.0) then                                  ! Symmetric case
        icolu=0
        do jevaw=1,nevaw
          elrhs(jevaw)=0.0d0
          do ievaw=jevaw,nevaw
            elrhs(jevaw)=elrhs(jevaw)+elmat(icolu+ievaw)
     .        *elvec(ievaw)
          end do
          icolu=icolu+nevaw-jevaw
        end do
        do ievaw=2,nevaw
          icolu=0
          do jevaw=1,ievaw-1
            elrhs(ievaw)=elrhs(ievaw)+elmat(icolu+ievaw)
     .        *elvec(jevaw)
            icolu=icolu+nevaw-jevaw
          end do
        end do
      else                                                  ! Non symmetric
        do ievaw=1,nevaw                                    ! case
          icolu=0
          elrhs(ievaw)=0.0
          do jevaw=1,nevaw
            elrhs(ievaw)=elrhs(ievaw)+elmat(icolu+ievaw)
     .        *elvec(jevaw)
            icolu=icolu+nevaw
          end do
        end do
      end if

      end
      subroutine elmodi(iffix,lnods,elmat,preel,rhsid,forel,
     .  nnode,nnodw,ndofn,nevat,npoin,ntotv,
     .  ksymm,zeroc,iflag)
c*****************************************************************************
c
c**** This routine computes the contribution to the RHS of the Dirichlet
c**** conditions and modifies the element matrices to account for them      
c
c*****************************************************************************
      implicit none
      integer    nnode,nnodw,ndofn,nevat,npoin,ntotv,ksymm,
     .  ievat,kevat,inode,ipoin,itott,idofn,ipos1,
     .  ipos2,iflag
      integer    lnods(nnodw),   iffix(ntotv)
      real*8     zeroc,adiag,adiao
      real*8     elmat(nevat*nevat),   preel(nevat),
     .  rhsid(ntotv)      ,   forel(nevat)
c
c*** Multiply A^e x^e, where x^e contains the prescribed d.o.f.
c
      if(iflag.eq.1) call elmbve(elmat,preel,forel,nevat,ksymm)
c
c*** Prescribe all the Dirichlet conditions by modifying the element
c*** matrices and RHS
c
      ievat=0
      adiao=1.0d0
      do inode=1,nnode                                      ! Main nodes
        ipoin=lnods(inode)
        if(ipoin.ne.0) then
          itott=(ipoin-1)*ndofn
          do idofn=1,ndofn
            itott=itott+1
            ievat=ievat+1
            if(iffix(itott).eq.1) then
              if(ksymm.eq.0) then
                ipos1=(ievat-1)*nevat-(ievat-1)*(ievat-2)/2+1
                adiag=elmat(ipos1)
                if(abs(adiag).lt.zeroc) adiag=adiao
                adiao=adiag
                elmat(ipos1)=adiag
                do kevat=ievat+1,nevat
                  ipos1=ipos1+1
                  elmat(ipos1)=0.0d0
                end do
                do kevat=1,ievat-1
                  ipos1=(kevat-1)*nevat-(kevat-1)*(kevat-2)/2
     .              +ievat-kevat+1
                  elmat(ipos1)=0.0d0
                end do
              else if(ksymm.eq.1) then
                ipos1=(ievat-1)*nevat+ievat
                adiag=elmat(ipos1)
                if(abs(adiag).lt.zeroc) adiag=adiao
                adiao=adiag
                do kevat=1,nevat
                  ipos1=(kevat-1)*nevat+ievat
                  ipos2=(ievat-1)*nevat+kevat
                  elmat(ipos1)=0.0d0
                  elmat(ipos2)=0.0d0
                end do
                ipos1=(ievat-1)*nevat+ievat
                elmat(ipos1)=adiag
              end if
              if(iflag.eq.1) forel(ievat)=-adiag*preel(ievat)
            end if
          end do
        end if
      end do
      
      adiao=1.0d0
      do inode=nnode+1,nnodw                                ! Secondary nodes
        ipoin=lnods(inode)
        if(ipoin.ne.0) then
          itott=npoin*ndofn+ipoin-npoin
          ievat=nnode*ndofn+inode-nnode
          if(iffix(itott).eq.1) then
            if(ksymm.eq.0) then
              ipos1=(ievat-1)*nevat-(ievat-1)*(ievat-2)/2+1
              adiag=elmat(ipos1)
              if(abs(adiag).lt.zeroc) adiag=adiao
              adiao=adiag
              elmat(ipos1)=adiag
              do kevat=ievat+1,nevat
                ipos1=ipos1+1
                elmat(ipos1)=0.0d0
              end do
              do kevat=1,ievat-1
                ipos1=(kevat-1)*nevat-(kevat-1)*(kevat-2)/2
     .            +ievat-kevat+1
                elmat(ipos1)=0.0d0
              end do
            else if(ksymm.eq.1) then
              ipos1=(ievat-1)*nevat+ievat
              adiag=elmat(ipos1)
              if(abs(adiag).lt.zeroc) adiag=adiao
              adiao=adiag
              do kevat=1,nevat
                ipos1=(kevat-1)*nevat+ievat
                ipos2=(ievat-1)*nevat+kevat
                elmat(ipos1)=0.0d0
                elmat(ipos2)=0.0d0
              end do
              ipos1=(ievat-1)*nevat+ievat
              elmat(ipos1)=adiag
            end if
            if(iflag.eq.1) forel(ievat)=-adiag*preel(ievat)
          end if
        end if
      end do
      
      end
      subroutine gmrabx(
     .  eamat,lengp,x,diag,nequa,y,b,job,lnods,elpve,
     .  elwve,nevaw,nnode,nnodt,nnodw,ndofn,npoin,npois,
     .  nelem,ksymm)
c*************************************************************************
c
c**** This routine performs the following operations:
c      
c     JOB = 0 : compute Ax,
c         = 1 : compute b-Ax
c      
c**** The result is put in the vector Y
c     
c*************************************************************************   
      implicit none
      integer    eamat                                      ! pointer
      integer    nevaw,nnode,nnodt,nnodw,ndofn,npoin,
     .           npois,nelem
      integer    lnods(nnodt,nelem),job,nequa,ksymm,lengp
      integer    i
      real*8     diag(nequa), x(nequa), y(nequa), b(nequa),
     .           elpve(nevaw), elwve(nevaw)

      call mabvec(
     .  lnods,eamat,lengp,x,y,elpve,elwve,nequa,
     .  nevaw,nnode,nnodt,nnodw,ndofn,npoin,npois,
     .  nelem,ksymm)
      
      do i = 1, nequa
        y(i) = y(i) / diag(i)
      end do
      
      if (job.eq.1) then
        do i = 1, nequa
          y(i) = b(i) - y(i)
        end do
      end if

      end
      subroutine gmrcor(eamat,lengp,diag,work,work1,x, idim, k0, 
     .                  itmax, toler, yours, iprint, lnods,
     .                  elpve,elwve,nevaw,nnode,nnodt,nnodw,
     .                  ndofn,npoin,npois,nelem,ksymm,lures,
     .                  nwork,iffix)
c*************************************************************************
c
c**** Basic routine of the GMRES iterative solver
c      
c     DATA:
c      
c     idim   : dimension of the problem (the matrix is (idim,idim))
c     k0     : dimension of the krylov subspace
c     itmax  : maximal number of iterations
c     toler  : relative precision needed to stop the process
c     Ax     : external procedure wich computes Ax or b-Ax
c     iprint : = 0 nothing is printed
c              = 1 the initial residual and the final relative residual
c                  are printed
c              = 2 the inial residual and all the relative residuals are
c                  printed
c                                       
c    RESULT:
c      
c    x : solution of Ax=b
c                                                
c    WORK ARRAYS :
c      
c    work  : it's a work array I need for gmres
c            the dimension of this array must be, at least of
c            idim*(k0+1)+(k0*(k0+1))/2+ 4*k0+2
c    yours : RHS of the equation (i.e., vector b)
c                                                   
c*************************************************************************
      implicit real*8 (a-h,o-z)
      integer   eamat                                       ! pointer
      integer   nevaw,nnode,nnodt,nnodw,ndofn,npoin,npois,nelem,lengp
      integer   idim, k0, itmax, iprint
      integer   lu, luj, luipl1, lsin, lcos, le, ly, lwork
      integer   lh, lhi, lhip1i, lhji, iter
      integer   ksymm,lures,nwork,lnods(nnodt,nelem),iffix(idim)
      real*8    resrel,rinit,reslast
      real*8    x(idim), work(nwork), yours(idim)
      real*8    diag(idim), elpve(nevaw), elwve(nevaw), work1(idim,k0+1)

      do i=1,idim
        if (abs(diag(i)).lt.1.d-15) diag(i)=1.0d0
      end do
c
c*** Memory allocation for the basis of Krylov subspaces, for the Hessenberg
c*** matrix, for cos(teta) and sin(teta), for y and e
c      
      lu   = 1
      lh   = 1
      lcos = lh+(k0*(k0+1))/2+1
      lsin = lcos+k0
      le   = lsin+k0
      ly   = le+k0+1
      lwork= ly+k0
      iter = 0
c
c*** Reference norm to compare the residuals (norm of the RHS excluding
c*** the prescribed degrees of freedom) and preconditioning     
c
      rinit=0.0d0
      do i = 1, idim
        yours(i) = yours(i)/diag(i)
        if(iffix(i).ne.1) rinit=rinit+yours(i)*yours(i)
      end do
      rinit=sqrt(rinit)
      residu=rinit
c
c*** Label to define the external loop
c
   20 continue
c      
c*** Compute the first residual from zero or from the last solution obtained
c      
      call gmrabx(eamat,lengp,x,diag,idim,work1(1,lu),yours, 
     .                1,lnods,elpve,elwve,nevaw,nnode,nnodt,
     .            nnodw,ndofn,npoin,npois,nelem,ksymm)

      prod=svecdo(work1(1,lu),work1(1,lu),idim)
      prod   = dsqrt(prod)
      prodin = 1.0d0/prod
      call svecbs(idim,work1(1,lu),work1(1,lu),prodin)
      work(le) = prod
      if ((iprint.gt.0).and.(iter.eq.0)) then
        write(lures,100) prod/rinit        
      endif
c
c*** Initialisation of the error vector
c      
      do i = 1, k0
        work(le+i) = 0.d0
      end do

      lhi = lh
      i = 0
c
c*** Label to define the internal loop : minimisation of the residual
c*** on a Krylov subspace
c      
   35 i = i+1
      luipl1 = lu+i
      call gmrabx(eamat,lengp,work1(1,luipl1-1),diag,idim,
     .            work1(1,luipl1),yours,0,
     .            lnods,elpve,elwve,nevaw,nnode,nnodt,nnodw,ndofn,
     .            npoin,npois,nelem,ksymm)
c      
c*** Gram-Schmidt orthogonalisation and construction of the column number
c*** i of the Hessenberg matrix :
c
      do j = 1, i
        luj  = lu+(j-1)
        lhji = lhi+j-1
        work(lhji)=svecdo(work1(1,luipl1),work1(1,luj),idim)
        coef = -work(lhji)
        call svecco(1.d0,work1(1,luipl1),coef,work1(1,luj),
     .              work1(1,luipl1),idim)
      end do

      lhip1i = lhi+i  
      work(lhip1i)=svecdo(work1(1,luipl1),work1(1,luipl1),idim)
      work(lhip1i) = dsqrt(work(lhip1i))
      
      coef = 1.d0/work(lhip1i)
      call svecbs(idim,work1(1,luipl1),work1(1,luipl1),coef)
c                                     
c*** QR factorisation of the Hessenberg matrix
c
      do j = 1, i-1
        lhji         = lhi+j-1
        lhjp1i       = lhi+j
        hji          = work(lhji)   
        hjp1i        = work(lhjp1i) 
        work(lhji)   = work(lcos+j-1)*hji+work(lsin+j-1)*hjp1i
        work(lhjp1i) = work(lcos+j-1)*hjp1i-work(lsin+j-1)*hji
      end do          
      hii            = work(lhi+i-1)
      hip1i          = work(lhi+i)  
      sqrhi          = dsqrt(hii*hii+hip1i*hip1i)
      work(lcos+i-1) = hii/sqrhi
      work(lsin+i-1) = hip1i/sqrhi
      work(lhi+i-1)  = sqrhi
      ei             = work(le+i-1)
      work(le+i-1)   = work(lcos+i-1)*ei
      work(le+i)     = -work(lsin+i-1)*ei
      reslast        = residu
      residu         = abs(work(le+i))/rinit
      iter = iter+1
      lhi = lhi+i
c
c*** End of the internal loop
c      
      if (((residu.gt.toler).and.(iter.lt.itmax)).and.(i.lt.k0)) goto 35
c           
c*** Resolution of the upper triangular system
c
      call gmrsol(work(lh),work(ly),work(le),i)
c                        
c*** Computation of the new estimation of the solution
c
      do j = 1, i
        luj= lu+(j-1)
        call svecco(1.d0,x,work(ly+j-1),work1(1,luj),x,idim)
      end do
c                              
c*** Test for reinitialization (if positive, repeat external loop)
c      
      if ((residu.gt.toler).and.(iter.lt.itmax)) goto 20
      resrel=residu/reslast
      if (iprint.gt.0) write(lures,110) residu,resrel,iter
      
  100 format(11x,'*** INITIAL RESIDUAL NORM  : ',e12.5)
  110 format(11x,'*** FINAL RESIDUAL NORM    : ',e12.5,/,
     .       11x,'*** LAST RESIDUAL FACTOR   : ',e12.5,/,
     .       11x,'*** NUMBER OF ITERATIONS   : ',i12)
      
      end
      subroutine gmress(eamat,nequa,bvect,xvect,toler,kmaxi,nkryl,
     .                  lnods,nnode,nnodt,nnodw,nelem,ndofn,npoin,
     .                  ksymm,lengp,iffix,lures,svome,svomt,kprec,
     .                  outso)
c*************************************************************************
c      
c**** Solves non-symmetric linear systems of equations by the precon- 
c**** ditioned GMRES method
c      
c*************************************************************************  
      implicit none
      integer    a,eamat,p_mat,p_wo1                        ! pointer
      integer    nnode,nnodt,nnodw,nevaw,ndofn,npoin,npois,
     .           ksymm,nelem,lures,nequa,kmaxi,nkryl,lengp, 
     .           ielem,adiag,elpve,elwve,work2,lbyts,ndimt,
     .           nwor2,svome,kprec,svomt,iprin,lbypw
      integer    lnods(nnodt,nelem), iffix(nequa)
      real*8     bvect(nequa),  xvect(nequa)
      real*8     toler
      logical*1  outso
      integer    lbyta,lbypa
      data       lbyta,lbypa/0,0/
      save       lbyta,lbypa,a,p_wo1
c
c*** Initializations
c
      iprin=0
      if(outso) iprin=2
c
c*** Auxiliar memory pointers
c      
      npois=nequa-npoin*ndofn
      nevaw=nnode*ndofn+nnodw-nnode
      nwor2=(nkryl*(nkryl+1))/2+4*nkryl+2
      adiag=0
      elpve=adiag+nequa*8
      elwve=elpve+nevaw*8
      work2=elwve+nevaw*8
c
c*** Allocate memory ,if it is neccesary, and initialize auxiliar memory
c
      lbyts=work2+nwor2*8
      lbypw=(nkryl+1)*nequa*8
      if (lbyta.eq.0) then                                  !First pass
        call sadmem(0,a,lbyts,svomt)                        !Allocate memory
        lbyta=lbyts
        call sadmem(0,p_wo1,lbypw,svomt)
        lbypa=lbypw
      else
        if (lbyts.gt.lbyta) then
          svome=max(svome,svomt)
          call sadmem(2,a,lbyta,svomt)                      !Deallocate memory
          call sadmem(0,a,lbyts,svomt)                      !Alloc. new memory
          lbyta=lbyts
        end if
        if (lbypw.gt.lbypa) then
          svome=max(svome,svomt)
          call sadmem(2,p_wo1,lbypa,svomt)
          call sadmem(0,p_wo1,lbypw,svomt)
          lbypa=lbypw
        end if
      end if
      ndimt=lbyts/8
      call svecze(ndimt,%val(a))
      call svecze(lbypw/8,%val(p_wo1))
c
c*** Assembly of the diagonal, if needed
c      
      if(kprec.eq.-1) then                                  ! No precond.
        call svecpu(1.0,nequa,%val(a+adiag))
      else if(kprec.eq.-2) then                             ! Diagonal precond.
        do ielem=1,nelem
          call spoias(%val(eamat+(ielem-1)*lengp),p_mat)
          call assdia(lnods(1,ielem),%val(p_mat),          
     .      %val(a+adiag),nnode,nnodt,nnodw,nevaw,
     .      ndofn,npoin,npois,nequa,ksymm)
        end do
      end if
c
c*** Solve the algebraic system using the preconditionned GMRES method
c
      call gmrcor(eamat,lengp,%val(a+adiag),%val(a+work2),
     .            %val(p_wo1),xvect,nequa,nkryl,kmaxi,
     .            toler,bvect,iprin,lnods,%val(a+elpve), 
     .            %val(a+elwve),nevaw,nnode,nnodt,nnodw,
     .            ndofn,npoin,npois,nelem,ksymm,lures,
     .            nwor2,iffix)

      end
      subroutine gmrsol(h,y,e,idim)
c*************************************************************************
c
c**** This routine solves an upper linear system of the form H Y = E. The
c**** (idim x idim) matrix H and the vector E are data. The result is Y
c
c*************************************************************************
      implicit none
      integer  idim, i, j, lhij, lhii, ly
      real*8   h(*), y(*), e(*)

      do i=1,idim
        ly=idim-i+1
        y(ly)=e(ly)
        do j=ly+1,idim
          lhij=((j-1)*j)/2+ly
          y(ly)=y(ly)-h(lhij)*y(j)
        end do
        lhii=((ly-1)*ly)/2+ly
        y(ly)=y(ly)/h(lhii)
      end do

      end         
      subroutine ite_ee(
     .  iffix,lnods,unkno,rhsid,nnode,nnodt,
     .  nnodw,nelem,ndofn,ntotv,npoin,nsist,
     .  iprob,kfact,outso)
c****************************************************************************
c
c**** This routine drives the library to solve a linear system
c**** of equations by EE iterative solvers 
c
c****************************************************************************
      implicit none
      include  'MatMan.h'
      integer   nnode,nnodt,nnodw,nelem,ndofn,npoin,ntotv,
     .          nsist,nevat,iprob,kfact,isist,itotv
      integer   iffix(ntotv),        lnods(nnodt,nelem)
      real*8    unkno(ntotv,nsist),  rhsid(ntotv,nsist)
      logical*1 outso
      integer   p_mat,p_rhs,p_for,p_pre,iwoso               ! pointer
      integer   lbya1,lbya2,lbyt1,lbyt2,svomt
      data      lbya1,lbya2,svomt/0,0,0/
      save      lbya1,lbya2,iwoso,p_for,p_pre
c
c*** Initializations
c
      nevat=ndofn*nnode+nnodw-nnode
c
c*** Store the original RHS and drive auxiliar memory
c
      lbyt1=8*ntotv*nsist
      if (lbya1.eq.0) then
        lbya1=lbyt1
        call sadmem(0,iwoso,lbya1,svomt)
      else
        if (lbyt1.gt.lbya1) then
          svome(iprob)=max(svome(iprob),svomt)
          call sadmem(2,iwoso,lbya1,svomt)
          lbya1=lbyt1
          call sadmem(0,iwoso,lbya1,svomt)
        end if
      end if

      lbyt2=8*nevat
      if (lbya2.eq.0) then
        lbya2=lbyt2
        call sadmem(0,p_for,lbya2,svomt)
        call sadmem(0,p_pre,lbya2,svomt)
      else
        if (lbyt2.gt.lbya2) then
          svome(iprob)=max(svome(iprob),svomt)
          call sadmem(2,p_pre,lbya2,svomt)
          call sadmem(2,p_for,lbya2,svomt)
          lbya2=lbyt2
          call sadmem(0,p_for,lbya2,svomt)
          call sadmem(0,p_pre,lbya2,svomt)
        end if
      end if
        
      call slveca(ntotv*nsist,rhsid,%val(iwoso))
c
c*** If needed, prescribe Dirichlet boundary conditions
c
      call spoias(%val(irhsx+(iprob-1)*lengp),p_rhs)
      call spoias(%val(iamat+(iprob-1)*lengp),p_mat)
      if(p_rhs.ne.0) then                                   ! kimpo was not 0
        if(kfact.eq.1) then                                 ! A x_p computed
          call iterhs(iffix,lnods,rhsid,%val(p_rhs),
     .                unkno,%val(p_pre),%val(p_for),
     .                ntotv,nsist,nnode,nnodt,nnodw,
     .                ndofn,nevat,nelem,npoin,lengp,
     .                p_mat,ksymm(iprob),zeroc) 
        end if
        do isist=1,nsist
          do itotv=1,ntotv
            if(iffix(itotv).eq.1) rhsid(itotv,isist)=0.0
          end do
        end do
        call svecad(nsist*ntotv,rhsid,%val(p_rhs),rhsid)
      end if
c
c*** Conjugate gradient solver
c
      if (kites(iprob).eq.21) then
        do isist=1,nsist
          call conjgr(lnods,p_mat,unkno(1,isist),rhsid(1,isist),
     .                toler(iprob),itmax(iprob),nnode,nnodt,nnodw,
     .                nelem,nequa(iprob),ndofn,npoin,lengp,iffix,
     .                lusol,svome(iprob),svomt,kmsip(iprob),outso)
        end do
c
c*** GMRES solver
c
      else if (kites(iprob).eq.22) then
        do isist=1,nsist
          call gmress(p_mat,nequa(iprob),rhsid(1,isist),unkno(1,isist),
     .                toler(iprob),itmax(iprob),kkryl(iprob),lnods,
     .                nnode,nnodt,nnodw,nelem,ndofn,npoin,ksymm(iprob),
     .                lengp,iffix,lusol,svome(iprob),svomt,kmsip(iprob),
     .                outso)
        end do
      end if
c
c*** Restore the original RHS and deallocate the work space used for it
c
      call slveca(nsist*ntotv,%val(iwoso),rhsid)
c
c*** Write temporary memory requirements 
c
      svome(iprob)=max(svome(iprob),svomt)
      if(outso) write(lusol,120) svome(iprob)
c
c*** Formats
c
  120 format(11x,'*** VOLATILE MEMORY (BYTES): ',i12)
      
      end
      subroutine iterhs(iffix,lnods,rhsid,rhsax,unkno,preel,
     .                  forel,ntotv,nsist,nnode,nnodt,nnodw,
     .                  ndofn,nevat,nelem,npoin,lengp,eamat,
     .                  ksymm,zeroc) 
c****************************************************************************
c
c**** This routine prescribes the Dirichlet boundary conditions when itera-
c**** tive solvers are used. This is done by substracting to the original
c**** RHS the product A x_p and replacing the equations for the prescribed
c**** d.o.f. by the trivial relations a_ii x_i = a_ii x_i,p (x_i,p: pres-
c**** cribed value, a_ii: diagonal term of A^e)      
c
c****************************************************************************
      implicit none
      integer ntotv,nsist,nnode,nnodt,nnodw,ndofn,nevat,
     .        nelem,npoin,lengp,ksymm,iflag,isist,itotv,
     .        ielem,idofn,itott,ievat,ipoin,inode
      integer eamat,p_mat                                   ! pointer
      integer iffix(ntotv),         lnods(nnodt,nelem)
      real*8  rhsid(nsist,ntotv),   rhsax(nsist,ntotv),
     .        unkno(ntotv,nsist),   preel(nevat),
     .        forel(nevat) 
      real*8  xcomp,zeroc
c
c*** Initialization
c
      call svecze(nsist*ntotv,rhsax)
c
c*** Check if x_p is zero 
c
      do isist=1,nsist
        iflag=0
        do itotv=1,ntotv
          if (iffix(itotv).eq.1) then
            xcomp=unkno(itotv,isist)
            xcomp=xcomp*xcomp
            if(xcomp.gt.zeroc) iflag=1
          end if
        end do
c
c*** Compute the product A x_p and modify the element matrices and RHS
c
        do ielem=1,nelem
          if(iflag.eq.1) then
            call svecze(nevat,preel)
            call svecze(nevat,forel)
            ievat=0
c
c*** Gather of PREEL (prescribed values)
c            
            do inode=1,nnode         
              ipoin=lnods(inode,ielem)
              if(ipoin.ne.0) then
                itott=(ipoin-1)*ndofn
                do idofn=1,ndofn
                  ievat=ievat+1
                  itott=itott+1
                  if(iffix(itott).eq.1)
     .              preel(ievat)=unkno(itott,isist)
                end do
              end if
            end do
            do inode=nnode+1,nnodw
              ipoin=lnods(inode,ielem)
              if(ipoin.ne.0) then
                ievat=ievat+1
                itott=ndofn*npoin+ipoin-npoin
                if(iffix(itott).eq.1)
     .            preel(ievat)=unkno(itott,isist)
              end if
            end do
          end if
c
c*** Calculation of f^e = A^e x^e and modification of A^e
c            
          call spoias(%val(eamat+(ielem-1)*lengp),p_mat)
          call elmodi(iffix,lnods(1,ielem),%val(p_mat),
     .                preel,rhsid,forel,nnode,nnodw,ndofn,
     .                nevat,npoin,ntotv,ksymm,zeroc,iflag)
c
c*** Scatter of f^e
c
          if(iflag.eq.1) then
            ievat=0
            do inode=1,nnode
              ipoin=lnods(inode,ielem)
              if(ipoin.ne.0) then
                itott=(ipoin-1)*ndofn
                do idofn=1,ndofn
                  ievat=ievat+1
                  itott=itott+1
                  rhsax(itott,isist)=rhsax(itott,isist)
     .              -forel(ievat)
                end do
              end if
            end do
            do inode=nnode+1,nnodw
              ipoin=lnods(inode,ielem)
              if(ipoin.ne.0) then
                ievat=ievat+1
                itott=ndofn*npoin+ipoin-npoin
                rhsax(itott,isist)=rhsax(itott,isist)
     .            -forel(ievat)
              end if
            end do
          end if
        end do
      end do

      end
      subroutine mabvec(
     .  lnods,eamat,lengp,xvect,bvect,elxve,elbve,
     .  nequa,nevaw,nnode,nnodt,nnodw,ndofn,npoin,
     .  npois,nelem,ksymm)
c****************************************************************************
c
c**** This routine multiplies a non-assembled matrix by a vector A*x=b
c
c****************************************************************************
      implicit none
      integer   eamat,p_mat                                 ! pointer
      integer
     .  nnodt,nevaw,nnode,ndofn,nnodw,nequa,npoin,lengp,
     .  npois,ksymm,nelem,ielem,ievat,inode,itott,ntotm,
     .  idofn,iequa,ipoin
      integer   lnods(nnodt,nelem)
      real*8    elxve(nevaw), xvect(nequa), elbve(nevaw),
     .          bvect(nequa)

      ntotm=nequa-npois                                     ! = ndofn*npoin
      do iequa=1,nequa
        bvect(iequa)=0.0d0
      end do
      do ielem=1,nelem
        ievat=0                                             ! Gather of the
        do inode=1,nnode                                    ! XVECT->ELXVE
          ipoin=lnods(inode,ielem)
          if(ipoin.ne.0) then
            itott=(ipoin-1)*ndofn
            do idofn=1,ndofn
              ievat=ievat+1
              itott=itott+1
              elxve(ievat)=xvect(itott)
            end do
          end if
        end do
        do inode=nnode+1,nnodw
          ipoin=lnods(inode,ielem)
          if(ipoin.ne.0) then
            ievat=ievat+1
            itott=ntotm+ipoin-npoin
            elxve(ievat)=xvect(itott)
          end if
        end do
        call spoias(%val(eamat+(ielem-1)*lengp),p_mat)
        call elmbve(%val(p_mat),elxve,elbve,nevaw,          ! Multiply A.x=b
     .              ksymm)
        ievat=0                                             ! Scatter of the
        do inode=1,nnode                                    ! ELBVE->BVECT
          ipoin=lnods(inode,ielem)
          if(ipoin.ne.0) then
            itott=(ipoin-1)*ndofn
            do idofn=1,ndofn
              ievat=ievat+1
              itott=itott+1
              bvect(itott)=bvect(itott)+elbve(ievat)
            end do
          end if
        end do
        do inode=nnode+1,nnodw
          ipoin=lnods(inode,ielem)
          if(ipoin.ne.0) then
            ievat=ievat+1
            itott=ntotm+ipoin-npoin
            bvect(itott)=bvect(itott)+elbve(ievat)
          end if
        end do
      end do      

      end
      subroutine amux (n, x, y, a,ja,ia) 
c*****************************************************************************
c
c         A times a vector
c
c Multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c ---------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c ----------
c y     = real array of length n, containing the product y=Ax
c
c
c*****************************************************************************
      implicit none
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
      real*8 t
      integer i, k
      
      do i = 1,n
        t = 0.0d0
        do k=ia(i), ia(i+1)-1 
          t = t + a(k)*x(ja(k))
        end do
        y(i) = t
      end do

      end
      subroutine atmux (n, x, y, a, ja, ia)
c*****************************************************************************
c
c         transp( A ) times a vector
c
c Multiplies the transpose of a matrix by a vector when the original
c matrix is stored in compressed sparse row storage. Can also be
c viewed as the product of a matrix by a vector when the original
c matrix is stored in the compressed sparse column format.
c ----------------------------------------------------------------------
c
c on entry:
c ---------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c ----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c*****************************************************************************
      implicit none
      real*8 x(*), y(*), a(*) 
      integer n, ia(*), ja(*)
      integer i, k 

      do i=1,n
        y(i) = 0.0d0
      end do

      do i = 1,n
        do  k=ia(i), ia(i+1)-1 
          y(ja(k)) = y(ja(k)) + x(i)*a(k)
        end do
      end do
      
      end
      subroutine auxdri(nrows,ipara,fpara,iters,resit,ounit,amatr,
     .                  jacol,iarow,aumat,jauco,jurow,works,again)
c********************************************************************************
c
c**** This routine performs some auxiliar operations to ITDRIV
c      
c********************************************************************************
      implicit none
      integer   nrows,iters,ounit, 
     .          jacol(*),iarow(*),ipara(*),jauco(*), jurow(*)
      real*8    resit, fpara(*), works(*), amatr(*),aumat(*)
      logical*1 again
c
c*** Convergence can be written by using the following lines:
c      
c     if (ipara(7).ne.iters) then
c       write (ounit, *) iters, real(resit)
c       iters = ipara(7)
c     endif
      resit = fpara(5)
c
c*** Decide what to do depending on ipara(1).
c
      if (ipara(1).eq.1) then
        call amux  (nrows,works(ipara(8)),works(ipara(9)),
     .              amatr,jacol,iarow)
      else if (ipara(1).eq.2) then
        call atmux (nrows,works(ipara(8)),works(ipara(9)),
     .              amatr,jacol,iarow)
      else if (ipara(1).eq.3 .or. ipara(1).eq.5) then
        call lusol (nrows,works(ipara(8)),works(ipara(9)),
     .              aumat,jauco,jurow)
      else if (ipara(1).eq.4 .or. ipara(1).eq.6) then
        call lutsol(nrows,works(ipara(8)),works(ipara(9)),
     .              aumat,jauco,jurow)
      else if (ipara(1).le.0) then
        again = .false. 
        if (ipara(1).eq.0) then
          continue ! Iterative solver has satisfied convergence test
        else if (ipara(1).eq.-1) then
          continue
c-w       write(ounit,*)                                ! This warning may be
c-w  .      '     >>>> WARNING: ITERATIVE SOLVER HAS',  ! written, if wanted
c-w  .      ' ITERATED TOO MANY TIMES <<<< '
        else if (ipara(1).eq.-2) then
          call srunen(
     .      'AUXDRI: ITERATIVE SOLVER WAS NOT GIVEN ENOUGH SPACE.',
     .      ounit)
          call srunei(
     .      'AUXDRI: REQUIRED NUMBER OF ELEMENTS OF THIS SPACE',
     .      ipara(4),ounit)
        else if (ipara(1).eq.-3) then
          call srunen(
     .      'AUXDRI: ITERATIVE SOLVER IS FACING A BREAK-DOWN',ounit)
        else
          call srunei(
     .      'AUXDRI: ITERATIVE SOLVER TERMINATED, CODE = ',
     .      ipara(1),ounit)
        endif
      endif

      end
      subroutine bcg(n,rhs,sol,ipar,fpar,w)
c***************************************************************************
c      
c     BCG: Bi Conjugate Gradient method. Programmed with reverse
c     communication, see the header for detailed specifications
c     of the protocol.
c
c     in this routine, before successful return, the fpar's are
c     fpar(3) == initial residual norm
c     fpar(4) == target residual norm
c     fpar(5) == current residual norm
c     fpar(7) == current rho (rhok = <r, s>)
c     fpar(8) == previous rho (rhokm1)
c
c     w(:,1) -- r, the residual
c     w(:,2) -- s, the dual of the 'r'
c     w(:,3) -- p, the projection direction
c     w(:,4) -- q, the dual of the 'p'
c     w(:,5) -- v, a scratch vector to store A*p, or A*q.
c     w(:,6) -- a scratch vector to store intermediate results
c     w(:,7) -- changes in the solution
c      
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8  fpar(16), rhs(n), sol(n), w(n,*)
c
c***  External routines used
c
      real*8 vecdo
      logical stopbis,brkdn
      external vecdo, stopbis, brkdn

      real*8 one
      parameter(one=1.0D0)
c
c***  Local variables
c
      integer i
      real*8 alpha
      logical rp, lp
      save
c
c***  Status of the program
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
c
c***  Initialization, initial residual
c
      call bisinit(ipar,fpar,7*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c***  compute initial residual, request a matvecc
c
      ipar(1) = 1
      ipar(8) = 3*n+1
      ipar(9) = ipar(8) + n
      do i = 1, n
        w(i,4) = sol(i)
      enddo
      ipar(10) = 1
      return
   10 ipar(7) = ipar(7) + 1
      do i = 1, n
        w(i,1) = rhs(i) - w(i,5)
      enddo
      fpar(11) = fpar(11) + n
      if (lp) then
        ipar(1) = 3
        ipar(8) = 1
        ipar(9) = n+1
        ipar(10) = 2
        return
      endif

   20 if (lp) then
        do i = 1, n
          w(i,1) = w(i,2)
          w(i,3) = w(i,2)
          w(i,4) = w(i,2)
        enddo
      else
        do i = 1, n
          w(i,2) = w(i,1)
          w(i,3) = w(i,1)
          w(i,4) = w(i,1)
        enddo
      endif

      fpar(7) = vecdo(n,w,w)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(fpar(7))
      fpar(5) = fpar(3)
      fpar(8) = one
      if (abs(ipar(3)).eq.2) then
        fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
        fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
        fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
        fpar(6) = fpar(5)
        goto 900
      endif
c
c***  End of initialization, begin iteration, v = A p
c
   30 if (rp) then
        ipar(1) = 5
        ipar(8) = n + n + 1
        if (lp) then
          ipar(9) = 4*n + 1
        else
          ipar(9) = 5*n + 1
        endif
        ipar(10) = 3
        return
      endif

   40 ipar(1) = 1
      if (rp) then
        ipar(8) = ipar(9)
      else
        ipar(8) = n + n + 1
      endif
      if (lp) then
        ipar(9) = 5*n + 1
      else
        ipar(9) = 4*n + 1
      endif
      ipar(10) = 4
      return

   50 if (lp) then
        ipar(1) = 3
        ipar(8) = ipar(9)
        ipar(9) = 4*n + 1
        ipar(10) = 5
        return
      endif

   60 ipar(7) = ipar(7) + 1
      alpha = vecdo(n,w(1,4),w(1,5))
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha,ipar)) goto 900
      alpha = fpar(7) / alpha
      do i = 1, n
        w(i,7) = w(i,7) + alpha * w(i,3)
        w(i,1) = w(i,1) - alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + 4 * n
      if (ipar(3).eq.999) then
        ipar(1) = 10
        ipar(8) = 6*n + 1
        ipar(9) = 5*n + 1
        ipar(10) = 6
        return
      endif
   70 if (ipar(3).eq.999) then
        if (ipar(11).eq.1) goto 900
      else if (stopbis(n,ipar,1,fpar,w,w(1,3),alpha)) then
        goto 900
      endif
c
c***  A^t * x
c
      if (lp) then
        ipar(1) = 4
        ipar(8) = 3*n + 1
        if (rp) then
          ipar(9) = 4*n + 1
        else
          ipar(9) = 5*n + 1
        endif
        ipar(10) = 7
        return
      endif

   80 ipar(1) = 2
      if (lp) then
        ipar(8) = ipar(9)
      else
        ipar(8) = 3*n + 1
      endif
      if (rp) then
        ipar(9) = 5*n + 1
      else
        ipar(9) = 4*n + 1
      endif
      ipar(10) = 8
      return

   90 if (rp) then
        ipar(1) = 6
        ipar(8) = ipar(9)
        ipar(9) = 4*n + 1
        ipar(10) = 9
        return
      endif

  100 ipar(7) = ipar(7) + 1
      do i = 1, n
        w(i,2) = w(i,2) - alpha * w(i,5)
      enddo
      fpar(8) = fpar(7)
      fpar(7) = vecdo(n,w,w(1,2))
      fpar(11) = fpar(11) + 4 * n
      if (brkdn(fpar(7), ipar)) return
      alpha = fpar(7) / fpar(8)
      do i = 1, n
        w(i,3) = w(i,1) + alpha * w(i,3)
        w(i,4) = w(i,2) + alpha * w(i,4)
      enddo
      fpar(11) = fpar(11) + 4 * n
c
c***  End of the iterations
c
      goto 30
c
c***  Some clean up job to do
c
  900 if (rp) then
        ipar(1) = 5
        ipar(8) = 6*n + 1
        ipar(9) = ipar(8) - n
        ipar(10) = 10
        return
      endif
  110 if (rp) then
        call tidycg(n,ipar,fpar,sol,w(1,6))
      else
        call tidycg(n,ipar,fpar,sol,w(1,7))
      endif
      return
      end
      subroutine bcgstab(n, rhs, sol, ipar, fpar, w)
c***************************************************************************
c
c     BCGSTAB --- Bi Conjugate Gradient stabilized (BCGSTAB)
c     This is an improved BCG routine. (1) no matrix transpose is
c     involved. (2) the convergence is smoother.
c
c
c     Algorithm:
c     Initialization - r = b - A x, r0 = r, p = r, rho = (r0, r),
c     Iterate -
c     (1) v = A p
c     (2) alpha = rho / (r0, v)
c     (3) s = r - alpha v
c     (4) t = A s
c     (5) omega = (t, s) / (t, t)
c     (6) x = x + alpha * p + omega * s
c     (7) r = s - omega * t
c     convergence test goes here
c     (8) beta = rho, rho = (r0, r), beta = rho * alpha / (beta * omega)
c         p = r + beta * (p - omega * v)
c
c     in this routine, before successful return, the fpar's are
c     fpar(3) == initial (preconditionied-)residual norm
c     fpar(4) == target (preconditionied-)residual norm
c     fpar(5) == current (preconditionied-)residual norm
c     fpar(6) == current residual norm or error
c     fpar(7) == current rho (rhok = <r, r0>)
c     fpar(8) == alpha
c     fpar(9) == omega
c
c     Usage of the work space W
c     w(:, 1) = r0, the initial residual vector
c     w(:, 2) = r, current residual vector
c     w(:, 3) = s
c     w(:, 4) = t
c     w(:, 5) = v
c     w(:, 6) = p
c     w(:, 7) = tmp, used in preconditioning, etc.
c     w(:, 8) = delta x, the correction to the answer is accumulated
c               here, so that the right-preconditioning may be applied
c               at the end
c      
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,8)
      
      real*8 vecdo
      logical stopbis, brkdn
      external vecdo, stopbis, brkdn
      real*8 one
      parameter(one=1.0D0)

      integer i
      real*8 alpha,beta,rho,omega
      logical lp, rp
      save lp, rp
c
c***  Where to go
c
      if (ipar(1).gt.0) then
        goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110) ipar(10)
      else if (ipar(1).lt.0) then
        goto 900
      endif
c
c***  Call the initialization routine
c
      call bisinit(ipar,fpar,8*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c***  Perform a matvec to compute the initial residual
c
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      do i = 1, n
        w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
   10 ipar(7) = ipar(7) + 1
      do i = 1, n
        w(i,1) = rhs(i) - w(i,2)
      enddo
      fpar(11) = fpar(11) + n
      if (lp) then
        ipar(1) = 3
        ipar(10) = 2
        return
      endif

   20 if (lp) then
        do i = 1, n
          w(i,1) = w(i,2)
          w(i,6) = w(i,2)
        enddo
      else
        do i = 1, n
          w(i,2) = w(i,1)
          w(i,6) = w(i,1)
        enddo
      endif

      fpar(7) = vecdo(n,w,w)
      fpar(11) = fpar(11) + 2 * n
      fpar(5) = sqrt(fpar(7))
      fpar(3) = fpar(5)
      if (abs(ipar(3)).eq.2) then
        fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
        fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
        fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
      if (ipar(3).ge.0) fpar(6) = fpar(5)
      if (ipar(3).ge.0 .and. fpar(5).le.fpar(4) .and.
     +  ipar(3).ne.999) then
        goto 900
      endif
c
c***  Beginning of the iterations
c
c***  Step (1), v = A p
c      
   30 if (rp) then
        ipar(1) = 5
        ipar(8) = 5*n+1
        if (lp) then
          ipar(9) = 4*n + 1
        else
          ipar(9) = 6*n + 1
        endif
        ipar(10) = 3
        return
      endif
c
   40 ipar(1) = 1
      if (rp) then
        ipar(8) = ipar(9)
      else
        ipar(8) = 5*n+1
      endif
      if (lp) then
        ipar(9) = 6*n + 1
      else
        ipar(9) = 4*n + 1
      endif
      ipar(10) = 4
      return
   50 if (lp) then
        ipar(1) = 3
        ipar(8) = ipar(9)
        ipar(9) = 4*n + 1
        ipar(10) = 5
        return
      endif
c
   60 ipar(7) = ipar(7) + 1
c
c***  Step (2)
c      
      alpha = vecdo(n,w(1,1),w(1,5))
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha, ipar)) goto 900
      alpha = fpar(7) / alpha
      fpar(8) = alpha
c
c***  Step (3)
c      
      do i = 1, n
        w(i,3) = w(i,2) - alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + 2 * n
c
c***  Step (4): the second matvec -- t = A s
c
      if (rp) then
        ipar(1) = 5
        ipar(8) = n+n+1
        if (lp) then
          ipar(9) = ipar(8)+n
        else
          ipar(9) = 6*n + 1
        endif
        ipar(10) = 6
        return
      endif

   70 ipar(1) = 1
      if (rp) then
        ipar(8) = ipar(9)
      else
        ipar(8) = n+n+1
      endif
      if (lp) then
        ipar(9) = 6*n + 1
      else
        ipar(9) = 3*n + 1
      endif
      ipar(10) = 7
      return
   80 if (lp) then
        ipar(1) = 3
        ipar(8) = ipar(9)
        ipar(9) = 3*n + 1
        ipar(10) = 8
        return
      endif
   90 ipar(7) = ipar(7) + 1
c
c***  step (5)
c      
      omega = vecdo(n,w(1,4),w(1,4))
      fpar(11) = fpar(11) + n + n
      if (brkdn(omega,ipar)) goto 900
      omega = vecdo(n,w(1,4),w(1,3)) / omega
      fpar(11) = fpar(11) + n + n
      if (brkdn(omega,ipar)) goto 900
      fpar(9) = omega
      alpha = fpar(8)
c
c***  Steps (6) and (7)
c      
      do i = 1, n
        w(i,7) = alpha * w(i,6) + omega * w(i,3)
        w(i,8) = w(i,8) + w(i,7)
        w(i,2) = w(i,3) - omega * w(i,4)
      enddo
      fpar(11) = fpar(11) + 6 * n + 1
c
c***  Convergence test
c      
      if (ipar(3).eq.999) then
        ipar(1) = 10
        ipar(8) = 7*n + 1
        ipar(9) = 6*n + 1
        ipar(10) = 9
        return
      endif
      if (stopbis(n,ipar,2,fpar,w(1,2),w(1,7),one))  goto 900
  100 if (ipar(3).eq.999.and.ipar(11).eq.1) goto 900
c
c***  Step (8): computing new p and rho
c      
      rho = fpar(7)
      fpar(7) = vecdo(n,w(1,2),w(1,1))
      omega = fpar(9)
      beta = fpar(7) * fpar(8) / (fpar(9) * rho)
      do i = 1, n
        w(i,6) = w(i,2) + beta * (w(i,6) - omega * w(i,5))
      enddo
      fpar(11) = fpar(11) + 6 * n + 3
      if (brkdn(fpar(7),ipar)) goto 900
c
c***  End of an iteration
c
      goto 30
c
c***  Some clean up job to do
c
  900 if (rp) then
        ipar(1) = 5
        ipar(8) = 7*n + 1
        ipar(9) = ipar(8) - n
        ipar(10) = 10
        return
      endif
  110 if (rp) then
        call tidycg(n,ipar,fpar,sol,w(1,7))
      else
        call tidycg(n,ipar,fpar,sol,w(1,8))
      endif

      return
      end

      subroutine bisinit(ipar,fpar,wksize,dsc,lp,rp,wk)
c***************************************************************************
c
c**** Some common initializations for the iterative solvers
c      
c***************************************************************************
      implicit none
      integer i,ipar(16),wksize,dsc
      logical lp,rp
      real*8  fpar(16),wk(*)
      real*8 zero, one
      parameter(zero=0.0D0, one=1.0D0)
c
c     ipar(1) = -2 indicate that there are not enough space in the work
c     array
c
      if (ipar(4).lt.wksize) then
         ipar(1) = -2
         ipar(4) = wksize
         return
      endif
c
      if (ipar(2).gt.2) then
         lp = .true.
         rp = .true.
      else if (ipar(2).eq.2) then
         lp = .false.
         rp = .true.
      else if (ipar(2).eq.1) then
         lp = .true.
         rp = .false.
      else
         lp = .false.
         rp = .false.
      endif
      if (ipar(3).eq.0) ipar(3) = dsc
c     .. clear the ipar elements used
      ipar(7) = 0
      ipar(8) = 0
      ipar(9) = 0
      ipar(10) = 0
      ipar(11) = 0
      ipar(12) = 0
c
c     fpar(1) must be between (0, 1), fpar(2) must be positive,
c     fpar(1) and fpar(2) can NOT both be zero
c     Normally return ipar(1) = -4 to indicate any of above error
c
      if (fpar(1).lt.zero .or. fpar(1).ge.one .or. fpar(2).lt.zero .or.
     &     (fpar(1).eq.zero .and. fpar(2).eq.zero)) then
         if (ipar(1).eq.0) then
            ipar(1) = -4
            return
         else
            fpar(1) = 1.0D-6
            fpar(2) = 1.0D-16
         endif
      endif
c     .. clear the fpar elements
      do i = 3, 10
         fpar(i) = zero
      enddo
      if (fpar(11).lt.zero) fpar(11) = zero
c     .. clear the used portion of the work array to zero
      do i = 1, wksize
         wk(i) = zero
      enddo
c
      return
      end
      logical function brkdn(alpha, ipar)
c***************************************************************************
c
c     test whether alpha is zero or an abnormal number, if yes,
c     this routine will return .true.
c
c     If alpha == 0, ipar(1) = -3,
c     if alpha is an abnormal number, ipar(1) = -9.
c
c***************************************************************************
      implicit none
      integer ipar(16)
      real*8 alpha, beta, zero, one
      parameter (zero=0.0D0, one=1.0D0)
      brkdn = .false.
      if (alpha.gt.zero) then
         beta = one / alpha
         if (.not. beta.gt.zero) then
            brkdn = .true.
            ipar(1) = -9
         endif
      else if (alpha.lt.zero) then
         beta = one / alpha
         if (.not. beta.lt.zero) then
            brkdn = .true.
            ipar(1) = -9
         endif
      else if (alpha.eq.zero) then
         brkdn = .true.
         ipar(1) = -3
      else
         brkdn = .true.
         ipar(1) = -9
      endif
      return
      end
      subroutine cg(n, rhs, sol, ipar, fpar, w)
c***************************************************************************
c
c     This is a implementation of the Conjugate Gradient (CG) method
c     for solving linear system.
c
c     NOTE: This is not the PCG algorithm. It is a regular CG algorithm.
c     To be consistent with the other solvers, the preconditioners are
c     applied by performing Ml^{-1} A Mr^{-1} P in place of A P in the
c     CG algorithm. The PCG uses its preconditioners very differently.
c
c     fpar(7) is used here internally to store <r, r>.
c     w(:,1) -- residual vector
c     w(:,2) -- P, the conjugate direction
c     w(:,3) -- A P, matrix multiply the conjugate direction
c     w(:,4) -- temporary storage for results of preconditioning
c     w(:,5) -- change in the solution (sol) is stored here until
c               termination of this solver
c
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
      real*8 vecdo
      logical stopbis, brkdn
      external vecdo, stopbis, brkdn, bisinit
      integer i
      real*8 alpha
      logical lp,rp
      save
c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80), ipar(10)
c
c     initialization
c
      call bisinit(ipar,fpar,5*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
      ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = ipar(8) + n
      ipar(10) = 1
      do i = 1, n
         w(i,2) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,2) = rhs(i) - w(i,3)
      enddo
      fpar(11) = fpar(11) + n
c
c     if left preconditioned
c
      if (lp) then
         ipar(1) = 3
         ipar(9) = 1
         ipar(10) = 2
         return
      endif
c
 20   if (lp) then
         do i = 1, n
            w(i,2) = w(i,1)
         enddo
      else
         do i = 1, n
            w(i,1) = w(i,2)
         enddo
      endif
c
      fpar(7) = vecdo(n,w,w)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(fpar(7))
      fpar(5) = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
c
c     before iteration can continue, we need to compute A * p, which
c     includes the preconditioning operations
c
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = n + 1
         if (lp) then
            ipar(9) = ipar(8) + n
         else
            ipar(9) = 3*n + 1
         endif
         ipar(10) = 3
         return
      endif
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + 1
      endif
      if (lp) then
         ipar(9) = 3*n+1
      else
         ipar(9) = n+n+1
      endif
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = n+n+1
         ipar(10) = 5
         return
      endif
c
c     continuing with the iterations
c
 60   ipar(7) = ipar(7) + 1
      alpha = vecdo(n,w(1,2),w(1,3))
      fpar(11) = fpar(11) + 2*n
      if (brkdn(alpha,ipar)) goto 900
      alpha = fpar(7) / alpha
      do i = 1, n
         w(i,5) = w(i,5) + alpha * w(i,2)
         w(i,1) = w(i,1) - alpha * w(i,3)
      enddo
      fpar(11) = fpar(11) + 4*n
c
c     are we ready to terminate ?
c
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 4*n + 1
         ipar(9) = 3*n + 1
         ipar(10) = 6
         return
      endif
 70   if (ipar(3).eq.999) then
         if (ipar(11).eq.1) goto 900
      else if (stopbis(n,ipar,1,fpar,w,w(1,2),alpha)) then
         goto 900
      endif
c
c     continue the iterations
c
      alpha = fpar(5)*fpar(5) / fpar(7)
      fpar(7) = fpar(5)*fpar(5)
      do i = 1, n
         w(i,2) = w(i,1) + alpha * w(i,2)
      enddo
      fpar(11) = fpar(11) + 2*n
      goto 30
c
c     clean up -- necessary to accommodate the right-preconditioning
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 7
         return
      endif
 80   if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,4))
      else
         call tidycg(n,ipar,fpar,sol,w(1,5))
      endif
c
      return
      end
      subroutine cgnr(n,rhs,sol,ipar,fpar,wk)
c***************************************************************************
c
c     CGNR -- Using CG algorithm solving A x = b by solving
c     Normal Residual equation: A^T A x = A^T b
c     As long as the matrix is not singular, A^T A is symmetric
c     positive definite, therefore CG (CGNR) will converge.
c
c     Usage of the work space:
c     wk(:,1) == residual vector R
c     wk(:,2) == the conjugate direction vector P
c     wk(:,3) == a scratch vector holds A P, or A^T R
c     wk(:,4) == a scratch vector holds intermediate results of the
c                preconditioning
c     wk(:,5) == a place to hold the modification to SOL
c
c     size of the work space WK is required = 5*n
c
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n),sol(n),fpar(16),wk(n,*)
      real*8 vecdo
      logical stopbis, brkdn
      external vecdo, stopbis, brkdn, bisinit
c
c     local variables
c
      integer i
      real*8 alpha, zz, zzm1
      logical lp, rp
      save
c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
c
c     initialization
c
      call bisinit(ipar,fpar,5*n,1,lp,rp,wk)
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      ipar(10) = 1
      do i = 1, n
         wk(i,1) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      do i = 1, n
         wk(i,1) = rhs(i) - wk(i,2)
      enddo
      fpar(11) = fpar(11) + n
c
c     if left preconditioned, precondition the initial residual
c
      if (lp) then
         ipar(1) = 3
         ipar(10) = 2
         return
      endif
c
 20   if (lp) then
         do i = 1, n
            wk(i,1) = wk(i,2)
         enddo
      endif
c
      zz = vecdo(n,wk,wk)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(zz)
      fpar(5) = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
c
c     normal iteration begins here, first half of the iteration
c     computes the conjugate direction
c
 30   continue
c
c     request the caller to perform a A^T r --> wk(:,3)
c
      if (lp) then
         ipar(1) = 4
         ipar(8) = 1
         if (rp) then
            ipar(9) = n + n + 1
         else
            ipar(9) = 3*n + 1
         endif
         ipar(10) = 3
         return
      endif
c
 40   ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 1
      endif
      if (rp) then
         ipar(9) = 3*n + 1
      else
         ipar(9) = n + n + 1
      endif
      ipar(10) = 4
      return
c
 50   if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = n + n + 1
         ipar(10) = 5
         return
      endif
c
 60   ipar(7) = ipar(7) + 1
      zzm1 = zz
      zz = vecdo(n,wk(1,3),wk(1,3))
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(zz,ipar)) goto 900
      if (ipar(7).gt.3) then
         alpha = zz / zzm1
         do i = 1, n
            wk(i,2) = wk(i,3) + alpha * wk(i,2)
         enddo
         fpar(11) = fpar(11) + 2 * n
      else
         do i = 1, n
            wk(i,2) = wk(i,3)
         enddo
      endif
c
c     before iteration can continue, we need to compute A * p
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = n + 1
         if (lp) then
            ipar(9) = ipar(8) + n
         else
            ipar(9) = 3*n + 1
         endif
         ipar(10) = 6
         return
      endif
c
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + 1
      endif
      if (lp) then
        ipar(9) = 3*n+1
      else
         ipar(9) = n+n+1
      endif
      ipar(10) = 7
      return
c
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = n+n+1
         ipar(10) = 8
         return
      endif
c
c     update the solution -- accumulate the changes in w(:,5)
c
 90   ipar(7) = ipar(7) + 1
      alpha = vecdo(n,wk(1,3),wk(1,3))
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha,ipar)) goto 900
      alpha = zz / alpha
      do i = 1, n
         wk(i,5) = wk(i,5) + alpha * wk(i,2)
         wk(i,1) = wk(i,1) - alpha * wk(i,3)
      enddo
      fpar(11) = fpar(11) + 4 * n
c
c     are we ready to terminate ?
c
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 4*n + 1
         ipar(9) = 3*n + 1
         ipar(10) = 9
         return
      endif
 100  if (ipar(3).eq.999) then
         if (ipar(11).eq.1) goto 900
      else if (stopbis(n,ipar,1,fpar,wk,wk(1,2),alpha)) then
         goto 900
      endif
c
c     continue the iterations
c
      goto 30
c
c     clean up -- necessary to accommodate the right-preconditioning
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call tidycg(n,ipar,fpar,sol,wk(1,4))
      else
         call tidycg(n,ipar,fpar,sol,wk(1,5))
      endif
      return
      end
      subroutine dbcg (n,rhs,sol,ipar,fpar,w)
c***************************************************************************
c
c Quasi GMRES method for solving a linear
c system of equations a * sol = y.  double precision version.
c this version is without restarting and without preconditioning.
c parameters :
c -----------
c n     = dimension of the problem
c
c y     = w(:,1) a temporary storage used for various operations
c z     = w(:,2) a work vector of length n.
c v     = w(:,3:4) size n x 2
c w     = w(:,5:6) size n x 2
c p     = w(:,7:9) work array of dimension n x 3
c del x = w(:,10)  accumulation of the changes in solution
c tmp   = w(:,11)  a temporary vector used to hold intermediate result of
c                  preconditioning, etc.
c
c sol   = the solution of the problem . at input sol must contain an
c         initial guess to the solution.
c    ***  note:   y is destroyed on return.
c
c-----------------------------------------------------------------------
c subroutines and functions called:
c 1) matrix vector multiplication and preconditioning through reverse
c     communication
c
c 2) implu, uppdir, vecdo
c-----------------------------------------------------------------------
c aug. 1983  version.    author youcef saad. yale university computer
c science dept. some  changes made july 3, 1986.
c references: siam j. sci. stat. comp., vol. 5, pp. 203-228 (1984)
c
c***************************************************************************
      implicit none
      integer n,ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
      real*8 one,zero
      parameter(one=1.0D0,zero=0.0D0)
      real*8 t,sqrt,vecdo,ss,res,beta,ss1,delta,x,zeta,umm
      integer k,j,i,i2,ip2,ju,lb,lbm1,np,indp
      logical lp,rp,full, perm(3)
      real*8 ypiv(3),u(3),usav(3)
      external tidycg
      save
c
c     where to go
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (110, 120, 130, 140, 150, 160, 170, 180, 190, 200) ipar(10)
c
c     initialization, parameter checking, clear the work arrays
c
      call bisinit(ipar,fpar,11*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
      perm(1) = .false.
      perm(2) = .false.
      perm(3) = .false.
      usav(1) = zero
      usav(2) = zero
      usav(3) = zero
      ypiv(1) = zero
      ypiv(2) = zero
      ypiv(3) = zero
c-----------------------------------------------------------------------
c     initialize constants for outer loop :
c-----------------------------------------------------------------------
      lb = 3
      lbm1 = 2
c
c     get initial residual vector and norm
c
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      do i = 1, n
         w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
 110  ipar(7) = ipar(7) + 1
      if (lp) then
         do i = 1, n
            w(i,1) = rhs(i) - w(i,2)
         enddo
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = n+n+1
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i,3) = rhs(i) - w(i,2)
         enddo
      endif
      fpar(11) = fpar(11) + n
c
 120  fpar(3) = sqrt(vecdo(n,w(1,3),w(1,3)))
      fpar(11) = fpar(11) + n + n
      fpar(5) = fpar(3)
      fpar(7) = fpar(3)
      zeta = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
         fpar(11) = fpar(11) + 2*n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * zeta + fpar(2)
      endif
      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
         fpar(6) = fpar(5)
         goto 900
      endif
c
c     normalize first arnoldi vector
c
      t = one/zeta
      do 22 k=1,n
         w(k,3) = w(k,3)*t
         w(k,5) = w(k,3)
 22   continue
      fpar(11) = fpar(11) + n
c
c     initialize constants for main loop
c
      beta = zero
      delta = zero
      i2 = 1
      indp = 0
      i = 0
c
c     main loop: i = index of the loop.
c
c-----------------------------------------------------------------------
 30   i = i + 1
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = (1+i2)*n+1
         if (lp) then
            ipar(9) = 1
         else
            ipar(9) = 10*n + 1
         endif
         ipar(10) = 3
         return
      endif
c
 130  ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (1+i2)*n + 1
      endif
      if (lp) then
         ipar(9) = 10*n + 1
      else
         ipar(9) = 1
      endif
      ipar(10) = 4
      return
c
 140  if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 1
         ipar(10) = 5
         return
      endif
c
c     A^t * x
c
 150  ipar(7) = ipar(7) + 1
      if (lp) then
         ipar(1) = 4
         ipar(8) = (3+i2)*n + 1
         if (rp) then
            ipar(9) = n + 1
         else
            ipar(9) = 10*n + 1
         endif
         ipar(10) = 6
         return
      endif
c
 160  ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (3+i2)*n + 1
      endif
      if (rp) then
         ipar(9) = 10*n + 1
      else
         ipar(9) = n + 1
      endif
      ipar(10) = 7
      return
c
 170  if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = n + 1
         ipar(10) = 8
         return
      endif
c-----------------------------------------------------------------------
c     orthogonalize current v against previous v's and
c     determine relevant part of i-th column of u(.,.) the
c     upper triangular matrix --
c-----------------------------------------------------------------------
 180  ipar(7) = ipar(7) + 1
      u(1) = zero
      ju = 1
      k = i2
      if (i .le. lbm1) ju = 0
      if (i .lt. lb) k = 0
 31   if (k .eq. lbm1) k=0
      k=k+1
c
      if (k .ne. i2) then
         ss  = delta
         ss1 = beta
         ju = ju + 1
         u(ju) = ss
      else
         ss = vecdo(n,w(1,1),w(1,4+k))
         fpar(11) = fpar(11) + 2*n
         ss1= ss
         ju = ju + 1
         u(ju) = ss
      endif
c
      do 32  j=1,n
         w(j,1) = w(j,1) - ss*w(j,k+2)
         w(j,2) = w(j,2) - ss1*w(j,k+4)
 32   continue
      fpar(11) = fpar(11) + 4*n
c
      if (k .ne. i2) goto 31
c
c     end of Mod. Gram. Schmidt loop
c
      t = vecdo(n,w(1,2),w(1,1))
c
      beta   = sqrt(abs(t))
      delta  = t/beta
c
      ss = one/beta
      ss1 = one/ delta
c
c     normalize and insert new vectors
c
      ip2 = i2
      if (i2 .eq. lbm1) i2=0
      i2=i2+1
c
      do 315 j=1,n
         w(j,i2+2)=w(j,1)*ss
         w(j,i2+4)=w(j,2)*ss1
 315  continue
      fpar(11) = fpar(11) + 4*n
c-----------------------------------------------------------------------
c     end of orthogonalization.
c     now compute the coefficients u(k) of the last
c     column of the  l . u  factorization of h .
c-----------------------------------------------------------------------
      np = min0(i,lb)
      full = (i .ge. lb)
      call implu(np, umm, beta, ypiv, u, perm, full)
c-----------------------------------------------------------------------
c     update conjugate directions and solution
c-----------------------------------------------------------------------
      do 33 k=1,n
         w(k,1) = w(k,ip2+2)
 33   continue
      call uppdir(n, w(1,7), np, lb, indp, w, u, usav, fpar(11))
c-----------------------------------------------------------------------
      if (i .eq. 1) goto 34
      j = np - 1
      if (full) j = j-1
      if (.not.perm(j)) zeta = -zeta*ypiv(j)
 34   x = zeta/u(np)
      if (perm(np))goto 36
      do 35 k=1,n
         w(k,10) = w(k,10) + x*w(k,1)
 35   continue
      fpar(11) = fpar(11) + 2 * n
c-----------------------------------------------------------------------
 36   if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 9*n + 1
         ipar(9) = 10*n + 1
         ipar(10) = 9
         return
      endif
      res = abs(beta*zeta/umm)
      fpar(5) = res * sqrt(vecdo(n, w(1,i2+2), w(1,i2+2)))
      fpar(11) = fpar(11) + 2 * n
      if (ipar(3).lt.0) then
         fpar(6) = x * sqrt(vecdo(n,w,w))
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).le.3) then
            fpar(3) = fpar(6)
            if (ipar(3).eq.-1) then
               fpar(4) = fpar(1) * sqrt(fpar(3)) + fpar(2)
            endif
         endif
      else
         fpar(6) = fpar(5)
      endif
c---- convergence test -----------------------------------------------
 190  if (ipar(3).eq.999.and.ipar(11).eq.0) then
         goto 30
      else if (fpar(6).gt.fpar(4) .and. (ipar(6).gt.ipar(7) .or.
     +        ipar(6).le.0)) then
         goto 30
      endif
c-----------------------------------------------------------------------
c     here the fact that the last step is different is accounted for.
c-----------------------------------------------------------------------
      if (.not. perm(np)) goto 900
      x = zeta/umm
      do 40 k = 1,n
         w(k,10) = w(k,10) + x*w(k,1)
 40   continue
      fpar(11) = fpar(11) + 2 * n
c
c     right preconditioning and clean-up jobs
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 9*n + 1
         ipar(9) = ipar(8) + n
         ipar(10) = 10
         return
      endif
 200  if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,11))
      else
         call tidycg(n,ipar,fpar,sol,w(1,10))
      endif
      return
      end
      subroutine diapre(
     .  diago,unkno,amatr,roffs,amcol,rhsid,nrows,
     .  ioutp,itask)
c****************************************************************************
c
c**** For ITASK = 1, this routine performs a diagonal scaling of a matrix
c**** A (stored in AMATR) as  A' = D^-1/2 A D^-1/2. The RHS and the initial
c**** guess are also modified. For ITASK = 2, the real unkonwn x = D^-1/2 x'
c**** is computed.
c
c****************************************************************************
      implicit none
      integer
     .  nrows,irows,izeco,jcomp,izecf,jzeco,ioutp,itask
      integer
     .  amcol(*), roffs(nrows)
      real*8
     .  diago(nrows), amatr(*), unkno(nrows), rhsid(nrows)
      real*8
     .  diagt,diaga
c
c***  Obtain the diagonal of the matrix
c
      if(itask.eq.1) then
        do irows = 1,nrows
          izeco = roffs(irows)
          jcomp = amcol(izeco)
          do while(jcomp.ne.irows)
            izeco=izeco+1
            jcomp=amcol(izeco)
          end do
          diagt = amatr(izeco)
          diaga = abs(diagt)
          if(diaga.lt.1.0e-10) then
            call srunei(
     .        'DIAPRE: ZERO DIAGONAL FOR EQUATION',irows,ioutp)
          else
            diago(irows) = 1.0/sqrt(diaga)
          end if
        end do
c
c***  Scaling of A and the RHS and modification of the intial guess
c
        do irows = 1,nrows
          izeco = roffs(irows)
          izecf = roffs(irows+1) - 1
          diagt = diago(irows)
          do jzeco = izeco,izecf
            jcomp = amcol(jzeco)
            amatr(jzeco) = amatr(jzeco)*diagt*diago(jcomp)
          end do
          rhsid(irows) = rhsid(irows)*diagt
          unkno(irows) = unkno(irows)/diagt
        end do
      end if
c
c***  Recovering of the unknown and restoring of A
c
      if(itask.eq.2) then
        do irows = 1,nrows
          izeco = roffs(irows)
          izecf = roffs(irows+1) - 1
          diagt = diago(irows)
          do jzeco = izeco,izecf
            jcomp = amcol(jzeco)
            amatr(jzeco) = amatr(jzeco)/(diagt*diago(jcomp))
          end do
          unkno(irows) = unkno(irows)*diagt
        end do
      end if

      end
      subroutine dqgmres(n, rhs, sol, ipar, fpar, w)
c***************************************************************************
c
c     DQGMRES -- Flexible Direct version of Quasi-General Minimum
c     Residual method. The right preconditioning can be varied from
c     step to step.
c
c     Work space used = n + lb * (2*n+4)
c     where lb = ipar(5) + 1 (default 16 if ipar(5) <= 1)
c
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
      real*8 one,zero,deps
      parameter(one=1.0D0,zero=0.0D0)
      parameter(deps=1.0D-33)
c
      integer i,ii,j,jp1,j0,k,ptrw,ptrv,iv,iw,ic,is,ihm,ihd,lb,ptr
      real*8 alpha,beta,psi,c,s,vecdo
      logical lp,rp,full
      external vecdo,bisinit
      save
c
c     where to go
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70) ipar(10)
c
c     locations of the work arrays. The arrangement is as follows:
c     w(1:n) -- temporary storage for the results of the preconditioning
c     w(iv+1:iw) -- the V's
c     w(iw+1:ic) -- the W's
c     w(ic+1:is) -- the COSINEs of the Givens rotations
c     w(is+1:ihm) -- the SINEs of the Givens rotations
c     w(ihm+1:ihd) -- the last column of the Hessenberg matrix
c     w(ihd+1:i) -- the inverse of the diagonals of the Hessenberg matrix
c
      if (ipar(5).le.1) then
         lb = 16
      else
         lb = ipar(5) + 1
      endif
      iv = n
      iw = iv + lb * n
      ic = iw + lb * n
      is = ic + lb
      ihm = is + lb
      ihd = ihm + lb
      i = ihd + lb
c
c     parameter check, initializations
c
      full = .false.
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
      ipar(1) = 1
      if (lp) then
         do ii = 1, n
            w(iv+ii) = sol(ii)
         enddo
         ipar(8) = iv+1
         ipar(9) = 1
      else
         do ii = 1, n
            w(ii) = sol(ii)
         enddo
         ipar(8) = 1
         ipar(9) = iv+1
      endif
      ipar(10) = 1
      return
c
 10   ipar(7) = ipar(7) + 1
      if (lp) then
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = iv+1
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(iv+i) = rhs(i) - w(iv+i)
         enddo
      endif
      fpar(11) = fpar(11) + n
c
 20   alpha = sqrt(vecdo(n, w(iv+1), w(iv+1)))
      fpar(11) = fpar(11) + (n + n)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
         fpar(11) = fpar(11) + 2*n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * alpha + fpar(2)
      endif
      fpar(3) = alpha
      fpar(5) = alpha
      psi = alpha
      if (alpha.le.fpar(4)) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 80
      endif
      alpha = one / alpha
      do i = 1, n
         w(iv+i) = w(iv+i) * alpha
      enddo
      fpar(11) = fpar(11) + n
      j = 0
c
c     iterations start here
c
 30   j = j + 1
      if (j.gt.lb) j = j - lb
      jp1 = j + 1
      if (jp1.gt.lb) jp1 = jp1 - lb
      ptrv = iv + (j-1)*n + 1
      ptrw = iv + (jp1-1)*n + 1
      if (.not.full) then
         if (j.gt.jp1) full = .true.
      endif
      if (full) then
         j0 = jp1+1
         if (j0.gt.lb) j0 = j0 - lb
      else
         j0 = 1
      endif
c
c     request the caller to perform matrix-vector multiplication and
c     preconditioning
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = ptrv
         ipar(9) = ptrv + iw - iv
         ipar(10) = 3
         return
      else
         do i = 0, n-1
            w(ptrv+iw-iv+i) = w(ptrv+i)
         enddo
      endif
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = ptrv
      endif
      if (lp) then
         ipar(9) = 1
      else
         ipar(9) = ptrw
      endif
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = ptrw
         ipar(10) = 5
         return
      endif
c
c     compute the last column of the Hessenberg matrix
c     modified Gram-schmidt procedure, orthogonalize against (lb-1)
c     previous vectors
c
 60   continue
      call mgsro(full,n,n,lb,jp1,fpar(11),w(iv+1),w(ihm+1),
     $     ipar(12))
      if (ipar(12).lt.0) then
         ipar(1) = -3
         goto 80
      endif
      beta = w(ihm+jp1)
c
c     incomplete factorization (QR factorization through Givens rotations)
c     (1) apply previous rotations [(lb-1) of them]
c     (2) generate a new rotation
c
      if (full) then
         w(ihm+jp1) = w(ihm+j0) * w(is+jp1)
         w(ihm+j0) = w(ihm+j0) * w(ic+jp1)
      endif
      i = j0
      do while (i.ne.j)
         k = i+1
         if (k.gt.lb) k = k - lb
         c = w(ic+i)
         s = w(is+i)
         alpha = w(ihm+i)
         w(ihm+i) = c * alpha + s * w(ihm+k)
         w(ihm+k) = c * w(ihm+k) - s * alpha
         i = k
      enddo
      call givens(w(ihm+j), beta, c, s)
      if (full) then
         fpar(11) = fpar(11) + 6 * lb
      else
         fpar(11) = fpar(11) + 6 * j
      endif
c
c     detect whether diagonal element of this column is zero
c
      if (abs(w(ihm+j)).lt.deps) then
         ipar(1) = -3
         goto 80
      endif
      w(ihd+j) = one / w(ihm+j)
      w(ic+j) = c
      w(is+j) = s
c
c     update the W's (the conjugate directions) -- essentially this is one
c     step of triangular solve.
c
      ptrw = iw+(j-1)*n + 1
      if (full) then
         do i = j+1, lb
            alpha = -w(ihm+i)*w(ihd+i)
            ptr = iw+(i-1)*n+1
            do ii = 0, n-1
               w(ptrw+ii) = w(ptrw+ii) + alpha * w(ptr+ii)
            enddo
         enddo
      endif
      do i = 1, j-1
         alpha = -w(ihm+i)*w(ihd+i)
         ptr = iw+(i-1)*n+1
         do ii = 0, n-1
            w(ptrw+ii) = w(ptrw+ii) + alpha * w(ptr+ii)
         enddo
      enddo
c
c     update the solution to the linear system
c
      alpha = psi * c * w(ihd+j)
      psi = - s * psi
      do i = 1, n
         sol(i) = sol(i) + alpha * w(ptrw-1+i)
      enddo
      if (full) then
         fpar(11) = fpar(11) + lb * (n+n)
      else
         fpar(11) = fpar(11) + j * (n+n)
      endif
c
c     determine whether to continue,
c     compute the desired error/residual norm
c
      ipar(7) = ipar(7) + 1
      fpar(5) = abs(psi)
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = 1
         ipar(10) = 6
         return
      endif
      if (ipar(3).lt.0) then
         alpha = abs(alpha)
         if (ipar(7).eq.2 .and. ipar(3).eq.-1) then
            fpar(3) = alpha*sqrt(vecdo(n, w(ptrw), w(ptrw)))
            fpar(4) = fpar(1) * fpar(3) + fpar(2)
            fpar(6) = fpar(3)
         else
            fpar(6) = alpha*sqrt(vecdo(n, w(ptrw), w(ptrw)))
         endif
         fpar(11) = fpar(11) + 2 * n
      else
         fpar(6) = fpar(5)
      endif
      if (ipar(1).ge.0 .and. fpar(6).gt.fpar(4) .and. (ipar(6).le.0
     +     .or. ipar(7).lt.ipar(6))) goto 30
 70   if (ipar(3).eq.999 .and. ipar(11).eq.0) goto 30
c
c     clean up the iterative solver
c
 80   fpar(7) = zero
      if (fpar(3).ne.zero .and. fpar(6).ne.zero)
     +     fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).ne.0) then
            ipar(1) = 0
         else if (fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
      return
      end
      subroutine fgmres(n, rhs, sol, ipar, fpar, w)
c***************************************************************************
c
c     This a version of FGMRES implemented with reverse communication.
c
c     ipar(5) == the dimension of the Krylov subspace
c
c     the space of the `w' is used as follows:
c     >> V: the bases for the Krylov subspace, size n*(m+1);
c     >> W: the above bases after (left-)multiplying with the
c     right-preconditioner inverse, size m*n;
c     >> a temporary vector of size n;
c     >> the Hessenberg matrix, only the upper triangular portion
c     of the matrix is stored, size (m+1)*m/2 + 1
c     >> three vectors, first two are of size m, they are the cosine
c     and sine of the Givens rotations, the third one holds the
c     residuals, it is of size m+1.
c
c     TOTAL SIZE REQUIRED == n*(2m+1) + (m+1)*m/2 + 3*m + 2
c     Note: m == ipar(5). The default value for this is 15 if
c     ipar(5) <= 1.
c
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
      real*8 vecdo
      external vecdo
c
      real*8 one, zero
      parameter(one=1.0D0, zero=0.0D0)
c
c     local variables, ptr and p2 are temporary pointers,
c     hess points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
c
      integer i,ii,idx,iz,k,m,ptr,p2,hess,vc,vs,vrn
      real*8 alpha, c, s
      logical lp, rp
      save
c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60) ipar(10)
c
c     initialization
c
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      endif
      idx = n * (m+1)
      iz = idx + n
      hess = iz + n*m
      vc = hess + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do ii = 1, n
         w(ii+n) = sol(ii)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      fpar(11) = fpar(11) + n
      if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
      endif
c
 20   alpha = sqrt(vecdo(n,w,w))
      fpar(11) = fpar(11) + n + n
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         endif
         fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      endif
      alpha = one / alpha
      do ii = 1, n
         w(ii) = w(ii) * alpha
      enddo
      fpar(11) = fpar(11) + n
c
c     request for (1) right preconditioning
c     (2) matrix vector multiplication
c     (3) left preconditioning
c
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         ipar(9) = iz + ipar(8)
         ipar(10) = 3
         return
      else
         do ii = 0, n-1
            w(iz+k*n-ii) = w(k*n-ii)
         enddo
      endif
c
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      endif
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      endif
      ipar(10) = 4
      return
c
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      endif
c
c     Modified Gram-Schmidt orthogonalization procedure
c     temporary pointer 'ptr' is pointing to the current column of the
c     Hessenberg matrix. 'p2' points to the new basis vector
c
 50   ptr = k * (k - 1) / 2 + hess
      p2 = ipar(9)
      ipar(7) = ipar(7) + 1
      call mgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),
     $     ipar(12))
      if (ipar(12).lt.0) goto 200
c
c     apply previous Givens rotations and generate a new one to eliminate
c     the subdiagonal element.
c
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      enddo
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
      fpar(11) = fpar(11) + 6 * k
c
c     end of one Arnoldi iteration, alpha will store the estimated
c     residual norm at current stage
c
      alpha = abs(alpha)
      fpar(5) = alpha
      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4))
     +      .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
c
c     update the approximate solution, first solve the upper triangular
c     system, temporary pointer ptr points to the Hessenberg matrix,
c     p2 points to the right-hand-side (also the solution) of the system.
c
 200  ptr = hess + k * (k + 1 ) / 2
      p2 = vrn + k
      if (w(ptr).eq.zero) then
c
c     if the diagonal elements of the last column is zero, reduce k by 1
c     so that a smaller trianguler system is solved [It should only
c     happen when the matrix is singular!]
c
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            goto 300
         endif
      endif
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
         enddo
         p2 = p2 - 1
         w(p2) = w(p2) / w(ptr)
      enddo
c
      do i = 0, k-1
         ptr = iz+i*n
         do ii = 1, n
            sol(ii) = sol(ii) + w(p2)*w(ptr+ii)
         enddo
         p2 = p2 + 1
      enddo
      fpar(11) = fpar(11) + 2*k*n + k*(k+1)
c
c     process the complete stopping criteria
c
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = idx + 1
         ipar(10) = 6
         return
      else if (ipar(3).lt.0) then
         if (ipar(7).le.m+1) then
            fpar(3) = abs(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         endif
         alpha = abs(w(vrn+k))
      else if (ipar(3).ne.999) then
         fpar(6) = alpha
      endif
c
c     do we need to restart ?
c
 60   if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0).and.
     +     ((ipar(3).eq.999.and.ipar(11).eq.0) .or.
     +     (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
c
c     termination, set error code, compute convergence rate
c
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. alpha.le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero) then
         fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
      else
         fpar(7) = zero
      endif
      return
      end
      subroutine fom(n, rhs, sol, ipar, fpar, w)
c***************************************************************************
c
c     This a version of The Full Orthogonalization Method (FOM) 
c     implemented with reverse communication. It is a simple restart 
c     version of the FOM algorithm and is implemented with plane 
c     rotations similarly to GMRES.
c
c  parameters:
c  ----------- 
c     ipar(5) == the dimension of the Krylov subspace
c     after every ipar(5) iterations, the FOM will restart with
c     the updated solution and recomputed residual vector.
c
c     the work space in `w' is used as follows:
c     (1) the basis for the Krylov subspace, size n*(m+1);
c     (2) the Hessenberg matrix, only the upper triangular
c     portion of the matrix is stored, size (m+1)*m/2 + 1
c     (3) three vectors, all are of size m, they are
c     the cosine and sine of the Givens rotations, the third one holds
c     the residuals, it is of size m+1.
c
c     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
c     Note: m == ipar(5). The default value for this is 15 if
c     ipar(5) <= 1.
c
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
      real*8 vecdo
      external vecdo
c
      real*8 one, zero
      parameter(one=1.0D0, zero=0.0D0)
c
c     local variables, ptr and p2 are temporary pointers,
c     hes points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
c
      integer i,ii,idx,k,m,ptr,p2,prs,hes,vc,vs,vrn
      real*8 alpha, c, s 
      logical lp, rp
      save
c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
c
c     initialization
c
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      endif
      idx = n * (m+1)
      hes = idx + n
      vc = hes + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do i = 1, n
         w(n+i) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
c
 20   alpha = sqrt(vecdo(n,w,w))
      fpar(11) = fpar(11) + 2*n + 1
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         endif
         fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      endif
      alpha = one / alpha
      do ii = 1, n
         w(ii) = alpha * w(ii)
      enddo
      fpar(11) = fpar(11) + n
c
c     request for (1) right preconditioning
c     (2) matrix vector multiplication
c     (3) left preconditioning
c
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         if (lp) then
            ipar(9) = k*n + 1
         else
            ipar(9) = idx + 1
         endif
         ipar(10) = 3
         return
      endif
c
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      endif
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      endif
      ipar(10) = 4
      return
c
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      endif
c
c     Modified Gram-Schmidt orthogonalization procedure
c     temporary pointer 'ptr' is pointing to the current column of the
c     Hessenberg matrix. 'p2' points to the new basis vector
c
 50   ipar(7) = ipar(7) + 1
      ptr = k * (k - 1) / 2 + hes
      p2 = ipar(9)
      call mgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),
     $     ipar(12))
      if (ipar(12).lt.0) goto 200
c
c     apply previous Givens rotations to column.
c
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      enddo
c
c     end of one Arnoldi iteration, alpha will store the estimated
c     residual norm at current stage
c
      fpar(11) = fpar(11) + 6*k

      prs = vrn+k
      alpha = fpar(5) 
      if (w(p2) .ne. zero) alpha = abs(w(p2+1)*w(prs)/w(p2)) 
      fpar(5) = alpha
c
      if (k.ge.m .or. (ipar(3).ge.0 .and. alpha.le.fpar(4))
     +     .or. (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)))
     +     goto 200
c
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      alpha = - s * w(prs)
      w(prs) = c * w(prs)
      w(prs+1) = alpha
c
      if (w(p2).ne.zero) goto 110
c
c     update the approximate solution, first solve the upper triangular
c     system, temporary pointer ptr points to the Hessenberg matrix,
c     prs points to the right-hand-side (also the solution) of the system.
c
 200  ptr = hes + k * (k + 1) / 2
      prs = vrn + k
      if (w(ptr).eq.zero) then
c
c     if the diagonal elements of the last column is zero, reduce k by 1
c     so that a smaller trianguler system is solved
c
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            goto 300
         endif
      endif
      w(prs) = w(prs) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(prs) * w(ptr+ii)
         enddo
         prs = prs - 1
         w(prs) = w(prs) / w(ptr)
      enddo
c
      do ii = 1, n
         w(ii) = w(ii) * w(prs)
      enddo
      do i = 1, k-1
         prs = prs + 1
         ptr = i*n
         do ii = 1, n
            w(ii) = w(ii) + w(prs) * w(ptr+ii)
         enddo
      enddo
      fpar(11) = fpar(11) + 2*(k-1)*n + n + k*(k+1)
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 1
         ipar(9) = idx + 1
         ipar(10) = 6
         return
      endif
c
 60   if (rp) then
         do i = 1, n
            sol(i) = sol(i) + w(idx+i)
         enddo
      else
         do i = 1, n
            sol(i) = sol(i) + w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
c
c     process the complete stopping criteria
c
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = idx + 1
         ipar(10) = 7
         return
      else if (ipar(3).lt.0) then
         if (ipar(7).le.m+1) then
            fpar(3) = abs(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         endif
         alpha = abs(w(vrn+k))
      endif
      fpar(6) = alpha
c
c     do we need to restart ?
c
 70   if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
      if (ipar(7).lt.ipar(6) .or. ipar(6).le.0) then
         if (ipar(3).ne.999) then
            if (fpar(6).gt.fpar(4)) goto 100
         else
            if (ipar(11).eq.0) goto 100
         endif
      endif
c
c     termination, set error code, compute convergence rate
c
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. alpha.le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero) then
         fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
      else
         fpar(7) = zero
      endif
      return
      end
      subroutine givens(x,y,c,s)
      real*8 x,y,c,s
c-----------------------------------------------------------------------
c     Given x and y, this subroutine generates a Givens' rotation c, s.
c     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
c     (See P 202 of "matrix computation" by Golub and van Loan.)
c-----------------------------------------------------------------------
      real*8 t,one,zero
      parameter (zero=0.0D0,one=1.0D0)
c
      if (x.eq.zero .and. y.eq.zero) then
         c = one
         s = zero
      else if (abs(y).gt.abs(x)) then
         t = x / y
         x = sqrt(one+t*t)
         s = sign(one / x, y)
         c = t*s
      else if (abs(y).le.abs(x)) then
         t = y / x
         y = sqrt(one+t*t)
         c = sign(one / y, x)
         s = t*c
      else
c
c     X or Y must be an invalid floating-point number, set both to zero
c
         x = zero
         y = zero
         c = one
         s = zero
      endif
      x = abs(x*y)
c
c     end of givens
c
      return
      end
      subroutine gmres(n, rhs, sol, ipar, fpar, w)
c***************************************************************************
c
c     This a version of GMRES implemented with reverse communication.
c     It is a simple restart version of the GMRES algorithm.
c
c     ipar(5) == the dimension of the Krylov subspace
c     after every ipar(5) iterations, the GMRES will restart with
c     the updated solution and recomputed residual vector.
c
c     the space of the `w' is used as follows:
c     (1) the basis for the Krylov subspace, size n*(m+1);
c     (2) the Hessenberg matrix, only the upper triangular
c     portion of the matrix is stored, size (m+1)*m/2 + 1
c     (3) three vectors, all are of size m, they are
c     the cosine and sine of the Givens rotations, the third one holds
c     the residuals, it is of size m+1.
c
c     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
c     Note: m == ipar(5). The default value for this is 15 if
c     ipar(5) <= 1.
c
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
      real*8 vecdo
      external vecdo
c
      real*8 one, zero
      parameter(one=1.0D0, zero=0.0D0)
c
c     local variables, ptr and p2 are temporary pointers,
c     hess points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
c
      integer i,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn
      real*8 alpha, c, s
      logical lp, rp
      save
c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
c
c     initialization
c
      if (ipar(5).le.1) then
        m = 15
      else
        m = ipar(5)
      endif
      idx = n * (m+1)
      hess = idx + n
      vc = hess + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
  100 ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do i = 1, n
        w(n+i) = sol(i)
      enddo
      return
   10 ipar(7) = ipar(7) + 1
      if (lp) then
        do i = 1, n
          w(n+i) = rhs(i) - w(i)
        enddo
        ipar(1) = 3
        ipar(10) = 2
        return
      else
        do i = 1, n
          w(i) = rhs(i) - w(i)
        enddo
      endif
      fpar(11) = fpar(11) + n
c
   20 alpha = sqrt(vecdo(n,w,w))
      fpar(11) = fpar(11) + 2*n
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
        if (abs(ipar(3)).eq.2) then
          fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
          fpar(11) = fpar(11) + 2*n
        else
          fpar(4) = fpar(1) * alpha + fpar(2)
        endif
        fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
        ipar(1) = 0
        fpar(6) = alpha
        goto 300
      endif
      alpha = one / alpha
      do ii = 1, n
        w(ii) = alpha * w(ii)
      enddo
      fpar(11) = fpar(11) + n
c
c     request for (1) right preconditioning
c     (2) matrix vector multiplication
c     (3) left preconditioning
c
  110 k = k + 1
      if (rp) then
        ipar(1) = 5
        ipar(8) = k*n - n + 1
        if (lp) then
          ipar(9) = k*n + 1
        else
          ipar(9) = idx + 1
        endif
        ipar(10) = 3
        return
      endif
c
   30 ipar(1) = 1
      if (rp) then
        ipar(8) = ipar(9)
      else
        ipar(8) = (k-1)*n + 1
      endif
      if (lp) then
        ipar(9) = idx + 1
      else
        ipar(9) = 1 + k*n
      endif
      ipar(10) = 4
      return
c
   40 if (lp) then
        ipar(1) = 3
        ipar(8) = ipar(9)
        ipar(9) = k*n + 1
        ipar(10) = 5
        return
      endif
c
c     Modified Gram-Schmidt orthogonalization procedure
c     temporary pointer 'ptr' is pointing to the current column of the
c     Hessenberg matrix. 'p2' points to the new basis vector
c
   50 ipar(7) = ipar(7) + 1
      ptr = k * (k - 1) / 2 + hess
      p2 = ipar(9)
      call mgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),
     $  ipar(12))
      if (ipar(12).lt.0) goto 200
c
c     apply previous Givens rotations and generate a new one to eliminate
c     the subdiagonal element.
c
      p2 = ptr + 1
      do i = 1, k-1
        ptr = p2
        p2 = p2 + 1
        alpha = w(ptr)
        c = w(vc+i)
        s = w(vs+i)
        w(ptr) = c * alpha + s * w(p2)
        w(p2) = c * w(p2) - s * alpha
      enddo
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
c
c     end of one Arnoldi iteration, alpha will store the estimated
c     residual norm at current stage
c
      fpar(11) = fpar(11) + 6*k + 2
      alpha = abs(alpha)
      fpar(5) = alpha
      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4))
     +  .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
c
c     update the approximate solution, first solve the upper triangular
c     system, temporary pointer ptr points to the Hessenberg matrix,
c     p2 points to the right-hand-side (also the solution) of the system.
c
  200 ptr = hess + k * (k + 1) / 2
      p2 = vrn + k
      if (w(ptr).eq.zero) then
c
c     if the diagonal elements of the last column is zero, reduce k by 1
c     so that a smaller trianguler system is solved [It should only
c     happen when the matrix is singular, and at most once!]
c
        k = k - 1
        if (k.gt.0) then
          goto 200
        else
          ipar(1) = -3
          ipar(12) = -4
          goto 300
        endif
      endif
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
        ptr = ptr - i - 1
        do ii = 1, i
          w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
        enddo
        p2 = p2 - 1
        w(p2) = w(p2) / w(ptr)
      enddo
c
      do ii = 1, n
        w(ii) = w(ii) * w(p2)
      enddo
      do i = 1, k-1
        ptr = i*n
        p2 = p2 + 1
        do ii = 1, n
          w(ii) = w(ii) + w(p2) * w(ptr+ii)
        enddo
      enddo
      fpar(11) = fpar(11) + 2*k*n - n + k*(k+1)
c
      if (rp) then
        ipar(1) = 5
        ipar(8) = 1
        ipar(9) = idx + 1
        ipar(10) = 6
        return
      endif
c
   60 if (rp) then
        do i = 1, n
          sol(i) = sol(i) + w(idx+i)
        enddo
      else
        do i = 1, n
          sol(i) = sol(i) + w(i)
        enddo
      endif
      fpar(11) = fpar(11) + n
c
c     process the complete stopping criteria
c
      if (ipar(3).eq.999) then
        ipar(1) = 10
        ipar(8) = -1
        ipar(9) = idx + 1
        ipar(10) = 7
        return
      else if (ipar(3).lt.0) then
        if (ipar(7).le.m+1) then
          fpar(3) = abs(w(vrn+1))
          if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
        endif
        alpha = abs(w(vrn+k))
      else
        fpar(6) = alpha
      endif
c
c     do we need to restart ?
c
   70 if (ipar(12).ne.0) then
        ipar(1) = -3
        goto 300
      endif
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and.
     +  ((ipar(3).eq.999.and.ipar(11).eq.0) .or.
     +  (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
c
c     termination, set error code, compute convergence rate
c
      if (ipar(1).gt.0) then
        if (ipar(3).eq.999 .and. ipar(11).eq.1) then
          ipar(1) = 0
        else if (ipar(3).ne.999 .and. alpha.le.fpar(4)) then
          ipar(1) = 0
        else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
          ipar(1) = -1
        else
          ipar(1) = -10
        endif
      endif
  300 if (fpar(3).ne.zero .and. fpar(6).ne.zero) then
        fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
      else
        fpar(7) = zero
      endif
      return
      end
      subroutine ilut (n,a,ja,ia,lfil,tol,alu,jlu,ju,iwk,
     .  wu,wl,jr,jwl,jwu,ierr)
c*************************************************************************      
c
c                      *** ILUT preconditioner ***                     
c                      ---------------------------                     
c      incomplete LU factorization with dual truncation mechanism      
c      VERSION 2 : sorting  done for both L and U.                     
c
c---- Dual drop-off strategy works as follows.                         
c                                                                      
c     1) Thresholding in L and U as set by tol. Any element whose size
c        is less than some tolerance (relative to the norm of current  
c        row in u) is dropped.                                         
c                                                                      
c     2) Keeping only the largest lfil+il(i) elements in the i-th row  
c        of L and the largest lfil+iu(i) elements in the i-th row of   
c        U where il(i), iu(i) are the original number of nonzero       
c        elements of the L-part and the U-part of the i-th row of A    
c                                                                      
c Flexibility: one can use tol=0 to get a strategy based on keeping the
c largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n
c will give the usual threshold strategy (however, fill-in is then     
c impredictible).                                                      
c                                                                      
c
c PARAMETERS
c ----------
c
c on entry:
c =========
c n       = integer. The dimension of the matrix A.
c
c a,ja,ia = matrix stored in Compressed Sparse Row format.
c
c lfil    = integer. The fill-in parameter. Each row of L and
c           each row of U will have a maximum of lfil elements
c           in addition to their original number of nonzero elements.
c           Thus storage can be determined beforehand.
c           lfil must be .ge. 0.
c
c iwk     = integer. The minimum length of arrays alu and jlu
c
c On return:
c ==========
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero row encountered.
c
c work arrays:
c ============
c jr,jwu,jwl 	  = integer work arrays of length n.
c wu, wl          = real work arrays of length n+1, and n resp.
c
c Notes:
c ------
c A must have all nonzero diagonal elements.
c
c*************************************************************************      
      implicit real*8 (a-h,o-z)
      real*8 a(*), alu(*), wu(n+1), wl(n), tol
      integer ja(*),ia(n+1),jlu(*),ju(n),jr(n), jwu(n),
     *  jwl(n), n, lfil, iwk, ierr
      if (lfil .lt. 0) goto 998
c-------------------------------
c initialize ju0 (points to next element to be added to alu,jlu)
c and pointer.
c
      ju0 = n+2
      jlu(1) = ju0
c
c  integer double pointer array.
c
      do 1 j=1, n
        jr(j)  = 0
    1 continue
c-----------------------------------------------------------------------
c  beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
        j1 = ia(ii)
        j2 = ia(ii+1) - 1
        tnorm = 0.0d0
        do 501 k=j1,j2
          tnorm = tnorm+abs(a(k))
  501   continue
        if (tnorm .eq. 0.0) goto 999
        tnorm = tnorm/real(j2-j1+1)
c
c--- unpack L-part and U-part of row of A in arrays wl, wu --
c
        lenu = 1
        lenl = 0
        jwu(1) = ii
        wu(1) = 0.0
        jr(ii) = 1
c
        do 170  j = j1, j2
          k = ja(j)
          t = a(j)
          if (abs(t) .lt. tol*tnorm .and. k .ne. ii) goto 170
          if (k .lt. ii) then
            lenl = lenl+1
            jwl(lenl) = k
            wl(lenl) = t
            jr(k) = lenl
          else if (k .eq. ii) then
            wu(1) = t
          else
            lenu = lenu+1
            jwu(lenu) = k
            wu(lenu) = t
            jr(k) = lenu
          endif
  170   continue
        tnorm = tnorm/real(j2-j1+1)
        lenl0 = lenl
        lenu0 = lenu
        jj = 0
        nl = 0
c-------------------------------------------------------------------
c---------------------- eliminate previous rows --------------------
c-------------------------------------------------------------------
  150   jj = jj+1
        if (jj .gt. lenl) goto 160
c-------------------------------------------------------------------
c in order to do the elimination in the correct order we need to
c exchange the current row number with the one that has
c smallest column number, among jj,jj+1,...,lenl.
c-------------------------------------------------------------------
        jrow = jwl(jj)
        k = jj
c
c determine smallest column index
c
        do 151 j=jj+1,lenl
          if (jwl(j) .lt. jrow) then
            jrow = jwl(j)
            k = j
          endif
  151   continue
c
c exchange in jwl
c
        if (k .ne. jj) then
          j = jwl(jj)
          jwl(jj) = jwl(k)
          jwl(k) = j
c
c exchange in jr
c
          jr(jrow) = jj
          jr(j) = k
c
c exchange in wl
c
          s = wl(jj)
          wl(jj) = wl(k)
          wl(k) = s
        endif
c
        if (jrow .ge. ii) goto 160
c---------get the multiplier for row to be eliminated: jrow
        fact = wl(jj)*alu(jrow)
c zero out element in row by setting jr(jrow) = 0
        jr(jrow) = 0
        if (abs(fact)*wu(n+2-jrow) .le. tol*tnorm) goto 150
c-------------------------------------------------------------------
c------------ combine current row and row jrow ---------------------
c-------------------------------------------------------------------
        do 203 k = ju(jrow), jlu(jrow+1)-1
          s = fact*alu(k)
          j = jlu(k)
          jpos = jr(j)
c
c if fill-in element is small then disregard:
c
          if (abs(s) .lt. tol*tnorm .and. jpos .eq. 0) goto 203
          if (j .ge. ii) then
c
c     dealing with upper part.
c
            if (jpos .eq. 0) then
c     this is a fill-in element
              lenu = lenu+1
              if (lenu .gt. n) goto 995
              jwu(lenu) = j
              jr(j) = lenu
              wu(lenu) = - s
            else
c     no fill-in element --
              wu(jpos) = wu(jpos) - s
            endif
          else
c
c     dealing with lower part.
c
            if (jpos .eq. 0) then
c     this is a fill-in element
              lenl = lenl+1
              if (lenl .gt. n) goto 995
              jwl(lenl) = j
              jr(j) = lenl
              wl(lenl) = - s
            else
c     no fill-in element --
              wl(jpos) = wl(jpos) - s
            endif
          endif
  203   continue
        nl = nl+1
        wl(nl) = fact
        jwl(nl)  = jrow
	goto 150
c----------------------------------------------------------
c------------ update l-matrix -----------------------------
c----------------------------------------------------------
  160   len = min0(nl,lenl0+lfil)
c 160    len = min0(nl,lfil)

  	call qsplit (wl,jwl,nl,len)
c
        do 204 k=1, len
          if (ju0 .gt. iwk) goto 996
          alu(ju0) =  wl(k)
          jlu(ju0) =  jwl(k)
          ju0 = ju0+1
  204   continue
c
c  save pointer to beginning of row ii of U
c
        ju(ii) = ju0
c
c  reset double-pointer jr to zero (L-part - except first
c  jj-1 elements which have already been reset)
	do 306 k= jj, lenl
          jr(jwl(k)) = 0
  306   continue
        len = min0(lenu,lenu0+lfil+1)                       ! +1 added 11/23/94
c        len = min0(lenu,lfil)
	call qsplit (wu(2), jwu(2), lenu-1,len)
c----------------------------------------------------------
c------------ update u-matrix -----------------------------
c----------------------------------------------------------
        t = abs(wu(1))
        if (len + ju0 .gt. iwk) goto 997
        do 302 k=2, len
          jlu(ju0) = jwu(k)
          alu(ju0) = wu(k)
          t = t + abs(wu(k) )
          ju0 = ju0+1
  302   continue
c
c     save norm in wu (backwards). Norm is in fact average abs value
c
        wu(n+2-ii) = t / real(len+1)
c
c     store inverse of diagonal element of u
c
        if (wu(1) .eq. 0.0) wu(1) = (0.0001 + tol)*tnorm
c
        alu(ii) = 1.0d0/ wu(1)
c
c     update pointer to beginning of next row of U.
c
	jlu(ii+1) = ju0
c
c     reset double-pointer jr to zero (U-part)
c
	do 308 k=1, lenu
          jr(jwu(k)) = 0
  308   continue
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
  500 continue
      ierr = 0
      return
c
c     zero pivot :
c
c 900    ierr = ii
c        return
c
c     incomprehensible error. Matrix must be wrong.
c
  995 ierr = -1
      return
c
c     insufficient storage in L.
c
  996 ierr = -2
      return
c
c     insufficient storage in U.
c
  997 ierr = -3
      return
c
c     illegal lfil entered.
c
  998 ierr = -4
      return
c
c     zero row encountered
c
  999 ierr = -5
      return

      end
      
c-----end-of-dbcg-------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine implu(np,umm,beta,ypiv,u,permut,full)
      real*8 umm,beta,ypiv(*),u(*),x, xpiv
      logical full, perm, permut(*)
      integer np,k,npm1
c-----------------------------------------------------------------------
c     performs implicitly one step of the lu factorization of a
c     banded hessenberg matrix.
c-----------------------------------------------------------------------
      if (np .le. 1) goto 12
      npm1 = np - 1
c
c     -- perform  previous step of the factorization-
c
      do 6 k=1,npm1
         if (.not. permut(k)) goto 5
         x=u(k)
         u(k) = u(k+1)
         u(k+1) = x
 5       u(k+1) = u(k+1) - ypiv(k)*u(k)
 6    continue
c-----------------------------------------------------------------------
c     now determine pivotal information to be used in the next call
c-----------------------------------------------------------------------
 12   umm = u(np)
      perm = (beta .gt. abs(umm))
      if (.not. perm) goto 4
      xpiv = umm / beta
      u(np) = beta
      goto 8
 4    xpiv = beta/umm
 8    permut(np) = perm
      ypiv(np) = xpiv
      if (.not. full) return
c     shift everything up if full...
      do 7 k=1,npm1
         ypiv(k) = ypiv(k+1)
         permut(k) = permut(k+1)
 7    continue
      return
c-----end-of-implu
      end
      subroutine itdriv(n,ns,rhs,sol,ipar,fpar,a,ja,ia,
     .  sotyp,iou,svomt,lenin,lenre,outso)
c******************************************************************************
c                                                                             *
c**** This routine drives the iterative solver to be used. The input is:      *
c                                                                             *
c     n:          size of the linear system                                   *
c     ns:         number of systems to be solved. Maximum = 10.               *
c     rhs(n,ns):  right-hand side of the linear system                        *
c     sol(n,ns):  solution to the linear system                               *
c     ipar(16):   integer parameter array for the reverse communication       *
c                 protocol. The components supplied on input are:             *
c                 ipar(2) < 0 No preconditioning                              *
c                         = 0 Diagonal preconditioning                        *
c                         = 1 Left preconditioning only                       *
c                         = 2 Right preconditioning only                      *
c                         = 3 Both left and right preconditioning             *
c                 ipar(5) = Size the Krylov subsapce for FOM, GMRES,          *
c                           FGMRES and DQGMRES                                *
c                 ipar(6) = maximum number of matrix-vector multiplies,       * 
c                           if not a positive number the iterative solver     *
c                           will run till convergence test is satisfied.      *
c                 ipar(13)= The fill-in parameter for the ILU precondi-       *
c                           tioner. Each row of L and each row of U           *
c                           will have a maximum of ipar(13) elements in       *
c                           addition to their original number of nonzero      *
c                           elements.                                         *
c     fpar(16):   floating-point parameter array storing information to       *
c                 and from the iterative solvers. The components supplied     *
c                 on input are:                                               *
c                 fpar(1) = relative convergence tolerance                    *
c                 fpar(13)= Thresholding tolerance in L and U. Any element    *
c                           whose size is less than this tolerance (relative  *
c                           to the norm of current row in U) is dropped (=0   *
c                           is allowed).                                      *
c     a           Matrix A of the system to be solved stored in CSR format    *
c     ja          Column indices of A                                         *
c     ia          Pointers to the beginning of each row of A                  *
c     sotyp       Integer with the solver type to be used (see above)         *
c     iou         Output unit                                                 *
c     svomt       Counter of volatile memory employed by the code (maximum)   *
c     lenin       Length of integers and pointers                             *
c     lenre       Length of reals                                             *
c     outso       Logical flag to decide whether to write the output of not   *
c                                                                             *
c******************************************************************************
      implicit none
      integer
     .  n,ns,sotyp,iou,svomt,oldme,newme,lenin,lenre,
     .  ipar(16), ia(n+1), ja(*)
      real*8
     .  fpar(16),rhs(n,ns),sol(n,ns),a(*)
      integer
     .  is, its, nwk, ierr, lbyt1, lbyt2, lbyt3,
     .  lbyt4,lbyt5,l1,l2,lby1o,lby2o,lby3o,lby4o,lby5o
      integer
     .  p_wk, p_au, p_jau, p_ju, p_di                       ! Pointers
      logical*1
     .  again,outso,tosol,dosys(10)
      real*8
     .  res,resini,rhsno(10),vecdo
      external
     .  vecdo
      data
     .  oldme,lby1o,lby2o,lby3o,lby4o,lby5o/0,0,0,0,0,0/
      save
     .  oldme,lby1o,lby2o,lby3o,lby4o,lby5o,its,res,
     .  p_wk,p_au, p_jau, p_ju, p_di
c
c***  Initializations:
c***  Convergence check: | residual | < rtol * | rhs | + atol
c***  where atol = fpar(2)
c
      newme=0
      ipar(3) = 2                                           ! Type of conv. test 
      fpar(2) = 1.0d-5*fpar(1)                              ! ATOL for conv.
      if(ns.gt.10) call srunen(
     .  'ITDRIV: TOO MANY SYSTEMS',iou)
c
c***  Memory allocation
c
      if((sotyp.eq.1).or.(sotyp.eq.2)) then
        ipar(4) = 5*n
      else if(sotyp.eq.3) then
        ipar(4) = 7*n
      else if((sotyp.eq.4).or.(sotyp.eq.6)) then
        ipar(4) = 11*n
      else if(sotyp.eq.5) then
        ipar(4) = 8*n
      else if((sotyp.eq.7).or.(sotyp.eq.8)) then
        ipar(4) = (n+3)*(ipar(5)+2) + (ipar(5)+1)*ipar(5)/2
      else if(sotyp.eq.9) then
        ipar(4) = n*(2*ipar(5)+1) + (ipar(5)+1)*ipar(5)/2
     .    + 3*ipar(5) + 2 
      else if(sotyp.eq.10) then
        ipar(4) =n + (ipar(5)+1)*(2*n+4)
      end if
      lbyt1=lenre*ipar(4)

      call maxmem(p_wk,lbyt1,lby1o,newme)
c
c***  Initial residual (for the scaled matrix, if needed)
c
      if(ipar(2).ge.0) then
        lbyt5=n*lenre
        call maxmem(p_di,lbyt5,lby5o,newme)
        call diapre(%val(p_di),sol,a,ia,ja,rhs,n,iou,1)     !A' = D^-1/2 A D^-1/2 
      end if
      tosol = .false.
      do is=1,ns
        resini = 0.0
        rhsno(is) = sqrt(vecdo(n,rhs(1,is),rhs(1,is)))
        if(rhsno(is).gt.1.0e-10) then
          call amux  (
     .      n,sol(1,is),%val(p_wk),a,ja,ia)                 ! Ax
          call sveccl(
     .      n,1.0,%val(p_wk),-1.0,rhs(1,is),%val(p_wk))     ! Ax - f
          resini=sqrt(
     .      vecdo(n,%val(p_wk),%val(p_wk)))/rhsno(is)
        end if
        if(resini.le.fpar(1)) then
          dosys(is)=.false.
          if(outso) write (iou,100) resini
        else
          dosys(is)=.true.
          tosol=.true.
        end if
      end do
      if(.not.tosol) then
        if(ipar(2).ge.0) then
          do is = 1,ns
            call diapre(
     .        %val(p_di),sol(1,is),a,ia,ja,rhs,n,iou,2)
          end do
        end if
        return
      end if
c
c***  Set-up the preconditioner (after a diagonal scaling)
c
      if(ipar(2).gt.0) then
        lbyt2 =  (n+1)*lenin
        call maxmem(p_ju,lbyt2,lby2o,newme)
        nwk = n*ipar(13) - n                                ! 1st estimate ( - n)
        ierr = -3
        do while((ierr.eq.-2).or.(ierr.eq.-3)) 
          nwk = nwk + n
          lbyt3 = lenre*nwk
          lbyt4 = lenin*nwk
          call maxmem(p_au ,lbyt3,lby3o,newme)
          call maxmem(p_jau,lbyt4,lby4o,newme)
          l1 = lenre*n
          l2 = lenin*n
          call ilut(n, a, ja, ia, ipar(13), fpar(13),
     .      %val(p_au), %val(p_jau), %val(p_ju), nwk,
     .      %val(p_wk), %val(p_wk +l1),
     .      %val(p_wk + 2*l1         ), 
     .      %val(p_wk + 2*l1 +  l2), 
     .      %val(p_wk + 2*l1 +2*l2), 
     .      ierr)
        end do
        if(ierr.ne.0)
     .    call srunei(
     .    'ITDRIV: ERROR IN THE ILU PREC., CODE = ',
     .    ierr,iou)
      end if

      do is = 1,ns
        if(dosys(is)) then
c
c***  Normal execution
c
          its = 0
          res = 0.0d0
          ipar(1) = 0
          if(sotyp.eq.1) then                               !CG method
            again=.true.
            do while(again)
              call cg    (
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.2) then                          !CGNR method
            again=.true.
            do while(again)
              call cgnr  (
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.3) then                          !BCG method
            again=.true.
            do while(again)
              call bcg   (
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.4) then                          !DBCG method
            again=.true.
            do while(again)
              call dbcg  (
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.5) then                          !BCGSTAB method
            again=.true.
            do while(again)
              call bcgstab(
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri (
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.6) then                          !TFQMR method
            again=.true.
            do while(again)
              call tfqmr (
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.7) then                          !FOM method
            again=.true.
            do while(again)
              call fom   (
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.8) then                          !GMRES method
            again=.true.
            do while(again)
              call gmres (
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.9) then                          !FGMRES method
            again=.true.
            do while(again)
              call fgmres(
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri(
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          else if(sotyp.eq.10) then                         !DQGMRES method
            again=.true.
            do while(again)
              call dqgmres(
     .          n,rhs(1,is),sol(1,is),ipar,fpar,%val(p_wk))
              call auxdri (
     .          n,ipar,fpar,its,res,iou,a,ja,ia,
     .          %val(p_au),%val(p_jau),%val(p_ju),
     .          %val(p_wk),again)
            end do
          end if
c
c*** Write to output file and recover of the unknown
c
          if(ipar(2).ge.0)                                  !x = D^-1/2 x' 
     .      call diapre(
     .      %val(p_di),sol(1,is),a,ia,ja,rhs,n,iou,2)
          if(outso) write (iou,110) fpar(3)/rhsno(is),
     .      fpar(6)/rhsno(is),fpar(7),ipar(7)
        else
          if(ipar(2).ge.0)
     .      call diapre(
     .      %val(p_di),sol(1,is),a,ia,ja,rhs,n,iou,2)
        end if
      end do

  100 format(11x,'*** CONVERGED GUESS, ERROR : ',e12.5)
  110 format(11x,'*** INITIAL RESIDUAL NORM  : ',e12.5,/,
     .       11x,'*** FINAL RESIDUAL NORM    : ',e12.5,/,
     .       11x,'*** CONVERGENCE RATE       : ',e12.5,/,
     .       11x,'*** NUMBER OF ITERATIONS   : ',i12)
c
c***  To deallocate memory, use the lines commented out
c
      if(oldme.lt.newme) then
        oldme = newme
        svomt = svomt + newme - oldme
      end if

c-m   call sadmem(2,p_wk ,lbyt1,svomt)
c-m   if(ipar(2).gt.0) then
c-m     call sadmem(2,p_ju ,lbyt2,svomt)
c-m     call sadmem(2,p_au ,lbyt3,svomt)
c-m     call sadmem(2,p_jau,lbyt4,svomt)
c-m     call sadmem(2,p_di ,lbyt5,svomt)
c-m   end if

      end
      subroutine ite_sk(
     .  rhsid,unkno,ntotv,nsist,kfact,iprob,outso)
c****************************************************************************
c
c**** This routine drives the library to solve a linear system of equations
c**** by SK iterative methods.
c
c****************************************************************************
      implicit none
      include   'MatMan.h'
      integer
     .  ntotv,nsist,iprob,kfact,iflag,svomt,
     .  lbyta,lbyts,sotyp,ipar(16)
      real*8
     .  rhsid(ntotv,nsist), unkno(ntotv,nsist)
      integer
     .  p_mat,iwoso,p_col,p_lpn,p_rof,p_rhs                 ! pointer
      logical*1   outso
      real*8      fpar(16)
      data        lbyta,svomt/0,0/
      save        lbyta,iwoso,svomt
c
c*** Initializations
c
      call svecze(16,fpar)
      call slvecz(16,ipar)
c
c*** Allocate a work space for the redefinition of RHSID and UNKNO (this
c*** work space is kept)      
c
      lbyts=8*ntotv*nsist
      if (lbyta.eq.0) then
        call sadmem(0,iwoso,lbyts,svomt)
        lbyta=lbyts
      else
        if (lbyts.gt.lbyta) then
          svome(iprob)=max(svome(iprob),svomt)
          call sadmem(2,iwoso,lbyta,svomt)
          call sadmem(0,iwoso,lbyts,svomt)
          lbyta=lbyts
        end if
      end if
      call spoias(%val(ialpn+(iprob-1)*lengp),p_lpn)
c
c*** Compute A_p x_p, where x_p are the prescribed values of the unknown
c*** (MODRHS in DIRMETH)      
c
      iflag=1
      call spoias(%val(iapma+(iprob-1)*lengp),p_mat)
      if(p_mat.eq.0) iflag=0
      call spoias(%val(irhsx+(iprob-1)*lengp),p_rhs)
      if((p_mat.ne.0).and.(kfact.le.1)) then
        call spoias(
     .    %val(iapii+(iprob-1)*lengp),p_rof)
        call spoias(
     .    %val(iapij+(iprob-1)*lengp),p_col)
        call modrhs(
     .    %val(p_lpn),%val(p_rof),%val(p_col),
     .    unkno,%val(p_rhs),%val(p_mat),ntotv,
     .    nzecp(iprob),nsist,nequa(iprob),iflag,
     .    zeroc)
      end if
c
c*** Redefine RHSID and UNKNO according to the renumbering strategy
c*** (REDEFR in DIRMETH). The prescribed values of UNKNO are put in 
c*** its nequa+1,...,ntotv components
c
      if(iflag.eq.0) p_rhs=iwoso ! not used 
      call reduit(
     .  unkno,%val(p_lpn),%val(iwoso),ntotv,nequa(iprob),
     .  nsist,1)
      call redefr(
     .  %val(iwoso),%val(p_lpn),%val(p_rhs),rhsid,ntotv,
     .  nequa(iprob),nsist,iflag)
c
c*** Call the driver for SK iterative solvers
c
      ipar( 2) = kmsip(iprob) - 1
      ipar( 5) = kkryl(iprob)
      ipar( 6) = itmax(iprob)
      ipar(13) = lfill(iprob)
      fpar( 1) = toler(iprob)
      fpar(13) = thres(iprob)
      call spoias(%val(iamat+(iprob-1)*lengp),p_mat)
      call spoias(%val(iafij+(iprob-1)*lengp),p_col)
      call spoias(%val(iafii+(iprob-1)*lengp),p_rof)
      sotyp = kites(iprob) - 30
      call itdriv(
     .  nequa(iprob),nsist,%val(iwoso),unkno,ipar,fpar,
     .  %val(p_mat),%val(p_col),%val(p_rof),sotyp,lusol,
     .  svomt,lengp,8,outso)
c
c*** Redefine UNKNO according to the renumbering strategy
c
      call reduit(
     .  unkno,%val(p_lpn),%val(iwoso),ntotv,nequa(iprob),
     .  nsist,2)
c
c*** Write temporary memory requirements 
c
      svome(iprob)=max(svome(iprob),svomt)
      if(outso) write(lusol,110) svome(iprob)
c
c*** Formats
c
  110 format(11x,'*** VOLATILE MEMORY (BYTES): ',i12)

      end
      subroutine lusol (n, y, x, alu, jlu, ju)
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)
c-----------------------------------------------------------------------
c
c performs a forward followed by a backward solve
c for LU matrix as produced by  ILUT
c
c-----------------------------------------------------------------------
c local variables
c
      integer i,k
c
c forward solve
c
      do 40 i = 1, n
        x(i) = y(i)
        do 41 k=jlu(i),ju(i)-1
          x(i) = x(i) - alu(k)* x(jlu(k))
   41   continue
   40 continue
c
c     backward solve.
c
      do 90 i = n, 1, -1
        do 91 k=ju(i),jlu(i+1)-1
          x(i) = x(i) - alu(k)*x(jlu(k))
   91   continue
        x(i) = alu(i)*x(i)
   90 continue

      return
      end
      subroutine lutsol (n, y, x, alu, jlu, ju) 
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)
c-----------------------------------------------------------------------
c
c Given a LU decomposition of a matrix stored in (alu, jlu, ju),
c this routine performs x = (LU)^{-T} y, i.e. LU transposed solve.
c 
c-----------------------------------------------------------------------
c local variables
c
      integer i,k
c
      do 10 i = 1, n
        x(i) = y(i)
   10 continue
c
c forward solve (with U^T)
c
      do 20 i = 1, n
        x(i) = x(i) * alu(i)
        do 30 k=ju(i),jlu(i+1)-1
          x(jlu(k)) = x(jlu(k)) - alu(k)* x(i)
   30   continue
   20 continue
c     
c     backward solve (with L^T)
c     
      do 40 i = n, 1, -1 
        do 50 k=jlu(i),ju(i)-1
          x(jlu(k)) = x(jlu(k)) - alu(k)*x(i)
   50   continue
   40 continue
c
      return
      end
      subroutine mgsro(full,lda,n,m,ind,ops,vec,hh,ierr)
c***************************************************************************
c
c     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
c               Orthogonalization
c     The ind'th vector of VEC is orthogonalized against the rest of
c     the vectors.
c
c     The test for performing re-orthogonalization is performed for
c     each indivadual vectors. If the cosine between the two vectors
c     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
c     performed. The norm of the 'new' vector is kept in variable NRM0,
c     and updated after operating with each vector.
c
c     full   -- .ture. if it is necessary to orthogonalize the ind'th
c               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
c               .false. only orthogonalize againt vec(:,1:ind-1)
c     lda    -- the leading dimension of VEC
c     n      -- length of the vector in VEC
c     m      -- number of vectors can be stored in VEC
c     ind    -- index to the vector to be changed
c     ops    -- operation counts
c     vec    -- vector of LDA X M storing the vectors
c     hh     -- coefficient of the orthogonalization
c     ierr   -- error code
c               0 : successful return
c               -1: zero input vector
c               -2: input vector contains abnormal numbers
c               -3: input vector is a linear combination of others
c
c***************************************************************************
      implicit none
      logical full
      integer lda,m,n,ind,ierr
      real*8  ops,hh(m),vec(lda,m)
      integer i,k
      real*8  nrm0, nrm1, fct, thr, vecdo, zero, one, reorth
      parameter (zero=0.0D0, one=1.0D0, reorth=0.98D0)
      external vecdo
c
c     compute the norm of the input vector
c
      nrm0 = vecdo(n,vec(1,ind),vec(1,ind))
      ops = ops + n + n
      thr = nrm0 * reorth
      if (nrm0.le.zero) then
         ierr = - 1
         return
      else if (nrm0.gt.zero .and. one/nrm0.gt.zero) then
         ierr = 0
      else
         ierr = -2
         return
      endif
c
c     Modified Gram-Schmidt loop
c
      if (full) then
         do 40 i = ind+1, m
            fct = vecdo(n,vec(1,ind),vec(1,i))
            hh(i) = fct
            do 20 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 20         continue
            ops = ops + 4 * n + 2
            if (fct*fct.gt.thr) then
               fct = vecdo(n,vec(1,ind),vec(1,i))
               hh(i) = hh(i) + fct
               do 30 k = 1, n
                  vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 30            continue
               ops = ops + 4*n + 1
            endif
            nrm0 = nrm0 - hh(i) * hh(i)
            if (nrm0.lt.zero) nrm0 = zero
            thr = nrm0 * reorth
 40      continue
      endif
c
      do 70 i = 1, ind-1
         fct = vecdo(n,vec(1,ind),vec(1,i))
         hh(i) = fct
         do 50 k = 1, n
            vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 50      continue
         ops = ops + 4 * n + 2
         if (fct*fct.gt.thr) then
            fct = vecdo(n,vec(1,ind),vec(1,i))
            hh(i) = hh(i) + fct
            do 60 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 60         continue
            ops = ops + 4*n + 1
         endif
         nrm0 = nrm0 - hh(i) * hh(i)
         if (nrm0.lt.zero) nrm0 = zero
         thr = nrm0 * reorth
 70   continue
c
c     test the resulting vector
c
      nrm1 = sqrt(vecdo(n,vec(1,ind),vec(1,ind)))
      ops = ops + n + n
 75   hh(ind) = nrm1
      if (nrm1.le.zero) then
         ierr = -3
         return
      endif
c
c     scale the resulting vector
c
      fct = one / nrm1
      do 80 k = 1, n
         vec(k,ind) = vec(k,ind) * fct
 80   continue
      ops = ops + n + 1
c
c     normal return
c
      ierr = 0
      return
c     end surbotine mgsro
      end
c----------------------------------------------------------------------- 
      subroutine qsplit  (a, ind, n, ncut)
      real*8 a(n)
      integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
      real*8 tmp, abskey
      integer itmp, first, last
c-----
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
    1 mid = first
      abskey = abs(a(mid))
      do 2 j=first+1, last
        if (abs(a(j)) .gt. abskey) then
          mid = mid+1
c     interchange
          tmp = a(mid)
          itmp = ind(mid)
          a(mid) = a(j)
          ind(mid) = ind(j)
          a(j)  = tmp
          ind(j) = itmp
        endif
    2 continue
c
c     interchange
c
      tmp = a(mid)
      a(mid) = a(first)
      a(first)  = tmp
c
      itmp = ind(mid)
      ind(mid) = ind(first)
      ind(first) = itmp
c
c     test for while loop
c
      if (mid .eq. ncut) return
      if (mid .gt. ncut) then
        last = mid-1
      else
        first = mid+1
      endif
      goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
      end
      subroutine reduit(unkno,lpntn,rwork,ntotv,nequa,nsist,itask)
c****************************************************************************
c
c**** This routine redefines the vector of unknowns according to the 
c**** renumbering strategy for SK iterative solvers
c
c****************************************************************************
      implicit none
      integer  ntotv,nequa,nsist,itask,isist,ipres,itotv,itotn
      real*8   unkno(ntotv,nsist), rwork(ntotv,nsist)
      integer  lpntn(ntotv)

      call srveca(ntotv*nsist,unkno,rwork) 
c
c*** Original UNKNO ---> Modified UNKNO
c      
      if(itask.eq.1) then
        do isist=1,nsist
          ipres=nequa
          do itotv=1,ntotv
            itotn=lpntn(itotv)
            if (itotn.gt.0) then
              unkno(itotn,isist)=rwork(itotv,isist)
            else if(itotn.le.0) then
              ipres=ipres+1
              unkno(ipres,isist)=rwork(itotv,isist)
            end if
          end do
        end do
c
c*** Modified UNKNO ---> Original UNKNO
c      
      else if(itask.eq.2) then
        do isist=1,nsist
          ipres=nequa
          do itotv=1,ntotv
            itotn=lpntn(itotv)
            if (itotn.gt.0) then
              unkno(itotv,isist)=rwork(itotn,isist)
            else if(itotn.le.0) then
              ipres=ipres+1
              unkno(itotv,isist)=rwork(ipres,isist)
            end if
          end do
        end do
      end if
      
      end
      logical function stopbis(n,ipar,mvpi,fpar,r,delx,sx)
c***************************************************************************
c
c     function for determining the stopping criteria. return value of
c     true if the stopbis criteria is satisfied.
c
c***************************************************************************
      implicit none
      integer n,mvpi,ipar(16)
      real*8 fpar(16), r(n), delx(n), sx, vecdo
      external vecdo

      if (ipar(11) .eq. 1) then
         stopbis = .true.
      else
         stopbis = .false.
      endif
      if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
         ipar(1) = -1
         stopbis = .true.
      endif
      if (stopbis) return
c
c     computes errors
c
      fpar(5) = sqrt(vecdo(n,r,r))
      fpar(11) = fpar(11) + 2 * n
      if (ipar(3).lt.0) then
c
c     compute the change in the solution vector
c
         fpar(6) = sx * sqrt(vecdo(n,delx,delx))
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).lt.mvpi+mvpi+1) then
c
c     if this is the end of the first iteration, set fpar(3:4)
c
            fpar(3) = fpar(6)
            if (ipar(3).eq.-1) then
               fpar(4) = fpar(1) * fpar(3) + fpar(2)
            endif
         endif
      else
         fpar(6) = fpar(5)
      endif
c
c     .. the test is struct this way so that when the value in fpar(6)
c       is not a valid number, STOPBIS is set to .true.
c
      if (fpar(6).gt.fpar(4)) then
         stopbis = .false.
         ipar(11) = 0
      else
         stopbis = .true.
         ipar(11) = 1
      endif
c
      return
      end
      subroutine tfqmr(n, rhs, sol, ipar, fpar, w)
c***************************************************************************
c
c     TFQMR --- transpose-free Quasi-Minimum Residual method
c     This is developed from BCG based on the principle of Quasi-Minimum
c     Residual, and it is transpose-free.
c
c     It uses approximate residual norm.
c
c     Internally, the fpar's are used as following:
c     fpar(3) --- initial residual norm squared
c     fpar(4) --- target residual norm squared
c     fpar(5) --- current residual norm squared
c
c     w(:,1) -- R, residual
c     w(:,2) -- R0, the initial residual
c     w(:,3) -- W
c     w(:,4) -- Y
c     w(:,5) -- Z
c     w(:,6) -- A * Y
c     w(:,7) -- A * Z
c     w(:,8) -- V
c     w(:,9) -- D
c     w(:,10) -- intermediate results of preconditioning
c     w(:,11) -- changes in the solution
c
c***************************************************************************
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
      real*8 vecdo
      logical stopbis, brkdn
      external stopbis, brkdn, vecdo
      real*8 one,zero
      parameter(one=1.0D0,zero=0.0D0)
c
c     local variables
c
      integer i
      logical lp, rp
      real*8 eta,sigma,theta,te,alpha,rho,tao
      save
c
c     status of the call (where to go)
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10,20,40,50,60,70,80,90,100,110), ipar(10)
c
c     initializations
c
      call bisinit(ipar,fpar,11*n,2,lp,rp,w)
      if (ipar(1).lt.0) return
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + 6*n
      do i = 1, n
         w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
 10   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,1) = rhs(i) - w(i,7)
         w(i,9) = zero
      enddo
      fpar(11) = fpar(11) + n
c
      if (lp) then
         ipar(1) = 3
         ipar(9) = n+1
         ipar(10) = 2
         return
      endif
 20   continue
      if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,3) = w(i,2)
         enddo
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,3) = w(i,1)
         enddo
      endif
c
      fpar(5) = sqrt(vecdo(n,w,w))
      fpar(3) = fpar(5)
      tao = fpar(5)
      fpar(11) = fpar(11) + n + n
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(vecdo(n,rhs,rhs)) + fpar(2)
         fpar(11) = fpar(11) + n + n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * tao + fpar(2)
      endif
      te = zero
      rho = zero
c
c     begin iteration
c
 30   sigma = rho
      rho = vecdo(n,w(1,2),w(1,3))
      fpar(11) = fpar(11) + n + n
      if (brkdn(rho,ipar)) goto 900
      if (ipar(7).eq.1) then
         alpha = zero
      else
         alpha = rho / sigma
      endif
      do i = 1, n
         w(i,4) = w(i,3) + alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + n + n
c
c     A * x -- with preconditioning
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 3*n + 1
         if (lp) then
            ipar(9) = 5*n + 1
         else
            ipar(9) = 9*n + 1
         endif
         ipar(10) = 3
         return
      endif
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 3*n + 1
      endif
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 5*n + 1
      endif
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 5*n + 1
         ipar(10) = 5
         return
      endif
 60   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,8) = w(i,6) + alpha * (w(i,7) + alpha * w(i,8))
      enddo
      sigma = vecdo(n,w(1,2),w(1,8))
      fpar(11) = fpar(11) + 6 * n
      if (brkdn(sigma,ipar)) goto 900
      alpha = rho / sigma
      do i = 1, n
         w(i,5) = w(i,4) - alpha * w(i,8)
      enddo
      fpar(11) = fpar(11) + 2*n
c
c     the second A * x
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         if (lp) then
            ipar(9) = 6*n + 1
         else
            ipar(9) = 9*n + 1
         endif
         ipar(10) = 6
         return
      endif
c
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 4*n + 1
      endif
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 6*n + 1
      endif
      ipar(10) = 7
      return
c
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 6*n + 1
         ipar(10) = 8
         return
      endif
 90   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,3) = w(i,3) - alpha * w(i,6)
      enddo
c
c     update I
c
      theta = vecdo(n,w(1,3),w(1,3)) / (tao*tao)
      sigma = one / (one + theta)
      tao = tao * sqrt(sigma * theta)
      fpar(11) = fpar(11) + 4*n + 6
      if (brkdn(tao,ipar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,4) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
         w(i,3) = w(i,3) - alpha * w(i,7)
      enddo
      fpar(11) = fpar(11) + 6 * n + 6
      if (ipar(7).eq.1) then
         if (ipar(3).eq.-1) then
            fpar(3) = eta * sqrt(vecdo(n,w(1,9),w(1,9)))
            fpar(4) = fpar(1)*fpar(3) + fpar(2)
            fpar(11) = fpar(11) + n + n + 4
         endif
      endif
c
c     update II
c
      theta = vecdo(n,w(1,3),w(1,3)) / (tao*tao)
      sigma = one / (one + theta)
      tao = tao * sqrt(sigma * theta)
      fpar(11) = fpar(11) + 8 + 2*n
      if (brkdn(tao,ipar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,5) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
      enddo
      fpar(11) = fpar(11) + 4*n + 3
c
c     this is the correct over-estimate
c      fpar(5) = sqrt(real(ipar(7)+1)) * tao
c     this is an approximation
      fpar(5) = tao
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 10*n + 1
         ipar(9) = 9*n + 1
         ipar(10) = 9
         return
      else if (ipar(3).lt.0) then
         fpar(6) = eta * sqrt(vecdo(n,w(1,9),w(1,9)))
         fpar(11) = fpar(11) + n + n + 2
      else
         fpar(6) = fpar(5)
      endif
      if (fpar(6).gt.fpar(4) .and. (ipar(7).lt.ipar(6)
     +     .or. ipar(6).le.0)) goto 30
 100  if (ipar(3).eq.999.and.ipar(11).eq.0) goto 30
c
c     clean up
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 10*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,10))
      else
         call tidycg(n,ipar,fpar,sol,w(1,11))
      endif
c
      return
      end
c-----end-of-stopbis
c-----------------------------------------------------------------------
      subroutine tidycg(n,ipar,fpar,sol,delx)
      implicit none
      integer i,n,ipar(16)
      real*8 fpar(16),sol(n),delx(n)
c-----------------------------------------------------------------------
c     Some common operations required before terminating the CG routines
c-----------------------------------------------------------------------
      real*8 zero
      parameter(zero=0.0D0)
c
      if (ipar(1).gt.0) then
         if ((ipar(3).eq.999 .and. ipar(11).eq.1) .or.
     +        fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
      if (fpar(3).gt.zero .and. fpar(6).gt.zero) then
         fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
      else
         fpar(7) = zero
      endif
      do i = 1, n
         sol(i) = sol(i) + delx(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine uppdir(n,p,np,lbp,indp,y,u,usav,flops)
      real*8 p(n,lbp), y(*), u(*), usav(*), x, flops
      integer k,np,n,npm1,j,ju,indp,lbp
c-----------------------------------------------------------------------
c     updates the conjugate directions p given the upper part of the
c     banded upper triangular matrix u.  u contains the non zero
c     elements of the column of the triangular matrix..
c-----------------------------------------------------------------------
      real*8 zero
      parameter(zero=0.0D0)
c
      npm1=np-1
      if (np .le. 1) goto 12
      j=indp
      ju = npm1
 10   if (j .le. 0) j=lbp
      x = u(ju) /usav(j)
      if (x .eq. zero) goto 115
      do 11 k=1,n
         y(k) = y(k) - x*p(k,j)
 11   continue
      flops = flops + 2*n
 115  j = j-1
      ju = ju -1
      if (ju .ge. 1) goto 10
 12   indp = indp + 1
      if (indp .gt. lbp) indp = 1
      usav(indp) = u(np)
      do 13 k=1,n
         p(k,indp) = y(k)
 13   continue
 208  return
c-----------------------------------------------------------------------
c-------end-of-uppdir---------------------------------------------------
      end
      function vecdo(n, a, b)
c***********************************************************************
c
c**** This function performs the vector dot product of the vectors
c**** a(n) and b(n)
c
c***********************************************************************
      implicit none
      integer n,i
      real*8  a(n), b(n), vecdo

      vecdo=0.0d0
      do i =  1, n
        vecdo = vecdo + a(i)*b(i)
      end do

      end
      subroutine arinde(nlist,nnode,nnodt,nelem,npoin,ndofn,ntotv,
     .                  ncoef,nzeco,liste,lnods,lpntn,touch,roffs,
     .                  apcol,izeco,itotv,nnodw,iblok)
c************************************************************************
c                 
c**** This routine constructs the arrays of indexes for a matrix A.
c**** These are organized as follows (CSR format):
c      
c**** ROFFS(ITOTV) = coefficient of A where eqn. ITOTV starts,
c**** APCOL(IZECO) = column of the IZECO coefficient of matrix A.
c                 
c************************************************************************
      implicit none
      integer   nlist,nnode,nnodt,nelem,npoin,ndofn,ntotv,
     .          ncoef,nzeco,ilist,ielem,jnode,jpoin,jposi,
     .          jdofn,jtotv,izeco,itotv,iblok,jtotn,klist,
     .          knode,kelem,kpoin,kposi,nnodw
      integer   liste(nlist),   lnods(nnodt,nelem), lpntn(ntotv),
     .          roffs(ntotv+1), apcol(nzeco)
      logical*1 touch(ncoef)

      roffs(itotv)=izeco
      do ilist=1,nlist
        ielem=liste(ilist)
c
c***  Columns for main degrees of freedom
c        
        do jnode=1,nnode
          jpoin=lnods(jnode,ielem)
          if(jpoin.gt.0) then
            jposi=(ilist-1)*nnodw+jnode
            if(.not.touch(jposi)) then
              do klist=1,nlist
                kelem=liste(klist)
                do knode=1,nnode
                  kpoin=lnods(knode,kelem)
                  if(kpoin.eq.jpoin) then
                    kposi=(klist-1)*nnodw+knode
                    touch(kposi)=.true.
                  end if
                end do
              end do
              do jdofn=1,ndofn
                jtotv=(jpoin-1)*ndofn+jdofn
                jtotn=lpntn(jtotv)
                if(((iblok.eq.1).and.(jtotn.gt.0))) then
                  apcol(izeco)=jtotn
                  izeco=izeco+1
                else if (((iblok.eq.2).and.(jtotn.lt.0))) then
                  apcol(izeco)=jtotv
                  izeco=izeco+1
                end if
              end do
            end if
          end if
        end do
c
c*** Columns for the Lagrange multiplier (it can not appear in A_p)
c
        if((iblok.eq.1).and.(nnodw.gt.nnode)) then
          do jnode=nnode+1,nnodw
            jpoin=lnods(jnode,ielem)
            if(jpoin.gt.0) then
              jposi=(ilist-1)*nnodw+jnode
              if(.not.touch(jposi)) then
                do klist=1,nlist
                  kelem=liste(klist)
                  do knode=nnode+1,nnodw
                    kpoin=lnods(knode,kelem)
                    if(kpoin.eq.jpoin) then
                      kposi=(klist-1)*nnodw+knode
                      touch(kposi)=.true.
                    end if
                  end do
                end do
                jtotv=(ndofn-1)*npoin+jpoin
                jtotn=lpntn(jtotv)
                if(jtotn.gt.0) then
                  apcol(izeco)=jtotn
                  izeco=izeco+1
                end if
              end if
            end if
          end do
        end if
      end do

      roffs(ntotv+1)=nzeco+1

      end
      subroutine delmem(iprob,nelem)
c**********************************************************************
c                 
c**** This routine deallocates the permanent memory for problem IPROB
c                 
c**********************************************************************
      implicit none
      include  'MatMan.h'
      integer   iprob,nelem,eamat,ielem
c
c***  Delete memory for element matrices when SOFAM = 2
c      
      if(sofam(iprob).eq.2) then
        call spoias(%val(iamat+(iprob-1)*lengp),eamat)
        do ielem = 1,nelem
          call delpoi(%val(eamat+(ielem-1)*lengp))
        end do
      end if
c
c***  Delete the rest of permanent memory
c      
      call delpoi(%val(iamat+(iprob-1)*lengp))
      call delpoi(%val(incol+(iprob-1)*lengp))
      call delpoi(%val(inrow+(iprob-1)*lengp))
      call delpoi(%val(ialpn+(iprob-1)*lengp))
      call delpoi(%val(ialpo+(iprob-1)*lengp))
      call delpoi(%val(iapma+(iprob-1)*lengp))
      call delpoi(%val(iapii+(iprob-1)*lengp))
      call delpoi(%val(iapij+(iprob-1)*lengp))
      call delpoi(%val(iafij+(iprob-1)*lengp))
      call delpoi(%val(irhsx+(iprob-1)*lengp))
c
c***  Re-initialize parameters
c      
      call slvecz(mprob,nband)
      call slvecz(mprob,nonze)
      call slvecz(mprob,nmult)
      call slvecz(mprob,svome) 
      call slvecz(mprob,speme)
      call slvecz(mprob,nzecp)
      call slvecz(mprob,nzecf)
      call slvecz(mprob,sizms)
      call slvecz(mprob,nequa)
      call slvecz(mprob,kmadi)

      end

      subroutine delpoi(point)
      integer point                                         ! pointer
      integer idumy
      if(point.gt.0) then
        call sadmem(2,point,idumy,idumy)
        point = 0
      end if

      end

      subroutine elepoi(
     .  nelem,nnodw,nnodt,npoit,lengp,lnods,svomt,
     .  p_nepoi,p_lelpo,itask)
c**********************************************************************
c                 
c**** This routine constructs the list of elements to which the nodal
c**** points belong      
c                 
c**********************************************************************
      implicit none
      integer nelem,nnodt,npoit,lengp,lbyts,ielem,ipoin,
     .        inode,nlist,svomt,nnodw,itask
      integer lnods(nnodt,nelem)
      integer p_lelpo,p_nepoi,iadd1,iadd2                   ! pointer
c
c*** Construct the array NEPOI (number of elements per node) and allocate
c*** memory for the list
c
      lbyts=lengp*npoit
      if(itask.eq.1) then                                   ! Allocate
        call sadmem(0,p_nepoi,lbyts,svomt)                  ! N. of elements
        call sadmem(0,p_lelpo,lbyts,svomt)                  ! Adresses 
        call nelepo(p_lelpo,%val(p_nepoi),lnods,nelem,
     .              nnodw,nnodt,npoit,lengp,svomt,itask)
      else if(itask.eq.2) then                              ! Deallocate
        call nelepo(p_lelpo,%val(p_nepoi),lnods,nelem,
     .              nnodw,nnodt,npoit,lengp,svomt,itask)
        call sadmem(2,p_nepoi,lbyts,svomt)   
        call sadmem(2,p_lelpo,lbyts,svomt)   
      end if
      if(itask.eq.2) return
c
c*** Construct the list
c
      do ielem=1,nelem
        do inode=1,nnodw
          ipoin=lnods(inode,ielem)
          if(ipoin.gt.0) then
            call idenum(%val(p_nepoi),npoit,ipoin,nlist)
            iadd1=p_lelpo+(ipoin-1)*lengp
            call spoias(%val(iadd1),iadd2)
            call filnum(%val(iadd2),nlist,ielem)
          end if
        end do
      end do

      end
      subroutine filnum(liste,nlist,ielem)
c**********************************************************************
c                 
c**** This routine writes an integer in a given memory adress
c                 
c**********************************************************************
      implicit none
      integer  nlist,ielem,ilist
      integer  liste(nlist)

      ilist=1
      do while(liste(ilist).ne.0)
        ilist=ilist+1
      end do
      liste(ilist)=ielem

      end
      subroutine idenum(nepoi,npoin,ipoin,nlist)
c**********************************************************************
c                 
c**** This routine identifies a number in a vector whose address
c**** is known
c                 
c**********************************************************************
      implicit none
      integer  npoin,ipoin,nlist
      integer  nepoi(npoin)

      nlist=nepoi(ipoin)

      end
      subroutine inispa(ndimt,amatr,lncol,lnrow)
c**********************************************************************
c                 
c**** This routine initializes the stiffness matrix
c                 
c**********************************************************************
      implicit   real*8(a-h,o-z)
      real*8     amatr(ndimt)
      integer    lncol(ndimt),lnrow(ndimt)

      do i=1,ndimt
        amatr(i)=0.0d0
        lncol(i)=0
        lnrow(i)=0
      end do

      end
      subroutine mameba(
     .  lnods,lpntn,ndofn,nnode,nnodt,nnodw,
     .  nelem,npoin,ntotv,nequa,ksymm,nband,
     .  iamat,iprob,lengp,ndimt,speme,sizms,
     .  lusol)
c**********************************************************************
c                 
c**** This routine allocates memory for the system matrix when it is
c**** stored in band
c                 
c**********************************************************************
      implicit   none
      integer    iamat                                      ! pointer
      integer
     .  ndofn,nnode,nnodt,nnodw,nelem,npoin,ntotv,
     .  nequa,ksymm,nband,iprob,lengp,ndimt,speme,
     .  sizms,lusol
      integer
     .  lnods(nnodt,nelem), lpntn(ntotv)
      integer
     .  ielem,inode,ipoin,idofn,itotv,jnode,jpoin,
     .  jdofn,jtotv,iband,lbyts
c
c*** Half-bandwidth
c
      nband=0
      do ielem=1,nelem
        do inode=1,nnodw
          ipoin=lnods(inode,ielem)
          if(ipoin.gt.0) then
            idofn=0
            do while(idofn.lt.ndofn)
              if(ipoin.le.npoin) then
                idofn=idofn+1
                itotv=(ipoin-1)*ndofn+idofn
              else
                itotv=npoin*(ndofn-1)+ipoin
                idofn=ndofn
              end if
              if (lpntn(itotv).gt.0) then
                do jnode=inode+1,nnodw
                  jpoin=lnods(jnode,ielem)
                  if(jpoin.gt.0) then
                    jdofn=0
                    do while(jdofn.lt.ndofn)
                      if(jpoin.le.npoin) then
                        jdofn=jdofn+1
                        jtotv=(jpoin-1)*ndofn+jdofn
                      else
                        jtotv=npoin*(ndofn-1)+jpoin
                        jdofn=ndofn
                      end if
                      if(lpntn(jtotv).gt.0)
     .                  iband=abs(lpntn(itotv)-lpntn(jtotv))
                      if(iband.gt.nband) nband=iband
                    end do
                  end if
                end do
              end if
            end do
          end if
        end do
      end do
c
c*** Memory allocation
c
      ndimt=(nband*(1+ksymm)+1)*nequa
      lbyts=ndimt*8
      call matpoi(iprob,iamat,lengp,sizms,speme,lbyts)
c
c*** Write solver information
c
      write(lusol,900) nband,ndimt 
      
  900 format(11x,'*** HALF BANDWIDTH              :',i10,/
     .       11x,'*** SIZE OF THE MATRIX          :',i10,' (R*8)')
      
      end
      subroutine mamelm(eamat,ksymm,nnode,nnodw,ndofn,nelem,
     .                  lures,lengp,speme)
c*************************************************************************
c                 
c**** This routine allocates memory for the elemental matrices when they
c**** are stored for iterative solvers
c                 
c*************************************************************************
      implicit none
      integer       ksymm,nnode,nnodw,ndofn,nelem,nevat,
     .              ndimt,ielem,lures,lengp,speme
      integer       eamat,p_mat                             ! pointer
      
      nevat=(nnode*ndofn+nnodw-nnode)
      ndimt=(ksymm+1)*(nevat*nevat-nevat)/2+nevat
      do ielem=1,nelem
        call sadmem(0,p_mat,ndimt*8,speme)
        call spoias(p_mat,%val(eamat+(ielem-1)*lengp))
        call svecze(ndimt,%val(p_mat))
      end do
      p_mat=0
      call spoias(p_mat,%val(eamat+nelem*lengp))

      write(lures,900) nelem,ndimt

  900 format(11x,'*** NUMBER OF ELEMENT MATRICES  :',i10,/
     .       11x,'*** SIZE OF EACH ELEMENT MATRIX :',i10,' (R*8)')

      end
      subroutine mamepr(lpont,nequa,ntotv,ksymm,ndimt,iamat,
     .                  iprob,lengp,speme,sizms,lusol)
c**********************************************************************
c                 
c**** This routine allocates memory for the system matrix when it is
c**** stored in profile
c                 
c**********************************************************************
      implicit   none
      integer    speme,nlast,nequa,ntotv,meanh,nsizz,
     .           ksymm,ndimt,iprob,lengp,lbyts,lusol,
     .           sizms
      integer    iamat                                      ! pointer
      real*8     hmean
      integer    lpont(ntotv)
c
c*** Write solver information
c
      hmean=0.0
      nlast=lpont(nequa)
      if(nequa.ne.0) hmean=float(nlast+nequa)/float(nequa)
      meanh=nint(hmean)
      nsizz=nequa+nlast*(1+ksymm)
      write(lusol,900) meanh,nsizz
c
c*** Allocate memory
c      
      ndimt=nequa+(1+ksymm)*nlast
      lbyts=8*ndimt
      call matpoi(iprob,iamat,lengp,sizms,speme,lbyts)
c
c*** Formats
c
  900 format(11x,'*** MEAN HALF-WIDTH             :',i10,/
     .       11x,'*** SIZE OF THE MATRIX          :',i10,' (R*8)')

      end
      subroutine mamesp(
     .  nequa,lnods,lpntn,nnodt,nelem,nnodw,
     .  npoin,ntotv,ndofn,nonze,ndimt,iamat,
     .  incol,inrow,lengp,speme)
c**********************************************************************
c                 
c**** This routine allocates memory for the system matrix when it is
c**** stored in sparse
c                 
c**********************************************************************
      implicit   real*8(a-h,o-z)
      integer    lnods(nnodt,nelem), lpntn(ntotv)
      integer    speme,lengp

      nonze=0
      do ielem=1,nelem
        do inode=1,nnodw
          ipoin=lnods(inode,ielem)
          if(ipoin.gt.0) then
            idofn=0
            do while(idofn.lt.ndofn)
              if(ipoin.le.npoin) then
                idofn=idofn+1
                itotv=(ipoin-1)*ndofn+idofn
              else
                itotv=npoin*(ndofn-1)+ipoin
                idofn=ndofn
              end if
              if (lpntn(itotv).gt.0) then
                do jnode=inode,nnodw
                  jpoin=lnods(jnode,ielem)
                  if(jpoin.gt.0) then
                    jdofn=0
                    do while(jdofn.lt.ndofn)
                      if(jpoin.le.npoin) then
                        jdofn=jdofn+1
                        jtotv=(jpoin-1)*ndofn+jdofn
                      else
                        jtotv=npoin*(ndofn-1)+jpoin
                        jdofn=ndofn
                      end if
                      if(lpntn(jtotv).gt.0)
     .                  nonze=nonze+1
                    end do
                  end if
                end do
              end if
            end do
          end if
        end do
      end do
      nonze=nonze*2
c
c*** Memory allocation
c
      ndimt=nonze*3
      lbyts=ndimt*8
      call sadmem(0,iamat,lbyts,speme)
      lbyts=ndimt*lengp
      call sadmem(0,incol,lbyts,speme)
      lbyts=ndimt*lengp
      call sadmem(0,inrow,lbyts,speme)
      nonze=0
      
      end
      subroutine matpoi(iprob,iamat,lengp,sizms,speme,lbyts)
c**********************************************************************
c                 
c**** This routine allocates memory for the system matrix if more than
c**** what is available is needed      
c                 
c**********************************************************************
      implicit none
      integer  iprob,iamat,lengp,sizms,speme,lbyts,p_mat

      call spoias(%val(iamat+(iprob-1)*lengp),p_mat)
      if(lbyts.gt.sizms) then
        if(p_mat.ne.0)
     .    call sadmem(2,%val(iamat+(iprob-1)*lengp),sizms,speme)
        call sadmem(0,%val(iamat+(iprob-1)*lengp),lbyts,speme)
        sizms = lbyts
      end if

      end
      subroutine maxmem(addre,newby,oldby,memor)
c******************************************************************************
c
c**** This is an auxiliary routine to allocate the maximum memory
c
c******************************************************************************
      implicit none
      integer newby,oldby,memor
      integer addre                                         ! pointer

      if(oldby.le.0) then 
        call sadmem(0,addre,newby,memor)
        oldby = newby
      else
        if(newby.gt.oldby) then
          call sadmem(2,addre,oldby,memor)
          call sadmem(0,addre,newby,memor)
          oldby = newby
        end if
      end if

      end
      
      subroutine memdof(nnode,nnodt,nelem,npoin,ndofn,ntotv,lpntn,
     .                  lnods,svomt,svome,speme,sizms,lengp,iblok,
     .                  nzeco,p_amatr,p_roffs,p_acolu) 
c************************************************************************
c                 
c**** This routine allocates memory for the coefficients of a matrix A
c**** stored in compressed sparse row (CSR) format. This matrix is
c**** assumed to come from a problem of the form     
c
c     A_f x_f = F - A_p x_p
c
c**** where x_f is unknown and x_p are prescribed values. Moreover, the
c**** vector x is allowed to be of the form x = ( y , z ) (when it is 
c**** properly ordered), where y are referred to as the `main' degrees
c**** of freedom and z is thought of as a Lagrange multiplier.
c
c**** INPUT:
c      
c     nnode:   Number of nodes per element of y
c     nnodt:   Total number of nodes per element (of y and z)
c     nelem:   Number of elements            
c     npoin:   Number of nodal points of y
c     ndofn:   Number of degrees of freedom per node of y (ONLY 1 for z) 
c     ntotv:   Total number of degrees of freedom (y and z)            
c     lpntn:   Array of ordering of d.o.f. to define prescriptions and
c              to improve the efficiency (if wanted). It is assumed that
c              lpntn(i) > 0 if i is free and < 0 if i is prescribed
c     lnods:   Matrix of nodal connections (of both y and z. It has
c              dimension (nnodt,nelem))
c     svomt:   Counter for volatile memory
c     svome:   Maximum of volatile memory
c     speme:   Counter for permanent memory
c     sizms:   Size previously needed to store A and the arrays of indexes
c     lengp:   Length of integers (4 or 8)
c     iblok:   Code for the contruction of matrices, 1 for A_f, 2 for A_p
c
c**** OUTPUT:
c
c     nzeco:   Number of nonzero coefficients of the matrix constructed
c     p_amatr: Pointer to matrix A (amatr(nzeco))
c     p_roffs: Pointer to the array of rows roffs(ntotv+1) (roffs(i) = 
c              coefficient of A where eqn. i starts)
c     p_acolu: Pointer to the array of columns acolu(nzeco) (acolu (i)
c              = column of the i coefficient of matrix A)
c                 
c************************************************************************
      implicit none
      integer   nnode,nnodt,nelem,npoin,ndofn,ntotv,svomt,
     .          svome,speme,sizms,lengp,iblok,nzeco
      integer   p_amatr,p_roffs,p_acolu                     ! pointers
      integer   lpntn(ntotv),lnods(nnodt,nelem) 
      integer   nlist,ncoef,itotv,itotn,ipoin,kflag,idofn,
     .          lbyts,izeco,npoit,nnodw
      integer   p_nepoi,p_lelpo,p_touch,iadd1,iadd2         ! pointers
c
c*** Initializations
c
      npoit=npoin+(ntotv-ndofn*npoin)
      nnodw=nnode
      if(npoit.gt.npoin) nnodw=nnodt
c
c*** List of elements connected to the nodal points
c      
      call elepoi(nelem,nnodw,nnodt,npoit,lengp,lnods,svomt,
     .            p_nepoi,p_lelpo,1)
c
c*** Total number of nonzero coefficients
c
      nzeco=0
      itotv=0
      call sadmem(0,p_touch,nelem*nnodw,svomt)
      do ipoin=1,npoin                                      ! Main points
        kflag=0
        do idofn=1,ndofn
          itotv=itotv+1
          itotn=lpntn(itotv)                           
          if(itotn.gt.0) then                               ! Row of unknowns
            if(kflag.eq.0) then
              kflag=1
              call idenum(%val(p_nepoi),npoit,ipoin,nlist)
              ncoef=nlist*nnodw
            end if
            call slogze(ncoef,%val(p_touch))
            iadd1=p_lelpo+(ipoin-1)*lengp
            call spoias(%val(iadd1),iadd2)
            call nzecoe(nlist,nnode,nnodt,nelem,npoin,ndofn,
     .                  ntotv,ncoef,nzeco,%val(iadd2),
     .                  lnods,lpntn,%val(p_touch),nnodw,iblok)
          end if
        end do
      end do
      itotv=ndofn*npoin
      do ipoin=npoin+1,npoit                                ! Secondary points
        itotv=itotv+1
        itotn=lpntn(itotv) 
        if(itotn.gt.0) then                                 ! Row of unknowns
          call idenum(%val(p_nepoi),npoit,ipoin,nlist)
          ncoef=nlist*nnodw
          call slogze(ncoef,%val(p_touch))
          iadd1=p_lelpo+(ipoin-1)*lengp
          call spoias(%val(iadd1),iadd2)
          call nzecoe(nlist,nnode,nnodt,nelem,npoin,ndofn,
     .                ntotv,ncoef,nzeco,%val(iadd2),
     .                lnods,lpntn,%val(p_touch),nnodw,iblok)
        end if
      end do
c
c*** Allocate memory for A and the arrays of indexes. These are organized
c*** as follows: ROFFS(ITOTV) = coefficient of A where eqn. ITOTV starts,
c*** ACOLU(IZECO) = column of the IZECO coefficient of matrix A.
c
      lbyts=8*nzeco
      call sadmem(0,p_amatr,lbyts,speme)
      lbyts=lengp*(ntotv+1)
      call sadmem(0,p_roffs,lbyts,speme)
      call slvecz(ntotv+1,%val(p_roffs))
      lbyts=lengp*nzeco
      call sadmem(0,p_acolu,lbyts,speme)
      call slvecz(nzeco,%val(p_acolu))
c
c*** Construct the array of indexes
c
      itotv=0
      izeco=1
      do ipoin=1,npoin                                      ! Main points
        kflag=0
        do idofn=1,ndofn
          itotv=itotv+1
          itotn=lpntn(itotv) 
          if(itotn.gt.0) then                               ! Row of unknowns
            if(kflag.eq.0) then
              kflag=1
              call idenum(%val(p_nepoi),npoit,ipoin,nlist)
              ncoef=nlist*nnodw
            end if
            call slogze(ncoef,%val(p_touch))
            iadd1=p_lelpo+(ipoin-1)*lengp
            call spoias(%val(iadd1),iadd2)
            call arinde(nlist,nnode,nnodt,nelem,npoin,ndofn,
     .                  ntotv,ncoef,nzeco,%val(iadd2),
     .                  lnods,lpntn,%val(p_touch),
     .                  %val(p_roffs),%val(p_acolu),
     .                  izeco,itotv,nnodw,iblok)
          end if
        end do
      end do
      itotv=ndofn*npoin
      do ipoin=npoin+1,npoit                                ! Secondary points
        itotv=itotv+1
        itotn=lpntn(itotv) 
        if(itotn.gt.0) then                                 ! Row of unknowns
          call idenum(%val(p_nepoi),npoit,ipoin,nlist)
          ncoef=nlist*nnodw
          call slogze(ncoef,%val(p_touch))
          iadd1=p_lelpo+(ipoin-1)*lengp
          call spoias(%val(iadd1),iadd2)
          call arinde(nlist,nnode,nnodt,nelem,npoin,ndofn,
     .                ntotv,ncoef,nzeco,%val(iadd2),
     .                lnods,lpntn,%val(p_touch),
     .                %val(p_roffs),%val(p_acolu),
     .                izeco,itotv,nnodw,iblok)
        end if
      end do
c
c*** Deallocate volatile memory 
c
      svome = max(svome,svomt)
      call sadmem(2,p_touch,nelem*nnodw,svomt)
      call elepoi(nelem,nnodw,nnodt,npoit,lengp,lnods,svomt,
     .            p_nepoi,p_lelpo,2)
      
      end
      subroutine memmat(
     .  lnods,ndofn,nnode,nnodt,nnodw,nelem,
     .  npoin,ntotv,nsist,kimpo,iprob)
c**********************************************************************
c                 
c**** This routine allocates memory for the matrix of problem IPROB
c                 
c**********************************************************************
      implicit none
      include 'MatMan.h'
      integer
     .  p_mat,p_lpn,p_lpo,p_rhs,p_row,p_col,p_tem           ! pointer
      integer
     .  ndofn,nnode,nnodt,nnodw,nelem,npoin,ntotv,
     .  iprob,ndimt,kimpo,nsist,svomt,inode,ielem,
     .  lbyts,losiz
      integer      lnods(nnodt,nelem)
      character*50 solty
      data         svomt/0/
      save         svomt
c
c*** Write permanent memory information
c
      if(sofam(iprob).eq.1) then
        if(kites(iprob).eq.11) then
          if(kmsip(iprob).eq.0) then
            solty='  LDU (BANDED STORAGE)'
          else if(kmsip(iprob).eq.1) then
            solty=' LDU (SKYLINE STORAGE)'
          else
            call srunen(
     .        'MEMMAT: PIVOTING NOT AVAILABLE',lusol)
          end if
        else if(kites(iprob).eq.12) then
          if(kmsip(iprob).eq.2) then
            solty='  LDU (SPARSE STORAGE)'
          else
            call srunen(
     .        'MEMMAT: ONLY PIVOTING AVAILABLE',lusol)
          end if
        else
          call srunen('MEMMAT: WRONG SOLVER TYPE',lusol)
        end if
      else if(sofam(iprob).eq.2) then
        nequa(iprob)=ntotv
        if(kites(iprob).eq.21) then
          solty=' EE_CG METHOD'
        else if(kites(iprob).eq.22) then
          solty=' EE_GMRES METHOD'
        else
          call srunen('MEMMAT: WRONG SOLVER TYPE',lusol)
        end if
        if(kmsip(iprob).eq.-2) then
          solty = solty//' (DIAG. PREC.)'
        end if
      else if(sofam(iprob).eq.3) then
        if(kites(iprob).eq.31) then
          solty=' SK_CG METHOD'
        else if(kites(iprob).eq.32) then
          solty=' SK_CGNR METHOD'
        else if(kites(iprob).eq.33) then
          solty=' SK_BCG METHOD'
        else if(kites(iprob).eq.34) then
          solty=' SK_DBCG METHOD'
        else if(kites(iprob).eq.35) then
          solty=' SK_BCGSTAB METHOD'
        else if(kites(iprob).eq.36) then
          solty=' SK_TFQMR METHOD'
        else if(kites(iprob).eq.37) then
          solty=' SK_FOM METHOD'
        else if(kites(iprob).eq.38) then
          solty=' SK_GMRES METHOD'
        else if(kites(iprob).eq.39) then
          solty=' SK_FGMRES METHOD'
        else if(kites(iprob).eq.40) then
          solty=' SK_DQGMRES METHOD'
        else
          call srunen('MEMMAT: WRONG SOLVER TYPE',lusol)
        end if
        if(kmsip(iprob).eq.1) then
          solty = solty//' (DIAG. PREC.)'
        else if(kmsip(iprob).eq.2) then
          solty = solty//' (L_ILUT PREC.)'
        else if(kmsip(iprob).eq.3) then
          solty = solty//' (R_ILUT PREC.)'
        else if(kmsip(iprob).eq.4) then
          solty = solty//' (L_R_ILUT PREC.)'
        end if
      else
        call srunen('MEMMAT: WRONG SOLVER TYPE',lusol)
      end if
      
      write(lusol,900)
     .  wosol(iprob),solty,ntotv,nequa(iprob)
      if(kkryl(iprob).gt.0) write(lusol,910) kkryl(iprob)
      call spoias(%val(ialpn+(iprob-1)*lengp),p_lpn)

      if(sofam(iprob).eq.1) then
c
c*** Band matrix storage for direct solvers
c
        if(kmsip(iprob).eq.0) then
          call mameba(
     .      lnods,%val(p_lpn),ndofn,nnode,nnodt,
     .      nnodw,nelem,npoin,ntotv,nequa(iprob),
     .      ksymm(iprob),nband(iprob),iamat,iprob,
     .      lengp,ndimt,speme(iprob),sizms(iprob),
     .      lusol)
          kmadi(iprob)=ndimt
c
c*** Profile matrix storage for direct solvers
c
        else if(kmsip(iprob).eq.1) then
          call spoias(%val(ialpo+(iprob-1)*lengp),p_lpo)
          call mamepr(
     .      %val(p_lpo),nequa(iprob),ntotv,ksymm(iprob),
     .      ndimt,iamat,iprob,lengp,speme(iprob),
     .      sizms(iprob),lusol)
          kmadi(iprob)=ndimt
c
c*** Sparse matrix storage for direct solvers. The same routine as for the
c*** case of a banded storage is employed
c
        else if(kmsip(iprob).eq.2) then
          call mameba(
     .      lnods,%val(p_lpn),ndofn,nnode,nnodt,
     .      nnodw,nelem,npoin,ntotv,nequa(iprob),
     .      ksymm(iprob),nband(iprob),iamat,iprob,
     .      lengp,ndimt,speme(iprob),sizms(iprob),
     .      lusol)
          kmadi(iprob)=ndimt
        end if
c
c*** Store the vector of pointers to elemental matrices for EE iterative solvers
c
      else if(sofam(iprob).eq.2) then
        lbyts = lengp*(nelem+1)
        call sadmem(0,p_mat,lbyts,speme(iprob))
        call spoias(p_mat,%val(iamat+(iprob-1)*lengp))
        call mamelm(
     .    p_mat,ksymm(iprob),nnode,nnodw,ndofn,nelem,
     .    lusol,lengp,speme(iprob))
        kmadi(iprob)=1
      else if(sofam(iprob).eq.3) then
c
c*** Memory allocation for the sparse matrix used in SK iterative solvers
c
        call memdof(
     .    nnode,nnodt,nelem,npoin,ndofn,ntotv,
     .    %val(p_lpn),lnods,svomt,svome(iprob),
     .    speme(iprob),sizms(iprob),lengp,   1,
     .    nzecf(iprob),p_mat,p_row,p_col)
        lbyts = lengp*ntotv
        call sadmem(0,p_tem,lbyts,svomt)
        call ordspa(%val(p_row),%val(p_lpn),%val(p_tem),ntotv,
     .              nequa(iprob),nzecf(iprob))
        svome(iprob) = max(svome(iprob),svomt)
        call sadmem(2,p_tem,lbyts,svomt)
        call spoias(p_mat,%val(iamat+(iprob-1)*lengp))
        call spoias(p_row,%val(iafii+(iprob-1)*lengp))
        call spoias(p_col,%val(iafij+(iprob-1)*lengp))
        kmadi(iprob)=1
        write(lusol,911) nzecf(iprob)
      end if
c
c*** Allocate memory for the terms of the matrix corresponding to the
c*** prescribed degrees of freedom (only for direct and SK iterative solvers)
c
      if(((sofam(iprob).eq.1).or.(sofam(iprob).eq.3))
     .  .and.(kimpo.ne.0)) then
        if(nequa(iprob).lt.ntotv) then
          lbyts = 8*nsist*ntotv
          losiz = 0
          call sadmem(0,p_rhs,lbyts,speme(iprob))
          call spoias(p_rhs,%val(irhsx+(iprob-1)*lengp))
          call memdof(
     .      nnode,nnodt,nelem,npoin,ndofn,ntotv,
     .      %val(p_lpn),lnods,svomt,svome(iprob),
     .      speme(iprob),losiz,lengp,   2,
     .      nzecp(iprob),p_mat,p_row,p_col)
          call spoias(p_mat,%val(iapma+(iprob-1)*lengp))
          call spoias(p_row,%val(iapii+(iprob-1)*lengp))
          call spoias(p_col,%val(iapij+(iprob-1)*lengp))
        end if
c
c*** For EE iterative solvers, allocate memory only for the product A x_p
c*** (A is not explicitly built up)        
c
      else if((sofam(iprob).eq.2).and.(kimpo.ne.0)) then
        lbyts = 8*nsist*ntotv
        call sadmem(0,p_rhs,lbyts,speme(iprob))
        call spoias(p_rhs,%val(irhsx+(iprob-1)*lengp))
      end if
c
c*** Update volatile memory and write permanent memory 
c
      write(lusol,920) speme(iprob)
c
c*** If there are negative values in LNODS, change their sign
c
      if(conde(iprob)) then
        do ielem=1,nelem
          do inode=1,nnodt
            lnods(inode,ielem)=abs(lnods(inode,ielem))
          end do
        end do
      end if
c
c*** Formats
c      
  900 format(//,5x,67('#'),//,5x,
     .  '>>>> PERMANENT MEMORY INFORMATION FOR:  ',a15,/,5x,
     .  '     --------------------------------   ',//,
     .       11x,'*** SOLVER TYPE                 :',a35,/,     
     .       11x,'*** NUMBER OF D.O.F.            :',i10,/,
     .       11x,'*** NUMBER OF EQUATIONS         :',i10)
  910 format(11x,'*** KRYLOV DIMENSION            :',i10)
  911 format(11x,'*** NONZERO MATRIX COEFFICIENTS :',i10)
  920 format(11x,'*** TOTAL PERMANENT MEMORY      :',i10,' (BYTES)')

      end
      subroutine nelepo(p_lelpo,nepoi,lnods,nelem,nnodw,nnodt,
     .                  npoit  ,lengp,svomt,itask)
c**********************************************************************
c                 
c**** This routine computes the number of elements to which the nodal
c**** points belong and allocates or deallocates memory for the list 
c                 
c**********************************************************************
      implicit none
      integer nelem,nnodt,npoit,ielem,inode,ipoin,lbyts,
     .        lengp,svomt,nnodw,itask
      integer nepoi(npoit), lnods(nnodt,nelem)
      integer p_lelpo,p_ele                                 ! pointer
c
c*** Number of elements per nodal point
c
      if(itask.eq.1) then
        call slvecz(npoit,nepoi)
        do ielem=1,nelem
          do inode=1,nnodw
            ipoin=lnods(inode,ielem)
            if(ipoin.gt.0)
     .        nepoi(ipoin)=nepoi(ipoin)+1
          end do
        end do
c
c*** Allocate memory for the list of elements
c
        do ipoin=1,npoit
          lbyts=lengp*nepoi(ipoin)
          call sadmem(0,p_ele,lbyts,svomt)
          call spoias(p_ele,%val(p_lelpo+(ipoin-1)*lengp))
          call slvecz(nepoi(ipoin),%val(p_ele))
        end do
      else if(itask.eq.2) then
c
c*** Deallocate memory for the list of elements
c
        do ipoin=1,npoit
          lbyts=lengp*nepoi(ipoin)
          call spoias(%val(p_lelpo+(ipoin-1)*lengp),p_ele)
          call sadmem(2,p_ele,lbyts,svomt)
        end do
      end if
        
      end
      subroutine nzecoe(nlist,nnode,nnodt,nelem,npoin,ndofn,
     .                  ntotv,ncoef,nzeco,liste,lnods,lpntn,
     .                  touch,nnodw,iblok)
c************************************************************************
c                 
c**** This routine computes the number of non-zero coefficients of a
c**** matrix A stored in compressed sparse row (CSR) format 
c                 
c************************************************************************
      implicit none
      integer   nlist,nnode,nnodt,nelem,npoin,ndofn,ntotv,
     .          ncoef,nzeco,iblok,ilist,ielem,jnode,jpoin,
     .          jposi,jdofn,jtotv,jtotn,klist,kelem,knode,
     .          kpoin,kposi,nnodw
      integer   liste(nlist), lnods(nnodt,nelem), lpntn(ntotv)
      logical*1 touch(ncoef)
      
      do ilist=1,nlist
        ielem=liste(ilist)
c
c***  Columns for main degrees of freedom
c        
        do jnode=1,nnode
          jpoin=lnods(jnode,ielem)
          if(jpoin.gt.0) then
            jposi=(ilist-1)*nnodw+jnode
            if(.not.touch(jposi)) then
              do klist=1,nlist
                kelem=liste(klist)
                do knode=1,nnode
                  kpoin=lnods(knode,kelem)
                  if(kpoin.eq.jpoin) then
                    kposi=(klist-1)*nnodw+knode
                    touch(kposi)=.true.
                  end if 
                end do
              end do
              do jdofn=1,ndofn
                jtotv=(jpoin-1)*ndofn+jdofn
                jtotn=lpntn(jtotv)           
                if(  ((iblok.eq.1).and.(jtotn.gt.0)).       ! Free column
     .            or.((iblok.eq.2).and.(jtotn.lt.0)))       ! or prescribed
     .            nzeco=nzeco+1 
              end do
            end if
          end if
        end do
c
c*** Columns for the Lagrange multiplier (it can not appear in A_p)
c        
        if((iblok.eq.1).and.(nnodw.gt.nnode)) then
          do jnode=nnode+1,nnodw
            jpoin=lnods(jnode,ielem)
            if(jpoin.gt.0) then
              jposi=(ilist-1)*nnodw+jnode
              if(.not.touch(jposi)) then
                do klist=1,nlist
                  kelem=liste(klist)
                  do knode=nnode+1,nnodw
                    kpoin=lnods(knode,kelem)
                    if(kpoin.eq.jpoin) then
                      kposi=(klist-1)*nnodw+knode
                      touch(kposi)=.true.
                    end if 
                  end do
                end do
                jtotv=(ndofn-1)*npoin+jpoin
                jtotn=lpntn(jtotv)
                if(jtotn.gt.0) nzeco = nzeco+1 
              end if
            end if
          end do
        end if
      end do

      end
      subroutine ordspa(roffs,lpntn,newof,ntotv,nequa,nzeco)
c*************************************************************************
c                 
c**** This routine orders the array ROFFS used for sparse matrix storage 
c**** to account for the prescribed degrees of freedom. Recall that, 
c**** originally,
c      
c**** ROFFS(ITOTV) = coefficient of A where eqn. ITOTV starts
c
c**** WARNING: Renumbering is NOT possible using a profile matrix storage
c****          due to the fact that the arrays ROFFS and ACOLU in MEMDOF
c****          are computed using the consecutive order of the original
c****          numbering. This is used in particular to compute in diffe-
c****          rent routines the first and last elements of a row IROWS,
c****          which are ROFFS(IROWS) and ROFFS(IROWS+1)-1.
c                 
c*************************************************************************
      implicit none
      integer  ntotv,itotv,itotn,nequa,nzeco
      integer  roffs(ntotv+1),lpntn(ntotv),newof(ntotv)

      do itotv=1,ntotv
        newof(itotv)=0
      end do
      do itotv=1,ntotv
        itotn=lpntn(itotv)
        if(itotn.gt.0) newof(itotn) = roffs(itotv)
      end do
      do itotv=1,ntotv
        roffs(itotv)=newof(itotv)
      end do

      roffs(nequa+1)=nzeco+1

      end
      subroutine almepo
c**********************************************************************
c                 
c**** This routine allocates memory for pointers in MatMan.h
c                 
c**********************************************************************
      implicit none
      include 'MatMan.h'
      integer  idumy

      call sadmem(0,iamat,mprob*lengp,idumy)
      call sadmem(0,incol,mprob*lengp,idumy)
      call sadmem(0,inrow,mprob*lengp,idumy)
      call sadmem(0,ialpn,mprob*lengp,idumy)
      call sadmem(0,ialpo,mprob*lengp,idumy)
      call sadmem(0,iapma,mprob*lengp,idumy)
      call sadmem(0,iapii,mprob*lengp,idumy)
      call sadmem(0,iapij,mprob*lengp,idumy)
      call sadmem(0,iafii,mprob*lengp,idumy)
      call sadmem(0,iafij,mprob*lengp,idumy)
      call sadmem(0,irhsx,mprob*lengp,idumy)      
      call slvecz(mprob,%val(iamat))
      call slvecz(mprob,%val(incol))
      call slvecz(mprob,%val(inrow))
      call slvecz(mprob,%val(ialpn))
      call slvecz(mprob,%val(ialpo))
      call slvecz(mprob,%val(iapma))
      call slvecz(mprob,%val(iapii))
      call slvecz(mprob,%val(iapij))
      call slvecz(mprob,%val(iafii))
      call slvecz(mprob,%val(iafij))
      call slvecz(mprob,%val(irhsx))

      end
      subroutine checpa
c**********************************************************************
c                 
c**** This routine checks some parameters
c                 
c**********************************************************************
      implicit none
      include   'MatMan.h'
      integer ipass
      data ipass/0/
      save ipass

      if(ipass.eq.0) then
        ipass=1
      end if

      end
      subroutine cheren(lpntn,touch,lusol,ntotv)
c**********************************************************************
c                 
c**** This routine checks the array lpntn
c                 
c**********************************************************************
      implicit none
      integer   ntotv,lusol,itotv,itotn
      integer   lpntn(ntotv)
      logical*1 touch(ntotv)

      call slogze(ntotv,touch)
      do itotv = 1,ntotv
        itotn = lpntn(itotv)
        if(itotn.gt.0) then
          if(touch(itotn)) then
            call srunen('CHEREN: THE RENUMBERING HAS FAILED',lusol)
          else
            touch(itotn)=.true.
          end if
        end if
      end do

      end
      subroutine condco(lpoin,lwork,lpote,npoit,npowt)
c*************************************************************************
c                 
c**** This routine performs the composition LWORK^-1(LPOIN(*)), where 
c**** LWORK is the first renumbering done to avoid condensed nodes in 
c**** the list of nodal points      
c
c*************************************************************************
      implicit none
      integer   npoit,npowt,ipoin,jpoin
      integer   lpote(npoit),lwork(npoit),lpoin(npoit)
c
c*** Inverse of LWORK stored in LPOTE
c
      do ipoin=1,npoit
        jpoin=lwork(ipoin)
        if(jpoin.gt.0) then
          lpote( jpoin)= ipoin
        else if(jpoin.lt.0) then
          lpote(-jpoin)=-ipoin
        end if
      end do
c
c*** Composition LPOTE(LPOIN(*)), using LWORK as bridge array
c
      do ipoin=1,npowt
        lwork(ipoin)=lpote(lpoin(ipoin))
      end do
      do ipoin=npowt+1,npoit
        lwork(ipoin)=lpote(ipoin)
      end do
      do ipoin=1,npoit
        lpoin(ipoin)=lwork(ipoin)
      end do

      end
      subroutine condre(nnode,nnodt,nnodw,nelem,npoit,npoiw,
     .                  npows,lnods,lnodw,lwork)
c*************************************************************************
c                 
c**** This routine modifies array LNODS to deal with bubble condensation.
c**** It is assumed that nodes with condensed degrees of freedom are
c**** given with a negative value in the original array LNODS, although
c**** the output is again positive. 
c
c*************************************************************************
      implicit none
      integer   nnode,nnodt,nnodw,nelem,npoit,npoiw,npows,
     .          ielem,inode,jnode,ipoin,jpoin,npoib,icoun
      integer   lnods(nnodt,nelem), lnodw(nnodt,nelem),
     .          lwork(npoit)
c
c*** Store LNODS in LNODW and put negative nodes of LNODS at the end with
c*** value 0
c
      call slvecz(npoit,lwork)
      do ielem=1,nelem
        do inode=1,nnodw
          lnodw(inode,ielem)=lnods(inode,ielem)
        end do
        inode=0
        do while(inode.lt.nnodw)
          inode=inode+1
          ipoin=lnods(inode,ielem)
          if(ipoin.lt.0) then
            do jnode=inode+1,nnodw
              lnods(jnode-1,ielem)=lnods(jnode,ielem)
            end do
            lwork(-ipoin)=ipoin
            lnods(nnodw,ielem)=0
            inode=inode-1
          end if
        end do
      end do
c
c*** Renumber of non-negative nodes 
c
      npoiw=0
      npows=0
      icoun=0
      do ielem=1,nelem
        do inode=1,nnodw
          ipoin=lnods(inode,ielem)
          if(ipoin.gt.0) then
            if(lwork(ipoin).eq.0) then
              if(inode.le.nnode) then
                npoiw=npoiw+1
              else if(inode.gt.nnode) then
                npows=npows+1
              end if
              icoun=icoun+1
              lwork(ipoin)=icoun
            end if
            lnods(inode,ielem)=lwork(ipoin)
          end if
        end do
      end do
c
c*** Give to negative (condensed) nodes the highest numbers in LWORK
c
      npoib=npoiw+npows
      do ipoin=1,npoit
        jpoin=lwork(ipoin)
        if(jpoin.lt.0) then
          npoib=npoib+1
          lwork(ipoin)=-npoib
        end if
      end do
      
      end
      subroutine notodf(lpntn,lpoin,iffix,ndofn,npoin,npois,
     .                  nequa)
c**********************************************************************
c                 
c**** This routine sets up array LPNTN from the array LPOIN, which 
c**** only contains the renumber of the nodes, not of the degrees of 
c**** freedom. LPOIN(I) contains the old node I e.g. LPOIN(3)=5 means
c**** that the new node 3 was the old node 5.
c                 
c**********************************************************************
      implicit  none
      integer   ndofn,npoin,npois,nequa,ntota,ipoin,joldf,ioldf,
     .          npres,kpoin
      integer   lpntn(npoin*ndofn+npois), lpoin(npoin+npois),
     .          iffix(npoin*ndofn+npois)
c
c*** Initialization
c      
      nequa=0
      npres=0
      ntota=ndofn*npoin
      do ipoin=1,npoin+npois
        kpoin=lpoin(ipoin)
c
c*** Condensed secondary degrees of freedom
c        
        if (kpoin.lt.-npoin) then
          ioldf=ntota-kpoin-npoin
          lpntn(ioldf)=0
c
c*** Condensed degrees of freedom associated to main nodes
c        
        else if(kpoin.lt.0) then
          joldf=(-kpoin-1)*ndofn
          do ioldf=joldf+1,joldf+ndofn
            lpntn(ioldf)=0
          end do
c
c*** Renumbering of degrees of freedom associated to main nodes
c          
        else if (lpoin(ipoin).le.npoin) then 
          joldf=(lpoin(ipoin)-1)*ndofn
          do ioldf=joldf+1,joldf+ndofn
            if (iffix(ioldf).ne.1) then
              nequa=nequa+1
              lpntn(ioldf)=nequa
            else
              npres=npres-1
              lpntn(ioldf)=npres
            end if
          end do
c
c*** Renumbering of secondary degrees of freedom 
c          
        else if (lpoin(ipoin).gt.npoin) then
          ioldf=ntota+lpoin(ipoin)-npoin
          if (iffix(ioldf).ne.1) then
            nequa=nequa+1
            lpntn(ioldf)=nequa
          else
            npres=npres-1
            lpntn(ioldf)=npres
          end if
        end if
      end do

      end
      
      subroutine optban(lpoin,lcone,numco,lolne,lneol,nnodt,
     .                  nnode,nelem,npoin,mxcon,nband)
c***************************************************************************   
c
c*** This routine optimizes the bandwidth.  LPOIN(I) contains the new node
c*** number for the old node I
c
c**************************************************************************** 
      integer   lpoin(npoin),numco(npoin),lcone(npoin,mxcon),
     .          lolne(npoin),lneol(npoin)

      do ipoin=1,npoin                                      ! Initializes LPOIN
        lpoin(ipoin)=ipoin
      end do

      do ipoin=1,npoin                                      ! Start the renum.
        do kpoin=1,npoin                                    ! with each nodal
          lolne(kpoin)=0                                    ! point. Initia-
          lneol(kpoin)=0                                    ! lizes LOLNE LNEOL
        end do
        iband=0
        lneol(1)=ipoin
        lolne(ipoin)=1
        nasig=1
        ineol=0
        do while(nasig.ne.npoin)
          ineol=ineol+1
          icone=0
          do while(icone.ne.numco(lneol(ineol)))
            icone=icone+1
            jpoin=lcone(lneol(ineol),icone)
            if (lolne(jpoin).eq.0) then
              nasig=nasig+1
              lolne(jpoin)=nasig
              lneol(nasig)=jpoin
              iband=max(iband,abs(ineol-nasig))
              if (iband.gt.nband) then
                icone=numco(lneol(ineol))
                nasig=npoin
              end if
            end if
          end do
        end do
        if (iband.lt.nband) then
          nband=iband
          do kpoin=1,npoin
            lpoin(kpoin)=lneol(kpoin)
          end do
        end if
      end do
      
      end
      subroutine renban(lnods,lpntn,lwork,iffix,nnodw,nnodt,
     .                  npoin,npoiw,npows,nelem,ndofn,npois,
     .                  ntotv,nequa,lures,krenu,lengp,conde,
     .                  svome,svomt)
c***************************************************************************
c
c*** This routine optimizes the bandwidth. LPOIN(I) contains the new node
c*** number for the old node I
c
c***************************************************************************
      implicit none
      integer   nnodw,nnodt,npoin,npoiw,npows,nelem,ndofn,
     .          npois,ntotv,nequa,lures,krenu,lengp,svome,
     .          svomt,npoit,ipoin,ielem,inode,mxele,mxcon,
     .          itota,nband,idumy
      integer   lnods(nnodt,nelem), lpntn(ntotv),
     .          lwork(npoin),       iffix(ntotv)
      integer   ipadd,iadr1,iadd1,iadd2,iadd3,iadd4,iadd5   ! pointer
      logical*1 conde
c
c*** If no renumbering is needed, set LPNTN = identity (with possibly
c*** negative values if there are condensed nodes)      
c      
      if (krenu.eq.0) then
        npoit=npoin+npois
        do ielem=1,nelem
          do inode=1,nnodw
            ipoin=lnods(inode,ielem)
            lpntn(abs(ipoin))=ipoin
          end do
        end do
        call sadmem(0,iadr1,lengp*npoit,svomt)              ! Memory for LPOIN
        call slveca(npoit,lpntn,%val(iadr1))
        call notodf(lpntn,%val(iadr1),iffix,ndofn,
     .              npoin,npois,nequa)
        svome = max(svome,svomt)
        call sadmem(2,iadr1,lengp*npoit,svomt)              ! Deallocate LPOIN
c
c*** Renumbering 
c
      else
        npoit=npoiw+npows
        call slvecz(npoit,lpntn)              ! Use LPNTN as an auxiliar vector
        do ielem=1,nelem                      ! to compute max. elements that
          do inode=1,nnodw                    ! a node belongs mxele
            ipoin=lnods(inode,ielem)
            if(ipoin.ne.0)
     .        lpntn(ipoin)=lpntn(ipoin)+1
          end do
        end do
        mxele=lpntn(1)
        do ipoin=2,npoiw+npows
          if (lpntn(ipoin).gt.mxele) mxele=lpntn(ipoin)
        end do
        mxcon=mxele*(nnodw-1)                 ! Max.of conectivities/node
        iadd1=0                               ! LPOIN(npoin+npois)
        iadd2=iadd1+(npoin+npois)*lengp       ! LCONE(npoiw+npows,mxcon)
        iadd3=iadd2+(npoiw+npows)*mxcon*lengp ! NUMCO(npoiw+npows)
        iadd4=iadd3+(npoiw+npows)*lengp       ! LOLNE(npoiw+npows)
        iadd5=iadd4+(npoiw+npows)*lengp       ! LNEOL(npoiw+npows)
        itota=iadd5+(npoiw+npows)*lengp
        call sadmem(0,ipadd,itota,svomt)      ! Allocate auxiliar memory
        call slvecz(itota/lengp,%val(ipadd))
        call setcon(lnods,%val(ipadd+iadd2),  ! Set conectivities NUMCO & LCONE
     .              %val(ipadd+iadd3),
     .              npoiw+npows,nnodt,nnodw,
     .              nelem,mxcon,nband)
        call optban(%val(ipadd+iadd1),        ! Optimizes bandwidth
     .              %val(ipadd+iadd2),
     .              %val(ipadd+iadd3),
     .              %val(ipadd+iadd4),
     .              %val(ipadd+iadd5),
     .              nnodt,nnodw,nelem,
     .              npoiw+npows,mxcon,nband)
        if(conde)
     .    call condco(%val(ipadd+iadd1),lwork,
     .                %val(ipadd+iadd2),npoin+npois,
     .                npoiw+npows)
        call notodf(lpntn,%val(ipadd+iadd1),
     .              iffix,ndofn,npoin,npois,
     .              nequa)
        svome = max(svome,svomt)
        call sadmem(2,ipadd,itota,svomt)      ! Deallocate auxiliar memory
      end if
      
      end
      subroutine renum0(lnods,lpntn,lword,nelem,nnode,nnodt,npoin,
     .                  lures,iwork)
c**********************************************************************
c                 
c**** This routine sets up array LPNTN if node renumbering for
c**** profile minimization is desired
c                 
c**********************************************************************
      implicit real*8(a-h,o-z)
      integer   lnods(nnodt,nelem), lpntn(npoin), iwork(*)
      logical   stamp
c
c*** Computes the connections between degrees of freedom
c
      n1 = 1
      call renum1(lnods,iwork(n1),nnode,nnodt,npoin,
     .            nelem,nposi)
c
c*** Auxiliar memory for renumbering
c
      n2 = n1 + nposi          !  nodal connections
      n3 = n2 + npoin          !  nstart
      n4 = n3 + npoin - n1     !  level
      if(n4.gt.lword) then
        write(lures,900) n4, lword
        call srunen('RENUM0: INSUFFICIENT WORK SPACE',lures)
      end if
c
c*** Renumbering routines
c
      call renum2(iwork(n1),iwork(n2),iwork(n3),lpntn,npoin,ns)
      call renum3(iwork(n1),iwork(n2),iwork(n3),lpntn,npoin,ns)
c
c*** Obtain the inverse of array lpntn
c
      do ipoin=1,npoin
        iwork(lpntn(ipoin))=ipoin
      end do
c
c*** Print renumbered nodes
c
      stamp=.false.
      if(stamp) then
        write(lures,901)
        do ipoin=1,npoin
          write(lures,902) ipoin,lpntn(ipoin),ipoin,iwork(ipoin)
        end do
      end if

      do ipoin=1,npoin
        lpntn(ipoin)=iwork(ipoin)
      end do
c
c*** Formats
c
  900 format(//,5x,67('#'),//,5x,
     .  '>>>> RENUMBERING MODULE REQUIRES MORE WORKING SPACE:',/,
     .  11x,'*** NUMBER OF REQUIRED  INTEGERS :',i10,/,
     .  11x,'*** NUMBER OF ALLOCATED INTEGERS :',i10   )
  901 format(//,5x,67('#'),//,5x,
     .  '>>>> RENUMBERED NODES :',//,
     .  11x,'OLD',10x,'NEW',10x,'NEW',10x,'OLD',/)
  902 format( 11x,   i10,3x,   i10,3x,   i10,3x,i10)

      end
      subroutine renum1(lnods,nodad,nnode,nnodt,npoin,
     .                  nelem,nposi)
c**********************************************************************
c
c*** Computes the connections between degrees of freedom (stored as 
c*** a linked list)
c
c**********************************************************************
      implicit real*8 (a-h,o-z)
      integer   lnods(nnodt,nelem), nodad(*)

      do ipoin = 1,npoin
         nodad(ipoin) = 0
      end do
      nfree = npoin+1

      do ielem = 1,nelem
        do inode = 1,nnode
          ipoin = lnods(inode,ielem)
          if(ipoin.ne.0) then
            do jnode = 1,nnode
              if(jnode.ne.inode) then
                jpoin = lnods(jnode,ielem)
                if(jpoin.ne.0) then
                   nadre = nodad(ipoin)
                   nadro = ipoin
                   do while (nadre.gt.0)
                     kpoin = nodad(nadre)
                     if(kpoin.eq.jpoin) go to 10
                     nadro = nadre + 1
                     nadre = nodad(nadro)
                   end do
                   nodad(nadro) = nfree
                   nodad(nfree) = jpoin
                   nodad(nfree+1) = 0
                   nfree = nfree + 2
   10              continue
                 end if
              end if
            end do
          end if
        end do
      end do
      nposi = nfree

      end


      subroutine renum2(nadj,nstart,lev,naux,nodes,ns)
c**********************************************************************
c
c**** Computes a set of posibles startings nodes
c
c**********************************************************************
      implicit real*8 (a-h,o-z)
      logical better
      dimension naux(*),nadj(*),lev(*),nstart(*)
c
c*** Begin iteration
c*** Select initial root node arbitrarily and generate its level
c*** structure
c
      iroot = 1
      better = .true.
      do while (better)

        call renum4(nstart,lev,idepth,nadj,iwidth,nodes,iroot,lhw)
c
c*** Create a list of nodes which are at maximum distance from root
c*** node and store the root 
c
        ns = iroot
c
c*** Loop over nodes at maximum distance from root node
c*** Generate level structure for each node
c*** Set switch if a level structure of greater depth occurs
c
        better = .false.

        do i = 1,lhw
          nnoded = nstart(i)
          call renum4(naux,lev,ndepth,nadj,nwidth,nodes,nnoded,lr)
          if(ndepth.ge.idepth) then
            if((ndepth.ne.idepth).or.(nwidth.lt.iwidth)) then
              iroot  = nnoded
              idepth = ndepth
              iwidth = nwidth
              better = .true.
            end if
          end if
        end do
      end do

      end

      subroutine renum3(nadj,nact,noda,newnn,nodes,i)
c**********************************************************************
c
c**** Resequence nodes for minimum profile
c
c**********************************************************************
      implicit real*8 (a-h,o-z)
      integer   nadj(*),newnn(*),nact(*),noda(*)
 
      large = 5**5
c
c*** King's scheme
c
      do j = 1,nodes
        newnn(j) = 0
        noda(j) = 0
      end do
      newnn(i) = 1
      nac = 0
c
c*** Negate all ndeg entries for nodes which are
c*** adjacent to starting node i
c
      maxfrt = 0
      nad = nadj(i)
      do while (nad.gt.0)
        maxfrt = maxfrt + 1
        npj = nadj(nad)
        if(noda(npj).eq.0) then
          nac = nac + 1
          noda(npj) = nac
          nact(nac) = npj
        end if
        nad = nadj(nad+1)
      end do
      noda(i) = large
c
c*** Loop over nodes to be renumbered
c
      do k = 2,nodes
        minnew = large
        lmin = large
c
c*** Loop over active nodes
c*** Skip to next node if old node is already renumbered
c
        do iact = 1,nac
          j = nact(iact)
          if(newnn(j).le.0) then
            new = -1
            min = large
c
c*** Compute the increment in active nodes for each node j
c*** Compute when this node was first active by checking for renumbered
c*** neighbours with lowest numbers
c
            nad = nadj(j)
            do while (nad.gt.0)
              n = nadj(nad)
              if(noda(n).eq.0) new = new + 1
              if(newnn(n).ne.0) then
                if(newnn(n).lt.min)min = newnn(n)
              end if
              nad = nadj(nad+1)
            end do
c
c*** Select node with smallest increment in active nodes
c*** in the case of a tie, select node which has been longest active
c
            if(new.le.minnew) then
              if((new.ne.minnew).or.(min.lt.lmin)) then
                minnew = new
                lmin = min
                next = j
              end if
            end if
          end if
        end do
c
c*** Renumber node and compute number of active nodes
c
        newnn(next) = k
        nif = nif+minnew
        if(nif.gt.maxfrt) maxfrt = nif
c
c*** Set nodes which are adjacent to the node just renumbered
c*** as actives nodes, deactivate next
c
        npos = noda(next)
        ilast = nact(nac)
        noda(ilast) = npos
        nact(npos) = ilast
        nac = nac - 1

        if(minnew .ne. -1) then
          nad = iabs(nadj(next))
          do while (nad.gt.0)
            n = nadj(nad)
            if(noda(n).eq.0) then
              nac = nac + 1
              noda(n) = nac
              nact(nac) = n
            end if
            nad = nadj(nad+1)
          end do
        end if
      end do

      end
      subroutine renum4(ndeg,lev,lsd,nadj,mlw,nodes,nroot,nloc)
c**********************************************************************
c
c**** Compute level structure rooted at nroot
c
c**********************************************************************
      implicit real*8 (a-h,o-z)
      logical*1 back
      integer   lev(*),ndeg(*),nadj(*)
c
c*** Initialization
c
      do i = 1,nodes
        lev(i) = 0
      end do
      lev(nroot) = 1    
      ndeg(1) = nroot
      back = .false.
      lsd = 1
      kount = 1
      nloc = 1
      nlocn = 0
      mlw = 1
c
c*** Assign levels
c
      do while(kount.lt.nodes)
        do il = 1,nloc
          if(.not.back) ip = ndeg(il)
          if(back) ip = ndeg(nodes+1-il)
          nad = nadj(ip)
          do while (nad.gt.0)
            jp = nadj(nad)
            if(lev(jp).eq.0) then
              nlocn = nlocn + 1
              lev(jp) = lsd + 1
              kount = kount + 1
              if(back) ndeg(nlocn) = jp
              if(.not.back) ndeg(nodes+1-nlocn) = jp
            end if
            nad = nadj(nad+1)
          end do
        end do
        nloc = nlocn
        nlocn = 0
        if(nloc.gt.mlw) mlw = nloc
        lsd = lsd + 1
        back = .not.back
      end do

      if(back) then
        nhalf = nodes/2
        do i = 1,nhalf
          nn = ndeg(i)
          ndeg(i) = ndeg(nodes+1-i)
          ndeg(nodes+1-i) = nn
        end do
      end if

      end


      subroutine renumb(
     .  lnods,iffix,nnode,nnodt,npoin,nelem,ndofn,
     .  ntotv,iprob)
c****************************************************************************
c
c**** This routine defines the address of the array LPNTN (ialpn) for
c**** the problem we are considering and sets up the components of this
c**** array according to the renumbering strategy used. The basic tree is
c
c     RENUMB                Main routine
c           CONDRE          Modification of LNODS (for condensed d.o.f.)
c           RENBAN          Renumbering for band storage
c                 SETCON    Set connectivities
c                 OPTBAN    Optimize bandwidth
c                 CONDCO    Undo the original modification for condensation
c           RENUMN          Renumbering for profile storage
c                 RENUM0    Optimize profile
c                 CONDCO    Undo the original modification for condensation
c           RENBAN          Renumbering for band storage (used for sparse,
c                           only to minimize the fill-in)      
c           SKYINI          Construct LPONT (profile storage)
c
c****************************************************************************
      implicit none
      include   'MatMan.h'
      integer    iwoso,p_lpn,p_lpo,p_lnw,p_tou              ! pointer
      integer    nnode,nnodt,npoin,nelem,ndofn,ntotv,iprob,
     .           npois,nnodw,nevat,npoiw,npows,ldumy,lbyts,
     .           inode,ielem,npoit,svomt,ipass
      integer    lnods(nnodt,nelem), iffix(ntotv)
      data       svomt,ipass/0,0/
      save       svomt,ipass
c
c*** Check parameters, allocate memory for pointers and initialize flags
c
      if(ipass.eq.0) then
        call almepo
        call setdat
        call slvecz(mprob,svome) 
        call slvecz(mprob,speme)
        call slvecz(mprob,sizms)
        ipass = 1
      end if
c
c*** Calculation of some integers
c
      npois=ntotv-ndofn*npoin                     ! Points with 1 d.o.f.
      npoiw=npoin                                 ! Main working points
      npows=npois                                 ! Secondary working points
      nnodw=nnode                                 ! Working nodes
      nevat=ndofn*nnode                           ! Working d.o.f.
      npoit=npoiw+npows                           ! Total number of points
      if (npois.gt.0) then
        nnodw=nnodt
        nevat=nevat+nnodt-nnode
      end if
c
c*** Memory allocation for LPNTN and LPONT
c
      if ((sofam(iprob).eq.1).or
     .  .(sofam(iprob).eq.3)) then
        lbyts=lengp*ntotv
        call sadmem(0,p_lpn,lbyts,speme(iprob))
        call sadmem(0,p_lpo,lbyts,speme(iprob))
        call slvecz(ntotv,%val(p_lpn))
        call slvecz(ntotv,%val(p_lpo))
        call spoias(
     .    p_lpn,%val(ialpn+(iprob-1)*lengp))
        call spoias(
     .    p_lpo,%val(ialpo+(iprob-1)*lengp))
      end if
c
c*** Modification of LNODS if some degrees of freedom have been
c*** condensed
c      
      if (conde(iprob).and.
     .  (krenu(iprob).gt.0)) then
        lbyts=lengp*nelem*nnodt
        call sadmem(0,p_lnw,lbyts,svomt)                    ! Memory for LNODW
        call condre(
     .    nnode,nnodt,nnodw,nelem,npoit,                    ! Use LPONT as 
     .    npoiw,npows,lnods,%val(p_lnw),                    ! auxiliar array
     .    %val(p_lpo))
      end if
c
c*** Renumber strategy according to the type of matrix storage for direct
c*** solvers (SOFAM = 1). For sk iterative solvers (SOFAM = 3) the routine
c*** RENBAN is used to initialize LPNTN, which is employed for identifying
c*** prescribed and free degrees-of-freedom      
c
      if(sofam(iprob).eq.1) then
        if (kmsip(iprob).eq.0) then                         ! BAND STORAGE
          call renban(
     .      lnods,%val(p_lpn),%val(p_lpo),
     .      iffix,nnodw,nnodt,npoin,npoiw,npows,
     .      nelem,ndofn,npois,ntotv,nequa(iprob),
     .      lusol,krenu(iprob),lengp,conde(iprob),
     .      svome(iprob),svomt)
        else if (kmsip(iprob).eq.1) then                    ! PROFILE
          call renumn(
     .      lnods,%val(p_lpn),%val(p_lpo), 
     .      iffix,nnodw,nnodt,npoin,npoiw,npows,
     .      nelem,ndofn,npois,ntotv,nequa(iprob),
     .      lusol,krenu(iprob),lengp,conde(iprob),
     .      svome(iprob),svomt)
        else if (kmsip(iprob).eq.2) then                    ! SPARSE STORAGE
          call renban(
     .      lnods,%val(p_lpn),%val(p_lpo),
     .      iffix,nnodw,nnodt,npoin,npoiw,npows,
     .      nelem,ndofn,npois,ntotv,nequa(iprob),
     .      lusol,krenu(iprob),lengp,conde(iprob),
     .      svome(iprob),svomt)
        end if
      else if(sofam(iprob).eq.2) then
        nequa(iprob)=ntotv
      else if(sofam(iprob).eq.3) then
        call renban(
     .    lnods,%val(p_lpn),%val(p_lpo),
     .    iffix,nnodw,nnodt,npoin,npoiw,npows,
     .    nelem,ndofn,npois,ntotv,nequa(iprob),
     .    lusol,    0,lengp,conde(iprob),
     .    svome(iprob),svomt)
      end if
c
c*** Recover LNODS and deallocate memory for LNODW, if needed
c
      if (conde(iprob).and.(krenu(iprob).gt.0)) then
        call slveca(nnodt*nelem,%val(p_lnw),lnods)
        svome(iprob) = max(svome(iprob),svomt)
        call sadmem(2,p_lnw,lengp*nnodt*nelem,svomt)
      end if
c
c*** Construct LPONT in the case of a slyline matrix storage
c      
      if ((sofam(iprob).eq.1).and.(kmsip(iprob).eq.1)) then
        call sadmem(0,iwoso,lengp*nevat,ldumy)
        call skyini(iprob,lnods,ndofn,nelem,        
     .              nequa(iprob),nevat,nnodt,nnodw, 
     .              npoin,ntotv,%val(p_lpn),
     .              %val(p_lpo),%val(iwoso))
        call sadmem(2,iwoso,lengp*nevat,ldumy)
      end if
c
c*** This routine can be called to check the renumber strategy
c
c-k   if(krenu(iprob).gt.0) then
c-k     call sadmem(0,p_tou,ntotv,svomt)
c-k     call cheren(%val(p_lpn),%val(p_tou),lusol,ntotv)
c-k     svome(iprob) = max(svome(iprob),svomt)
c-k     call sadmem(2,p_tou,ntotv,svomt)
c-k   end if
c
c*** Put all values of LNODS positive, if needed
c
      if(conde(iprob)) then
        do ielem=1,nelem
          do inode=1,nnodt
            lnods(inode,ielem)=abs(lnods(inode,ielem))
          end do
        end do
      end if

      end
      subroutine renumn(
     .  lnods,lpntn,lwork,iffix,nnode,nnodt,npoin,
     .  npoiw,npows,nelem,ndofn,npois,ntotv,nequa,
     .  lures,krenu,lengp,conde,svome,svomt)
c**********************************************************************
c                 
c**** This routine sets up array LPNTN. In this routine only the nodes
c**** are renumbered, so that the array LPNTN has the NDOFN components
c**** equal.
c                 
c**********************************************************************
      implicit  none
      integer   iadr1,iadr2                                 ! pointer
      integer
     .  nnode,nnodt,npoin,npoiw,npows,nelem,ndofn,
     .  npois,ntotv,nequa,lures,krenu,lengp,svome,
     .  svomt,lword,npoit,ipoin,lbytr,inode,ielem
      integer   lnods(nnodt,nelem), lpntn(ntotv),
     .          lwork(npoin),       iffix(ntotv)
      logical*1 conde
c
c*** If no renumbering is needed, set LPNTN = identity (with possibly
c*** negative values if there are condensed nodes)      
c      
      npoit=npoin+npois
      call sadmem(0,iadr1,lengp*npoit,svomt)                ! Memory for LPOIN
      if (krenu.eq.0) then
        if(conde) then
          npoit=npoin+npois
          do ielem=1,nelem
            do inode=1,nnodt
              ipoin=lnods(inode,ielem)
              if(ipoin.ne.0) lpntn(abs(ipoin))=ipoin
            end do
          end do
        else
          do ipoin = 1,npoin
            lpntn(ipoin) = ipoin
          end do
        end if
        call slveca(npoit,lpntn,%val(iadr1))
      end if

      if(krenu.eq.1) then
        npoit=npoiw+npows
c
c*** Node renumbering for profile minimization, only nodes are renumbered
c
        lword=(nnode*(nnode-1))*nelem*2+3*npoit+1
        lbytr=lengp*lword                                   ! Memory in bytes
        call sadmem(0,iadr2,lbytr,svomt)                    ! Allocate memory
        call renum0(lnods,%val(iadr1),lword,nelem,
     .              nnode,nnodt,npoit,lures,
     .              %val(iadr2))
        if(conde)
     .    call condco(%val(iadr1),lwork,%val(iadr2),
     .    npoin+npois,npoiw+npows)
        svome = max(svome,svomt)
        call sadmem(2,iadr2,lbytr,svomt)                    ! Deallocate memory
      end if
c
c*** Constructs the rest of components of LPNTN 
c
      call notodf(lpntn,%val(iadr1),iffix,ndofn,npoin,
     .            npois,nequa)
      svome = max(svome,svomt)
      call sadmem(2,iadr1,lengp*npoit,svomt)                ! Deallocate memory

      end
      subroutine setcon(lnods,lcone,numco,npoin,nnodt,nnode,nelem,
     .                  mxcon,nband)
c*************************************************************************
c
c*** This routine builds the conectivity arrays NUMCO and LCONE,
c*** and computes the initial bandwidth. The arrays are:
c*** NUMCO(I)   contains the number of nodes related with node I
c*** LCONE(I,J) contains the identities of those nodes
c      
c*************************************************************************
      logical   assig
      integer   lnods(nnodt,nelem),numco(npoin),lcone(npoin,mxcon)

      nband=0
      do ipoin=1,npoin
        numco(ipoin)=0
      end do
      
      do ielem=1,nelem
        do inode=1,nnode
          ipoin=lnods(inode,ielem)
          if(ipoin.ne.0) then
            do jnode=1,nnode
              jpoin=lnods(jnode,ielem)
              if ((ipoin.ne.jpoin).and.(jpoin.ne.0)) then
                assig=.true.
                icone=1
                do while((assig).and.(icone.le.numco(ipoin)))
                  if (lcone(ipoin,icone).eq.jpoin) assig=.false.
                  icone=icone+1
                end do
                if (assig) then
                  numco(ipoin)=icone
                  lcone(ipoin,icone)=jpoin
                  iband=abs(ipoin-jpoin)
                  if (iband.gt.nband) nband=iband
                end if
              end if
            end do
          end if
        end do
      end do

      do ipoin=1,npoin
        do icone=1,numco(ipoin)-1
          do jcone=icone+1,numco(ipoin)
            if (numco(lcone(ipoin,icone)).gt.
     .        numco(lcone(ipoin,jcone))) then
              itemp=lcone(ipoin,icone)
              lcone(ipoin,icone)=lcone(ipoin,jcone)
              lcone(ipoin,jcone)=itemp
            end if
          end do
        end do
      end do

      end
      subroutine setdat
c**********************************************************************
c                 
c**** This routine initializes some data
c                 
c**********************************************************************
      implicit none
      include   'MatMan.h'
      integer kprob
c
c*** Type of solver family
c      
      do kprob=1,mprob
        if (     (kites(kprob).gt.10)
     .    .and.  (kites(kprob).le.20)) then                 ! Direct solvers
          sofam(kprob) = 1
        else if( (kites(kprob).gt.20)
     .      .and.(kites(kprob).le.30)) then                 ! EE Iterative solvers
          sofam(kprob) = 2
        else if( (kites(kprob).gt.30)
     .      .and.(kites(kprob).le.40)) then                 ! SK Iterative solvers
          sofam(kprob) = 3
        else
          if(kites(kprob).ne.0) 
     .      call srunei(
     .      'SETDAT: WRONG TYPE OF SOLVER FOR PROBLEM',
     .      kprob,lusol)
        end if
c
c*** Initialization of some flags
c
        if((sofam(kprob).eq.2).and.(kmsip(kprob).gt.0))
     .    kmsip(kprob) = -1                                 ! No prec. by default
      end do
      
      end
      subroutine setldo(lnods,ldofn,nnodt,nnode,ndoft,ndofn,
     .                  nelem,npoin)
c**********************************************************************
c                 
c*** Build LDOFN
c                 
c**********************************************************************
      implicit real*8(a-h,o-z)
      integer   lnods(nnodt,nelem), ldofn(ndoft,nelem)

      do ielem=1,nelem
        idofn=0
        do inode=1,nnode
          ipoin=lnods(inode,ielem)
          if(ipoin.ne.0) then
            do kdofn=1,ndofn
              idofn=idofn+1
              ldofn(idofn,ielem)=(ipoin-1)*ndofn+kdofn
            end do
          end if
        end do
      end do
      
      if (ndofn*nnode.ne.ndoft) then
        ntova=npoin*ndofn
        do ielem=1,nelem
          idofn=nnode*ndofn
          do inode=nnode+1,nnodt
            ipoin=lnods(inode,ielem)
            if(ipoin.ne.0) then
              idofn=idofn+1
              ldofn(idofn,ielem)=ntova+ipoin-npoin
            end if
          end do
        end do
      end if
      
      end
      
      subroutine skyini(
     .  iprob,lnods,ndofn,nelem,neqns,nevac,
     .  nnodt,nnodw,npoin,ntotv,lnueq,lpont,
     .  leqns)
c***********************************************************************
c
c**** This routine performs the skyline solver initialization when
c**** renumbering is employed
c
c.... Input parameters
c
c     lnods(nnodt,nelem)  -  element connectivity
c     lnueq(ntotv)        -  equation numbers
c     ndofn               -  no. of d.o.f per node
c     nelem               -  no. of elements
c     neqns               -  number of equations
c     nnodt               -  number of element nodes
c     npoin               -  number of points
c     ntotv               -  number of variables
c     nnodw               -  number of nodes to work
c
c.... Output parameters
c
c     lpont(neqns)        -  position of the term just above the diagonal
c                            in the arrays gstup & gstlo
c
c.... Auxiliary working parameters
c
c     leqns(nevac)        -  auxiliary array to store the equation
c                            numbers of the current element
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      include 'MatMan.h'
      integer
     .  lpont(neqns), leqns(nevac),
     .  lnueq(ntotv), lnods(nnodt,nelem)
c
c*** Find the position of the diagonal terms
c
      do ieqns=1,neqns
        lpont(ieqns)=0
      end do

      maxpr = 0
      do ielem = 1,nelem
        mneqn = 0
        nevac = 0
        do inode = 1,nnodw
          ipont = lnods(inode,ielem)
          if (ipont.gt.0) then
            if (ipont.le.npoin) then
              itotv=(ipont-1)*ndofn
              do idofn = 1,ndofn
                itotv=itotv+1
                ieqns = lnueq(itotv)
                if(ieqns.gt.0) then
                  if(mneqn.eq.0) mneqn = ieqns
                  mneqn = min0(mneqn,ieqns)
                  nevac = nevac + 1
                  leqns(nevac) = ieqns
                end if
              end do
            else if (ipont.ne.0) then
              itotv=npoin*ndofn+(ipont-npoin)
              ieqns = lnueq(itotv)
              if(ieqns.gt.0) then
                if(mneqn.eq.0) mneqn = ieqns
                mneqn = min0(mneqn,ieqns)
                nevac = nevac + 1
                leqns(nevac) = ieqns
              end if
            end if
          end if
        end do

        if(nevac.gt.0) then
          do ievac = 1,nevac
            ieqns = leqns(ievac)
            npont = max0(lpont(ieqns),ieqns-mneqn)
            lpont(ieqns) = npont
            maxpr = max0(maxpr,lpont(ieqns))
          end do
        end if
      end do ! ielem = 1,nelem
c
c*** Compute diagonal pointers for profile
c
      lpont(1) = 0
      if(neqns.gt.1) then
        do ieqns = 2,neqns
          lpont(ieqns) = lpont(ieqns) + lpont(ieqns-1)
        end do
      end if

      end
      subroutine MatMan(
     .  lnods,iffix,rhsid,unkno,ndofn,nnode,nnodt,
     .  nnodw,nelem,npoin,ntotv,nsist,lures,kfact,
     .  kimpo,iprob,inte1,inte2,real1,itask)
c****************************************************************************
c
c**** This routine drives the library to deal with a linear system
c**** of equations. The tasks are:
c
c     itask = 1  Renumber of the equations
c     itask = 2  Memory allocation for the system matrix of problem IPROB      
c     itask = 3  Construction of the matrix
c     itask = 4  Solution of a linear system
c     itask = 5  Deallocate memory for problem IPROB
c
c     Ramon Codina
c     Last revision: March 1997 
c      
c****************************************************************************
      implicit   none
      include   'MatMan.h'
      integer
     .  ndofn,nnode,nnodt,nnodw,nelem,npoin,ntotv,
     .  lures,iprob,inte1,nsist,itask,kfact,kimpo
      real*8     timei,timef
      integer
     .  lnods(nnodt,nelem),   iffix(ntotv),
     .  inte2(*)
      real*8
     .  rhsid(ntotv,nsist),   unkno(ntotv,nsist),
     .  real1(*)
c
c***  Renumber
c      
      if(itask.eq.1) then
        lusol=lures                 
        call scputi(timei) 
        call renumb(
     .    lnods,iffix,nnode,nnodt,npoin,nelem,ndofn,
     .    ntotv,iprob)
        call scputi(timef)
        real1(1) = timef-timei 
c
c***  Memory allocation
c      
      else if(itask.eq.2) then
        if(nequa(iprob).gt.0) then
          lusol=lures                            
          call memmat(
     .      lnods,ndofn,nnode,nnodt,nnodw,nelem,
     .      npoin,ntotv,nsist,kimpo,iprob)
        end if
c
c***  Construction of matrices
c
      else if(itask.eq.3) then
        if(nequa(iprob).gt.0) then
          call builma(
     .      iprob,real1,inte1,inte2,nelem,ndofn,
     .      nnode,nnodw)
        end if
c
c***  Solution of the linear system
c
      else if(itask.eq.4) then
        if(nequa(iprob).gt.0) then
          lusol=lures                            
          call scputi(timei)
          call solsys(
     .      iffix,lnods,rhsid,unkno,ntotv,nelem,nnode,
     .      nnodt,nnodw,ndofn,npoin,nsist,kfact,iprob,
     .      inte1,inte2,timei)
          call scputi(timef)
          real1(1) = timef-timei
        end if                  
c
c***  Deallocate memory
c
      else if(itask.eq.5) then
        call delmem(iprob,nelem)
      else
        call srunen('MATMAN: WRONG TASK NUMBER',lures)
      end if

      end
      subroutine sadmem(ioptn,iaddr,lbyts,memor)
c***********************************************************************
c*
c*    This routine allocates and releases virtual memory. The length of 
c*    the pointers must be that of the machine. This must be declared 
c*    at the moment of compiling.
c*
c***********************************************************************
      implicit none
      include 'MatMan.h'
      integer  iaddr                                        ! pointer
      integer  ioptn,lbyts,memor
      integer  smemal
      external smemal
      
      if (ioptn.eq.0) then                                  ! Memory allocation
        memor = memor+lbyts
        iaddr = smemal(lbyts)                                   
        if (iaddr.eq.0. or .lbyts.lt.0)                      
     .    call srunen(
     .    'SADMEM: ERROR WHEN CALLING MALLOC  ',
     .    lusol)                                             
      else if (ioptn.eq.1) then                             ! To be used for 
        call srunen('SADMEM: TASK NOT READY ',lusol)        ! re-allocation
      else if (ioptn.eq.2) then                             ! Release of memory
        memor = memor-lbyts
        call smemde(iaddr)                                 
      end if                                                

      end
      subroutine scputi(rtime)
c***********************************************************************
c*
c**** This routine finds out the CPU time in seconds
c*
c***********************************************************************
      implicit none
      real*8 rtime
      real*4 stime(2)
      real*4 ETIME
      external ETIME

      rtime = ETIME(stime)

      end
      subroutine slogze(n,v)
c*****************************************************************************
c
c**** This routine initialises a logical vector to false
c
c*****************************************************************************
      implicit none
      integer    n,i
      logical*1  v(n)
      
      do i=1,n
        v(i)=.false.
      end do
      
      end
      subroutine slveca(n,v1,v2)
c*****************************************************************************
c
c***  Vector assign:    v2(i) = v1(i)   i=1..n
c
c*****************************************************************************
      implicit none
      integer   n,i
      integer   v1(n),v2(n)
      
      do i=1,n
        v2(i)=v1(i)
      end do

      end
      subroutine slvecz(n,v)
c*****************************************************************************
c
c**** This routine initialises a vector to zero
c
c*****************************************************************************
      implicit none
      integer   n,i
      integer   v(n)
      
      do i=1,n
        v(i)=0
      end do
      
      end
      subroutine spoias(p1,p2)
c***********************************************************************
c
c**** This routine assigns pointers P2 <- P1
c
c***********************************************************************
      implicit none
      integer  p1,p2                                        ! pointer

      p2=p1

      end
      subroutine srunei(message,integ,lures)
c*************************************************************************
c
c**** This routine stops the run and writes a message with an integer
c
c*************************************************************************
      implicit none
      integer       integ,lures
      character     message*(*)

      write(lures,100) message,integ
      stop ' *** SOLVER INTERRUPTED *** '
      
  100 format(//,5x,67('#'),//,25x,'AN ERROR HAS BEEN',
     .  ' DETECTED BY MATMAN',//,5x,a,1x,i5,//,5x,67('#'),//)
      
      end
      subroutine srunen(message,lures)
c*************************************************************************
c
c**** This routine stops the run 
c
c*************************************************************************
      implicit none
      integer       lures
      character     message*(*)

      write(lures,100) message
      stop ' *** SOLVER INTERRUPTED *** '
      
  100 format(//,5x,67('#'),//,25x,'AN ERROR HAS BEEN',
     .  ' DETECTED BY MATMAN',//,5x,a,//,5x,67('#'),//)
      
      end
      subroutine sruner(message,realn,lures)
c*****************************************************************************
c
c**** This routine stops the run and writes a message with a floating number
c
c*****************************************************************************
      implicit none
      integer       lures
      real*8        realn
      character     message*(*)

      write(lures,100) message,realn
      stop ' *** SOLVER INTERRUPTED *** '
      
  100 format(//,5x,67('#'),//,25x,'AN ERROR HAS BEEN',
     .  ' DETECTED BY MATMAN',//,5x,a,1x,f10.5//,5x,67('#'),//)
      
      end
      subroutine srveca(n,v1,v2)
c*****************************************************************************
c
c***  Vector assign:    v2(i) = v1(i)   i=1..n
c
c*****************************************************************************
      implicit none
      integer   n,i
      real*8    v1(n),v2(n)
      
      do i=1,n
        v2(i)=v1(i)
      end do

      end
      subroutine svecad(n,v1,v2,v3)
c*****************************************************************************
c
c**** Vector addition:    v3(i) = v2(i) + v1(i)   i=1..n
c
c*****************************************************************************
      implicit none
      integer   n,i
      real*8    v1(n),v2(n),v3(n)
      
      do i=1,n
	 v3(i) = v2(i) + v1(i)
      end do

      end
      subroutine svecbs(n,v,w,x)
c*****************************************************************************
c
c**** This routine multiplies a vector V by a scalar X
c
c*****************************************************************************
      implicit none
      integer    n,i
      real*8     v(n),w(n),x
      
      do i=1,n
        w(i)=x*v(i)
      end do
      
      end
      subroutine sveccl(n,x1,v1,x2,v2,v)
c*****************************************************************************
c
c**** Linear combination of two vectors: v = x1 * v1 + x2 * v2
c
c*****************************************************************************
      implicit none
      integer   n,i
      real*8    x1,x2,v1(n),v2(n),v(n)
      
      do i=1,n
        v(i) = x1*v1(i) + x2*v2(i)
      end do

      end
      subroutine svecco(a,x,b,y,z,n)
c*****************************************************************************
c
c**** This routine computes the linear combination:
c
c     z = a*x + b*y      
c
c*****************************************************************************
      implicit none
      integer    n,i
      real*8     x(n),y(n),z(n),a,b
      
      do i=1,n
        z(i)=a*x(i)+b*y(i)
      end do
      
      end
      real*8 function svecdo(avect,bvect,nterm)
c***********************************************************************
c
c**** This function performs the vector dot product of the vectors
c**** avect(nterm) and bvect(nterm)
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      real*8 avect(nterm), bvect(nterm)

      svecdo=0.0
      do iterm=1,nterm
        svecdo=svecdo+avect(iterm)*bvect(iterm)
      end do

      end
      subroutine svecpu(x,n,v)
c*****************************************************************************
c
c**** This routine initialises a vector to a value x
c
c*****************************************************************************
      implicit none
      integer    n,i
      real*8     x,v(n)
      
      do i=1,n
        v(i)=x
      end do
      
      end
      subroutine svecze(n,v)
c*****************************************************************************
c
c**** This routine initialises a vector to zero
c
c*****************************************************************************
      implicit none
      integer    n,i
      real*8     v(n)
      
      do i=1,n
        v(i)=0.0
      end do
      
      end
