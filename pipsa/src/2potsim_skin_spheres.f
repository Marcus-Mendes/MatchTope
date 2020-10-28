c  program PIPSA    
c  PIPSA version 4.0.2 (2019)
c
c  Authors: R.R.Gabdoulline, Stefan Richter, R.C.Wade
c  Copyright (c) 2004,2019
c  HITS gGmbH, Schloss-Wolfsbrunnenweg 33, 69118 Heidelberg,
c  Germany, and
c  European Molecular Biology Laboratory, Heidelberg, Germany
c
c  This copyright notice must be attached to all copies, or extracts,
c  of this software.
c
c  Please obtain a license before using this software from the PIPSA
c  website or "mcmsoft@h-its.org".
c  Conditions of usage are as specified in the license agreement given 
c  as LICENSE.TXT in this distribution
c
c  Queries about problems with this software should be addressed to
c  "mcmsoft@h-its.org".
c
c  References:
c  N. Blomberg, R.R. Gabdoulline, M. Nilges and R. C. Wade.
c   Classification of protein sequences by homology modeling and
c   quantitative analysis of electrostatic similarity.
c   Proteins: Str., Function and Genetics, 37:379-387 (1999)
c  Wade,R.C., Gabdoulline, R.R. and De Rienzo, F.
c   Protein Interaction Property Similarity Analysis
c   Intl. J. Quant. Chem., 83:122-127 (2001).
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c
      include "maxdim.inc"
CccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
c
c Computes similarity of 2 potential grids on the molecular skin
c and around the specified point
c
c Input:
c - command line arg -g1 - the name of the file with potential grid for protein 1 
c   in UHBD format, default grd1.grd
c - command line arg -g2 - the name of the file with potential grid for protein 2
c   in UHBD format, default grd2.grd
c - command line arg -p1 - the name of the file with atom coordinates of protein 1
c   in PDB format, default pdb1.pdb
c - command line arg -p2 - the name of the file with atom coordinates of protein 2
c   in PDB format, default pdb2.pdb
c  - command line arg -pa - the name of the file with the list of directions
c   default name "spheres" - a list of (up to 99) sphere center coordinates and radii
c   to define the region where similarity is calculated
c   format of the file "spheres" is: xsc,radius where xsc is 3 coordinates 
c   (in &Aring;) of the sphere center and radius is the value of its radius
c - command line arg -pr - probe radius, default is 3 A
c - command line arg -sk - skin thickness, default is 4 A
c + default values are used if no input given
c
c The program will construct 2 skins (for protein 1 and 2) having thickness
c "skin" and at distance "probes" from van der Waals surface of the proteins.
c (i.e. from "probes" to "probes+skin" distance), using the points of the
c potential grids.  The potential values outside this skin will never be used.
c The followings will be computed:
c np1  = no of points of the skin of the protein 1
c np2  = no of points of the skin of the protein 2
c npoi = no of points of the intersecion of 2 skins
c aa0  = square of the norm of the grid 1 on its skin
c aa   = square of the norm of the grid 1 on intersection of the skins
c bb0  = square of the norm of the grid 2 on its skin
c bb   = square of the norm of the grid 2 on intersection of the skins
c ab   = scalar product of 2 potentials (on intersection of skins)
c amb  = average difference of the electrostatic potentials a and b (a-b)/npoi
c ambl = log(sum[exp(a)]/sum[exp(b)])
c ambm = log(sum[exp(-a)]/sum[exp(-b)])
c
c  Output:
c - fort.66 - some info about constructed skins
c - fort.6 - standard output has following data in one line:
c
c   si_hodgkin = 2*ab/(aa+bb)
c   si_carbo   = ab/sqrt(aa*bb)
c   aa
c   bb
c   ab
c   aa0
c   bb0
c   si_hodgkin_shape = 2.*float(npoi)/float(np1+np2)
c   si_carbo_shape   = float(npoi)/sqrt(float(np1*np2))
c   np1
c   np2
c   npoi
c   amb
c   ambl
c   ambm
c + this info is printed for every sphere defined in the file "spheres"
c
CccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
c
c 20.01.2003 - corrected wrong calculation of np2 and bb0
c
CccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
      character*256 fpdb1,fpdb2,fgrd1,fgrd2,arg,fspheres
      integer iform
      integer im, jm, km
      dimension xsc(3,999),ar(999),rsc(3)
      dimension grid1(im3_max)
      dimension grid2(im3_max)
      integer*1 isk1(im_max,im_max,im_max)
      integer*1 isk2(im_max,im_max,im_max)
      integer*1 isk1a(im_max,im_max,im_max)
      integer*1 isk2a(im_max,im_max,im_max)
	  character*54 str
      common/ius/iulog
	  common/probc/probes,skin
c
	  iupdb=8
	  iugrd=10
      iuparam=88
      iusphere=77
      iulog=66
c
c. hard coded here change if need
      fpdb1='pdb1.pdb'
      fpdb2='pdb2.pdb'
      iform=0
      fgrd1='grd1.grd'
      fgrd2='grd2.grd'
      fspheres='spheres'
      probes=3.0
      skin=4.0
c
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: 2potsim_skin_spheres -g1 grd1.grd -g2 grd2.grd ',
     ,'-p1 pdb1.pdb -p2 pdb2.pdb -pa spheres -pr probe -sk skin')
      call exit
      end if
c
      nargs=iargc()
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-g1') then
      call getarg(ia+1,fgrd1)
      end if
      if(arg(1:3).eq.'-g2') then
      call getarg(ia+1,fgrd2)
      end if
      if(arg(1:3).eq.'-p1') then
      call getarg(ia+1,fpdb1)
      end if
      if(arg(1:3).eq.'-p2') then
      call getarg(ia+1,fpdb2)
      end if
      if(arg(1:3).eq.'-pa') then
      call getarg(ia+1,fspheres)
      end if
      if(arg(1:3).eq.'-pr') then
      call getarg(ia+1,arg)
      read(arg,*) probes
      end if
      if(arg(1:3).eq.'-sk') then
      call getarg(ia+1,arg)
      read(arg,*) skin
      end if
      end do
c
c begin code
c
       open(iusphere,file=fspheres)
       nsphere=0
 71    continue
       read(iusphere,*,end=72) rsc,arr
       nsphere=nsphere+1
       ar(nsphere)=arr
        do m=1,3
        xsc(m,nsphere)=rsc(m)
        end do
        goto 71
 72     continue
        close(iusphere)
        write(iulog,*) 'will work on',nsphere,'  spheres'
c
c read the 3D grid
c
	  write(iulog,*)
	  write(iulog,*) 'Grid 1 file reading now... wait...'
      call inpgrd(fgrd1,iugrd,iform,im,jm,km,
     &             h,xo,yo,zo,grid1)
	  write(iulog,611) im,jm,km,h,xo,yo,zo
c
      call  mkskin(ifpdb,fpdb1,im,jm,km,h,xo,yo,zo,np1,isk1)
c
	  write(iulog,*)
	  write(iulog,*) 'Grid 2 file reading now... wait...'
      call inpgrd(fgrd2,iugrd,iform,im,jm,km,
     &             h,xo,yo,zo,grid2)
	  write(iulog,611) im,jm,km,h,xo,yo,zo
c
      call  mkskin(ifpdb,fpdb2,im,jm,km,h,xo,yo,zo,np2,isk2)
c
ccc
                DO isphere=1,nsphere
      do 421 i=1,im   
      do 422 j=1,jm   
      do 423 k=1,km   
      isk1a(i,j,k)=isk1(i,j,k)
      isk2a(i,j,k)=isk2(i,j,k)
 423  continue
 422  continue
 421  continue
                 do m=1,3
                 rsc(m)=xsc(m,isphere)
                 end do
                 scarr  = ar(isphere)**2
                 write(iulog,'(7f8.3)') rsc,scarr
                call rebuild(im,jm,km,h,xo,yo,zo,
     ,          rsc,scarr,isk1a,isk2a)
      call rearra(grid1,grid2,im,jm,km,h,
     &                  xo,yo,zo,isk1a,isk2a,np1a,np2a,
     &                  npoi,aa,bb,ab,aa0,bb0,amb,ambl,ambm)
	  write(iulog,611) im,jm,km,h,xo,yo,zo
	  write(iulog,*) 'No of points on skin:',npoi
 611  format('grid dimensions:  ',i3,2('.',i3),'  spacing:',f8.3,/,
     ,'grid origin:',3f8.3)
       pnp12=float(np1a)*float(np2a)
       write(6,612) 2*ab/(aa+bb),ab/sqrt(aa*bb),aa,bb,ab,aa0,bb0
     ,,2.*float(npoi)/float(np1a+np2a)
     ,,float(npoi)/sqrt(pnp12),np1a,np2a,npoi,amb,ambl,ambm
                 END DO
 612  format(2f7.3,5e15.7,2f7.3,3i7,5e15.7)
c
 1000  continue
      stop
   79 format(6e15.5)
   80 format(e15.5)
   81 format(3f8.2)
   83 format(10f8.3)
c
      end
      subroutine rearra(grid1,grid2,im,jm,km,h,
     &                  xo,yo,zo,isk1,isk2,np1a,np2a,
     &                  npoi,aa,bb,ab,aa0,bb0,amb,ambl,ambm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include "maxdim.inc"
      integer im, jm, km
	  integer*1 isk1(im_max,im_max,im_max)
	  integer*1 isk2(im_max,im_max,im_max)
      dimension grid1(1)
      dimension grid2(1)
      double precision daa,dbb,dab
      double precision daa0,dbb0,damb,dambla,damblb
      double precision dambma,dambmb
      common/ius/iulog
c
              write(iulog,81) xo, yo, zo
   81 format('rearrangement routine, grid origin:',3f8.3)
c
c ijk transform to n as n=i+(j-1)*im+(k-1)*im*jm
      imjm=im*jm
	  npoi=0
	  np1a=0
	  np2a=0
	  g1max=-1.e+20
	  g1min= 1.e+20
	  g2max=-1.e+20
	  g2min= 1.e+20
       daa=0.d0
       dbb=0.d0
       dab=0.d0
       daa0=0.d0
       dbb0=0.d0
       damb=0.d0
       dambla=0.d0
       damblb=0.d0
       dambma=0.d0
       dambmb=0.d0
      do 421 i=1,im   
      do 422 j=1,jm   
      do 423 k=1,km   
      nn=i+(j-1)*im+(k-1)*imjm
      if(isk1(i,j,k).ne.0) then
        daa0=daa0+grid1(nn)**2
        np1a=np1a+1
      end if
      if(isk2(i,j,k).ne.0) then
         dbb0=dbb0+grid2(nn)**2
         np2a=np2a+1
      end if
      if(isk1(i,j,k).eq.0) goto 410
      if(isk2(i,j,k).eq.0) goto 410
      npoi=npoi+1
      daa=daa+grid1(nn)**2 
      dbb=dbb+grid2(nn)**2 
      dab=dab+grid1(nn)*grid2(nn)
      damb=damb+grid1(nn)-grid2(nn)
      dambla=dambla+exp(grid1(nn))
      damblb=damblb+exp(grid2(nn))
      dambma=dambma+exp(-grid1(nn))
      dambmb=dambmb+exp(-grid2(nn))
	  if(grid1(npoi).gt.g1max) g1max=grid1(nn)
	  if(grid1(npoi).lt.g1min) g1min=grid1(nn)
	  if(grid2(npoi).gt.g2max) g2max=grid2(nn)
	  if(grid2(npoi).lt.g2min) g2min=grid2(nn)
 410  continue
 423  continue
 422  continue
 421  continue
       aa=daa 
       bb=dbb 
       ab=dab 
       aa0=daa0
       bb0=dbb0
       amb=damb/dfloat(npoi)
       ambl=dlog(dambla/damblb)
       ambm=dlog(dambma/dambmb)
              write(iulog,82) g1min, g1max
              write(iulog,82) g2min, g2max
   82 format('rearrangement routine, grid minmax:',2e11.3)
      return
      end
      subroutine inpgrd( filnam,iounit,iform,im,jm,km,
     &                   h,ox,oy,oz,grid) 
      include "maxdim.inc"
c
c. This is essentially a copy from UHBD version 4.1
      character*250 filnam
      integer iounit, iform
      character*72 title
      integer im, jm, km
	  dimension grid(1)
      integer grdflg
      integer stderr
      common/ius/iulog
c
      real dum2, dum3, dum4, dum5, dum6, dum7, dum8
      integer idum1, idum2, idum3, idum4, idum5, idum6, idum7
      integer i, j, k, l, kx
      real xbuf(im_max*im_max*2)
	  stderr=6
c
      if ( iform .eq. 0 ) then
         open(iounit, file=filnam, form='unformatted', status='old' )
         read(iounit) title, scale, dum2, idum1, idum2,
     &                idum3, idum4, idum5, im, jm, km,
     &                h, ox, oy, oz,
     &                dum3, dum4, dum5, dum6, dum7, dum8, idum6, idum7
         if ( im*jm .gt. im_max*im_max*2 ) then
            write(stderr,99001)
            stop
c???
c???
         endif
         if ( im*jm*km .gt. im3_max ) then
            write(stderr,99002)
            stop
         endif
         do  10  k = 1 , km
            read(iounit) kx, im, jm
            read(iounit) (xbuf(l), l = 1 , im*jm)
            do  9  l = 1 , im*jm
c              grid((k-1)*im*jm+l) = xbuf(l) / scale
               grid((k-1)*im*jm+l) = xbuf(l)        
    9       continue
   10    continue
         close(iounit)
      else if ( iform .eq. 1 ) then
         open( iounit, file=filnam, form='formatted', status='old' )
         read(iounit,77) title, scale, dum2, idum1, idum2,
     &                   idum3, idum4, idum5, im, jm, km,
     &                   h, ox, oy, oz,
     &                   dum3, dum4, dum5, dum6, dum7, dum8,
     &                   idum6, idum7
         if ( im*jm*km .gt. im3_max ) then
            write(stderr,99002)
            stop
         endif
         do  20  k = 1 , km
            read(iounit,78) kx, im, jm
            read(iounit,79) ( grid((k-1)*im*jm+l), l=1,im*jm )
   20    continue
         close(iounit)
      endif
      return
c
c..........
c
c               Format statements.
c
   77 format(a72/2e12.6,5i7/3i7,4e12.6/4e12.6/2e12.6,2i7)
   78 format(3i7)
   79 format( (6(1x,e12.6)) )
99001 format(1x,'FATAL: inpgrd(001) buffer too small.' )
99002 format(1x,'FATAL: inpgrd(002) grid too large. ',
     .          'recompile with larger im_max value.' )
      end
	  subroutine mkskin(ifpdb,fpdb,im,jm,km,h,xo,yo,zo,npp,isk)
      include "maxdim.inc"
      character*80 string
	  character*250 fpdb
	  character*4 atom
	  character*1 an
	  dimension xv(3)
	  character*1 iexsk(im_max),chex,chsk,chsp
c. this array is an internal
	  integer*1 iex(im_max,im_max,im_max)
c. this is supposed to be used later
	  integer*1 isk(im_max,im_max,im_max)
      common/ius/iulog
	  common/probc/probes,skin
	  data atom /'ATOM'/,chex,chsk,chsp/' ','x',' '/
c---
c
      do 301 i=1,im_max
      do 302 j=1,im_max
      do 303 k=1,im_max
      iex(i,j,k)=0
      isk(i,j,k)=0
 303  continue
 302  continue
 301  continue
c
	  write(iulog,*)
	  write(iulog,*) 'Constructing fitting region...'
      open(ifpdb,file=fpdb,form='formatted',status='old')
      dexcl=h
	  rexmax=1.9d0+probes
      hdexcl=dexcl/2.d0
      klim=(rexmax+dexcl*1.5d0)/dexcl
c		xexh=xo-hdexcl
c		yexh=yo-hdexcl
c		zexh=zo-hdexcl
		xexh=xo
		yexh=yo
		zexh=zo
c---find the center of the protein 
      natoms=0
	  npoi=0
 1    continue
	  read(ifpdb,500,end=99) string
	  if(string(1:4).ne.atom) goto 1
      read(string(1:54),501) an,xv
 500  format(a)
 501  format(13X,A1,16X,3F8.3)
      natoms=natoms+1
        call van(an,probes,rexcl2)
        rexcl=sqrt(rexcl2)
        x1=xv(1)
        x1a=x1+hdexcl  
        x2=xv(2)
        x2a=x2+hdexcl  
        x3=xv(3)
        x3a=x3+hdexcl  
C          FIND HOME GRID CELL IH,JH,KH
        ih=(x1a-xo)/dexcl
		jh=(x2a-yo)/dexcl
		kh=(x3a-zo)/dexcl
C
        DO 423 KK=-KLIM,KLIM
        K=KK+KH
        IF(k.le.0.or.k.gt.km)GOTO 423
c..here and further we correct the shifting by hdexcl:
c        ZK=dfloat(K)*dexcl-hdexcl+zo
c..here is a concerted change with change in exchk.f
        ZK=dfloat(K)*dexcl+zo
        ZTEMP = ZK - X3
        IF(ZTEMP.GT.rexcl)GOTO 423
                DO 422 JJ=-KLIM,KLIM
                J=JJ+JH
                IF(j.le.0.or.j.gt.jm)GOTO 422
                YJ=dFLOAT(J)*dexcl +yo
                YTEMP = YJ - X2
                IF(YTEMP.GT.rexcl)GOTO 422
                        DO 421 II=-KLIM,KLIM
                        I=II+IH
                        IF(i.le.0.or.i.gt.im)GOTO 421
                        IF(IEX(I,J,K).EQ.1)GOTO 421
                        XI=dFLOAT(I)*dexcl +xo
                        XTEMP = XI - X1
                        DIST2 =
     &                  XTEMP * XTEMP + YTEMP * YTEMP + ZTEMP * ZTEMP
                        IF(DIST2.GT.rexcl2)GOTO 421
                        IEX(I,J,K)=1
						npoi=npoi+1
 421                    CONTINUE
 422            CONTINUE
 423    CONTINUE
      goto 1
 99   continue
CCC      CLOSE(UNIT=iusas)
	  dismx1=sqrt(dismx1)
	  write(iulog,602) natoms,npoi
 602  format('No of atoms used for excl grid       :',i6,/,
     ,       'No of points assigned value 1        :',i6)
CCC
CCC__version 2.23__ repeat now the same with increased by skin probe
      dexcl=h   
	  probsk=probes+skin
c..the first number is the maximal atom radius appearing in pdb
	  rexmax=1.9d0+probsk        
      hdexcl=dexcl/2.d0
      klim=(rexmax+dexcl*1.5d0)/dexcl
CCC..skip the center finding 
c---read pdb file
      rewind (ifpdb)
      natosk=0
	  npoisk=0
 11   continue
	  read(ifpdb,500,end=991) string
	  if(string(1:4).ne.atom) goto 11
      read(string(1:54),501) an,xv
      natosk=natosk+1
        call van(an,probsk,rexcl2)
        rexcl=sqrt(rexcl2)
        x1=xv(1)
        x1a=x1+hdexcl  
        x2=xv(2)
        x2a=x2+hdexcl  
        x3=xv(3)
        x3a=x3+hdexcl  
C          FIND HOME GRID CELL IH,JH,KH
        ih=(x1a-xo)/dexcl
		jh=(x2a-yo)/dexcl
		kh=(x3a-zo)/dexcl
C
        DO 4231 KK=-KLIM,KLIM
        K=KK+KH
        IF(k.le.0.or.k.gt.km)GOTO 4231
c..here and further we correct the shifting by hdexcl:
c        ZK=dfloat(K)*dexcl-hdexcl+oez
c..here is a concerted change with change in exchk.f
        ZK=dfloat(K)*dexcl+ zo
        ZTEMP = ZK - X3
        IF(ZTEMP.GT.rexcl)GOTO 4231
                DO 4221 JJ=-KLIM,KLIM
                J=JJ+JH
                IF(j.le.0.or.j.gt.jm)GOTO 4221
                YJ=dFLOAT(J)*dexcl + yo
                YTEMP = YJ - X2
                IF(YTEMP.GT.rexcl)GOTO 4221
                        DO 4211 II=-KLIM,KLIM
                        I=II+IH
                        IF(i.le.0.or.i.gt.im)GOTO 4211
                        IF(IEX(I,J,K).EQ.1)GOTO 4211
                        IF(ISK(I,J,K).EQ.1)GOTO 4211
                        XI=dFLOAT(I)*dexcl + xo
                        XTEMP = XI - X1
                        DIST2 =
     &                  XTEMP * XTEMP + YTEMP * YTEMP + ZTEMP * ZTEMP
                        IF(DIST2.GT.rexcl2)GOTO 4211
                        ISK(I,J,K)=1
						npoisk=npoisk+1
 4211                   CONTINUE
 4221           CONTINUE
 4231   CONTINUE
      goto 11
 991  continue
      CLOSE(UNIT=ifpdb)

	  write(iulog,6021) natosk,npoisk
      npp=npoisk
 6021 format('No of atoms used for prot 1 excl+skin grid:',i6,/,
     ,       'No of points assigned value 1             :',i6)
CCC__version 2.23__ repeat ends                                     
CCC
CCC...Now write the prot1 ptofile 
      WRITE(iulog,603)
 603  format('A cross section of excl. grid of pr. #1, xOy plane',/)
      K=km/2
	  jskip=(jm-78)/2
	  if(jskip.lt.0) jskip=0 
	  iskip=(im-78)/2
	  if(iskip.lt.0) iskip=0 
      DO J=jskip+1,jm-jskip
		do 12 i=1,im
        iexsk(i)=chsp
		if(isk(i,j,k).eq.1) iexsk(i)=chsk
		if(iex(i,j,k).eq.1) iexsk(i)=chex
 12     continue
      WRITE(iulog,13)(iexsk(i),i=iskip+1,im-iskip)
      CONTINUE
 13   FORMAT(80a1)  
      END DO
	  write(iulog,*)
      RETURN
      END
      subroutine van(a,srad,v)
      character*1 a
c..These are OPLS parameters used in UHBD
      V=1.9  !DEFAULT VDW RADIUS
      IF(A.EQ.'H')V=1.20
      IF(A.EQ.'E')V=1.200
      IF(A.EQ.'C')V=1.90
      IF(A.EQ.'N')V=1.625
      IF(A.EQ.'O')V=1.48
      IF(A.EQ.'P')V=1.90
      IF(A.EQ.'S')V=1.85
      IF(A.EQ.'F')V=0.64  !ASSUMES FE
      IF(A.EQ.'Z')V=1.38  !ASSUMES ZN
      IF(A.EQ.'B')V=1.95  !ASSUMES BR
      IF(A.EQ.'I')V=2.15
      V=V+SRAD  !ADD PROBE RADIUS
      V=V*V
      RETURN
      END
      subroutine rebuild(im,jm,km,h,xo,yo,zo,
     ,          xsc,scarr,isk1a,isk2a)
      include "maxdim.inc"
      integer im, jm, km
      integer*1 isk1a(im_max,im_max,im_max)
      integer*1 isk2a(im_max,im_max,im_max)
      dimension xsc(3)
      common/ius/iulog
      do 421 i=1,im   
      do 422 j=1,jm   
      do 423 k=1,km   
          if((isk1a(i,j,k).eq.0).and.
     ,       (isk2a(i,j,k).eq.0)) goto 410
          xx=i*h+xo
          yy=j*h+yo
          zz=k*h+zo
          scp=
     +(xx-xsc(1))**2+(yy-xsc(2))**2+(zz-xsc(3))**2
          if(scp.gt.scarr) then
          isk1a(i,j,k)=0
          isk2a(i,j,k)=0
          end if
 410  continue
 423  continue
 422  continue
 421  continue
      return
      end

