c------------------------------------------------------------------------
c  program grid2i      
c  version 4.23 (February 2000) 
c
c  Authors: R.R.Gabdoulline, R.C.Wade
c  Copyright (c) 2000
c  European Molecular Biology Laboratory
c  Meyerhofstr. 1, Postfach 10.2209
c  D-69012, Heidelberg, Germany
c
c  Please send your contact address to get information on updates and
c  new features to "gabdoull@EMBL-Heidelberg.DE". Questions will be 
c  answered as soon as possible.
c
c  Please obtain a license before using this software from 
c  "wade@EMBL-Heidelberg.DE".  Conditions of usage are as specified 
c  in the license.
c
c  References:
c     R.R. Gabdoulline, R.C. Wade.     
c  Brownian dynamics simulation of protein-protein diffusional encounter.
c  (1998) Methods, 14, 329-341. 
c------------------------------------------------------------------------
c
      integer im, jm, km, grdflg, mxgrd, one, kk
      parameter (mxgrd=110)
      integer iunit, iform
      character*72 title
      character*80 argum
      REAL h, ox, oy, oz, grid(mxgrd,mxgrd,mxgrd), scale
      REAL x0, y0, z0, x1, y1, z1, xmax, ymax, zmax, 
     &     xstart, ystart, zstart, xend, yend, zend
      REAL dum2, dum3, dum4, dum5, dum6, dum7, dum8
      integer idum2, idum3, idum4
      integer i, j, k, ounit
      character*80 newfle
      character*80 flnm
      character*4 stem
      character*1 string
c
      nargum=iargc()
      if(nargum.eq.0) then
      write(6,501)
      goto 999
      end if
 501  format(1x,'grid2i - convert map from GRID to InsightII format',/,
     ,       1x,'usage: grid2i A_or_B GRID_file Insight_file',/,
     ,       1x,'where: A_or_B is A for ascii and B for binary grid')
      call getarg(1,argum)      
      if(argum(1:1).eq.'A'.or.argum(1:1).eq.'a') then
        iform = 1
      else
        iform = 0
      endif
      call getarg(2,argum)
      flnm=argum
      call getarg(3,argum)
      newfle=argum
c
c Statement functions
c begin code
c
      write(6,100) 
100   format(1x,'CONVERT MAP FROM GRID TO INSIGHTII FORMAT')
c
c      write(6,110) 
c110   format(1x,'Is map in ascii (A) or binary (B) format?'
c     &            ,/,'(Answer: A or B)')
c      read(5,120) string
c120   format(a1)
c      if (string.eq.'A'.or.string.eq.'a') then
c        iform = 1
c      else if (string.eq.'B'.or.string.eq.'b') then
c        iform = 0
c      endif
c      write(6,130) 
c130   format(1x,'Enter grid map file name:')
c      read(5,135) flnm
c135   format(a)
      write(6,140) flnm
140   format(1x,'Grid input file :',a)
      iunit=11
      if ( iform .eq. 0 ) then
        open(unit=11, file = flnm, status = 'old', form='unformatted')
      else
        open(unit=11, file = flnm, status = 'old')
      endif
c
c      write(6,145) 
c145   format(1x,'Enter output INSIGHTII map file name:')
c      read(5,135) newfle
      write(6,148) newfle
148   format(1x,'Grid output file :',a)
c
       ounit=12
         open(12, file = newfle, status='new', form='unformatted')
c            call uhopen(12, newfle, 0)
c
c
	if ( iform .eq. 0 ) then
	  read (iunit) title, scale, dum2, grdflg, idum2,
     &              km, one, km, im, jm, km,
     &                   h ,      ox ,      oy ,      oz ,
     &          dum3, dum4, dum5, dum6, dum7, dum8, idum3, idum4
      elseif ( iform .eq. 1 ) then
        read (iunit,77) title, scale, dum2, grdflg, idum2,
     &                    km, one, km, im, jm, km,
     &                         h ,      ox ,      oy ,      oz ,
     &                    dum3, dum4, dum5, dum6, dum7,
     &                    dum8, idum3, idum4
       endif
c
c     do conversions (see grid2i)
      X0=ox+h
      Y0=oy+h
      Z0=oz+h
C**
      X1=ox+h*im
      Y1=oy+h*jm
      Z1=oz+h*km
C**
C**   FIND THE CORNER OF THE GRID WHICH IS FURTHEST FROM THE
C**   ORIGIN OF THE CARTESIAN COORDINATES.
C**
      XMAX = X0
      YMAX = Y0
      ZMAX = Z0
C**
      IF (ABS(X1).GT.ABS(X0)) XMAX=X1
      IF (ABS(Y1).GT.ABS(Y0)) YMAX=Y1
      IF (ABS(Z1).GT.ABS(Z0)) ZMAX=Z1
C**
      XSTART=X0/(ABS(XMAX))
      YSTART=Y0/(ABS(YMAX))
      ZSTART=Z0/(ABS(ZMAX))
C**
      XEND=X1/(ABS(XMAX))
      YEND=Y1/(ABS(YMAX))
      ZEND=Z1/(ABS(ZMAX))
C**
      write(ounit) title
      write(ounit) 0, 4, 0, 
     &            real(abs(xmax)), real(abs(ymax)), real(abs(zmax)),
     &            90.0, 90.0, 90.0,  
     &            xstart, xend, ystart, yend, zstart,zend,
     &            im -1, jm - 1, km - 1 
c
      if ( iform .eq. 0 ) then
        do  10  k = 1 , km
          read(iunit) kk, im, jm
          read(iunit) ( (grid(i,j,k), i=1,im ), j=1,jm )
   10   continue
      elseif ( iform .eq. 1 ) then
        do  20  k = 1 , km
          read(iunit,78) kk, im, jm
          read(iunit,79) ( ( grid(i,j,k), i=1,im), j=1,jm )
   20   continue
      endif
c
      gmx= grid(1,1,1)
      gmn= grid(1,1,1)
      do 51 i=1,im
      do 52 j=1,jm
      do 53 k=1,km
      if(grid(i,j,k).gt.gmx) gmx=grid(i,j,k)
      if(grid(i,j,k).lt.gmn) gmn=grid(i,j,k)
 53   continue
 52   continue
 51   continue
c
      write(6,502) im,jm,km,h,ox,oy,oz
 502  format(1x,'Read map, dimensions:',3i4,' spacing',f8.3,/,
     .       1x,'origin:',3f8.3)
      write(6,*) 'min-max:',gmn,gmx
c
      do  30  k = 1 , km
        do  30  j = 1 , jm
          write(ounit) (grid(i,j,k), i=1, im )
   30 continue
c
      close(iunit)
      close(ounit)
c
c
c               Format statements.
c
   77 format(a72/2e12.6,5i7/3i7,4e12.6/4e12.6/2e12.6,2i7)
   78 format(3i7)
   79 format( (6(1x,e12.6)) )
 999  continue
      end
