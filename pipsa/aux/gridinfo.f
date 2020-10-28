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
      if(nargum.lt.2) then
      write(6,501)
      goto 999
      end if
 501  format(1x,'gridinfo - print min-max of a grid',/,
     ,       1x,'usage: grid2i A_or_B GRID_file ',/,
     ,       1x,'where: A_or_B is A for ascii and B for binary grid')
      call getarg(1,argum)      
      if(argum(1:1).eq.'A'.or.argum(1:1).eq.'a') then
        iform = 1
      else
        iform = 0
      endif
      call getarg(2,argum)
      flnm=argum
c      call getarg(3,argum)
c      newfle=argum
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
c      write(6,148) newfle
148   format(1x,'Grid output file :',a)
c
c       ounit=12
c         open(12, file = newfle, status='new', form='unformatted')
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
c      X0=ox+h
c      Y0=oy+h
c      Z0=oz+h
C**
c      X1=ox+h*im
c      Y1=oy+h*jm
c      Z1=oz+h*km
C**
C**   FIND THE CORNER OF THE GRID WHICH IS FURTHEST FROM THE
C**   ORIGIN OF THE CARTESIAN COORDINATES.
C**
c      XMAX = X0
c      YMAX = Y0
c      ZMAX = Z0
C**
c      IF (ABS(X1).GT.ABS(X0)) XMAX=X1
c      IF (ABS(Y1).GT.ABS(Y0)) YMAX=Y1
c      IF (ABS(Z1).GT.ABS(Z0)) ZMAX=Z1
C**
c      XSTART=X0/(ABS(XMAX))
c      YSTART=Y0/(ABS(YMAX))
c      ZSTART=Z0/(ABS(ZMAX))
C**
c      XEND=X1/(ABS(XMAX))
c      YEND=Y1/(ABS(YMAX))
c      ZEND=Z1/(ABS(ZMAX))
C**
c      write(ounit) title
c      write(ounit) 0, 4, 0, 
c     &            real(abs(xmax)), real(abs(ymax)), real(abs(zmax)),
c     &            90.0, 90.0, 90.0,  
c     &            xstart, xend, ystart, yend, zstart,zend,
c     &            im -1, jm - 1, km - 1 
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
c      do  30  k = 1 , km
c        do  30  j = 1 , jm
c          write(ounit) (grid(i,j,k), i=1, im )
c   30 continue
cc
      close(iunit)
c      close(ounit)
c
c
c               Format statements.
c
   77 format(a72/2e12.6,5i7/3i7,4e12.6/4e12.6/2e12.6,2i7)
   78 format(3i7)
   79 format( (6(1x,e12.6)) )
 999  continue
      end
