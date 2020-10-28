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
c Computes similarity of 2 potential grids on all points
c
c Input:
c - command line arg -g1 - the name of the file with potential grid for protein 1 
c   in UHBD format, default grd1.grd
c - command line arg -g2 - the name of the file with potential grid for protein 2 
c   in UHBD format, default grd2.grd
c
c The program will use the values of electrostatic potentials at 
c each point of the grids and derive the similarity index.
c The followings will be computed:
c aa   = square of the norm of the grid 1
c bb   = square of the norm of the grid 2
c ab   = scalar product of 2 potentials
c
c  Output:
c - fort.66 - some info about constructed skins
c - standard output (fort.6) has following data in one line:
c
c   si_hodgkin = 2*ab/(aa+bb)
c   si_carbo   = ab/sqrt(aa*bb)
c   aa
c   bb
c   ab
c
CccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
      integer iform
      character*256 fgrd1,fgrd2,arg
      integer im, jm, km
      dimension grd1(im3_max)
      dimension grd2(im3_max)
c
      anull=0.0
c
      fgrd1='grd1.grd'
      fgrd2='grd2.grd'
c
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: 2potsim_noskin -g1 grd1.grd -g2 grd2.grd')
      call exit
      end if
      nargs=iargc()
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-g1') then
      call getarg(ia+1,fgrd1)
      end if
      if(arg(1:3).eq.'-g2') then
      call getarg(ia+1,fgrd2)
      end if
      end do
c
      ilog=13
      iform=0
c. Above are the hardcoded conventions
      iugrd1=11
      iugrd2=12
c
c begin code
c
c read the 3D grid
c
	  write(ilog,*)
	  write(ilog,*) 'Grid # 1 file reading now... wait...'
      call inpgrd(fgrd1,iugrd2,iform,im,jm,km,
     &             h,xo,yo,zo,grd1)
	  write(ilog,611) im,jm,km,h,xo,yo,zo
	  write(ilog,*)
	  write(ilog,*) 'Grid # 2 file reading now... wait...'
      call inpgrd(fgrd2,iugrd2,iform,im,jm,km,
     &             h,xo,yo,zo,grd2)
	  write(ilog,611) im,jm,km,h,xo,yo,zo
c
 611  format('grid dimensions:  ',i3,2('.',i3),'  spacing:',f8.3,/,
     ,'grid origin:',3f8.3)
      imjmkm=im*jm*km
	  call over(grd1,grd1,an11,imjmkm)
	  call over(grd2,grd2,an22,imjmkm)
	  call over(grd1,grd2,an12,imjmkm)

	  write(6,612) 2.*an12/(an11+an22),an12/sqrt(an11*an22),
     ,an11,an22,an12,anull,anull
c 612  format(2f8.4,3e12.4)
 612  format(2f7.3,5e15.7,2f7.3,3i7)

 1000  continue
      stop
c
      end
      subroutine over(grida,gridb,ov,npoi)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dimension grida(npoi)
      dimension gridb(npoi)
c
      integer i,npoi
	  real*8 ovd
c
c begin code
c
c calculate the overlapping      
c
	  ovd=0.0d0
      do 130 i=1,npoi     
	  ovd=ovd+grida(i)*gridb(i)
  130 continue
      ov=ovd
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
