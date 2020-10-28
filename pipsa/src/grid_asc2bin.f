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
c Converts grid from GRID ASCII format to UHBD binary format
c GRID binary GRIDKONT can be converted to this ASCII format by k2a program
c
      include "maxdim.inc"
      integer im, jm, km, grdflg, mxgrd, one, kk
      integer iunit, iform
      character*72 title
      real h, ox, oy, oz, grid(im_max,im_max,im_max), scale
      real x0, y0, z0, x1, y1, z1, xmax, ymax, zmax, 
     &     xstart, ystart, zstart, xend, yend, zend
      real dum2, dum3, dum4, dum5, dum6, dum7, dum8
      integer idum2, idum3, idum4
      integer i, j, k, ounit
      character*4 stem
      character*1 string
      character*256 flnm,newfle,s
c
c begin code
c
      if(iargc().lt.2) then
      write(6,601)
 601  format('Usage: grid_asc2bin asci_grid_name bin_grid_name')
      call exit
      end if
c
      call getarg(1,flnm)
      call getarg(2,newfle)
      open(11, file = flnm, status = 'old')
      open(12, file = newfle, status='new', form='unformatted')
c
c        read (11,77) title, scale, dum2, grdflg, idum2,
c     &                    km, one, km, im, jm, km,
c     &                         h ,      ox ,      oy ,      oz ,
c     &                    dum3, dum4, dum5, dum6, dum7,
c     &                    dum8, idum3, idum4
      scale=1.0
      dum2=0.0
      grdflg=-1
      idum2=0
      one=1
      dum3=0.0
      dum4=0.0
      dum5=0.0
      dum6=0.0
      dum7=0.0
      dum8=0.0
      idum3=0
      idum4=0
 111  continue
      read(11,'(a)',end=99) s
      if(s(1:8).ne.' HEADER ') goto 111
      read(11,'(a)',end=99) s
      read(11,'(a)',end=99) title
c
 112  continue
      read(11,'(a)',end=99) s
      if(s(1:10).ne.' GRID size') goto 112
      read(11,'(a)',end=99) s
      read(11,'(a)',end=99) s
c  NNX:  65             NNY:  65             NNZ:  65
      do i=1,256
      if(s(i:i+3).eq.'NNX:') read(s(i+4:i+8),*) im
      if(s(i:i+3).eq.'NNY:') read(s(i+4:i+8),*) jm
      if(s(i:i+3).eq.'NNZ:') read(s(i+4:i+8),*) km
      end do
 113  continue
      read(11,'(a)',end=99) s
      if(s(1:10).ne.' GRID spac') goto 113
      read(11,'(a)',end=99) s
      read(11,'(a)',end=99) s
c   RA:   1.000          RX: -11.651          RY: -19.853          RZ: -21.540
            do i=1,256
      if(s(i:i+2).eq.'RA:') read(s(i+3:i+11),*) h 
      if(s(i:i+2).eq.'RX:') read(s(i+3:i+11),*) ox
      if(s(i:i+2).eq.'RY:') read(s(i+3:i+11),*) oy
      if(s(i:i+2).eq.'RZ:') read(s(i+3:i+11),*) oz
      end do
 114  continue
      read(11,'(a)',end=99) s
      if(s(1:10).ne.' Results f') goto 114
c
        write (12) title, scale, dum2, grdflg, idum2,
     &                    km, one, km, im, jm, km,
     &                         h ,      ox ,      oy ,      oz ,
     &                    dum3, dum4, dum5, dum6, dum7,
     &                    dum8, idum3, idum4
c
c        do  20  k = 1 , km
c          read(11,78) kk, im, jm
c          read(11,79) ( ( grid(i,j,k), i=1,im), j=1,jm )
c   20   continue
c
      do k=1,km
c      write(6,*) k
      do i=1,im
      read(11,*)
      read(11,501) (grid(i,j,k),j=1,jm)
      end do
 115  continue
      read(11,'(a)',end=98) s
      if(s(1:10).ne.' Results f') goto 115
      end do
c
c   1    1.223  1.236  1.250  1.264  1.278  1.293  1.307  1.321  1.335  1.349
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
 501  format(6x,10f7.3)
c
 98   continue
      do  30  k = 1 , km
          write(12) kk, im, jm
          write(12) ( ( grid(i,j,k), i=1,im), j=1,jm )
   30 continue
c
      close(11)   
      close(12)   
 99   continue
c
c
c               Format statements.
c
   77 format(a72/2e12.6,5i7/3i7,4e12.6/4e12.6/2e12.6,2i7)
   78 format(3i7)
   79 format( (6(1x,e12.6)) )
      end
