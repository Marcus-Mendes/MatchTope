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
c By Dave Sept, modified from UHBD utility
c Modified for pipsa 2.0 April 2004
c 
c Converts grid from UHBD ASCII format to UHBD binary format
c
c     iunit        i     unit number to use for writing.
c     iform        i     0=unformatted file, 1=formatted file.
c     title        i     title for file.
c     im           i     numbers of points in x - direction.
c     jm           i     numbers of points in y - direction.
c     km           i     numbers of points in z - direction.
c     h            i     distance between grid points.
c     ox \             / coordinate of grid point (0,0,0)
c     oy  >        i  <  such that the coordinates of (i,j,k) are
c     oz /             \ ( ox+h*i, oy+h*j, oz+h*k).
c     grdflg       i     flag indicating the type of grid.
c     grid         i     grid of values to be printed out.
c     scale        i     factor multiplied times the grid that is output.
c
      include "maxdim.inc"
      integer im, jm, km, grdflg, one, kk
      integer iunit, iform
      character*72 title
      real h, ox, oy, oz, grid(im_max,im_max,im_max), scale
      real x0, y0, z0, x1, y1, z1, xmax, ymax, zmax, 
     &     xstart, ystart, zstart, xend, yend, zend
  
c     dum#         dummy real variables
c     idum#        dummy integer variables
c
      real dum2, dum3, dum4, dum5, dum6, dum7, dum8
      integer idum2, idum3, idum4
c
      integer i, j, k, ounit
c
c     flnm          name of trajectory file.
c     newfle        output DCD file name (stemunit.sub.DCD)
c     stem          output DCD file name stem
c     string
c
c      character*20 newfle
c      character*80 flnm
      character*4 stem
      character*1 string
      character*256 flnm,newfle
c
c Statement functions
c begin code
c
      if(iargc().lt.2) then
      write(6,601)
 601  format('Usage: uhbd_asc2bin asci_grid_name bin_grid_name')
      call exit
      end if
c
      call getarg(1,flnm)
      call getarg(2,newfle)
c
c      write(6,100) 
c100   format(1x,'Convert UHBD ascii grid to binary')
c
c      write(6,130) 
c130   format(1x,'Enter grid map file name:')
c      read(5,135) flnm
c135   format(a)
c      write(6,140) flnm
c140   format(1x,'Grid input file :',a)
      iunit=11
      open(unit=11, file = flnm, status = 'old')
c
c      write(6,145) 
c145   format(1x,'Enter output file name:')
c      read(5,135) newfle
c      write(6,148) newfle
c148   format(1x,'Grid output file :',a)
c
       ounit=12
         open(12, file = newfle, status='new', form='unformatted')
c            call uhopen(12, newfle, 0)
c
c
        read (iunit,77) title, scale, dum2, grdflg, idum2,
     &                    km, one, km, im, jm, km,
     &                         h ,      ox ,      oy ,      oz ,
     &                    dum3, dum4, dum5, dum6, dum7,
     &                    dum8, idum3, idum4

        write (ounit) title, scale, dum2, grdflg, idum2,
     &                    km, one, km, im, jm, km,
     &                         h ,      ox ,      oy ,      oz ,
     &                    dum3, dum4, dum5, dum6, dum7,
     &                    dum8, idum3, idum4
c
        do  20  k = 1 , km
          read(iunit,78) kk, im, jm
          read(iunit,79) ( ( grid(i,j,k), i=1,im), j=1,jm )
   20   continue
c
      do  30  k = 1 , km
          write(ounit) k, im, jm
          write(ounit) ( ( grid(i,j,k), i=1,im), j=1,jm )
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
      end
