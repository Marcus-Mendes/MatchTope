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
c program reads max namx atoms
c it will not crash, but write a note if protein has more than nmax atoms
c
      include "maxdim.inc"
	  character*54 string
	  character*4 atom
	  data atom/'ATOM'/
      dimension xx(namx),yy(namx),zz(namx)
c
c      open(55,file='4ccenter.pdb')
      natom=0
	  xc=0.0
	  yc=0.0
	  zc=0.0
 1    continue
	  read(5,500,end=99) string
	  if(string(1:4).ne.atom) goto 1 
      read(string(31:54),501) x,y,z
      natom=natom+1
      if(natom.gt.namx) goto 999
      xc=xc+x
      yc=yc+y
      zc=zc+z
      xx(natom)=x 
      yy(natom)=y 
      zz(natom)=z 
       goto 1
 99   continue
c      close(55)
      xc=xc/float(natom)
      yc=yc/float(natom)
      zc=zc/float(natom)
 500  format(a)
 501  format(3f8.3)
 601  format('#',3f8.3,8x,f8.3,2i8)
c
	  d2max=0.0
	  iamax=0
      do 98 ia=1,natom
      x=xx(ia)
      y=yy(ia)
      z=zz(ia)
	  xr=x-xc
	  yr=y-yc
	  zr=z-zc
	  d2=xr**2+yr**2+zr**2
	  if(d2.gt.d2max) then
	  d2max=d2
	  iamax=ia
	  end if
 98   continue
      write(6,601) xc,yc,zc,sqrt(d2max),iamax,natom
c      open(66,file='ccenter.dat')
c      write(66,601) xc,yc,zc,sqrt(d2max),iamax,natom
c      close(66)
c
      call exit
 999  continue
      write(6,*) 'ERROR: no of atoms exceeds', namx
      call exit
      end
