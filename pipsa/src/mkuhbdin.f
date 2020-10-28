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
c Rewrites uhbd input file changing center of grids 
c and ionic strength value
c
c Input:
c - command line arg -ft - location of template uhbd input file
c   default "uhbd.in_tmpl", this means that the template file is 
c   copied from pipsa directory to the current . 
c - command line arg -is - ionic strength in mM
c   default 50 mM
c - command line arg -fu - input file to be used for UHBD calculations
c   default "uhbd.in"
c - standard input - output of ccenter.f - a table of centers of 
c   proteins in a set
c
c Output:
c   Information about the spread of centers and extents of proteins

      include "maxdim.inc"
      character*256 arg,s,fuhbdt,fuhbdin
      character*1 key
      dimension x(nprmx),y(nprmx),z(nprmx),d(nprmx)
      data key/'#'/
c
      fuhbdt='uhbd.in_tmpl'
      fuhbdin='uhbd.in'
      aios=50.0
c
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: mkuhbdin -ft uhbd.in_tmpl -fu uhbd.in ',
     ,'-is ionic_str_in_mM < center.log')
      call exit
      end if
c
      nargs=iargc()
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-ft') then
      call getarg(ia+1,fuhbdt)
      end if
      if(arg(1:3).eq.'-fu') then
      call getarg(ia+1,fuhbdin)
      end if
      if(arg(1:3).eq.'-is') then
      call getarg(ia+1,arg)
      read(arg,*) aios
      end if
      end do
c
      xa=0.0
      ya=0.0
      za=0.0
      da=0.0
      ic=0
 1    continue
      read(5,'(a)',end=99) s
      if(s(1:1).ne.key) goto 1
      ic=ic+1 
      read(s(1:41),601) 
     r     x(ic),y(ic),z(ic),d(ic)
      xa=xa+x(ic)
      ya=ya+y(ic)
      za=za+z(ic)
      da=da+d(ic)
      goto 1
 99   continue
      nc=ic
      xa=xa/float(nc)
      ya=ya/float(nc)
      za=za/float(nc)
      da=da/float(nc)
c
      dx=0.0
      dd=0.0
      do ic=1,nc
      dx=dx+(xa-x(ic))**2+
     +(ya-y(ic))**2+(za-z(ic))**2
      dd=dd+(da-d(ic))**2
      end do
      dx=sqrt(dx/float(nc))
      dd=sqrt(dd/float(nc))
c
      open(81,file=fuhbdt)
      open(91,file=fuhbdin)
 81   continue
      read(81,'(a)',end=82) s
      call slen(s,ns)
      write(91,'(a)') s(1:ns)
      if(s(1:24).eq.'! The end of assignments') then
      write(91,602) xa,ya,za
      write(91,603) aios
      end if
      goto 81
 82   continue
      close(81)
      close(91)
c
 602  format('assign xc = ',f8.3,' end',/,
     ,       'assign yc = ',f8.3,' end',/,
     ,       'assign zc = ',f8.3,' end')
 603  format('assign ios = ',f8.3,' end')
c
      write(6,604) xa,ya,za,nc,dx,da
 604  format(/,'Ave center:',3f8.3,'  has been computed',/
     ,'rms of ',i3,' centers is',f8.3,'  compared to radii',f8.3)
c
 601  format(1x,3f8.3,8x,f8.3)
c
      call exit
      end
      subroutine slen(s,ns)
      character*256 s
      do i=256,1,-1
      ns=i
      if(s(i:i).ne.' ') return
      end do
      return
      end
