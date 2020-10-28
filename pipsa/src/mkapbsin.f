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
c This is for apbs version 0.3.2 and newer, taking into account 
c corrected uhbd grid origin writing
c
c Rewrites apbs input file changing center of grids 
c and ionic strength value
c
c Input:
c - command line arg -ft - location of template apbs input file
c   default "apbs.in_tmpl", this means that the template file is 
c   copied from pipsa directory to the current . 
c - command line arg -is - ionic strength in mM
c   default 50 mM
c - command line arg -fu - input file to be used for UHBD calculations
c   default "apbs.in"
c - standard input - output of ccenter.f - a table of centers of 
c   proteins in a set
c
c Output:
c   Information about the spread of centers and extents of proteins

      include "maxdim.inc"
      character*256 arg,s,fapbst,fapbsin,se
      character*1 key
      dimension x(nprmx),y(nprmx),z(nprmx),d(nprmx)
      data key/'#'/
c
      do ii=1,256
      se(ii:ii)=' '
      end do
c
      fapbst='apbs.in_tmpl'
      fapbsin='apbs.in'
      aios=50.0
c
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: mkapbsin -ft apbs.in_tmpl -fu apbs.in ',
     ,'-is ionic_str_in_mM < center.log')
      call exit
      end if
c
      nargs=iargc()
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-ft') then
      call getarg(ia+1,fapbst)
      end if
      if(arg(1:3).eq.'-fu') then
      call getarg(ia+1,fapbsin)
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
      open(81,file=fapbst)
c. get spacing from dime and glen to shift center 
 810  continue
      read(81,'(a)',end=820) s 
      if(s(1:4).eq.'dime') then
      read(s(5:256),*) napbsx,napbsy,napbsz
      end if
      if(s(1:4).eq.'glen') then
      read(s(5:256),*) eapbsx,eapbsy,eapbsz
      end if
      goto 810
 820  continue
c      spx=eapbsx/float(napbsx-1)
c      xa=xa-spx
c      spy=eapbsy/float(napbsy-1)
c      ya=ya-spy
c      spz=eapbsz/float(napbsz-1)
c      za=za-spz
      rewind(81)
c
      open(91,file=fapbsin)
 81   continue
      read(81,'(a)',end=82) s
      call slen(s,ns)
c
      if(s(1:3).eq.'ion') then
      read(s(4:ns),*) iq,aiosi,radi
      s=se
      write(s,'(3Hion,i3,2f6.3)') iq,aios/1000.,radi
      call slen(s,ns)
      end if 
      if(s(1:5).eq.'gcent') then
      s=se
      write(s,'(5Hgcent,3f8.3)') xa,ya,za
      call slen(s,ns)
      end if
c
      write(91,'(a)') s(1:ns)
      goto 81
 82   continue
      close(81)
      close(91)
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
