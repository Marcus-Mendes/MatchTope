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
c: Command File for Programme GRID
c: Channel Numbers and Output File names:
c LONT      6
c LONT  gridlont.dat                                                            
c KONT     20
c KONT  gridkont.dat                                                            
c: Channel Numbers and Input File names:
c INPT     10
c INPT  grinkout.dat                                                            
c: Control Parameters:
c CLER      5.000 
c DEEP      5.000 
c DPRO      4.000 
c DWAT     80.000 
c EACH      5.000 
c EMAX      5.000 
c FARH      5.000 
c FARR      8.000 
c KWIK      0     
c LEAU      0     
c LENG     12     
c LEVL      0     
c LIST      1     
c MOVE      0     
c NETA      1    
c NPLA      1     
c NUMB      1  
c VALU      0.000 
c DRY           
c TOPX     76.467 
c TOPY     66.882
c TOPZ     66.678
c BOTX    -33.533
c BOTY    -43.118
c BOTZ    -43.322
c IEND
c GRIG 1 A spacing, size - from top and botss
c       0    1
c:PROG:  /wade/wade3/owl2/local/gri17/grid
c:DEFA:                                                                         
c. reads the file center.log and writes grid.in
      character*41 string
      character*1 key
      character*256 fgridin,arg,probe
      data key/'#'/
c
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: mkgridin -fg grid.in -pr GRID_PROBE_NAME',
     ,' < center.log')
      call exit
      end if
c
      fgridin='grid.in'
      probe='PO4' 
c
      nargs=iargc()
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-fg') then
      call getarg(ia+1,fgridin)
      end if
      if(arg(1:3).eq.'-pr') then
      call getarg(ia+1,probe)
      end if
      end do
      call slen(probe,npr)
c
      xa=0.0
      ya=0.0
      za=0.0
      da=0.0
      dmax=0.0
      ic=0
 1    continue
      read(5,'(a)',end=99) string
      if(string(1:1).ne.key) goto 1
      ic=ic+1 
      read(string(1:41),601) 
     r     xic,yic,zic,dic
      xa=xa+xic
      ya=ya+yic
      za=za+zic
      da=da+dic
      if(dic.gt.dmax) dmax=dic
      goto 1
 99   continue
      nc=ic
      xa=xa/float(nc)
      ya=ya/float(nc)
      za=za/float(nc)
      da=da/float(nc)
c
c. grid size is 65^3, spacing 1 A, no checks in this version
c
      imax=65
      half=float(imax-1)/2.0
      tx=xa+half
      ty=ya+half
      tz=za+half
      bx=xa-half
      by=ya-half
      bz=za-half
      open(77,file=fgridin)
      write(77,602) probe(1:npr),tx,ty,tz,bx,by,bz
      close(77)
c
 601  format(1x,3f8.3,8x,f8.3)
 602  format(
     ,': Command File for Programme GRID',/,
     ,': Channel Numbers and Output File names:'/,
     ,' LONT       6',/,' LONT  gridlont.dat'/,
     ,' KONT      20',/,' KONT  gridkont.dat',/,
     ,': Channel Numbers and Input File names:',/,
     ,' INPT      10',/,' INPT  grinkout.dat',/,
     ,': Control Parameters:',/,
     ,' CLER     5.000',/,' DEEP     5.000',/,
     ,' DPRO     4.000',/,' DWAT    80.000',/,
     ,' EACH     5.000',/,' EMAX     5.000',/,
     ,' FARH     5.000',/,' FARR     8.000',/,
     ,' KWIK     0'    ,/,' LEAU     0',/,
     ,' LENG    12'    ,/,' LEVL     0',/,
     ,' LIST     1'    ,/,' MOVE     0'    ,/,
     ,' NETA     1'    ,/,' NPLA     1'    ,/,
     ,' NUMB     1'    ,/,' VALU     0.000',/,
     ,' ',a,/,
     ,' TOPX ',f8.3,/,' TOPY ',f8.3,/,' TOPZ ',f8.3,/,
     ,' BOTX ',f8.3,/,' BOTY ',f8.3,/,' BOTZ ',f8.3,/,
     ,' IEND',/,
     ,' GRID 1 A spacing, size - from top and bots',/,
     ,'       0    1',/,
     ,':PROG:  /wade/wade3/owl2/local/gri17/grid',/,
     ,':DEFA:'
     ,)
c
 999  continue
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
