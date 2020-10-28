c----------------------------------------------------------------------
c  program npotsim     
c  PIPSA version 2.0 beta (December 2004)
c
c  Authors: R.R.Gabdoulline, R.C.Wade
c  Copyright (c) 2004
c  EML Research gGmbH, Schloss-Wolfsbrunnenweg 33, 69118 Heidelberg,
c  Germany, and
c  European Molecular Biology Laboratory, Heidelberg, Germany
c
c  This copyright notice must be attached to all copies, or extracts,
c  of this software.
c
c  Please obtain a license before using this software from the PIPSA
c  website or "rebecca.wade@eml-r.villa-bosch.de".
c  Conditions of usage are as specified in the license agreement.
c
c  Queries about problems with this software should be addressed to
c  "razif.gabdoulline@eml-r.villa-bosch.de".
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
c
c drives similarity index calculations with 2potsim*
c
c Input:
c - command line arg -pg - similarity calculation program name
c   2potsim_noskin, 2potsim_skin or 2potsim_skin_parts, 
c   default is ../bin/2potsim_noskin
c - command line arg -fp - the directory where pdb files are located, 
c   default is ../pdbs
c - command line arg -fn - the name of the file with the names of 
c   proteins, each of which should have corresponding PDB file in the 
c   directory from -fp (../pdbs/), and corresponding potential file in ./
c   default is "names"
c - command line arg -lg - the name of the similarity matrix file
c   default is "sims.log"
c - command line arg -pr - the value of probe radius
c   default is 3 A
c - command line arg -sk - the value of skin thickness
c   default is 4 A
c - command line arg -pa - the name of the file with the list of directions
c   default name "parts" - a list of directions and angles to define conical parts
c
c Note that the program implies that it is executed in the directory, where 
c grid files are located and expects pdb files to be located in ../pdbs/ subdirectory
c all grid files corresponding to protein names in the file "names" must have 
c extension .grd and all pdb files must have extension .pdb
c
      include "maxdim.inc"
      character*256 pname,fname,lname,fparts,fpdbs,arg
      character*512 command,dum
      character*1 s
      character*256 ns(nprmx)
      character*512 target
      character*20 sgrpcount
      integer nmkcmds=10000
      integer icmdcount=int(0)
      integer igrpcnt=int(0)
      integer ismake=int(0)
C      character*512 target(maxgrps)
c
      do i=1,512
      dum(i:i)=' '
      end do
c
      pname='../bin/2potsim_noskin'
      fpdbs='../pdbs'
      fname='names'
      lname='sims.log'
      probes=3.0
      skin=4.0
      fparts='parts'
      nargs=iargc()
c
      call getarg(1,arg) 
      if((arg(1:2).eq.'').or.(arg(1:2).eq.'-h')) then
      write(6,499)
 499  format('usage: npotsim -pr program -fn fname -lg lname ',
     ,'-pr probe -sk skin -pa parts [-ma size_of_commandgroups ',
     ,'(to generate makefile)]')
      call exit
      end if
      ismake=0
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-pg') then
      call getarg(ia+1,pname)
      end if
      if(arg(1:3).eq.'-fp') then
      call getarg(ia+1,fpdbs)
      end if
      if(arg(1:3).eq.'-fn') then
      call getarg(ia+1,fname)
      end if
      if(arg(1:3).eq.'-lg') then
      call getarg(ia+1,lname)
      end if
      if(arg(1:3).eq.'-pr') then
      call getarg(ia+1,arg)
      read(arg,*) probes
      end if
      if(arg(1:3).eq.'-sk') then
      call getarg(ia+1,arg)
      read(arg,*) skin
      end if
      if(arg(1:3).eq.'-pa') then
      call getarg(ia+1,fparts)
      end if
      if(arg(1:3).eq.'-ma') then
      ismake=int(1)
      call getarg(ia+1,arg)
      read(arg,*) nmkcmds
      end if
      end do
c
      open(11,file=fname)
      np=0
 1    continue
      read(11,'(a)',end=91) arg
      if(arg(1:6).eq.'      ') goto 1
      np=np+1
      ns(np)=arg
      goto 1
 91   continue
      close(11)
c
      if(ismake.eq.1) then
      open(10,file="Makefile")
      open(9,file="Makefile.alltargets")
      agrpcnt=0
      write(10,506) "all: collect"
      write(10,506) "clean: "
      write(10,503) char(9),"rm *.log *.part collect"
      write(9,502) "collect:"
      else
      open(12,file=lname)
      write(12,'(1H#,i6)') np
      close(12)
      end if
      do i=1,np
      do j=1,i
c
      call slen(ns(i),li)
      call slen(ns(j),lj)
      if(ismake.eq.int(0)) then
      open(12,file=lname)
 121  continue
      read(12,'(a)',end=122) s
      goto 121
 122  continue
      write(12,'(a,1x,a)') ns(i)(1:li),ns(j)(1:lj)
      close(12)
      end if
c
      command=dum
      call slen(pname,lpg)
      call slen(lname,llg)
      call slen(fparts,lpa)
      call slen(fpdbs,lpdbs)
      if(ismake.eq.0)then
      write(command, 500) pname(1:lpg),ns(i)(1:li),ns(j)(1:lj),
     ,fpdbs(1:lpdbs),ns(i)(1:li),fpdbs(1:lpdbs),ns(j)(1:lj),
     ,fparts(1:lpa),probes,skin,lname(1:llg)
      call system(command)
      else
      if (mod(icmdcount,nmkcmds).eq.0) then
      write(sgrpcount,505) igrpcnt+1
      call llen20(sgrpcount,lgrp)
      write(target,501) sgrpcount(lgrp:20),".",lname(1:llg)
      call slen512(target,lta)
      write(10,503) target(1:lta),":"
      write(9,502) target(1:lta)
      igrpcnt = igrpcnt+1
      endif
      write(command, 510) pname(1:lpg),ns(i)(1:li),ns(j)(1:lj),
     ,fpdbs(1:lpdbs),ns(i)(1:li),fpdbs(1:lpdbs),ns(j)(1:lj),
     ,fparts(1:lpa),probes,skin,sgrpcount(lgrp:20),".",lname(1:llg),
     ,'.part'
      call slen512(command,lc)
      if (mod(icmdcount,nmkcmds).eq.0) then
      write(10,507) char(9),'@echo "',ns(i)(1:li),ns(j)(1:lj),'">',
     ,sgrpcount(lgrp:20),".",lname(1:llg),'.part'
      else
      write(10,507) char(9),'@echo "',ns(i)(1:li),ns(j)(1:lj),'">>',
     ,sgrpcount(lgrp:20),".",lname(1:llg),'.part'      
      end if
      write(10,503) char(9), command(1:lc)
      if ((mod(icmdcount+1,nmkcmds).eq.0).or.(j.eq.np)) then 
      write(10,504) char(9), "touch ",target
      end if
      icmdcount = icmdcount+1
      end if
      call slen512(command,lc)
      if(ismake.eq.0) then
      write(6,'(a)') command(1:lc)
      else
      incp=int(np*(np+1)/2)
      if ((mod(icmdcount,10000).eq.0).or.(icmdcount.eq.incp))then
      write(6,'(a,i9,a,i9,a)') 'wrote ',icmdcount,' of ',incp,
     ,' commands'
      end if
      end if
c
      end do
      end do
      if(ismake.eq.1) then
      write(10, 506) 'include Makefile.alltargets'
      close (10)
      write(9,*)
      write(9,508) char(9),'@echo "#',np,'" > ',lname(1:llg)
      write(9,501) char(9),
     ,'ls *.log.part|sort -n -t .|xargs cat >> ',lname(1:llg)
      write(9,503) char(9),"touch collect"
      close(9)
      end if
c
 500  format(a,' -g1 ',a,'.grd -g2 ',a,'.grd -p1 ',a,'/',a,
     ,'.pdb -p2 ',a,'/',a,'.pdb -pa ',a,' -pr ',f7.3,
     ,' -sk ',f7.3,' >> ',a)
 501  format(a,a,a)
 502  format(a,' ',$)
 503  format(a,a)
 504  format(a,a,a)
 505  format(i20)
 506  format(a)
 507  format(a,a,a,1x,a,a,a,a,a,a)
 508  format(a,a,i6,a,a)
 510  format(a,' -g1 ',a,'.grd -g2 ',a,'.grd -p1 ',a,'/',a,
     ,'.pdb -p2 ',a,'/',a,'.pdb -pa ',a,' -pr ',f7.3,
     ,' -sk ',f7.3,' >> ',a,a,a,a)
c
      call exit
      end  
c
      subroutine slen(arg,na)
      character*256 arg
      na=1
      do ia=2,253
      na=ia-1
      if(arg(ia:ia+3).eq.'    ') return
      end do 
      return
      end

      subroutine slen512(arg,na)
      character*512 arg
      do i=512,1,-1
      na=i
      if(arg(i:i).ne.' ') return
      end do 
      return
      end
      
      subroutine llen20(arg,na)
      character*20 arg
      do i=1,512
      na=i
      if(arg(i:i).ne.' ') return
      end do
      return
      end
