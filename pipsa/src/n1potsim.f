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
c drives similarity index calculations with 2potsim* when 1 extra protein 
c to be added to the set of originally processed proteins
c
c Input:
c - command line arg -pg - similarity calculation program name
c   2potsim_noskin, 2potsim_skin or 2potsim_skin_parts, 
c   default is ../bin/2potsim_noskin
c - command line arg -fp - the directory where pdb files are located, 
c   default is ../pdbs
c - command line arg -fn - the name of the file with the names of 
c   original set of proteins, each of which should have corresponding 
c   PDB file in the directory ../pdbs/, and corresponding potential file in ./
c   default is "names"
c - command line arg -p1 - the name of the protein to be added, 
c   should have PDB file in ../pdbs/, and the potential file in ./
c - command line arg -lg - the name of the similarity matrix file
c   default is "sims.log"
c - command line arg -pr - the value of probe radius
c   default is 3 A
c - command line arg -sk - the value of skin thickness
c   default is 4 A
c - command line arg -pa - the name of the file with the list of directions
c   default name "parts" - a list of directions and angles to define conical parts
c
c Note that 
c + the program implies that it is executed in the directory, where 
c   grid files are located and expects pdb files to be located in ../pdbs/ subdirectory
c + all grid files corresponding to protein names in the file "names" must have 
c   extension .grd and all pdb files must have extension .pdb
c + the file "names" after execution of this program will be renamed to "names-old" 
c   and a new file "names" will be created, which includes newly added protein name
c + the old similarity matrix file will be renamed to "sims.log-old" and new file 
c    "sims.log" (or any other name given after -lg) sill be created which has indices for added protein
c 
c
      include "maxdim.inc"
      character*256 pname,fname,lname,fparts,fpdbs,arg
      character*512 command,dum
      character*256 fnameo,lnameo,nap1
      character*256 ns(nprmx)
      character*1 s
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
      nap1='      '
      nargs=iargc()
c
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: n1potsim -pr program -fn fname -lg lname ',
     ,'-pr probe -sk skin -pa parts -p1 protein_name')
      call exit
      end if
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
      if(arg(1:3).eq.'-p1') then
      call getarg(ia+1,nap1)
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
      open(12,file=lname)
      read(12,'(1x,i6)') npo
      close(12)
      if(npo.ne.np) then
      write(6,*) 'file with names says there were',np,'  proteins'
      write(6,*) 'file with sim matrix says there were',npo,'  proteins'
      write(6,*) 'please resolve the inconsistency first'
      call exit
      end if
      if(nap1.eq.'      ') then
      write(6,*) 'the name of the protein to add to the existing set',
     ,' must be given'
      write(6,499)
      call exit
      end if
c
      call slen(fname,lfn)
      write(fnameo,701) fname(1:lfn)
      command=dum
      write(command,702)  fname(1:lfn), fname(1:lfn)
      call slen512(command,lc)
      write(6,'(a)') command(1:lc)
      call system(command)
      call slen(lname,lln)
      write(lnameo,701) lname(1:lln)
      write(command,702)  lname(1:lln), lname(1:lln)
      call slen512(command,lc)
      write(6,'(a)') command(1:lc)
      call system(command)
 701  format(a,'-old')
 702  format('mv -f ',a,' ',a,'-old')
c
      open(13,file=lname)
      write(13,'(1H#,i6)') np+1
c
      open(12,file=lnameo)
      read(12,'(a)',end=95) arg
 112  continue
      read(12,'(a)',end=95) arg
      call slen(arg,la)
      write(13,'(a)') arg(1:la)
      goto 112
 95   continue
      close(12)
      close(13)
c
      open(11,file=fname)
      do j=1,np
c
      call slen(nap1,lp1)
      call slen(ns(j),lj)
      open(13,file=lname)
 131  continue
      read(13,'(a)',end=132) s
      goto 131
 132  continue
      write(13,'(a,1x,a)') nap1(1:lp1),ns(j)(1:lj)
      close(13)
c
      command=dum
      call slen(pname,lpg)
      call slen(fpdbs,lpdbs)
      call slen(lname,llg)
      call slen(fparts,lpa)
      write(command, 500) pname(1:lpg),nap1(1:lp1),ns(j)(1:lj),
     ,fpdbs(1:lpdbs),nap1(1:lp1),fpdbs(1:lpdbs),ns(j)(1:lj),
     ,fparts(1:lpa),probes,skin,lname(1:llg)
      call system(command)
      call slen512(command,lc)
      write(6,'(a)') command(1:lc)
c
      write(11,'(a)') ns(j)(1:lj)
      end do
      open(13,file=lname)
 133  continue
      read(13,'(a)',end=134) s
      goto 133
 134  continue
      write(13,'(a,1x,a)') nap1(1:lp1),nap1(1:lp1)
      close(13)
      write(13,'(a,1x,a)') nap1(1:lp1),nap1(1:lp1)
      write(command, 500) pname(1:lpg),nap1(1:lp1),nap1(1:lp1),
     ,fpdbs(1:lpdbs),nap1(1:lp1),fpdbs(1:lpdbs),nap1(1:lp1),
     ,fparts(1:lpa),probes,skin,lname(1:llg)
      call system(command)
      call slen512(command,lc)
      write(6,'(a)') command(1:lc)
c
      write(11,'(a)') nap1(1:lp1)
      close(11)
      close(13)
c
 500  format(a,' -g1 ',a,'.grd -g2 ',a,'.grd -p1 ',a,'/',a,
     ,'.pdb -p2 ',a,'/',a,'.pdb -pa ',a,' -pr ',f7.3,
     ,' -sk ',f7.3,' >> ',a)

c
      call exit
      end  
c
      subroutine slen(arg,na)
      character*256 arg
      na=1
      do ia=2,249
      na=ia-1
      if(arg(ia:ia+7).eq.'        ') return
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

