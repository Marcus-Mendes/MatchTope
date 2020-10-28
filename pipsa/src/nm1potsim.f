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
c Simple operation: removes a given protein from pipsa analysis.
c This is done by reading the protein's name, removing it from the 
c list (file "names") and from related to it similarity matrix ("sims.log")
c
c Input:
c - command line arg -fn - the name of the file with the names of 
c   original set of proteins, each of which should have corresponding 
c   PDB file in the directory ../pdbs/, and corresponding potential file in ./
c   default is "names"
c - command line arg -p1 - the name of the protein to be removed,
c   its potential in ./ and pdb file in ../pdbs are not removed
c   and not used.
c - command line arg -p  - the same as above 
c - command line arg -lg - the name of the similarity matrix file
c   default is "sims.log" - will be modified
c
c Note that 
c + the program implies that it is executed in the directory, where 
c   grid files are located and expects pdb files to be located in ../pdbs/ subdirectory
c + the file "names" after execution of this program will be renamed to "names-old" 
c   and a new file "names" will be created, which does not have the removed protein name
c + the old similarity matrix file will be renamed to "sims.log-old" and new file 
c    "sims.log" (or any other name given after -lg) sill be created which does not have 
c    all entries related to the removed protein anymore
c + note that the grid and pdb files will not be removed.  These need to
c    be either removed separately, or replaced by a new versions, if
c    subsequently a new version of a protein is supposed to be added
c
      include "maxdim.inc"
      character*256 pname,fname,lname,fparts,fpdbs,arg
      character*256 pnam1 pnam2
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
 499  format('usage: m1potsim -fn fname -lg lname ',
     ,'-p1 protein_name')
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
      write(6,*) 'the name of the protein to remove from the set',
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
 7021 format('mv -f ',a,'-old ',a)
c
      open(10,file=fnameo)
      open(11,file=fname)
      call slen(nap1,lnap1)
      nprem=0
 111  continue
      read(10,'(a)',end=94) arg
      call slen(arg,larg)
      if(nap1(1:lnap1).eq.arg(1:larg)) then
      nprem=nprem+1
      goto 111
      end if
      write(11,'(a)') arg(1:larg)
      goto 111
  94  continue
      close(10)
      close(11)
c
      if(nprem.eq.0) then
      command=dum
      write(command,7021)  fname(1:lfn), fname(1:lfn)
      call slen512(command,lc)
      write(6,'(a)') command(1:lc)
      call system(command)
      command=dum
      write(command,7021)  lname(1:lln), lname(1:lln)
      call slen512(command,lc)
      write(6,'(a)') command(1:lc)
      call system(command)
      write(6,*) 'a given protein is not found in ', 
     ,'the original set, therefore nothing changed'
      write(6,'(a)') nap1(1:lnap1) 
      goto 999
      end if
c
      open(13,file=lname)
      write(13,'(1H#,i6)') np-1
      open(12,file=lnameo)
      read(12,'(a)',end=95) arg
 112  continue
      read(12,'(a)',end=95) arg
      call get2str(arg,n1b,n1e,n2b,n2e)
c remove if the name is given
      if(arg(n1b:n1e).eq.nap1(1:lnap1)) then
      read(12,'(a)',end=95) arg
      goto 112
      end if
      if(arg(n2b:n2e).eq.nap1(1:lnap1)) then
      read(12,'(a)',end=95) arg
      goto 112
      end if
c leave as is if not given
      call slen(arg,la)
      write(13,'(a)') arg(1:la)
      read(12,'(a)',end=95) arg
      call slen(arg,la)
      write(13,'(a)') arg(1:la)
      goto 112
 95   continue
      close(12)
      close(13)
c
 999  continue
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

      subroutine get2str(arg,n1b,n1e,n2b,n2e)
      character*256 arg
      do i=1,256
      n1b=i
      if(arg(i:i).ne.' ') goto 11
      end do
 11   continue
      do i=n1b,256
      n1e=i-1
      if(arg(i:i).eq.' ') goto 12
      end do
 12   continue
      do i=n1e+1,256
      n2b=i
      if(arg(i:i).ne.' ') goto 21
      end do
 21   continue
      do i=n2b,256
      n2e=i-1
      if(arg(i:i).eq.' ') goto 22
      end do
 22   continue
      return
      end

