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
CccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
c
c Computes pairwise electrostatic similarity of a list of proteins, 
c based on their monopole and dipole moments 
c
c Input: 
c - command line arg -fn - the name of the file with the names of proteins, 
c   each of which should have corresponding PDB file in the current directory
c   default is "pdbnames"
c - command line arg -fd - the name of the file where dipole moment 
c   information will be written
c   default name "dipoles" 
c - command line arg -r - the size of the proteins, approx average gyration radius
c   default value 9.815 A, valid for PH domains
c
c The program will assign formal charges to all charged residues 
c (+0.5 for NHX of Arg, +1 for NZ of Lys, -0.5 for OEX of Glu and ODX of Asp, 
c compute the monopole and dipole moments of protein.  
c The similarity index is computed then following the analytical formula 
c (6) from the Proteins paper, i.e. comparing monopole+dipole potentials 
c at the sphere of some radius R. Hardcoded value of R is 9.815, i.e. 
c parameter alpha (coded as scf) is 17 A**-2.   
c
c The following quantities are computed:
c aa   = square of the norm of the dipole+monopole potential of protein 1
c bb   = square of the norm of the dipole+monopole potential of protein 2
c ab   = scalar product of 2 dipole+monopole potentials
c
c  Output:
c - file "dipoles" where the following information about proteins are printed:
c   the name of pdb-file
c   total charge
c   norm of the dipole moment
c   3 (x,y,z) components of the dipole moment 
c   number of charge sites    
c - fort.6 - standard output has following data in one line:
c
c   si_hodgkin = 2*ab/(aa+bb)
c   si_carbo   = ab/sqrt(aa*bb)
c   aa
c   bb
c   ab
c   aa
c   bb
c   npoi
c
CccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
      include "maxdim.inc"
      dimension d(3,nprmx),c(3,nprmx),q(nprmx),cc(3),cv(3)
      dimension d1(3),d2(3)
      integer nq(nprmx),lfn(nprmx)
      character*256 fname(nprmx),fnamet,ffname,fdipoles,arg
c
      ffname='pdbnames'
      fdipoles='dipoles'
      rprot=9.815
      nargs=iargc()
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: qdipsim -fn pdbnames -fd dipoles -r rprot')
      call exit
      end if
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-fn') then
      call getarg(ia+1,ffname)
      end if
      if(arg(1:3).eq.'-fd') then
      call getarg(ia+1,fdipoles)
      end if
      if(arg(1:2).eq.'-r') then
      call getarg(ia+1,arg)
      read(arg,*) rprot
      end if
      end do
c
      open(12,file=ffname)
      np=0
 12   continue
      read(12,'(a)',end=92) fnamet
      if(fnamet(1:6).eq.'      ') goto 92
      call pdbfl(fnamet,lfnt)
      np=np+1
      call dipole1(fnamet,cv,q1,d1,nq1)
      nq(np)=nq1
      q(np)=q1
        do m=1,3
      d(m,np)=d1(m)
      c(m,np)=cv(m)
        end do
      fname(np)=fnamet
      lfn(np)=lfnt
      goto 12
 92   continue
      close(12)
      write(6,'(1H#,i6)') np
c
       do m=1,3
      cc(m)=0.0
      do i=1,np
      cc(m)=cc(m)+c(m,i)
      end do
      cc(m)=cc(m)/float(np)
       end do
      open(11,file='center')
      write(11,*) cc
      close(11)
c
      open(13,file=fdipoles)
      do i=1,np
       do m=1,3
       d(m,i)=d(m,i)-cc(m)*q(i)
       end do
       dlen=sqrt(d(1,i)**2+d(2,i)**2+d(3,i)**2)
      lfnt=lfn(i)
      write(13,'(a)') fname(i)(1:lfnt)
      write(13,502) q(i),dlen,(d(k,i),k=1,3),nq(i)
      end do
      close(13)
c
      scf=rprot**(-2)/3.0
c
      do i=1,np
      do j=1,i
      li=lfn(i)
      lj=lfn(j)
      write(6,'(a,x,a)') fname(i)(1:li),fname(j)(1:lj)
      d1(1)=d(1,i)
      d1(2)=d(2,i)
      d1(3)=d(3,i)
      d2(1)=d(1,j)
      d2(2)=d(2,j)
      d2(3)=d(3,j)
      call abc(q(i),d1,q(j),d2,scf,aa,bb,ab)
       write(6,612) 2*ab/(aa+bb),ab/sqrt(aa*bb),aa,bb,ab,aa,bb
c     ,,2.*float(npoi)/float(np1+np2),float(npoi)/sqrt(float(np1*np2)),
c     ,np1,np2,npoi
      end do
      end do
 612  format(2f7.3,5e15.7,2f7.3,3i7)
 502  format(f10.3,3x,f10.3,3x,3f10.3,i6)
c
      call exit
      end
      subroutine dipole1(fname,c,chnet,d,ncha)
      dimension xv(3),d(3),c(3)
c. c here - computed !!
      character*54 string
      character*40 fname 
	  character*4 atom
      character*3 aname,rname
      character*1 space 
      logical uncharged
	  data atom/'ATOM'/,space/' '/
c
        do m=1,3
      d(m)=0.
      c(m)=0.
        end do
      chnet=0.0
      achnet=0.0
      ncha=0
      na=0
      open(1,file=fname)
  1   continue
	  read(1,500,end=99) string
	  if(string(1:4).ne.atom) goto 1
	  read(string(31:54),501) xv
         na=na+1
         do m=1,3
         c(m)=c(m)+xv(m)
         end do
       aname=string(14:16)
       if(string(14:14).eq.space) then
       aname=string(15:17)
       if(string(15:15).eq.space) then
       aname=string(16:18)
       aname(3:3)=space
       end if
       end if
       rname=string(18:20)
      call chdef(aname,rname,uncharged,cha) 
      if(uncharged) goto 1
c      cha=cha*78.
c      write(6,'(2a3,f8.3)') aname,rname,cha
	  chnet=chnet+cha
	  achnet=achnet+abs(cha)
      d(1)=d(1)+xv(1)*cha
      d(2)=d(2)+xv(2)*cha
      d(3)=d(3)+xv(3)*cha
      ncha=ncha+1
	  goto 1
c
  99  continue
      close(1)
         do m=1,3
         c(m)=c(m)/float(na)
         end do
 500  format(a)
 501  format(3f8.3)
 502  format(f10.3,3x,f10.3,3x,3f10.3,i6)
c
	  return   
      end
	   subroutine chdef(aname,rname,uncharged,charge)
	   character*3 aname,rname
	   logical uncharged
       charge=0.00
	   uncharged=.true.
c
	   if(aname.eq.'N'.and.rname.eq.'NTR') uncharged=.false.
	   if(aname.eq.'N'.and.rname.eq.'NTR') charge=1.00  
	   if(aname.eq.'NZ'.and.rname.eq.'LYS') uncharged=.false.
	   if(aname.eq.'NZ'.and.rname.eq.'LYS') charge=1.00  
	   if(aname.eq.'NH1'.and.rname.eq.'ARG') uncharged=.false.
	   if(aname.eq.'NH1'.and.rname.eq.'ARG') charge=0.50  
	   if(aname.eq.'NH2'.and.rname.eq.'ARG') uncharged=.false.
	   if(aname.eq.'NH2'.and.rname.eq.'ARG') charge=0.50  
	   if(aname.eq.'OD1'.and.rname.eq.'ASP') uncharged=.false.
	   if(aname.eq.'OD1'.and.rname.eq.'ASP') charge=-0.50  
	   if(aname.eq.'OD2'.and.rname.eq.'ASP') uncharged=.false.
	   if(aname.eq.'OD2'.and.rname.eq.'ASP') charge=-0.50  
	   if(aname.eq.'OE1'.and.rname.eq.'GLU') uncharged=.false.
	   if(aname.eq.'OE1'.and.rname.eq.'GLU') charge=-0.50  
	   if(aname.eq.'OE2'.and.rname.eq.'GLU') uncharged=.false.
	   if(aname.eq.'OE2'.and.rname.eq.'GLU') charge=-0.50  
	   if(aname.eq.'OT1'.and.rname.eq.'CTR') uncharged=.false.
	   if(aname.eq.'OT1'.and.rname.eq.'CTR') charge=-0.50  
	   if(aname.eq.'OT2'.and.rname.eq.'CTR') uncharged=.false.
	   if(aname.eq.'OT2'.and.rname.eq.'CTR') charge=-0.50  
	   if(aname.eq.'O1'.and.rname.eq.'CTR') uncharged=.false.
	   if(aname.eq.'O1'.and.rname.eq.'CTR') charge=-0.50  
	   if(aname.eq.'O2'.and.rname.eq.'CTR') uncharged=.false.
	   if(aname.eq.'O2'.and.rname.eq.'CTR') charge=-0.50  
	   if(aname.eq.'O'.and.rname.eq.'CTR') uncharged=.false.
	   if(aname.eq.'O'.and.rname.eq.'CTR') charge=-0.50  
	   if(aname.eq.'OXT'.and.rname.eq.'CTR') uncharged=.false.
	   if(aname.eq.'OXT'.and.rname.eq.'CTR') charge=-0.50  
c
	   return
	   end
       subroutine abc(q1,d1,q2,d2,scf,aa,bb,ab)
       dimension d1(3),d2(3)
       aa=q1*q1+scf*(d1(1)*d1(1)+d1(2)*d1(2)+d1(3)*d1(3))
       bb=q2*q2+scf*(d2(1)*d2(1)+d2(2)*d2(2)+d2(3)*d2(3))
       ab=q1*q2+scf*(d1(1)*d2(1)+d1(2)*d2(2)+d1(3)*d2(3))
       return
       end
       subroutine pdbfl(name,l)
       character*256 name
       character*4 pdbe,sp4
       data pdbe/'.pdb'/,sp4/'    '/
c
       do 10 i=1,253
       l=i
       i3=i+3
       if(name(i:i3).eq.pdbe) goto 9
       if(name(i:i3).eq.sp4) goto 9
 10    continue
  9    continue
       l=l-1
       return
       end
