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
      implicit real*8 (a-h,o-z)
c. 
c. this version reads one of 14 results
c..
c.. reads from a single file :
c.. pair of names, 
c.. the matrix of ab,aa,bb, converts this to d(ab)
c.. convert d(ab)'s to 3D coordinates and writes in
c.. kinemage format
c..
      include "maxdim.inc"
      dimension dab1(nprmx2),xx(3,nprmx)
      character*40 argument
      character*30 pname1(nprmx2),pnamer
      character*30 pnames(nprmx,nprmx),pnamew
      dimension dab(nprmx,nprmx)
      dimension xav(3)
      real*8 dum(2)
      real*8 ans1,ans
c        
      call getarg(1,argument)
      read(argument(1:40),'(i40)') inom
      call getarg(2,argument)
      read(argument(1:40),'(i40)') nnom
c      write(6,*) ' analysing part #', inom,' of',nnom
c
      read(5,'(1x,i6)',end=99) np
      is1=0
 1    continue
      read(5,'(a)',end=99) pnamer
        do iom=1,inom-1
        read(5,*)
        end do
      read(5,*,end=99) si_h,si_c,aa,bb,ab,aa0,bb0
        do iom=inom+1,nnom
        read(5,*)
        end do
      is1=is1+1
      pname1(is1)=pnamer
c. distance based on norms at the intersection of skins
c      dab1(is1)=dsqrt(aa+bb-2.d0*ab)
c. distance based on norms on the whole skins
c      dab1(is1)=dsqrt(aa0+bb0-2.d0*ab)
c. distance based on similarity index
      dab1(is1)=dsqrt(2.d0-2.d0*si_c)
      goto 1
 99   continue
c
c       write(6,*) ' matrix read'
c
       ns1=is1
       ns1c=(np*np+np)/2 
       ns=np
       if(ns1c.ne.ns1) then
       write(6,*) ' bad matrix: expected',ns1c,'  points'
       write(6,*) 'got:', ns1,'  - please correct'
       goto 999
       end if
c
         fmax=0.0d0
         i1max=1
         i2max=2
       is1=0
       do i1=1,ns
       do i2=1,i1
       is1=is1+1
c       is1=i1+(i2-1)*ns
       dab(i1,i2)=dab1(is1)
       dab(i2,i1)=dab1(is1)
         if(dab(i1,i2).gt.fmax) then
         fmax=dab(i1,i2)
         i1max=i1
         i2max=i2
         end if
       pnames(i1,i2)=pname1(is1)
       end do 
       end do 
c
      do i=1,ns
      write(6,611) (dab(i,j),j=1,ns)
      end do
 611  format(999f8.3)
 999  continue
c
      call exit
      end 
c
      subroutine namelength(name,lpn)
      implicit real*8 (a-h,o-z)
      parameter (n30=30)
      character*30 name
      character*1 sp
      data sp/' '/
c
      ke=n30+1
      do k=1,n30
      if(name(k:k).eq.sp) then
      ke=k
      goto 41
      end if
      end do
 41   continue
      lpn=ke-1
      return
      end
      subroutine nameclear(name)
      implicit real*8 (a-h,o-z)
      parameter (n30=30)
      character*30 name
      character*1 sp
      data sp/' '/
c
      ke=n30+1
      do k=1,n30
      if(name(k:k).eq.sp) then
      ke=k
      goto 41
      end if
      end do
 41   continue
      do k=ke,n30
      name(k:k)=sp
      end do
      return
      end
