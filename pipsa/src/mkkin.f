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
      dimension dabn(nprmx,nprmx)
      dimension xav(3)
      real*8 dum2(2)
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
c       is1=i1+(i2-1)*ns
       is1=is1+1
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
c       do i=1,10
c       write(6,'(10f7.2)') (dab(i,j),j=1,10)
c       end do
c
          do 1111 i1=1,ns   
          do 2111 i2=1,ns   
          dab(i1,i2)=dab(i1,i2)/fmax
 2111 continue
 1111 continue
c
          dmax=0.
      do 112 i3=1,ns   
          d13=dab(i1max,i3)
          d23=dab(i2max,i3)
          aa=0.5*(1+d13**2-d23**2)
          bb=dsqrt(d13**2-aa**2)
          d2=bb
          if(d2.gt.dmax) then
          dmax=d2
          i3max=i3
          a=aa
          b=bb
          end if
 112  continue
c
          dmax2=0.
      do 1121 i4=1,ns   
      d14=dab(i1max,i4)
      d24=dab(i2max,i4)
      d34=dab(i3max,i4)
          cc=0.5*(1+d14**2-d24**2)
      dd=0.5/b*(d14**2-d34**2-a*(2.*cc-a)+b**2)
          ee=dsqrt(d14**2-cc**2-dd**2)
          d2=ee
          if(d2.gt.dmax2) then
          dmax2=d2
          i4max=i4
          c=cc
          d=dd
          e=ee
          end if
 1121  continue
c
c. largest tetraedr found and a,b,c,d,e are known
          d12=dab(i1max,i2max)
          d13=dab(i1max,i3max)
          d14=dab(i1max,i4max)
          d23=dab(i2max,i3max)
          d24=dab(i2max,i4max)
          d34=dab(i3max,i4max)
c      write(6,500) i1max,i2max,i3max,i4max,
c     ,d12,d13,d14,d23,d24,d34
c      write(6,500) i1max,i2max,i3max,i4max,a,b,c,d,e
c 500  format('#',4i5,6f8.3)
c
c. all others are projected:
c
          xav(1)=0.
          xav(2)=0.
          xav(3)=0.
      do 50 ii=1,ns   
          d1=dab(i1max,ii)
          d2=dab(i2max,ii)
          d3=dab(i3max,ii)
          d4=dab(i4max,ii)
          x=(1.+d1**2-d2**2)/2.
          y=(b**2+d1**2-d3**2-a*(2*x-a))/2./b
          z=(d3**2-d4**2-(a-x)**2+(c-x)**2-
     -(b-y)**2+(d-y)**2+e**2)/2./e
          xx(1,ii)=x*fmax
          xx(2,ii)=y*fmax
          xx(3,ii)=z*fmax
          xav(1)=xav(1)+x*fmax
          xav(2)=xav(2)+y*fmax
          xav(3)=xav(3)+z*fmax
 50   continue
      xav(1)=xav(1)/float(ns)   
      xav(2)=xav(2)/float(ns)   
      xav(3)=xav(3)/float(ns)   
c      write(6,501) xav
c 501  format('#',6f8.3)
 502  format(i5,3f11.6) 
c
 999  continue
c
      lpnmx=0
      do i=1,ns
      pnamew=pnames(i,i)
      call namelength(pnamew,lpn)
      if(lpn.gt.lpnmx) lpnmx=lpn
      end do
c
      dmoo=0.d0
      dmnn=0.d0
      dmon=0.d0
      do i1=1,ns
      do i2=1,ns
      dmoo=dmoo+dab(i1,i2)*dab(i1,i2)
      d12=sqrt((xx(1,i1)-xx(1,i2))**2+
     +         (xx(2,i1)-xx(2,i2))**2+
     +         (xx(3,i1)-xx(3,i2))**2)
      dmnn=dmnn+d12*d12
      dmon=dmon+d12*dab(i1,i2) 
      end do
      end do
      accuracy=2.d0*dmon/(dmoo+dmnn)
c
      call wheader(accuracy,xav)
c
      do i=1,ns
      pnamew=pnames(i,i)
      call nameclear(pnamew) 
      write(6,601) i,pnamew(1:lpnmx),(xx(m,i),m=1,3)
      end do
c
 601  format('{',i3,1x,a,' }',3f8.3)
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
      subroutine wheader(accuracy,xav)
      implicit real*8 (a-h,o-z)
      dimension xav(3)
c
      write(6,501) accuracy,xav
 501  format(
     1'@text',/,
     1'    Representation of the proteins as a points in 3D',/,
     1'based on similarities of their interaction fields.',/,
     1'    The distance between proteins i and j is defined as',/,
     1'D(i,j) = sqrt(2-2*SI(i,j)), where SI(i,j) is the',/, 
     1'similarity index of their interaction fields. Four',/, 
     1'points defining the 3D space were chosen as follows:',/,
     1'    a) The first two points (1 and 2) are chosen as',/, 
     1'those for which the distance D(1,2), is largest. This',/, 
     1'D(1,2) is then used for normalization of all other',/, 
     1'distances.',/,
     1'    b) Point 3 is chosen so that the area of the',/, 
     1'triangle formed by points 1,2 and 3 is the largest.',/,
     1'    c) Point 4 is chosen so that the volume of the',/,
     1'tetrahedron formed by points 1,2,3 and 4 is the largest.',/,
     1'    The 3D coordinates assigned to these 4 points are',/,
     1'(0,0,0), (1,0,0), (a,b,0), (c,d,e), where a, b, c, d',/, 
     1'and e are computed unambiguously from the 6 pairwise',/, 
     1'distances.  The other compounds are projected onto this',/, 
     1'defined 3D space with coordinates determined from their',/, 
     1'SI(i,j) values.',/,
     1'@kinemage 1',/,
     1'@caption',/,
     1'    Representation of the proteins as a points in 3D',/,
     1'based on similarities of their interaction fields.',/,
     1'(see MAGE Text)',/,
     1'    This representation is approximate and distorts',/, 
     1'original distance matrix.  Accuracy may be measured',/, 
     1'by similarity of the original matrix to the matrix',/, 
     1'derived from 3D distribution; Hodgkin index of this',/,
     1'similarity is : ',f8.3,/,  
     1'(more accurate means this index close to 1.000)',/,
     1'@onewidth',/,
     1'@matrix',/,
     1'1 0 0  0 1 0  0 0 1',/,
     1'@viewid {x y}',/,
     1'@2zoom 1.0',/,
     1'@2zslab 2',/,
     1'@2center',3f8.3,/,
     1'@2matrix',/,
     1'1 0 0   0 0 1  0 -1 0',/,
     1'@2viewid {x z}',/,
     1'@3zslab 2',/,
     1'@3matrix',/,
     1'0 0 1  0 1 0  -1 0 0',/,
     1'@3viewid {z y}',/,
     1'@group {ALL}',/,
     1'@balllist {ALL} color= white radius= 0.008'
     1)
c
      return
      end
