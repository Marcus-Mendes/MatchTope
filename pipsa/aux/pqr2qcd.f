c.
c. usage pqr2qcd filename-pqr filename-qcd
c.
c pqr
cATOM      3  C   TYR     1         6.194    14.920    54.642   0.550   1.70           C
cATOM      9 HB2  CYS     1        39.813     3.622     0.014   0.000   0.00           H
c123456789012345678901234567890123456789012345678901234567890123456789012345
c         1         2         3         4         5         6         7         8         9
c qcd
cATOM     3 ASN  O      18.000  41.580  20.651  -0.550   1.800
c1234567890123456789012345678901234567890123456789012345678901
c         1         2         3         4         5         6
c(a4,i6,1x,a4,1x,a4,  1x,5f8.3)
c.
c. converts pqr format file to qcd format
c.
c. All atoms with 0 charge and 0 radius will be deleted
c.
c. A simple filter, it will read first 30 characters as 
c. if they are exactly in pdb and extract ATOM,na,aname,rname,chname,nr
c  The coordinates will be read in free format
c. 
c. To be improved to
c. 1. avoid failure while reading coordinates (f8.3 as in pdb 
c.                                             sometimes leaves no spaces)
c. 2. avoid failure while reading atom name/number part
c.
      character*80 s
      character*80 fpqr,fqcd
      character*4 aname,rname
c
      if(iargc().lt.2) then
      write(6,*) 'usage: pqr2qcd filename-pqr filename-qcd'
      call exit
      end if
c
      call getarg(1,fpqr)
      call getarg(2,fqcd)
c
      open(11,file=fpqr)
      open(31,file=fqcd)
 1    continue
      read(11,'(a)',end=91) s
      if((s(1:4).ne.'ATOM').and.(s(1:4).ne.'HETA')) goto 1
c
      read(s(31:80),*) x,y,z,q,r
      if((q.eq.0).and.(r.eq.0)) goto 1
c
      read(s(7:11),'(i5)') ia
      read(s(13:16),'(a)') aname
      if(s(13:13).eq.' ') read(s(14:17),'(a)') aname
      read(s(18:21),'(a)') rname
c
      write(31,503) ia,rname,aname,x,y,z,q,r
      goto 1
 91   continue
      close(11)
      close(31)
c
 503  format('ATOM',i6,1x,a4,1x,a4,1x,5f8.3)
 504  format(a,3f8.3,'  1.00  5.00')
c
      call exit
      end

