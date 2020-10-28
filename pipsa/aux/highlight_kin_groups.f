c Reads kinemage and the filenames wich contain 
c protein names to be in different groups 
c and writes a kinemage with groups colured differently
c
      character*128 s,se,sa(999),namg(99),namp
      character*128 namgt,nampt
      character*9 col(19),colt
      data col/'magenta  ','cyan     ','green   ','orange    ',
     ,         'purple   ','yellow   ','red     ','yellowtint',
     ,         'brown    ','sky      ','gray    ','gold      ',
     ,         'sea      ','hotpink  ','bluetint' ,'greentint',
     ,         'pink     ','pinktint ','white   '/
c
      do i=1,128
      se(i:i)=' '
      end do
c
      ng=iargc()
      if(ng.lt.1) then
      write(6,*) 'usage: highlight_kin_groups fname1 fname2 ... ',
     ,'< in.kin > out.kin' 
      call exit
      end if
      do ig=1,ng
      call getarg(ig,namg(ig))
      end do
c
      np=0
 1    continue
      read(5,'(a)',end=99) s
      if(s(1:6).eq.'@balll') goto 91
      call slen(s,ns)
      write(6,'(a)') s(1:ns)
      goto 1
 91   continue
c. read coordinates
      np=0
 2    continue
      read(5,'(a)',end=92) s
      np=np+1
      sa(np)=s
      goto 2
 92   continue
c
      write(6,502)
      do ip=1,np
      s=sa(ip)
        call slen(s,ns)
        write(6,'(a)') s(1:ns)
      end do
 502  format('@balllist {ALL} color =white radius = 0.008')
c
      do ig=1,ng
      namgt=namg(ig)
      open(1,file=namgt)
      call slen(namgt,lnamgt)
      igc=ig
      if(igc.gt.19) igc=igc-19
      colt=col(ig)
      call slen(colt,lcolt)
      write(6,501) namgt(1:lnamgt),colt(1:lcolt)
 111  continue
      read(1,'(a)',end=112) namp
      call slen(namp,lnamp)
        do ip=1,np
        s=sa(ip)
        ll=lnamp+5
        if(s(6:ll).eq.namp(1:lnamp)) then
        call slen(s,ns)
        write(6,'(a)') s(1:ns)
        end if
        end do
      goto 111
 112  continue
      close(1)
      end do
 501  format('@balllist {',a,'} color =',a,' radius = 0.010')
c
 99   continue
      call exit
      end
      subroutine slen(s,ns)
      character*128 s
      do i=128,1,-1
      ns=i
      if(s(i:i).ne.' ') return
      end do
      return
      end
