c Reads kinemage and highlights in red the proteins 
c which names are given in the command line 
      character*128 s,se,namh(999),sa(999),namt
      dimension x(3,999)
c
      do i=1,128
      se(i:i)=' '
      end do
c
      nh=iargc()
      if(nh.lt.1) then
      write(6,*) 'usage: highlight_kin_points pname1 pname2 ...',
     ,' < in.kin > out.kin' 
      call exit
      end if
c
      do ih=1,nh
      call getarg(ih,namh(ih))
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
      do ih=1,nh
      namt=namh(ih)
      call slen(namt,nnamt)
      write(6,501) namt(1:nnamt)
        do ip=1,np
        s=sa(ip)
        ll=nnamt+5
        if(s(6:ll).eq.namt(1:nnamt)) then
        call slen(s,ns)
        write(6,'(a)') s(1:ns)
        end if
        end do
      end do
 501  format('@balllist {',a,'} color =red radius = 0.016')
      write(6,502)
      do ip=1,np
      s=sa(ip)
        call slen(s,ns)
        write(6,'(a)') s(1:ns)
      end do
 502  format('@balllist {ALL} color =white radius = 0.008')
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
