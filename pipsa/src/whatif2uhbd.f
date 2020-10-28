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
      program whatif2uhbd
      implicit real*4 (a-h,o-z)
c
c originally written by V. Lounnas
c additions by RRG
c will work only for a single chain proteins !!!
c waters are however permitted
c
      character*4 wrkres,wrkres0,wrkres1,wrkat,wrkat0,wrkat1
      character*4 atom,hid,hie,hip,hit
      character*5 wresi
      character*66 line,linet
      character*66 lines(17)
      logical endfil
      data atom/'ATOM'/
      data hid/'HID '/,hie/'HIE '/,hip/'HIP'/
c
      open(12,file='tmp.pdb',form='formatted',status='unknown')
c
      wrkat1 = ' ** '
      wresi  = ' ** '
      ifirst=0
10    read(5,100,end=999) line 
      if(line(1:4).ne.atom) then
      write(12,100) line
      goto 10
      end if
      if(ifirst.eq.0) then
      ifirst=1
      read(line(22:26),'(1x,i4)') irmin
      end if
      wrkres = line(18:21)
      wrkat  = line(13:16) 
      if(wrkat(2:3).eq.'O2') then
      write(line(22:26),'(a)') wresi
      read(wresi(1:5),'(1x,i4)') irmax
      else
      if(line(18:21).ne.'HOH ') wresi  = line(22:26)
      end if
c
      if(wrkat(2:2).ne.'C'.and.wrkat(2:2).ne.'N'.and.
     +   wrkat(2:2).ne.'O'.and.wrkat(2:2).ne.'S') then 
c      wrkat0 = ' H  '
c
      wrkat0 = ' H'//wrkat1(3:3)
c      if(wrkat1(2:2).eq.'N') then
         if(wrkat1(4:4).ne.' ') then 
         wrkat0 = ' H'//wrkat1(3:4)
         endif
c      endif
c
      if(wrkat1(2:2).eq.'N'.and.wrkat1(3:3).ne.' '.and.
     +   wrkat1(4:4).eq.' ') then 
      wrkat0 = wrkat0(1:3)//'1'
      endif
c
      if(wrkat1(2:3).eq.'OW') then 
      wrkat0 = wrkat0(1:3)//'1'
      endif
c
      if(wrkat1(2:2).eq.'H'.and.wrkat1(4:4).eq.'1') then
      wrkat0 = wrkat0(1:3)//'2'
      endif
c
      if(wrkat1(2:2).eq.'H'.and.wrkat1(4:4).eq.'2') then
      wrkat0 = wrkat0(1:3)//'3'
      endif
c
c
      if(wrkat1(2:2).eq.'N'.and.wrkat1(3:3).ne.' '.and.
     +   wrkat1(4:4).ne.' ') then    
      wrkat0(1:4) = wrkat0(2:4)//'1'
      endif
c
      if(wrkat1(1:1).eq.'H'.and.wrkat1(4:4).eq.'1') then
      wrkat0(1:4) = wrkat1(1:3)//'2'
      endif
c
c particular cases:
c      
      if(wrkat1(2:4).eq.'NE1') then 
      wrkat0 = ' HE1'
      endif
c
      if(wrkat1(2:4).eq.'OE2') then 
      wrkat0 = ' HE '
      endif
c
      if(wrkat1(2:4).eq.'OD2') then 
      wrkat0 = ' HD '
      endif
c
      if(wrkat1(2:4).eq.'OD1'.and.wrkat0(2:3).eq.'HD') then
      wrkat0 = ' HD '
      endif
c
      if(wrkat1(2:4).eq.'OE1'.and.wrkat0(2:3).eq.'HE') then
      wrkat0 = ' HE '
      endif
c
      if(wrkat1(2:4).eq.'NE ') then
      wrkat0 = ' HE '
      endif
c
      if(wrkres(1:3).eq.'HIS'.and.wrkat1(2:4).eq.'NE2') then
      wrkat0 = ' HE2'
      endif
c
      if(wrkres(1:3).eq.'HIS'.and.wrkat1(2:4).eq.'ND1') then
      wrkat0 = ' HD1'
      endif
c
      if(wrkres(1:3).eq.'ASP'.and.wrkat1(2:4).eq.'OD1') then
      wrkat0 = ' HD1'
      endif
c
      if(wrkres(1:3).eq.'ASP'.and.wrkat1(2:4).eq.'OD2') then
      wrkat0 = ' HD2'
      endif
c
      else
      wrkat0 = wrkat
      endif
c
      write(12,101) line(1:12),wrkat0,' ',wrkres,line(22:66)
c      write(6,101) line(1:12),wrkat0,' ',wrkres,line(22:66)
c
      wrkat1  = wrkat0
      goto 10
c
999   continue 
c      write(6,*) irmin,irmax
c.. post processing (12) "tmp.pdb" - rename HIS to HID or HIE and NTRM and CTRM
      rewind (12) 
 11   continue
      read(12,100,end=99) line
 111  continue
      if(line(1:4).ne.atom) then
      write(6,100) line
      goto 11
      end if
c. no tracing for waters
      if(line(18:21).eq.'HOH ') then
      write(6,100) line
      goto 11
      end if 
c
c
c. trace HID or HIE
      if(line(18:21).ne.'HIS ') goto 15
c. here read anything
      endfil=.false.
      nn=1
      lines(nn)=line
      read(line(22:26),'(1x,i4)') irthis
 151  continue
      read(12,100,end=152) line
      read(line(22:26),'(1x,i4)') ir
      if(ir.ne.irthis) goto 153
      nn=nn+1 
      lines(nn)=line
      goto 151
 152  continue 
      endfil=.true.
 153  continue 
c
      itishd=0
      itishe=0
      do 161 ii=1,nn
      linet=lines(ii)
      if(linet(13:16).eq.' HE2') itishe=1
      if(linet(13:16).eq.' HD1') itishd=1
 161  continue
      hit=hid
      if(itishe.eq.1) then
      if(itishd.ne.1) hit=hie 
      if(itishd.eq.1) hit=hip
      end if
      do 171 ii=1,nn
      linet=lines(ii)
      write(linet(18:21),100) hit
      call term(linet,irmin,irmax)
      write(6,100) linet
 171  continue
c      
      if(endfil) goto 99
      goto 111
 15   continue
c. nothing special found
c. check for a terminal though
      call term(line,irmin,irmax)
      write(6,100) line
      goto 11
 99   continue
      close (12)
      stop
100   format(a)
101   format(a12,a4,a1,a4,a45)
200   format(4x,i7,1x,a5,a4,1x,i4,4x,3f8.3)
      end
      subroutine term(line,irmin,irmax)
      character*66 line
      character*4 atom,ntr,ntrm,ntrp,ctr,ctrm,ctrp
      data atom/'ATOM'/,ntrm/'NTRM'/,ctrm/'CTRM'/
      data ntrp/'NTRP'/,ctrp/'CTRP'/
c
c. trace the NTRM
      read(line(22:26),'(1x,i4)') ir
      if(ir.ne.irmin) goto 12
      ntr=ntrm
      if(line(18:21).eq.'PRO ') then
      ntr=ntrp
      if(line(14:15).eq.'CD') write(line(18:21),100) ntr
      end if
      if(line(14:15).eq.'H ') write(line(18:21),100) ntr
      if(line(14:15).eq.'N ') write(line(18:21),100) ntr
      if(line(14:15).eq.'CA') write(line(18:21),100) ntr
      if(line(14:15).eq.'C ') write(line(18:21),100) ntr
      if(line(14:15).eq.'O ') write(line(18:21),100) ntr
c      write(6,100) line
      goto 14
 12   continue
c. trace the CTRM
      if(line(14:15).ne.'O2') goto 13
      write(line(18:21),100) ctrm
c      write(6,100) line
      goto 14
 13   continue
      read(line(22:26),'(1x,i4)') ir
      if(ir.ne.irmax) goto 14
      ctr=ctrm
      if(line(18:21).eq.'PRO ') then
      ctr=ctrp
      if(line(14:15).eq.'CD') write(line(18:21),100) ctr
      end if
      if(line(14:15).eq.'N ') write(line(18:21),100) ctr
      if(line(14:15).eq.'H ') write(line(18:21),100) ctr
      if(line(14:15).eq.'CA') write(line(18:21),100) ctr
      if(line(14:15).eq.'C ') write(line(18:21),100) ctr
      if(line(14:15).eq.'O ') write(line(18:21),100) ctr
c      write(6,100) line
 14   continue
100   format(a)
      return
      end
