c program delphi2uhbd
c
c Converts Delphi binary output grid to uhbd ascii format grid, 
c which can be then converted to binary uhbd by uhbd_asc2bin program 
c and used for PIPSA similarity analysis (or other purposes).
c
c The flag IBIOS while writing Delphi grid must be false  
c (i.e. grasp style phi-map is assumed).
c
c When compiled under linux, care must be taken depending on 
c under which platform the Dephi grids are computed and written.
c If SGI Irix version of Delphi was used, then the compilation 
c option -byteswapio must be used in order to grid be read correctly
c This option is not supported by g77, so one have to use, for example, 
c pgi or intel compilers. Executable compiled with pgi is supplied here.
c
c Should be no such problems, when both Delphi and delphi2uhbd are 
c the same platform executables.
c
c R. Gabdoulline, gabdoulline@eml.org, December 2004
c
      character*80 delphifile,uhbdfile
c
      character*20 uplbl
c      character*10 nxtlbl,character*60 toplbl
      character*10 nxtlbl
      character*60 toplbl
      real phi(65,65,65)
      character*16 botlbl
      real scale,oldmid(3)
c
      integer im, jm, km, grdflg, one, kk
      character*72 title
      real h, ox, oy, oz, uscale
      real dum2, dum3, dum4, dum5, dum6, dum7, dum8
      integer idum2, idum3, idum4
      data dum2,dum3,dum4,dum5,dum6,dum7,dum8
     ,    / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
      data idum2,idum3,idum4/0,0,0/
c
      if(iargc().lt.2) then
      write(6,601)
 601  format('Usage: delphi2uhbd delphi_grid_name uhbd_grid_name')
      call exit
      end if
c
      call getarg(1,delphifile)
      call getarg(2,uhbdfile)
      open (unit=11,file=delphifile,form='unformatted')
      write(6,602) delphifile,uhbdfile
 602  format('Converting delphi binary grid to uhbd ascii:',/,a,/,a)
c
      read (11) uplbl
      write(6,'(a)') uplbl
      read (11) nxtlbl,toplbl
      write(6,'(a)') nxtlbl
      write(6,'(a)') toplbl
      read (11) phi
      read (11) botlbl
      write(6,'(a)') botlbl
      read (11) dscale,oldmid
      write(6,*) dscale,oldmid
c
      close(11)
c
      write(title,'(a,2x,a)') nxtlbl,toplbl
      uscale=1.0
      grdflg=1
      one=1.0
      im=65
      jm=65
      km=65
      h=1.0/dscale
      ox=oldmid(1)-33.0*h
      oy=oldmid(2)-33.0*h
      oz=oldmid(3)-33.0*h
c
      open (unit=12,file=uhbdfile,status='new')
      write (12,77) title, uscale, dum2, grdflg, idum2,
     &               km, one, km, im, jm, km,
     &               h, ox, oy, oz,
     &               dum3, dum4, dum5, dum6, dum7,
     &               dum8, idum3, idum4
      do  20  k = 1 , km
      write(12,78) k,im,jm
      write(12,79) ((phi(i,j,k),i=1,im),j=1,jm)
   20 continue
c
      close(12)
c
   77 format(a72/2e12.6,5i7/3i7,4e12.6/4e12.6/2e12.6,2i7)
   78 format(3i7)
   79 format( (6(1x,e12.6)) )

      end
