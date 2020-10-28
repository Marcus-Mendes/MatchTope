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
c reads sims.log and exp data, derives regression from 
c known kinetic parameters (rate ratio vs ep difference)
c predicts unknown kinetic parameters
c  
c the protein names must be arranged the same way as in 
c the file "names" used to compute similarity matrix 
c and known parameters should be written on the line after 
c the protein name
c if no number after the protein name - then this parameter
c is assumed to be unknown and predicted
c
c this version aasumes sum(a-b) proportional to log(Ka/Kb)
c proportionality constant derived from known cases
c and used for predicting unknown cases
c
c Input
c - command line arg -fl - the name of pipsa similarity log file
c   default is sims.log
c - command line arg -fe - the name of the file with experimental data 
c   where in one line the protein name is followed by experimental data
c   when available and with nothing or 0.0 when needs to be predicted
c - command line arg -fo - the name of the output file with 
c   correlation between log(k1/k2) and    (ep1-ep2)         
c   default is smNex2cor.out
c - command line arg -fp - the name of the output file with 
c   all predictions for each case
c   default is smNextopred.pre
c - command line arg -sn - the ep difference measure to be used:
c   1 - difference of average ep : av(ep1)-av(ep2)
c   2 - log of ratio of exp(ep)  : log (sum(exp(ep1))/sum(exp(ep2)))
c   3 - log of ratio of exp(-ep) : log (sum(exp(-ep1))/sum(exp(-ep2)))
c   default is 1
c - command line arg -rc - user-defined correlation coefficient
c   default is define it correlating known cases
c   correlation is Drate = rc*Dpotential, 
c   i.e. how much rate changes by 1 kcal/mole change in potential
c 
c Output
c - Standard out - prediction results with errors
c
c
      character*256 flog,fexp,fout,fpre
      character*256 s,arg,pnams(99)
      dimension rt(99),sm1(99,99),sm2(99,99),sm3(99,99)
      dimension sm4(99,99),sm5(99,99),sm6(99,99),sm7(99,99)
      dimension drt(9999),s1msi1(9999),s1msi2(9999)
      dimension rt1m(99),rt2m(99),rt1p(99),rt2p(99)
      dimension rt1(2,99),rt2(2,99),rtpre(99)
      dimension rt11(99),rt12(99),rt21(99),rt22(99)
      dimension amb(3)
      integer known(99)
c
      flog='sims.log'
      fexp='exp'
      fout='smNextopred.out'
      fpre='smNextopred.pre'
      isn=1
c
      call getarg(1,arg) 
      if(arg(1:2).eq.'-h') then
      write(6,499)
 499  format('usage: smNextopred -fl sims.log -fe exp -fo ',
     ,        'smNextopred.out -fp smNextopred.pre',
     ,        ' -sn sim_data_num -rc regr_coef')
      call exit
      end if
c
      aa=0.0
      nargs=iargc()
      do ia=1,nargs-1
      call getarg(ia,arg)
      if(arg(1:3).eq.'-fl') then
      call getarg(ia+1,flog)
      end if
      if(arg(1:3).eq.'-fe') then
      call getarg(ia+1,fexp)
      end if
      if(arg(1:3).eq.'-fo') then
      call getarg(ia+1,fout)
      end if
      if(arg(1:3).eq.'-fp') then
      call getarg(ia+1,fpre)
      end if
      if(arg(1:3).eq.'-sn') then
      call getarg(ia+1,arg)
      read(arg,*) isn
      end if
      if(arg(1:3).eq.'-rc') then
      call getarg(ia+1,arg)
      read(arg,*) rc
      end if
      end do
c
      open(1,file=flog)
      read(1,'(1x,i6)') np
      do ip=1,np
      do jp=1,ip
      read(1,*)
c      read(1,*) sm1(ip,jp),sm2(ip,jp)
      read(1,612) sm1p,sm2p,aa,bb,ab,aa0,bb0,
     ,           sm1s,sm2s,np1,np2,npp,amb
      sm1(ip,jp)=amb(isn)
      sm1(jp,ip)=-sm1(ip,jp)
      end do
      end do
      close(1)
 612  format(2f7.3,5e15.7,2f7.3,3i7,3e15.7)
c
c. read exp and assign known=0 for unknown cases
c
      open(2,file=fexp)
      do ip=1,np
      read(2,'(a)') s
      pnams(ip)=s
      call sfsp(s,nfsp)
      rtt=0.0
      read(s(nfsp:256),*,end=100) rtt
 100  continue
      if(rtt.eq.0.) then
      known(ip)=0
      rt(ip)=0.0
      else
      known(ip)=1
      rt(ip)=alog(rtt)
      end if
      end do
      close(2)
c
c. write correlation for known cases
c
      open(3,file=fout)
      write(3,501)
      nreg=0
      do ip=1,np
      if(known(ip).eq.0) goto 302
      write(3,*)
       do jp=1,np
       if(known(jp).eq.0) goto 301
       nreg=nreg+1
       drt(nreg)=rt(ip)-rt(jp)
       s1msi1(nreg)=sm1(ip,jp)
       write(3,'(3f8.3)') drt(nreg),s1msi1(nreg)
 301   continue
       end do 
 302  continue
      end do 
      close(3)
 501  format('# log(k1)-log(k2)  ',
     ,       '  phi1 - phi2      ')
 5011 format('# pred(i)')
c
c. construct regression model
c
      if(rc.ne.0.0) then
      d1=0.0
      a1=rc
      goto 333
      end if
c
      sxy=0.0
      sxx=0.0
      syy=0.0
      do ireg=1,nreg
      sxy=sxy+drt(ireg)*s1msi1(ireg)
      sxx=sxx+s1msi1(ireg)*s1msi1(ireg)
      syy=syy+drt(ireg)*drt(ireg)
      end do
      a1=sxy/sxx 
      d1=(a1**2*sxx-2.0*a1*sxy+syy)/syy
 333  continue
      write(6,*) 'Coeff=',a1,'  Disp=',d1
c
c. predict: all possibilities and av w err
c
      open(4,file=fpre)
      do ip=1,np
      if(known(ip).ne.0) goto 402
      rti=0.0
      npre=0
        do jp=1,np
        if(known(jp).eq.0) goto 401
        npre=npre+1
        si1=sm1(ip,jp)      
        drt1=si1*a1
        rt1(1,npre)=rt(jp)+drt1
 401    continue
        end do 
c
        do 404 lp=1,npre
        rtpre(lp)=rt1(1,lp)     
 404    continue
        rtpred=0.0
          do 405 lp=1,npre
          rtpred=rtpred+rtpre(lp)
 405      continue
          rtpred=rtpred/float(npre)
          drtpred=0.0
          do 406 lp=1,npre
          drtpred=drtpred+(exp(rtpre(lp))-exp(rtpred))**2
 406      continue
          drtpred=sqrt(drtpred/float(npre))
c.  
      s=pnams(ip)
      call sfsp(s,nfsp)
      write(6,505) s(1:nfsp),ip,exp(rtpred),
     ,drtpred
      write(4,505) s(1:nfsp),ip,exp(rtpred),
     ,drtpred
c
 402  continue
      end do 
      close(4)
 505  format(a,i4,2e13.4,99e13.4)

c
      call exit
      end
      subroutine sfsp(s,nfsp)
      character*256 s
      do i=1,256
      nfsp=i
      if(s(i:i).eq.' ') return
      end do
      return 
      end
      subroutine slen(s,ns)
      character*256 s
      do i=256,1,-1
      ns=i
      if(s(i:i).ne.' ') return
      end do
      return 
      end
