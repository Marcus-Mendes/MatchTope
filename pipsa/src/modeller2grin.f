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
c. reads output from MODELLER and changes the names of atoms 
c. to the names recognizable by GRIN routine of GRID
c
cATOM      2  CA  ASN     1      16.035  -8.224   5.143  1.00  3.50       1SG   3
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7
c(a4,2x,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2)
c atom name: (14:16)
c
      parameter(nchmax=99)
      character*80 s
      character*3 na_old(nchmax) 
      character*3 na_new(nchmax) 
c
c. edit these if more needed
      nch=2
      na_old(1)='OT1'
      na_new(1)='O  '
      na_old(2)='OT2'
      na_new(2)='OXT'
c
 1    continue
      read(5,'(a)',end=91) s
      if(s(1:4).eq.'ATOM') then
        do i=1,nch
        if(s(14:16).eq.na_old(i)) s(14:16)=na_new(i)
        end do 
      end if
      write(6,'(a)') s
      goto 1
 91   continue
      call exit
      end
