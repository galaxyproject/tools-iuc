*************************************************************************
*     This is a program for protein sequence alignment using the standard
*     Needleman-Wunsch dynamic programming. The mutation matrix is from 
*     PAM30 with gap openning penaly=-11 and gap extension panalty=-1. 
*     The program can be freely copied and modified provided the notices 
*     on the head are retained. Comments and bug report should be addressed 
*     to Yang Zhang (Email: zhng@umich.edu). Last update is in 2010/08/03.
*
*     Instructions:
*     1, the program can be compiled by 
*        >gfortran -static -O3 -ffast-math -lm -o align align.f
*     2, simply running the program will give a brief note on how to use it
*     3, You can run the program in following convenient ways:
*        >align F1.fasta F2.fasta (align two sequences in fasta file)
*        >align F1.pdb F2.pdb 1   (align two sequences in PDB file)
*        >align F1.fasta F2.pdb 2 (align Sequence 1 in fasta and 2 in pdb)
*        >align GKDGL EVADELVSE 3 (align sequences typed by keyboard)
*        >align GKDGL F.fasta 4   (align Seq-1 by keyboard and 2 in fasta)
*        >align GKDGL F.pdb 5     (align Seq-1 by keyboard and 2 in pdb)
*************************************************************************
      
      program compares
      PARAMETER(ndim=6000)
      parameter(naa=24) !number of amino acid
      common/dpc/score(ndim,ndim),gap_open,gap_extn,j2i(ndim)
     &     ,nseq1,nseq2
      common/matra/imut(naa,naa)     !b,z,x are additional

      integer seq1(ndim),seq2(ndim)
      character*10000 fnam1,fnam2,fnam3,fnam4
      character*10000 s
      character*3 aa(naa),aanam
      character seqw(naa),upper
      character*100 du,ad
      character sequenceA(ndim),sequenceB(ndim),sequenceM(ndim)

*---------------------- 24 amino acids ---------------------
      data aa/'ALA','ARG','ASN','ASP','CYS','GLN','GLU',
     &     'GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
     &     'THR','TRP','TYR','VAL','ASX','GLX','UNK','STOP'/
      data seqw/'A','R','N','D','C','Q','E','G','H','I','L','K',
     &     'M','F','P','S','T','W','Y','V','B','Z','X','*'/

      call getarg(1,fnam1)
      call getarg(2,fnam2)
      call getarg(3,fnam3)
      call getarg(4,fnam4)

      if(fnam1.eq.' ')then
         write(*,*)'align F1.fasta F2.fasta ',
     &        '(align two sequences in fasta file)'
         write(*,*)'align F1.pdb F2.pdb 1   ',
     &        '(align two sequences in PDB file)'
         write(*,*)'align F1.fasta F2.pdb 2 ',
     &        '(align Sequence 1 in fasta and 2 in pdb)'
         write(*,*)'align GKDGL EVADELVSE 3 ',
     &        '(align two sequences typed by keyboard)'
         write(*,*)'align GKDGL F.fasta 4   ',
     &        '(align Sequence 1 by keyboard and 2 in fasta)'
         write(*,*)'align GKDGL F.pdb 5     ',
     &        '(align Sequence 1 by keyboard and 2 in pdb)'
         goto 999
      endif
      
*1******* read sequences ------------------------->
      if(fnam3.eq.'5')then      !direct, 555555555555555555
***   read sequence1:
         i=0
         do k=1,10000
            fnam1(k:k)=upper(fnam1(k:k))
            do j=1,naa
               if(fnam1(k:k).eq.seqw(j))then
                  i=i+1
                  seq1(i)=j
                  goto 5
               endif
            enddo
            if(fnam1(k:k).ne.'-')goto 55 !same time
 5          continue
            if(i.ge.ndim)goto 55
         enddo
 55      continue
         nseq1=i
***   read sequence2:
         open(unit=10,file=fnam2,status='old')
         i=0
         do while (.true.)
            read(10,1,end=551) s
            if(i.gt.0.and.s(1:3).eq.'TER')goto 551
            if(s(1:3).eq.'ATO')then
               if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.
     &              or.s(13:16).eq.'  CA')then
                  i=i+1
                  read(s,111)du,aanam
                  do j=1,naa
                     if(aanam.eq.aa(j))seq2(i)=j
                  enddo
               endif
            endif
            if(i.ge.ndim)goto 551
         enddo
 551     continue
         close(10)
         nseq2=i
      elseif(fnam3.eq.'4')then  !direct, 444444444444444444444444444
***   read sequence1:
         i=0
         do k=1,10000
            fnam1(k:k)=upper(fnam1(k:k))
            do j=1,naa
               if(fnam1(k:k).eq.seqw(j))then
                  i=i+1
                  seq1(i)=j
                  goto 4
               endif
            enddo
            if(fnam1(k:k).ne.'-')goto 44
 4          continue
            if(i.ge.ndim)goto 44
         enddo
 44      continue
         nseq1=i
***   read sequence2:
         open(unit=10,file=fnam2,status='old')
         i=0
         do while(.true.)
            read(10,1,end=443)s
            if(s(1:1).eq.'>')goto 442
            do k=1,10000
               s(k:k)=upper(s(k:k))
               do j=1,naa
                  if(s(k:k).eq.seqw(j))then
                     i=i+1
                     seq2(i)=j
                     goto 441
                  endif
               enddo
               if(s(k:k).ne.'-')goto 442 !same time
 441           continue
            enddo
 442        continue
            if(i.ge.ndim)goto 443
         enddo
 443     continue
         close(10)
         nseq2=i


      elseif(fnam3.eq.'3')then  !direct, 33333333333333333333333333333333333
***   read sequence1:
         i=0
         do k=1,10000
            fnam1(k:k)=upper(fnam1(k:k))
            do j=1,naa
               if(fnam1(k:k).eq.seqw(j))then
                  i=i+1
                  seq1(i)=j
                  goto 3
               endif
            enddo
            if(fnam1(k:k).ne.'-')goto 33
 3          continue
            if(i.ge.ndim)goto 33
         enddo
 33      continue
         nseq1=i
***   read sequence2:
         i=0
         do k=1,10000
            fnam2(k:k)=upper(fnam2(k:k))
            do j=1,naa
               if(fnam2(k:k).eq.seqw(j))then
                  i=i+1
                  seq2(i)=j
                  goto 331
               endif
            enddo
            if(fnam2(k:k).ne.'-')goto 332
 331        continue
            if(i.ge.ndim)goto 332
         enddo
 332     continue
         nseq2=i
      elseif(fnam3.eq.'1')then  !pdb,pdb, 11111111111111111111111111111
***   read sequence1:
         open(unit=10,file=fnam1,status='old')
         i=0
         do while (.true.)
            read(10,1,end=11) s
            if(i.gt.0.and.s(1:3).eq.'TER')goto 11
            if(s(1:3).eq.'ATO')then
               if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.
     &              or.s(13:16).eq.'  CA')then
                  i=i+1
                  read(s,111)du,aanam
                  do j=1,naa
                     if(aanam.eq.aa(j))seq1(i)=j
                  enddo
               endif
            endif
            if(i.ge.ndim)goto 11
         enddo
 1       format(A10000)
 11      continue
 111     format(A17,A3)
         close(10)
         nseq1=i
***   read sequence2:
         open(unit=10,file=fnam2,status='old')
         i=0
         do while (.true.)
            read(10,1,end=112) s
            if(i.gt.0.and.s(1:3).eq.'TER')goto 112
            if(s(1:3).eq.'ATO')then
               if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.
     &              or.s(13:16).eq.'  CA')then
                  i=i+1
                  read(s,111)du,aanam
                  do j=1,naa
                     if(aanam.eq.aa(j))seq2(i)=j
                  enddo
               endif
            endif
            if(i.ge.ndim)goto 112
         enddo
 112     continue
         close(10)
         nseq2=i
      elseif(fnam3.eq.'2')then  !seq,pdb 2222222222222222222222222222222
***   read sequence1:
         open(unit=10,file=fnam1,status='old')
         i=0
         do while(.true.)
            read(10,1,end=221)s
            if(s(1:1).eq.'>')goto 22
            do k=1,10000
               s(k:k)=upper(s(k:k))
               do j=1,naa
                  if(s(k:k).eq.seqw(j))then
                     i=i+1
                     seq1(i)=j
                     goto 2
                  endif
               enddo
               if(s(k:k).ne.'-')goto 22
 2             continue
            enddo
 22         continue
            if(i.ge.ndim)goto 221
         enddo
 221     continue
         close(10)
         nseq1=i
***   read sequence2:
         open(unit=10,file=fnam2,status='old')
         i=0
         do while (.true.)
            read(10,1,end=222) s
            if(i.gt.0.and.s(1:3).eq.'TER')goto 222
            if(s(1:3).eq.'ATO')then
               if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.
     &              or.s(13:16).eq.'  CA')then
                  i=i+1
                  read(s,111)du,aanam
                  do j=1,naa
                     if(aanam.eq.aa(j))seq2(i)=j
                  enddo
               endif
            endif
            if(i.ge.ndim)goto 222
         enddo
 222     continue
         close(10)
         nseq2=i
      else                      !seq,seq 00000000000000000000000000000000
***   read sequence1:
         open(unit=10,file=fnam1,status='old')
         i=0
         do while(.true.)
            read(10,1,end=881)s
            if(s(1:1).eq.'>')goto 88
            do k=1,10000
               s(k:k)=upper(s(k:k))
               do j=1,naa
                  if(s(k:k).eq.seqw(j))then
                     i=i+1
                     seq1(i)=j
                     goto 8
                  endif
               enddo
               if(s(k:k).ne.'-')goto 88
 8             continue
            enddo
 88         continue
            if(i.ge.ndim)goto 881
         enddo
 881     continue
         close(10)
         nseq1=i
***   read sequence2:
         open(unit=10,file=fnam2,status='old')
         i=0
         do while(.true.)
            read(10,1,end=884)s
            if(s(1:1).eq.'>')goto 883
            do k=1,10000
               s(k:k)=upper(s(k:k))
               do j=1,naa
                  if(s(k:k).eq.seqw(j))then
                     i=i+1
                     seq2(i)=j
                     goto 882
                  endif
               enddo
               if(s(k:k).ne.'-')goto 883
 882           continue
            enddo
 883        continue
            if(i.ge.ndim)goto 884
         enddo
 884     continue
         close(10)
         nseq2=i
      endif
      
*2**   read mutation matrix ---------->
      call matrix               !take pam
***   set unit mutation matrix ---------->
c      do i=1,naa
c         do j=1,naa
c            imut(i,j)=0
c         enddo
c      enddo
c      do i=1,naa
c         imut(i,i)=1
c      enddo

*3**   score------------------>
      do i=1,nseq1
         do j=1,nseq2
            score(i,j)=imut(seq1(i),seq2(j))
         enddo
      enddo

*4*****************************************************************
*     dynamatic program:
******************************************************************
      gap_open=-11
      gap_extn=-1
      call DP(score0)           !W(k)=Go+Ge*k1+Go+Ge*k2, standard NW
c      call DPalt(score0)        !W(k)=Go+Ge*k1+Ge*k2, alternative NW
      
*5**   calculate sequence identity---------------------------->
      L_id=0
      L_ali=0
      do j=1,nseq2
         if(j2i(j).gt.0)then
            i=j2i(j)
            L_ali=L_ali+1
            if(seq1(i).eq.seq2(j))L_id=L_id+1
         endif
      enddo
      
      write(*,*)
      write(*,101)nseq1,fnam1
 101  format('Length of sequence 1: ',I4,' ->',A10)
      write(*,102)nseq2,fnam2
 102  format('Length of sequence 2: ',I4,' ->',A10)
      write(*,103)L_ali
 103  format('Aligned length: ',I4)
      write(*,104)L_id
 104  format('Identical length: ',I4)
      write(*,105)float(L_id)/(nseq2+0.00000001),L_id,nseq2
 105  format('Sequence identity: ',F8.3,' (=',I4,'/',I4,')')
      write(*,*)

*6******************************************************************
***   output aligned sequences
      k=0                       !final aligned order
      i=1                       !on sequence 1
      j=1                       !on sequence 2
 800  continue
      if(i.gt.nseq1.and.j.gt.nseq2)goto 802
      if(i.gt.nseq1.and.j.le.nseq2)then !unaligned C on 1
         k=k+1
         sequenceA(k)='-'
         sequenceB(k)=seqw(seq2(j))
         sequenceM(k)=' '
         j=j+1
         goto 800
      endif
      if(i.le.nseq1.and.j.gt.nseq2)then !unaligned C on 2
         k=k+1
         sequenceA(k)=seqw(seq1(i))
         sequenceB(k)='-'
         sequenceM(k)=' '
         i=i+1
         goto 800
      endif
      if(i.eq.j2i(j))then    !if aligned
         k=k+1
         sequenceA(k)=seqw(seq1(i))
         sequenceB(k)=seqw(seq2(j))
         if(seq1(i).eq.seq2(j))then !identical
            sequenceM(k)=':'
         else
            sequenceM(k)=' '
         endif
         i=i+1
         j=j+1
         goto 800
      elseif(j2i(j).lt.0)then !if gap on 1
         k=k+1
         sequenceA(k)='-'
         sequenceB(k)=seqw(seq2(j))
         sequenceM(k)=' '
         j=j+1
         goto 800
      elseif(j2i(j).gt.0)then !if gap on 2
         k=k+1
         sequenceA(k)=seqW(seq1(i))
         sequenceB(k)='-'
         sequenceM(k)=' '
         i=i+1
         goto 800
      endif
 802  continue

      write(*,601)(sequenceA(i),i=1,k)
      write(*,601)(sequenceM(i),i=1,k)
      write(*,601)(sequenceB(i),i=1,k)
      write(*,602)(mod(i,10),i=1,k)
 601  format(2000A1)
 602  format(2000I1)
      write(*,*)

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c      STOP
 999  END
      
********************************************************************
*     This is a standard Needleman-Wunsch dynamic program (by Y. Zhang 2005)
*     1. Count multiple-gap.
*     2. The gap penality W(k)=Go+Ge*k1+Go+Ge*k2 if gap open on both sequences
*     
*     Input: score(i,j), gap_open, gap_extn
*     Output: j2i(j)
*     idir(i,j)=1,2,3, from diagonal, horizontal, vertical
*     val(i,j) is the cumulative score of (i,j)
********************************************************************
      subroutine DP(score0)
      PARAMETER(ndim=6000)
      common/dpc/score(ndim,ndim),gap_open,gap_extn,j2i(ndim)
     &     ,nseq1,nseq2
      
      dimension val(0:ndim,0:ndim),idir(0:ndim,0:ndim)
      dimension jpV(0:ndim,0:ndim),jpH(0:ndim,0:ndim)
      dimension preV(0:ndim,0:ndim),preH(0:ndim,0:ndim)
      real D,V,H
      
ccc   initializations --------------->
      val(0,0)=0.0 
      do i=1,nseq1
         val(i,0)=gap_extn*i
         preV(i,0)=val(i,0) !not use preV at the beginning
         idir(i,0)=0            !useless
         jpV(i,0)=1             !useless
         jpH(i,0)=i             !useless
      enddo
      do j=1,nseq2
         val(0,j)=gap_extn*j
         preH(0,j)=val(0,j)
         idir(0,j)=0
         jpV(0,j)=j
         jpH(0,j)=1
      enddo
      
ccc   DP ------------------------------>
      do 111 j=1,nseq2
         do 222 i=1,nseq1
ccc   D=VAL(i-1,j-1)+SCORE(i,j)--------------->
            D=val(i-1,j-1)+score(i,j) !from diagonal, val(i,j) is val(i-1,j-1)
ccc   H=H+gap_open ------->
            jpH(i,j)=1
            val1=val(i-1,j)+gap_open !gap_open from both D and V
            val2=preH(i-1,j)+gap_extn !gap_extn from horizontal
            if(val1.gt.val2) then !last step from D or V
               H=val1
            else                !last step from H
               H=val2
               if(i.gt.1)jpH(i,j)=jpH(i-1,j)+1 !record long-gap
            endif
ccc   V=V+gap_open --------->
            jpV(i,j)=1
            val1=val(i,j-1)+gap_open
            val2=preV(i,j-1)+gap_extn
            if(val1.gt.val2) then
               V=val1
            else
               V=val2
               if(j.gt.1)jpV(i,j)=jpV(i,j-1)+1
            endif
            preH(i,j)=H         !unaccepted H
            preV(i,j)=V         !unaccepted V
            
            if(D.gt.H.and.D.gt.V)then
               idir(i,j)=1
               val(i,j)=D
            elseif(H.gt.V)then
               idir(i,j)=2
               val(i,j)=H
            else
               idir(i,j)=3
               val(i,j)=V
            endif
 222     continue
 111  continue
      score0=val(nseq1,nseq2)   !alignment score
      
c     tracing back the pathway:
      do j=1,nseq2
         j2i(j)=-1              !all are not aligned
      enddo
      i=nseq1
      j=nseq2
      do while(i.gt.0.and.j.gt.0)
         if(idir(i,j).eq.1)then !from diagonal
            j2i(j)=i
            i=i-1
            j=j-1
         elseif(idir(i,j).eq.2)then !from horizonal
            it=jpH(i,j)
            do me=1,it
               if(i.gt.0) then
                  i=i-1                       
               endif
            enddo
         else
            it=jpV(i,j)
            do me=1,it
               if(j.gt.0) then
                  j=j-1                       
               endif
            enddo
         endif
      enddo
      
*^^^^^^^^^^^DP finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      return
      end
      
********************************************************************
*     This is an alternative implementation of Needleman-Wunsch dynamic program 
*     (by Y. Zhang 2005)
*     1. Count two-layer iteration and multiple-gaps
*     2. The gap penality W(k)=Go+Ge*k1+Ge*k2 if gap open on both sequences
*     
*     Input: score(i,j), gap_open, gap_extn
*     Output: j2i(j)
*     idir(i,j)=1,2,3, from diagonal, horizontal, vertical
*     val(i,j) is the cumulative score of (i,j)
********************************************************************
      subroutine DPalt(score0)
      PARAMETER(ndim=6000)
      common/dpc/score(ndim,ndim),gap_open,gap_extn,j2i(ndim)
     &     ,nseq1,nseq2
      
      dimension val(0:ndim,0:ndim),idir(0:ndim,0:ndim)
      dimension preV(0:ndim,0:ndim),preH(0:ndim,0:ndim),
     &     preD(0:ndim,0:ndim)
      dimension idirH(0:ndim,0:ndim),idirV(0:ndim,0:ndim)
      
ccc   initializations --------------->
      val(0,0)=0.0
      do i=1,nseq1
        val(i,0)=0
        idir(i,0)=0
        preD(i,0)=0.0
        preH(i,0)=-1000.0
        preV(i,0)=-1000.0
      enddo
      do j=1,nseq2
        val(0,j)=0
        idir(0,j)=0
        preD(0,j)=0.0
        preH(0,j)=-1000.0
        preV(0,j)=-1000.0
      enddo
      
ccc   DP ------------------------------>
      do 111 j=1,nseq2
         do 222 i=1,nseq1
ccc   preD=VAL(i-1,j-1)+SCORE(i,j)--------------->
            preD(i,j)=val(i-1,j-1)+score(i,j)
ccc   preH: pre-accepted H----------------------->
            D=preD(i-1,j)+gap_open
            H=preH(i-1,j)+gap_extn
            V=preV(i-1,j)+gap_extn
            if(D.gt.H.and.D.gt.V)then
               preH(i,j)=D
               idirH(i-1,j)=1
            elseif(H.gt.V)then
               preH(i,j)=H
               idirH(i-1,j)=2
            else
               preH(i,j)=V
               idirH(i-1,j)=3
            endif
ccc   preV: pre-accepted V----------------------->
            D=preD(i,j-1)+gap_open
            H=preH(i,j-1)+gap_extn
            V=preV(i,j-1)+gap_extn
            if(D.gt.H.and.D.gt.V)then
               preV(i,j)=D
               idirV(i,j-1)=1
            elseif(H.gt.V)then
               preV(i,j)=H
               idirV(i,j-1)=2
            else
               preV(i,j)=V
               idirV(i,j-1)=3
            endif
            
ccc   decide idir(i,j)----------->
            if(preD(i,j).gt.preH(i,j).and.preD(i,j).gt.preV(i,j))then
               idir(i,j)=1
               val(i,j)=preD(i,j)
            elseif(preH(i,j).gt.preV(i,j))then
               idir(i,j)=2
               val(i,j)=preH(i,j)
            else
               idir(i,j)=3
               val(i,j)=preV(i,j)
            endif
 222     continue
 111  continue
      score0=val(nseq1,nseq2)   !alignment score
      
ccc   tracing back the pathway:
      do j=1,nseq2
        j2i(j)=-1              !all are not aligned
      enddo
      i=nseq1
      j=nseq2
      do while(i.gt.0.and.j.gt.0)
         if(idir(i,j).eq.1)then !from diagonal
            j2i(j)=i
            i=i-1
            j=j-1
         elseif(idir(i,j).eq.2)then
            i=i-1
            idir(i,j)=idirH(i,j)
         else
            j=j-1
            idir(i,j)=idirV(i,j)
         endif
      enddo
      
*^^^^^^^^^^^DP finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      return
      end
      
********************************************************************
*     read matrix
*     
      subroutine matrix
      parameter(naa=24) !number of amino acid
      common/matra/imut(naa,naa)     !b,z,x are additional

*     following from PAM30:      
      imut(1,1)=6
      imut(1,2)=-7
      imut(1,3)=-4
      imut(1,4)=-3
      imut(1,5)=-6
      imut(1,6)=-4
      imut(1,7)=-2
      imut(1,8)=-2
      imut(1,9)=-7
      imut(1,10)=-5
      imut(1,11)=-6
      imut(1,12)=-7
      imut(1,13)=-5
      imut(1,14)=-8
      imut(1,15)=-2
      imut(1,16)=0
      imut(1,17)=-1
      imut(1,18)=-13
      imut(1,19)=-8
      imut(1,20)=-2
      imut(1,21)=-3
      imut(1,22)=-3
      imut(1,23)=-1
      imut(1,24)=-17
      imut(2,1)=-7
      imut(2,2)=8
      imut(2,3)=-6
      imut(2,4)=-10
      imut(2,5)=-8
      imut(2,6)=-2
      imut(2,7)=-9
      imut(2,8)=-9
      imut(2,9)=-2
      imut(2,10)=-5
      imut(2,11)=-8
      imut(2,12)=0
      imut(2,13)=-4
      imut(2,14)=-9
      imut(2,15)=-4
      imut(2,16)=-3
      imut(2,17)=-6
      imut(2,18)=-2
      imut(2,19)=-10
      imut(2,20)=-8
      imut(2,21)=-7
      imut(2,22)=-4
      imut(2,23)=-1
      imut(2,24)=-17
      imut(3,1)=-4
      imut(3,2)=-6
      imut(3,3)=8
      imut(3,4)=2
      imut(3,5)=-11
      imut(3,6)=-3
      imut(3,7)=-2
      imut(3,8)=-3
      imut(3,9)=0
      imut(3,10)=-5
      imut(3,11)=-7
      imut(3,12)=-1
      imut(3,13)=-9
      imut(3,14)=-9
      imut(3,15)=-6
      imut(3,16)=0
      imut(3,17)=-2
      imut(3,18)=-8
      imut(3,19)=-4
      imut(3,20)=-8
      imut(3,21)=6
      imut(3,22)=-3
      imut(3,23)=-1
      imut(3,24)=-17
      imut(4,1)=-3
      imut(4,2)=-10
      imut(4,3)=2
      imut(4,4)=8
      imut(4,5)=-14
      imut(4,6)=-2
      imut(4,7)=2
      imut(4,8)=-3
      imut(4,9)=-4
      imut(4,10)=-7
      imut(4,11)=-12
      imut(4,12)=-4
      imut(4,13)=-11
      imut(4,14)=-15
      imut(4,15)=-8
      imut(4,16)=-4
      imut(4,17)=-5
      imut(4,18)=-15
      imut(4,19)=-11
      imut(4,20)=-8
      imut(4,21)=6
      imut(4,22)=1
      imut(4,23)=-1
      imut(4,24)=-17
      imut(5,1)=-6
      imut(5,2)=-8
      imut(5,3)=-11
      imut(5,4)=-14
      imut(5,5)=10
      imut(5,6)=-14
      imut(5,7)=-14
      imut(5,8)=-9
      imut(5,9)=-7
      imut(5,10)=-6
      imut(5,11)=-15
      imut(5,12)=-14
      imut(5,13)=-13
      imut(5,14)=-13
      imut(5,15)=-8
      imut(5,16)=-3
      imut(5,17)=-8
      imut(5,18)=-15
      imut(5,19)=-4
      imut(5,20)=-6
      imut(5,21)=-12
      imut(5,22)=-14
      imut(5,23)=-1
      imut(5,24)=-17
      imut(6,1)=-4
      imut(6,2)=-2
      imut(6,3)=-3
      imut(6,4)=-2
      imut(6,5)=-14
      imut(6,6)=8
      imut(6,7)=1
      imut(6,8)=-7
      imut(6,9)=1
      imut(6,10)=-8
      imut(6,11)=-5
      imut(6,12)=-3
      imut(6,13)=-4
      imut(6,14)=-13
      imut(6,15)=-3
      imut(6,16)=-5
      imut(6,17)=-5
      imut(6,18)=-13
      imut(6,19)=-12
      imut(6,20)=-7
      imut(6,21)=-3
      imut(6,22)=6
      imut(6,23)=-1
      imut(6,24)=-17
      imut(7,1)=-2
      imut(7,2)=-9
      imut(7,3)=-2
      imut(7,4)=2
      imut(7,5)=-14
      imut(7,6)=1
      imut(7,7)=8
      imut(7,8)=-4
      imut(7,9)=-5
      imut(7,10)=-5
      imut(7,11)=-9
      imut(7,12)=-4
      imut(7,13)=-7
      imut(7,14)=-14
      imut(7,15)=-5
      imut(7,16)=-4
      imut(7,17)=-6
      imut(7,18)=-17
      imut(7,19)=-8
      imut(7,20)=-6
      imut(7,21)=1
      imut(7,22)=6
      imut(7,23)=-1
      imut(7,24)=-17
      imut(8,1)=-2
      imut(8,2)=-9
      imut(8,3)=-3
      imut(8,4)=-3
      imut(8,5)=-9
      imut(8,6)=-7
      imut(8,7)=-4
      imut(8,8)=6
      imut(8,9)=-9
      imut(8,10)=-11
      imut(8,11)=-10
      imut(8,12)=-7
      imut(8,13)=-8
      imut(8,14)=-9
      imut(8,15)=-6
      imut(8,16)=-2
      imut(8,17)=-6
      imut(8,18)=-15
      imut(8,19)=-14
      imut(8,20)=-5
      imut(8,21)=-3
      imut(8,22)=-5
      imut(8,23)=-1
      imut(8,24)=-17
      imut(9,1)=-7
      imut(9,2)=-2
      imut(9,3)=0
      imut(9,4)=-4
      imut(9,5)=-7
      imut(9,6)=1
      imut(9,7)=-5
      imut(9,8)=-9
      imut(9,9)=9
      imut(9,10)=-9
      imut(9,11)=-6
      imut(9,12)=-6
      imut(9,13)=-10
      imut(9,14)=-6
      imut(9,15)=-4
      imut(9,16)=-6
      imut(9,17)=-7
      imut(9,18)=-7
      imut(9,19)=-3
      imut(9,20)=-6
      imut(9,21)=-1
      imut(9,22)=-1
      imut(9,23)=-1
      imut(9,24)=-17
      imut(10,1)=-5
      imut(10,2)=-5
      imut(10,3)=-5
      imut(10,4)=-7
      imut(10,5)=-6
      imut(10,6)=-8
      imut(10,7)=-5
      imut(10,8)=-11
      imut(10,9)=-9
      imut(10,10)=8
      imut(10,11)=-1
      imut(10,12)=-6
      imut(10,13)=-1
      imut(10,14)=-2
      imut(10,15)=-8
      imut(10,16)=-7
      imut(10,17)=-2
      imut(10,18)=-14
      imut(10,19)=-6
      imut(10,20)=2
      imut(10,21)=-6
      imut(10,22)=-6
      imut(10,23)=-1
      imut(10,24)=-17
      imut(11,1)=-6
      imut(11,2)=-8
      imut(11,3)=-7
      imut(11,4)=-12
      imut(11,5)=-15
      imut(11,6)=-5
      imut(11,7)=-9
      imut(11,8)=-10
      imut(11,9)=-6
      imut(11,10)=-1
      imut(11,11)=7
      imut(11,12)=-8
      imut(11,13)=1
      imut(11,14)=-3
      imut(11,15)=-7
      imut(11,16)=-8
      imut(11,17)=-7
      imut(11,18)=-6
      imut(11,19)=-7
      imut(11,20)=-2
      imut(11,21)=-9
      imut(11,22)=-7
      imut(11,23)=-1
      imut(11,24)=-17
      imut(12,1)=-7
      imut(12,2)=0
      imut(12,3)=-1
      imut(12,4)=-4
      imut(12,5)=-14
      imut(12,6)=-3
      imut(12,7)=-4
      imut(12,8)=-7
      imut(12,9)=-6
      imut(12,10)=-6
      imut(12,11)=-8
      imut(12,12)=7
      imut(12,13)=-2
      imut(12,14)=-14
      imut(12,15)=-6
      imut(12,16)=-4
      imut(12,17)=-3
      imut(12,18)=-12
      imut(12,19)=-9
      imut(12,20)=-9
      imut(12,21)=-2
      imut(12,22)=-4
      imut(12,23)=-1
      imut(12,24)=-17
      imut(13,1)=-5
      imut(13,2)=-4
      imut(13,3)=-9
      imut(13,4)=-11
      imut(13,5)=-13
      imut(13,6)=-4
      imut(13,7)=-7
      imut(13,8)=-8
      imut(13,9)=-10
      imut(13,10)=-1
      imut(13,11)=1
      imut(13,12)=-2
      imut(13,13)=11
      imut(13,14)=-4
      imut(13,15)=-8
      imut(13,16)=-5
      imut(13,17)=-4
      imut(13,18)=-13
      imut(13,19)=-11
      imut(13,20)=-1
      imut(13,21)=-10
      imut(13,22)=-5
      imut(13,23)=-1
      imut(13,24)=-17
      imut(14,1)=-8
      imut(14,2)=-9
      imut(14,3)=-9
      imut(14,4)=-15
      imut(14,5)=-13
      imut(14,6)=-13
      imut(14,7)=-14
      imut(14,8)=-9
      imut(14,9)=-6
      imut(14,10)=-2
      imut(14,11)=-3
      imut(14,12)=-14
      imut(14,13)=-4
      imut(14,14)=9
      imut(14,15)=-10
      imut(14,16)=-6
      imut(14,17)=-9
      imut(14,18)=-4
      imut(14,19)=2
      imut(14,20)=-8
      imut(14,21)=-10
      imut(14,22)=-13
      imut(14,23)=-1
      imut(14,24)=-17
      imut(15,1)=-2
      imut(15,2)=-4
      imut(15,3)=-6
      imut(15,4)=-8
      imut(15,5)=-8
      imut(15,6)=-3
      imut(15,7)=-5
      imut(15,8)=-6
      imut(15,9)=-4
      imut(15,10)=-8
      imut(15,11)=-7
      imut(15,12)=-6
      imut(15,13)=-8
      imut(15,14)=-10
      imut(15,15)=8
      imut(15,16)=-2
      imut(15,17)=-4
      imut(15,18)=-14
      imut(15,19)=-13
      imut(15,20)=-6
      imut(15,21)=-7
      imut(15,22)=-4
      imut(15,23)=-1
      imut(15,24)=-17
      imut(16,1)=0
      imut(16,2)=-3
      imut(16,3)=0
      imut(16,4)=-4
      imut(16,5)=-3
      imut(16,6)=-5
      imut(16,7)=-4
      imut(16,8)=-2
      imut(16,9)=-6
      imut(16,10)=-7
      imut(16,11)=-8
      imut(16,12)=-4
      imut(16,13)=-5
      imut(16,14)=-6
      imut(16,15)=-2
      imut(16,16)=6
      imut(16,17)=0
      imut(16,18)=-5
      imut(16,19)=-7
      imut(16,20)=-6
      imut(16,21)=-1
      imut(16,22)=-5
      imut(16,23)=-1
      imut(16,24)=-17
      imut(17,1)=-1
      imut(17,2)=-6
      imut(17,3)=-2
      imut(17,4)=-5
      imut(17,5)=-8
      imut(17,6)=-5
      imut(17,7)=-6
      imut(17,8)=-6
      imut(17,9)=-7
      imut(17,10)=-2
      imut(17,11)=-7
      imut(17,12)=-3
      imut(17,13)=-4
      imut(17,14)=-9
      imut(17,15)=-4
      imut(17,16)=0
      imut(17,17)=7
      imut(17,18)=-13
      imut(17,19)=-6
      imut(17,20)=-3
      imut(17,21)=-3
      imut(17,22)=-6
      imut(17,23)=-1
      imut(17,24)=-17
      imut(18,1)=-13
      imut(18,2)=-2
      imut(18,3)=-8
      imut(18,4)=-15
      imut(18,5)=-15
      imut(18,6)=-13
      imut(18,7)=-17
      imut(18,8)=-15
      imut(18,9)=-7
      imut(18,10)=-14
      imut(18,11)=-6
      imut(18,12)=-12
      imut(18,13)=-13
      imut(18,14)=-4
      imut(18,15)=-14
      imut(18,16)=-5
      imut(18,17)=-13
      imut(18,18)=13
      imut(18,19)=-5
      imut(18,20)=-15
      imut(18,21)=-10
      imut(18,22)=-14
      imut(18,23)=-1
      imut(18,24)=-17
      imut(19,1)=-8
      imut(19,2)=-10
      imut(19,3)=-4
      imut(19,4)=-11
      imut(19,5)=-4
      imut(19,6)=-12
      imut(19,7)=-8
      imut(19,8)=-14
      imut(19,9)=-3
      imut(19,10)=-6
      imut(19,11)=-7
      imut(19,12)=-9
      imut(19,13)=-11
      imut(19,14)=2
      imut(19,15)=-13
      imut(19,16)=-7
      imut(19,17)=-6
      imut(19,18)=-5
      imut(19,19)=10
      imut(19,20)=-7
      imut(19,21)=-6
      imut(19,22)=-9
      imut(19,23)=-1
      imut(19,24)=-17
      imut(20,1)=-2
      imut(20,2)=-8
      imut(20,3)=-8
      imut(20,4)=-8
      imut(20,5)=-6
      imut(20,6)=-7
      imut(20,7)=-6
      imut(20,8)=-5
      imut(20,9)=-6
      imut(20,10)=2
      imut(20,11)=-2
      imut(20,12)=-9
      imut(20,13)=-1
      imut(20,14)=-8
      imut(20,15)=-6
      imut(20,16)=-6
      imut(20,17)=-3
      imut(20,18)=-15
      imut(20,19)=-7
      imut(20,20)=7
      imut(20,21)=-8
      imut(20,22)=-6
      imut(20,23)=-1
      imut(20,24)=-17
      imut(21,1)=-3
      imut(21,2)=-7
      imut(21,3)=6
      imut(21,4)=6
      imut(21,5)=-12
      imut(21,6)=-3
      imut(21,7)=1
      imut(21,8)=-3
      imut(21,9)=-1
      imut(21,10)=-6
      imut(21,11)=-9
      imut(21,12)=-2
      imut(21,13)=-10
      imut(21,14)=-10
      imut(21,15)=-7
      imut(21,16)=-1
      imut(21,17)=-3
      imut(21,18)=-10
      imut(21,19)=-6
      imut(21,20)=-8
      imut(21,21)=6
      imut(21,22)=0
      imut(21,23)=-1
      imut(21,24)=-17
      imut(22,1)=-3
      imut(22,2)=-4
      imut(22,3)=-3
      imut(22,4)=1
      imut(22,5)=-14
      imut(22,6)=6
      imut(22,7)=6
      imut(22,8)=-5
      imut(22,9)=-1
      imut(22,10)=-6
      imut(22,11)=-7
      imut(22,12)=-4
      imut(22,13)=-5
      imut(22,14)=-13
      imut(22,15)=-4
      imut(22,16)=-5
      imut(22,17)=-6
      imut(22,18)=-14
      imut(22,19)=-9
      imut(22,20)=-6
      imut(22,21)=0
      imut(22,22)=6
      imut(22,23)=-1
      imut(22,24)=-17
      imut(23,1)=-1
      imut(23,2)=-1
      imut(23,3)=-1
      imut(23,4)=-1
      imut(23,5)=-1
      imut(23,6)=-1
      imut(23,7)=-1
      imut(23,8)=-1
      imut(23,9)=-1
      imut(23,10)=-1
      imut(23,11)=-1
      imut(23,12)=-1
      imut(23,13)=-1
      imut(23,14)=-1
      imut(23,15)=-1
      imut(23,16)=-1
      imut(23,17)=-1
      imut(23,18)=-1
      imut(23,19)=-1
      imut(23,20)=-1
      imut(23,21)=-1
      imut(23,22)=-1
      imut(23,23)=-1
      imut(23,24)=-17
      imut(24,1)=-17
      imut(24,2)=-17
      imut(24,3)=-17
      imut(24,4)=-17
      imut(24,5)=-17
      imut(24,6)=-17
      imut(24,7)=-17
      imut(24,8)=-17
      imut(24,9)=-17
      imut(24,10)=-17
      imut(24,11)=-17
      imut(24,12)=-17
      imut(24,13)=-17
      imut(24,14)=-17
      imut(24,15)=-17
      imut(24,16)=-17
      imut(24,17)=-17
      imut(24,18)=-17
      imut(24,19)=-17
      imut(24,20)=-17
      imut(24,21)=-17
      imut(24,22)=-17
      imut(24,23)=-17
      imut(24,24)=1
      
      return
      end
      
      function upper(A)
      CHARACTER A,upper
      IF(A.LE.'z'.and.A.GE.'a')then
         A=CHAR(ICHAR(A)-32)
      endif
      upper=A
      RETURN
      END
