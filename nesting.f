      subroutine nesting
c
c-----CAMx v4.02 030709
c
c     NESTING is the driver for grid nesting algorithm.  It does the 
c     following tasks in a recursive order:
c       1. determines boundary conditions for children grids
c       2. calls EMISTRNS for each grid
c       3. calls CHEMRXN for each grid 
c       4. aggregates concentrations on children grids to parent grid
c     Chemistry and transport are performed for each grid on their own
c     time step. Up to 4 generations of grid nesting are currently allowed.
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        none
c
c     Input arguments:
c        none
c
c     Output arguments:
c        none
c
c     Subroutines Called:
c        SETBC
c        EMISTRNS
c        CHEMRXN
c        AGGR00
c        FGAVRG
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "grid.com"
c
      dimension iparnt(20),igrd0(20)
c
c-----Entry point
c
c-----Computation for children grids
c
      do igen=1,20
        iparnt(igen)=0
        igrd0(igen)=0
      enddo
c
c-----Generation 2
c
      igen=1
      iparnt(igen+1)=1
      call setbc(iparnt(igen+1))
      mch2=nchdrn(iparnt(igen+1))
      do 100 ic2=1,mch2
        igen=2
        igrd0(igen)=idchdrn(ic2,iparnt(igen))
c
c-----Perform emissions and transport for generation 2
c
        do 99 it2=1,ntim(igrd0(igen))
          call fgavrg(igrd0(igen))
          call emistrns(igrd0(igen))
c
c-----Generation 3
c
          iparnt(igen+1)=igrd0(igen)
          call setbc(iparnt(igen+1))
          mch3=nchdrn(iparnt(igen+1))
          do 90 ic3=1,mch3
            igen=3
            igrd0(igen)=idchdrn(ic3,iparnt(igen))
c
c-----Perform emissions and transport for generation 3
c
            do 89 it3=1,ntim(igrd0(igen))
              call fgavrg(igrd0(igen))
              call emistrns(igrd0(igen))
c
c-----Generation 4
c
              iparnt(igen+1)=igrd0(igen)
              call setbc(iparnt(igen+1))
              mch4=nchdrn(iparnt(igen+1))
              do 80 ic4=1,mch4
                igen=4
                igrd0(igen)=idchdrn(ic4,iparnt(igen))
c
c-----Perform emissions and transport for generation 4
c
                do 79 it4=1,ntim(igrd0(igen))
                  call fgavrg(igrd0(igen))
                  call emistrns(igrd0(igen))
c
c-----Generation X: more generations would be added here
c
c-----Perform chemistry for generation 4
c
                  igen=4
                  call chemrxn(igrd0(igen))
                  call fgavrg(igrd0(igen))
  79            continue
                icode = 3
                if (mch4.eq.1) then
                  icode = 0
                elseif (ic4.eq.1) then
                  icode = 1
                elseif (ic4.eq.mch4) then
                  icode = 2
                endif
                call aggr00(igrd0(igen),iparnt(igen),icode)
  80          continue
c
c-----Perform chemistry for generation 3
c
              igen=3
              call chemrxn(igrd0(igen))
              call fgavrg(igrd0(igen))
  89        continue
            icode = 3
            if (mch3.eq.1) then
              icode = 0
            elseif (ic3.eq.1) then
              icode = 1
            elseif (ic3.eq.mch3) then
              icode = 2
            endif
            call aggr00(igrd0(igen),iparnt(igen),icode)
  90      continue
c
c-----Perform chemistry for generation 4
c
          igen=2
          call chemrxn(igrd0(igen))
          call fgavrg(igrd0(igen))
  99    continue
        icode = 3
        if (mch2.eq.1) then
          icode = 0
        elseif (ic2.eq.1) then
          icode = 1
        elseif (ic2.eq.mch2) then
          icode = 2
        endif
        call aggr00(igrd0(igen),iparnt(igen),icode)
 100  continue
c
      return
      end
