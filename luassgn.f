      subroutine luassgn(ncol,nrow,nlu,io,jo,nmesh,ncolf,nrowf,fsurf,
     &                   fsurff)
c
c-----CAMx v4.02 030709
c
c     LUASSGN assigns land use value for fine grid from that of its 
c     parent grid
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modfications:
c        none
c
c     Input arguments:
c        ncol                number of columns
c        nrow                number of rows
c        nlu                 number of landuse categories
c        io                  starting i index for the fine grid
c        jo                  starting j index for the fine grid
c        nmesh               mesh number
c        ncolf               number of columns in fine grid
c        nrowf               number of rows in fine grid
c        fsurf               cell centered variable on coarse grid
c
c     Output arguments:
c        fsurff              cell centered variable on coarse grid
c
c     Routines called:
c        none
c
c     Called by:
c        STARTUP
c
      dimension fsurf(ncol,nrow,nlu),fsurff(ncolf,nrowf,nlu)
c
      do 50 l = 1,nlu
        do 40 jfin = 2,nrowf-1
          j = (jfin - 2)/nmesh + jo
          do 30 ifin = 2,ncolf-1
            i = (ifin - 2)/nmesh + io
            fsurff(ifin,jfin,l) = fsurf(i,j,l)
  30      continue
  40    continue
  50  continue
c
      return
      end
