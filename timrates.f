      subroutine timrates(ncol,nrow,nlay,fnow,fnxt,pfpt)
c  
c-----CAMx v4.02 030709
c 
c     TIMRATES calculates local (Eulerian) time-rate of change of the 
c     input field for a given grid.
c  
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c  
c     Modifications:   
c        none
c   
c     Input arguments:   
c        ncol                number of columns  
c        nrow                number of rows  
c        nlay                number of layers  
c        fnow                current parameter field 
c        fnxt                future parameter field
c              
c     Output arguments:   
c        pfpt                time-rate change of parameter field
c              
c     Routines Called:   
c        none  
c              
c     Called by:   
c        READINP
c        INTRPDAT
c
      include 'camx.prm'
      include 'camx.com'
      include 'bndary.com'
c
      real fnow(ncol,nrow,nlay),fnxt(ncol,nrow,nlay),
     &     pfpt(ncol,nrow,nlay)
c
c-----Entry point
c 
      do 10 k = 1,nlay 
        do 15 j = 1,nrow 
          do 20 i = 1,ncol 
            pfpt(i,j,k) = (fnxt(i,j,k) - fnow(i,j,k))/(60.*dtinp)
 20       continue 
 15     continue 
 10   continue 
c
      return
      end
