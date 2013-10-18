      subroutine chmprep
c
c-----CAMx v4.02 030709
c
c     CHMPREP reads the chemistry parameter file, the photolysis rates
c     lookup table, and the albedo/haze/ozone file
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
c     Routines Called:
c        READCHM, READPHT, AHOPREP
c
c     Called by:
c        STARTUP
c
      include 'camx.prm'
      include 'flags.com'
      include 'tracer.com'
c
c
c-----Entry point
c
c-----Read chemistry parameters
c
      call readchm
c
c-----Read photolysis rates lookup table
c
      if (lchem) call readpht()
c             
c-----Read albedo/haze/ozone file header
c             
      if (lchem) call ahoprep()
c 
      return
      end
