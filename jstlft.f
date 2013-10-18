      subroutine jstlft( string )
c
c-----CAMx v4.02 030709
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c
c   Description:
c     Left justifies a string
c
c   Arguments:
c     Inputs/Outputs: (the string arguments serves as both input and
c                      output)
c       string   C   string to left justify
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     08/10/98  -gmw-  original development
c
c-----------------------------------------------------------------------
c   Argument declaration:
c-----------------------------------------------------------------------
c
      character*(*) string
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      integer ibeg, i
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c   ---- it may already be left-justified ---
c
      if( string(1:1) .NE. ' ' ) goto 9999
c
c   ---- find the first non-blank character ---
c
      do 10 i=1,LEN( string )
         if( string(i:i) .NE. ' ' ) then
             ibeg = i
             goto 111
         endif
   10 continue
c
c   --- no non-blanks found, it's a blank string, nothing to do ----
c
      goto 9999
c
c   ---- move the string over, 2 char at a time ---
c
  111 continue
      do 20 i=1,LEN( string )-ibeg+1
         string(i:i) = string(i+ibeg-1:i+ibeg-1)
         string(i+ibeg-1:i+ibeg-1) = ' '
   20 continue
      goto 9999
c
c-----------------------------------------------------------------------
c   Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
