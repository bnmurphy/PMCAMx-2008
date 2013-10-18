      subroutine ahoprep
c
c-----CAMx v4.02 030709
c
c     AHOPREP reads the header records of the albedo/haze/ozone file and
c     reads the time-invariant albedo codes
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
c     Routines called:
c        none
c
c     Called by:
c        CHMPREP
c
      include 'camx.prm'
      include 'ahomap.com'
      include 'filunit.com'
      include 'grid.com'
c
      dimension rdalb(NALB),rdhaz(NHAZE),rdozn(NOZN)
      character*10 title,name(3)
      logical lerror
c
      data name/'ALBEDO    ','HAZE      ','OZONE COL '/
c
c-----Entry point
c
      do i=1,ngrid
         lrdalb(i) = .FALSE.
      enddo
      lerror = .false.
      read(iaho,*)
c
c-----Read albedo class and check inputs
c
      read(iaho,'(a10,5f10.0)') title,(rdalb(i),i=1,NALB)
      if (title.ne.name(1)) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Expecting keyword: ',name(1)
        write(iout,*) 'Read from Albedo/haze/ozone file: ',title
        write(iout,*) 'Perhaps your file order in CAMx.in is incorrect.'
        write(iout,*) 'Check the .out file to make sure.'
        call camxerr()
      endif
      do i = 1,NALB
        if (rdalb(i).ne.albcl(i)) lerror = .true.
      enddo
      if (lerror) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Mismatch in albedo'
        write(iout,*) 'Photolysis rates  file: ',(albcl(i),i=1,NALB)
        write(iout,*) 'Albedo/haze/ozone file: ',(rdalb(i),i=1,NALB)
        call camxerr()
      endif
c
c-----Read haze class and check inputs
c
      read(iaho,'(a10,3f10.0)') title,(rdhaz(i),i=1,NHAZE)
      if (title.ne.name(2)) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Expecting keyword: ',name(2)
        write(iout,*) 'Read from Albedo/haze/ozone file: ',title
        call camxerr()
      endif
      do i = 1,NHAZE
        if (rdhaz(i).ne.hazcl(i)) lerror = .true.
      enddo
      if (lerror) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Mismatch in haze class'
        write(iout,*) 'Photolysis rates  file: ',(hazcl(i),i=1,NHAZE)
        write(iout,*) 'Albedo/haze/ozone file: ',(rdhaz(i),i=1,NHAZE)
        call camxerr()
      endif
c
c-----Read ozone class and check inputs
c
      read(iaho,'(a10,5f10.0)') title,(rdozn(i),i=1,NOZN)
      if (title.ne.name(3)) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Expecting keyword: ',name(3)
        write(iout,*) 'Read from Albedo/haze/ozone file: ',title
        call camxerr()
      endif
      do i = 1,NOZN
        if (rdozn(i).ne.ozcl(i)) lerror = .true.
      enddo
      if (lerror) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Mismatch in ozone class'
        write(iout,*) 'Photolysis rates  file: ',(ozcl(i),i=1,NOZN)
        write(iout,*) 'Albedo/haze/ozone file: ',(rdozn(i),i=1,NOZN)
        call camxerr()
      endif
c
c-----Read coarse grid albedo index
c
      read(iaho,'(a10,3i10)') title,idd,nxalb,nyalb
      if (title.ne.name(1)) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Expecting keyword: ',name(1)
        write(iout,*) 'Read from Albedo/haze/ozone file: ',title
        call camxerr()
      endif
      if (nxalb.ne.ncol(1) .or. nyalb.ne.nrow(1)) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Mismatch in number of coarse cells'
        write(iout,*) 'NXALB,NYALB,NCOL,NROW: ',
     &                 nxalb,nyalb,ncol(1),nrow(1)
        call camxerr()
      endif
      do j = nyalb,1,-1
        read(iaho,'(9999i1)') (icdalb(i+(j-1)*nxalb),i=1,nxalb)
      enddo
      lrdalb(1) = .TRUE.
c 
c-----Read fine grid albedo index (if any)
c              
 100  read(iaho,'(a10,3i10)') title,iff,nxalb,nyalb 
      if (title.ne.name(1)) then 
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Expecting keyword: ',name(1)
        write(iout,*) 'Read from Albedo/haze/ozone file: ',title
        call camxerr()
      endif 
      if (iff.eq.0) goto 999
      if (iff.gt.nnest) then
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Albedo codes are given for more fine grids'
        write(iout,*) 'than have been specified in the control file'
        write(iout,*) 'IFF, NNEST ',iff,nnest
        call camxerr()
      endif
      if (nxalb.ne.ncol(iff+1) .or. nyalb.ne.nrow(iff+1)) then 
        write(iout,'(//,a)') 'ERROR in AHOPREP:'
        write(iout,*) 'Mismatch in number of fine cells' 
        write(iout,*) 'IFINE,NXALB,NYALB,NCOL,NROW: ',iff,
     &                 nxalb,nyalb,ncol(iff+1),nrow(iff+1) 
        call camxerr()
      endif 
      do j = nyalb,1,-1   
        read(iaho,'(9999i1)') (icdalb(iptr2d(iff+1)-1+i+(j-1)*nxalb),
     &                                                      i=1,nxalb)
      enddo 
      lrdalb(iff+1) = .TRUE.
      goto 100
c
 999  return
      end
