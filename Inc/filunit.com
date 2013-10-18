c-----CAMx v4.02 030709
c 
c     FILUNIT.COM contains all model I/O unit numbers
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        8/30/02    Cloud file now contains rain record, rain unit removed.
c                   Cloud/rain and water vapor files can be read in for all
c                   nests
c        1/10/03    Added deposition file units
c
c-----------------------------------------------------------------------
c     File units for ASCII and Fortran binary I/O files:
c
c     iout   -- message output file
c     idiag  -- diagnostic output file
c     imass  -- mass summary output file
c     iconc  -- coarse grid instantaneous concentration output file
c     ifconc -- fine grid instantaneous concentration output file
c     iavg   -- coarse grid average concentration output file
c     ifavg  -- fine grid average concentration output file
c     idep   -- coarse grid deposition output file
c     ifdep  -- fine grid deposition output file
c     ipig   -- PiG output file
c     ichem  -- chemistry parameters input file
c     iphot  -- photolysis lookup input file
c     ih2o   -- water vapor concentration input file
c     icld   -- cloud/rain input file
c     iic    -- initial conditions input file
c     ibc    -- boundary conditions input file
c     itopc  -- top conditions input file
c     iaho   -- albedo/haze/ozone input file
c     iptem  -- point source input file
c     ihtp   -- layer height/pressure input file
c     iarem  -- area emission input file
c     isurf  -- landuse input file
c     iwind  -- wind input file
c     itemp  -- temperature input file
c     ikv    -- vertical diffusivity input file
c     irstc  -- coarse grid restart input file
c     irstf  -- fine grid restart input file
c     irstp  -- PiG restart file
c-----------------------------------------------------------------------

      common /funit/ iout,idiag,imass,iconc(2),ifconc(2),iavg,ifavg,
     &               idep,ifdep,ipig,ichem,iphot,ih2o(MXGRID),
     &               icld(MXGRID),iic,ibc,
     &               itopc,iaho,iptem,ihtp(MXGRID),iarem(3,MXGRID),
     &               isurf(MXGRID),iwind(MXGRID),itemp(MXGRID),
     &               ikv(MXGRID),irstc,irstf,irstp
c
c========================= Process Analysis Begin ==============================
c
c     ipr_unit  -- unit number of Integrated Process Rates output file
c     irr_unit  -- unit number of Integrated Reaction Rates output file
c
       common /paunit/ ipr_unit, irr_unit
c
c========================= Process Analysis End ==============================
c
c
c-----------------------------------------------------------------------
c     File units for netCDF files:
c
c     ncinst  -- instantaneous concentration output file
c     ncavg   -- average concentration output file
c     ncsurf  -- landuse input file
c     ncic    -- initial conditions input file
c     ncbc    -- boundary conditions input file
c     ncpig   -- PiG output file
c     ncpnt   -- point source input file
c     ncbio   -- biogenic area source input file
c     ncgenrl -- general area source input file
c     ncmobil -- mobile area source input file
c     ncmet   -- meteorological input file
c     ncrstc  -- restart concentration input file
c     ncrstp  -- restart PiG input file
c-----------------------------------------------------------------------
c
      common /funitnc/ ncinst,ncavg,ncsurf,ncic,ncbc,ncpig,ncpnt,
     &                 ncbio,ncgenrl,ncmobil,ncmet,ncrstc,ncrstp
