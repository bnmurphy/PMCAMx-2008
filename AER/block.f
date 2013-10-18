c=========================================================================+
c       03/09/03: bkoo
c                 - commented out parameter nn2
c                   (not used; nsoap is undefined here)
c       June 8 2000: bkoo                                                 +
c                    - add default dt & pres                              +
c       MAY 2000: modified by bkoo                                        +
c                 - aerm & ver                                            +
c=========================================================================+
	block data properties

c       This subroutine contains physical properties for the new generation
c       multicomponent aerosol dynamics model (MADM)
c *************************************************************************
c               WRITTEN BY DR. CHRISTODOULOS PILINIS
c                          December, 1996
c *************************************************************************
c

	include 'dynamic.inc'
	parameter (nn1=norg+ninert)
cbk	parameter (nn2=ngas-nsoap) ! bkoo (03/09/03)
	parameter (nn3=nsec*nsp)

	data dsec	/nsecp1*0.0d0/
	data ps		/nn3*0.0d0/
	data qt		/nsec  *0.0d0/
	data qsource	/naer  *0.0d0/
	data c		/nsxint*0.0d0/
	
	data dsecf	/nsecp1*0.0d0/
	data hi		/nsxgas*0.0d0/
	data qtt	/nsp   *0.0d0/
	data vel	/nsp  *0.0d0/
	data lamda	/nsp  *0.0d0/
	data qn		/nsec  *0.0d0/
	data gsource	/ngas  *0.0d0/
	data qtot0	/nexti *0.0d0/
c	
c	molecular weights of external inorganics and norg.
c       sulfates as H2SO4, nitrates as HNO3, ammonium as NH3 and Cl as HCl
c	dummy organics set to 100, corrected later
c       
      data emw/18.0d0,23.0d0,98.0d0,63.0d0,17.0d0,36.5d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               180.0d0, 180.0d0, 180.0d0, 180.0d0,
     $               150.0d0, 150.0d0, 150.0d0, 150.0d0,
     $               220.0d0, 100.0d0, 100.0d0 / ! tmg (04/15/02)

      data gmw/17.0d0,63.0d0,98.0d0,36.5d0,
c     $         nn2*100.0d0 /
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               250.0d0, 250.0d0, 250.0d0, 250.0d0,
     $               180.0d0, 180.0d0, 180.0d0, 180.0d0,
     $               150.0d0, 150.0d0, 150.0d0, 150.0d0 / ! tmg (04/15/02)

      data intmw/18.0d0,  1.0d0, 23.0d0, 18.0d0, 35.5d0,96.0d0,97.0d0
     $          ,62.0d0, 58.5d0,142.0d0, 85.0d0,132.0d0,80.0d0,53.5d0
     $          ,98.0d0,115.0d0,120.0d0,247.0d0, 17.0d0/ ! bkoo (02/14/02)

      data pi    / 3.14159265358979d0 /
      data rgas  / 8.3144d0           /	! gas constant in SI units
      data diffus/ nsp*1.18e-5       /	! diffusion coefficients (m**2 sec-1)
      data delta / nsp*1.0d-1        /	! accomodation coefficient
      data tinys / 1.0d-9             /	! minimumum non-zero concentration for
c					   diffun and diffund subroutines
      data aerm /'EQUI'/		! default aerosol module
c      data ver /'v1.0      (05/20/00)'/	! version info.
c      data ver /'v1.1      (06/08/00)'/	! version info.
      data ver /'v1.1a     (06/09/00)'/	! version info.
c
      data dt / 600.0d0 /		! default timestep (sec)
      data pres / 1.0d0 /		! default pressure (ATM)
      data rh / 0.9d0 /		! default relative humidity (0-1)
      data temp / 298.d0 /		! default temperature (K)

      data ims / nsec*0 /               ! deliquescent aerosol - bkoo (03/05/02)

	end
