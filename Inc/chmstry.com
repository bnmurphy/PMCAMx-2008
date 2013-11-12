c-----CAMx v4.02 030709
c  
c     CHMSTRY.COM contains all chemistry variables 
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c       4/4/00     Added aerosol deposition variables to /aerochm/
c       1/9/02     Aerosol size cut points and density now a function of
c                  species
c       8/20/02    Added minimum CWC to define presence of clouds
c      12/12/02    Expanded species list for Mechanism 4
c       1/10/03    Added array for deposition output species names
c 
c-----------------------------------------------------------------------
c     Parameters for some of the switches:
c
c     CDCMC  -- string for requesting the CMC (standard) chemistry solver
c     CDIEH  -- string for requesting the IEH chemistry solver
c     IDCMC  -- code for using the CMC (standard) chemistry solver
c     IDIEH  -- code for using the IEH chemistry solver
c-----------------------------------------------------------------------
c
      character*10 CDCMC
      character*10 CDIEH
      integer      IDCMC
      integer      IDIEH
c
      parameter( CDCMC = "CMC       " )
      parameter( CDIEH = "IEH       " )
      parameter( IDCMC = 1 )
      parameter( IDIEH = 2 )
c 
c-----------------------------------------------------------------------
c    Variables for the number of species in input files:
c
c    ngas   --  number of gas species being modeled
c    naero  --  number of aersol species being modeled
c    nspec  --  total number of modeled species
c    nrad   --  number of radical species being modeled
c    nreact --  number of chemical reactions
c    nspfst --  number of "fast" species -- handled by the fast solver
c    iessrad--  number of radicals in steady state for IEH solver
c    idmech --  the code which determines which chemical mechanism is used
c    idsolv --  the code which determines which chemstry solver to use
c    navspc --  number of species to write to output average file
c    nicspc --  number of species in the initial conditions file
c    nbcspc --  number of species in the boundary conditions file
c    nptspc --  number of species in the point source emissions file
c    narspc --  number of species in the surface emissions file
c-----------------------------------------------------------------------
c
       integer   ngas
       integer   naero
       integer   nspec
       integer   nrad
       integer   nreact
       integer   nspfst
       integer   iessrad
       integer   idmech
       integer   idsolv
       integer   navspc
       integer   nicspc
       integer   nbcspc
       integer   nptspc
       integer   narspc
c
      common /chm1/ ngas, naero, nspec, nrad, nreact, nspfst, iessrad,
     &              idmech, idsolv, navspc, nicspc, nbcspc, nptspc, 
     &              narspc
c
c-----------------------------------------------------------------------
c     Variables for keeping track of where chmistry is being performed:
c     NOTE:  Used for diagnostic and error messages.
c
c     igrdchm  --  grid number of current chemistry step
c     ichm     --  column for the current chemistry step
c     jchm     --  row for the current chemistry step
c     kchm     --  layer for the current chemistry step
c-----------------------------------------------------------------------
c
      integer   igrdchm
      integer   ichm
      integer   jchm
      integer   kchm
c
      common /ijkgrd/ igrdchm, ichm, jchm, kchm
c$omp threadprivate(/ijkgrd/)
c
c-----------------------------------------------------------------------
c     Variables for storing chemical reaction data:
c
c     rk     -- reaction rate constant (ppm/hr)
c     ltdep  -- flag to determine if rate constant is temperature dependent
c     lpdep  -- flag to determine if rate constant is pressure dependent
c     bdnl   -- lower vound value for each modeled species (ppm)
c     bdlrad -- lower bound value for each radical species (ppm)
c-----------------------------------------------------------------------
c
      real    rk(MXRXN)
      logical ltdep(MXRXN)
      logical lpdep(MXRXN)
      real    bdnl(MXSPEC)
      real    bdlrad
      real    nflag  ! turns hno3 production off at night (tmg,06/19/04)
c
      common /chmratep/ rk, nflag
c$omp threadprivate(/chmratep/)
      common /chmrate/ ltdep, lpdep, bdnl, bdlrad
c
c-----------------------------------------------------------------------
c     Variables for photolysis data:
c
c     nphot1   -- number of primary photolysis reactions
c     nphot2   -- number of secondary (scaled) photolysis reactions
c     idphot1  -- ID of primary photolysis reactions
c     idphot2  -- ID of secondary (scaled) photolysis reactions 
c     idphot3  -- ID of primary photolysis reaction to scale to obtain
c                 the secondary photolysis reaction
c     phtscl   -- photolysis reaction scaling factor
c-----------------------------------------------------------------------
c
      integer   nphot1
      integer   nphot2
      integer   idphot1(MXPHT1)
      integer   idphot2(MXPHT2)
      integer   idphot3(MXPHT2)
      real      phtscl(MXPHT2)
c
      common /photmap/ nphot1, nphot2, idphot1, idphot2, idphot3, phtscl
c 
c-----------------------------------------------------------------------
c     Variables for species names:
c
c     spname  --  name of each modeled species
c     spavg   --  name of each species to be written to the output file
c     nmrad   --  name of each radical species
c     depsp   --  name of each deposition species output to file
c-----------------------------------------------------------------------
c
      character*10 spname(MXSPEC+1)
      character*10 spavg(MXSPEC)
      character*10 nmrad(MXRADCL+1)
      character*10 depsp(4*MXSPEC)
c
      common /cname/ spname, spavg, nmrad, depsp
c 
c-----------------------------------------------------------------------
c     Variables for mapping input species to internal model order:
c
c     krad     -- mapping of radical species to specific mechanism order
c     kmap     -- mapping of species on chemistry parameters file to
c                 internal order
c     lbcmap   -- mapping of species in the boundary condition file
c     lavmap   -- mapping of species written to average file
c     lptmap   -- mapping of species in the point source emissions file
c     lptrdmap -- mapping of species in the NetCDF point source file
c     larmap   -- mapping of species in the surface emissions file
c     licmap   -- mapping of species in the initial conditions file
c     lbcmapn  -- mapping of species on north edge of NetCDF boundary file
c     lbcmaps  -- mapping of species on south edge of NetCDF boundary file
c     lbcmape  -- mapping of species on east edge of NetCDF boundary file
c     lbcmapw  -- mapping of species on west edge of NetCDF boundary file
c     lgenrmap -- mapping of species in NetCDF general area emission file
c     lbiomap  -- mapping of species in NetCDF biogenic area emission file
c     lmoblmap -- mapping of species in NetCDF mobile area emission file
c     lspmap   -- mapping of species in NetCDF instant concentration file
c     lavwrmap -- mapping of species in NetCDF average concentration file
c-----------------------------------------------------------------------
c
      integer   krad(NRADNM)
      integer   kmap(NSPNAM)
      integer   lbcmap(MXSPEC)
      integer   lavmap(MXSPEC)
      integer   lptmap(MXSPEC)
      integer   lptrdmap(MXSPEC)
      integer   larmap(MXSPEC,MXGRID) 
      integer   licmap(MXSPEC,MXGRID) 
c
      common /kname/ krad, kmap, lbcmap,
     &               lavmap, lptmap, lptrdmap,
     &               larmap, licmap,
     &               lbcmapn(MXSPEC),lbcmaps(MXSPEC),lbcmape(MXSPEC),
     &               lbcmapw(MXSPEC),lgenrmap(MXSPEC,MXGRID),
     &               lbiomap(MXSPEC,MXGRID),lmoblmap(MXSPEC,MXGRID),
     &               lspmap(MXSPEC,MXGRID),lavwrmap(MXSPEC,MXGRID)
c
c  BNM 07/21/09 CHANGED MANY OF THESE ARRAYS TO ADD NTSOA
      integer   kno   ,kno2  ,ko3  
      integer   kpan  ,kcres ,kpan2
      integer   kmpan ,kpbzn ,knphe
      integer   krno3 ,kdcb2 ,kdcb3
      integer   khno4 ,kacet ,kald2
      integer   kalk1 ,kalk2 ,kalk3
      integer   kalk4 ,kalk5 ,karo1
      integer   karo2 ,kbacl ,kbald
      integer   kbcl1 ,kbcl2 ,kbuta
      integer   kccho ,kccrs  
      integer   kcl2  ,kco   ,kco2h
      integer   kco3h ,kcooh ,kcprm
      integer   kdcb1 ,keth  ,kethe
      integer   ketoh ,kfcrs ,kfmcl
      integer   kform ,kfprm ,kgly 
      integer   kh2o2 ,khc2h ,khcho
      integer   khcl  ,khono ,khno3
      integer   kho2h ,khocl ,kicl1
      integer   kicl2 ,kisop ,kispd
      integer   kmek  ,kmeoh ,kmeth
      integer   kmgly ,kmvk  ,kna  
      integer   knh3  ,kntr  ,knxoy
      integer   kole  ,kole1 ,kole2
      integer   kbpin ,klimo ,kmono ,ksesq
      integer   kopen ,kpar  ,kpcl 
      integer   kpec  ,kphen ,kpna 
      integer   kpnh4 ,kpno3 ,kpoa 
      integer   kprod ,kpso4 ,krc2h
      integer   krc3h ,krcho ,krooh
      integer   kso2
      integer   ksulf
      integer   kterp ,ktol  ,kxn  ,kxyl 
      integer   kcpo1, kcpo2, kcpo3, kcpo4
      integer   kcpo5, kcpo6, kcpo7, kcpo8
      integer   kcoo1, kcoo2, kcoo3, kcoo4
      integer   kcoo5, kcoo6, kcoo7, kcoo8
      integer   kcbs1, kcbs2, kcbs3, kcbs4, kcbs5
      integer   kcas1, kcas2, kcas3, kcas4, kcas5
      integer   kcns1, kcns2, kcns3, kcns4
      integer   kcns5, kcns6, kcns7, kcns8
c
      integer   kapo1, kapo2, kapo3, kapo4
      integer   kapo5, kapo6, kapo7, kapo8
      integer   kaoo1, kaoo2, kaoo3, kaoo4
      integer   kaoo5, kaoo6, kaoo7, kaoo8
      integer   kabs1, kabs2, kabs3, kabs4, kabs5
      integer   kaas1, kaas2, kaas3, kaas4, kaas5
      integer   kans1, kans2, kans3, kans4
      integer   kans5, kans6, kans7, kans8
c
      integer   kapo1_1  ,kapo1_2  ,kapo1_3  ,kapo1_4  ,kapo1_5
      integer   kapo1_6  ,kapo1_7  ,kapo1_8  ,kapo1_9  ,kapo1_10
      integer   kapo2_1  ,kapo2_2  ,kapo2_3  ,kapo2_4  ,kapo2_5
      integer   kapo2_6  ,kapo2_7  ,kapo2_8  ,kapo2_9  ,kapo2_10
      integer   kapo3_1  ,kapo3_2  ,kapo3_3  ,kapo3_4  ,kapo3_5
      integer   kapo3_6  ,kapo3_7  ,kapo3_8  ,kapo3_9  ,kapo3_10
      integer   kapo4_1  ,kapo4_2  ,kapo4_3  ,kapo4_4  ,kapo4_5
      integer   kapo4_6  ,kapo4_7  ,kapo4_8  ,kapo4_9  ,kapo4_10
      integer   kapo5_1  ,kapo5_2  ,kapo5_3  ,kapo5_4  ,kapo5_5
      integer   kapo5_6  ,kapo5_7  ,kapo5_8  ,kapo5_9  ,kapo5_10
      integer   kapo6_1  ,kapo6_2  ,kapo6_3  ,kapo6_4  ,kapo6_5
      integer   kapo6_6  ,kapo6_7  ,kapo6_8  ,kapo6_9  ,kapo6_10
      integer   kapo7_1  ,kapo7_2  ,kapo7_3  ,kapo7_4  ,kapo7_5
      integer   kapo7_6  ,kapo7_7  ,kapo7_8  ,kapo7_9  ,kapo7_10
      integer   kapo8_1  ,kapo8_2  ,kapo8_3  ,kapo8_4  ,kapo8_5
      integer   kapo8_6  ,kapo8_7  ,kapo8_8  ,kapo8_9  ,kapo8_10

      integer   kaoo1_1  ,kaoo1_2  ,kaoo1_3  ,kaoo1_4  ,kaoo1_5
      integer   kaoo1_6  ,kaoo1_7  ,kaoo1_8  ,kaoo1_9  ,kaoo1_10
      integer   kaoo2_1  ,kaoo2_2  ,kaoo2_3  ,kaoo2_4  ,kaoo2_5
      integer   kaoo2_6  ,kaoo2_7  ,kaoo2_8  ,kaoo2_9  ,kaoo2_10
      integer   kaoo3_1  ,kaoo3_2  ,kaoo3_3  ,kaoo3_4  ,kaoo3_5
      integer   kaoo3_6  ,kaoo3_7  ,kaoo3_8  ,kaoo3_9  ,kaoo3_10
      integer   kaoo4_1  ,kaoo4_2  ,kaoo4_3  ,kaoo4_4  ,kaoo4_5
      integer   kaoo4_6  ,kaoo4_7  ,kaoo4_8  ,kaoo4_9  ,kaoo4_10
      integer   kaoo5_1  ,kaoo5_2  ,kaoo5_3  ,kaoo5_4  ,kaoo5_5
      integer   kaoo5_6  ,kaoo5_7  ,kaoo5_8  ,kaoo5_9  ,kaoo5_10
      integer   kaoo6_1  ,kaoo6_2  ,kaoo6_3  ,kaoo6_4  ,kaoo6_5
      integer   kaoo6_6  ,kaoo6_7  ,kaoo6_8  ,kaoo6_9  ,kaoo6_10
      integer   kaoo7_1  ,kaoo7_2  ,kaoo7_3  ,kaoo7_4  ,kaoo7_5
      integer   kaoo7_6  ,kaoo7_7  ,kaoo7_8  ,kaoo7_9  ,kaoo7_10
      integer   kaoo8_1  ,kaoo8_2  ,kaoo8_3  ,kaoo8_4  ,kaoo8_5
      integer   kaoo8_6  ,kaoo8_7  ,kaoo8_8  ,kaoo8_9  ,kaoo8_10

      integer   kabs1_1  ,kabs1_2  ,kabs1_3  ,kabs1_4  ,kabs1_5
      integer   kabs1_6  ,kabs1_7  ,kabs1_8  ,kabs1_9  ,kabs1_10
      integer   kabs2_1  ,kabs2_2  ,kabs2_3  ,kabs2_4  ,kabs2_5
      integer   kabs2_6  ,kabs2_7  ,kabs2_8  ,kabs2_9  ,kabs2_10
      integer   kabs3_1  ,kabs3_2  ,kabs3_3  ,kabs3_4  ,kabs3_5
      integer   kabs3_6  ,kabs3_7  ,kabs3_8  ,kabs3_9  ,kabs3_10
      integer   kabs4_1  ,kabs4_2  ,kabs4_3  ,kabs4_4  ,kabs4_5
      integer   kabs4_6  ,kabs4_7  ,kabs4_8  ,kabs4_9  ,kabs4_10
      integer   kabs5_1  ,kabs5_2  ,kabs5_3  ,kabs5_4  ,kabs5_5
      integer   kabs5_6  ,kabs5_7  ,kabs5_8  ,kabs5_9  ,kabs5_10

      integer   kaas1_1  ,kaas1_2  ,kaas1_3  ,kaas1_4  ,kaas1_5
      integer   kaas1_6  ,kaas1_7  ,kaas1_8  ,kaas1_9  ,kaas1_10
      integer   kaas2_1  ,kaas2_2  ,kaas2_3  ,kaas2_4  ,kaas2_5
      integer   kaas2_6  ,kaas2_7  ,kaas2_8  ,kaas2_9  ,kaas2_10
      integer   kaas3_1  ,kaas3_2  ,kaas3_3  ,kaas3_4  ,kaas3_5
      integer   kaas3_6  ,kaas3_7  ,kaas3_8  ,kaas3_9  ,kaas3_10
      integer   kaas4_1  ,kaas4_2  ,kaas4_3  ,kaas4_4  ,kaas4_5
      integer   kaas4_6  ,kaas4_7  ,kaas4_8  ,kaas4_9  ,kaas4_10
      integer   kaas5_1  ,kaas5_2  ,kaas5_3  ,kaas5_4  ,kaas5_5
      integer   kaas5_6  ,kaas5_7  ,kaas5_8  ,kaas5_9  ,kaas5_10

      integer   kans1_1  ,kans1_2  ,kans1_3  ,kans1_4  ,kans1_5
      integer   kans1_6  ,kans1_7  ,kans1_8  ,kans1_9  ,kans1_10
      integer   kans2_1  ,kans2_2  ,kans2_3  ,kans2_4  ,kans2_5
      integer   kans2_6  ,kans2_7  ,kans2_8  ,kans2_9  ,kans2_10
      integer   kans3_1  ,kans3_2  ,kans3_3  ,kans3_4  ,kans3_5
      integer   kans3_6  ,kans3_7  ,kans3_8  ,kans3_9  ,kans3_10
      integer   kans4_1  ,kans4_2  ,kans4_3  ,kans4_4  ,kans4_5
      integer   kans4_6  ,kans4_7  ,kans4_8  ,kans4_9  ,kans4_10
      integer   kans5_1  ,kans5_2  ,kans5_3  ,kans5_4  ,kans5_5
      integer   kans5_6  ,kans5_7  ,kans5_8  ,kans5_9  ,kans5_10
      integer   kans6_1  ,kans6_2  ,kans6_3  ,kans6_4  ,kans6_5
      integer   kans6_6  ,kans6_7  ,kans6_8  ,kans6_9  ,kans6_10
      integer   kans7_1  ,kans7_2  ,kans7_3  ,kans7_4  ,kans7_5
      integer   kans7_6  ,kans7_7  ,kans7_8  ,kans7_9  ,kans7_10
      integer   kans8_1  ,kans8_2  ,kans8_3  ,kans8_4  ,kans8_5
      integer   kans8_6  ,kans8_7  ,kans8_8  ,kans8_9  ,kans8_10
c
      integer   kpoc_1
      integer   kpoc_2   ,kpoc_3   ,kpoc_4
      integer   kpoc_5   ,kpoc_6   ,kpoc_7
      integer   kpoc_8   ,kpoc_9   ,kpoc_10
      integer   kpec_1   ,kpec_2   ,kpec_3
      integer   kpec_4   ,kpec_5   ,kpec_6
      integer   kpec_7   ,kpec_8   ,kpec_9
      integer   kpec_10  ,kcrust_1 ,kcrust_2
      integer   kcrust_3 ,kcrust_4 ,kcrust_5
      integer   kcrust_6 ,kcrust_7 ,kcrust_8
      integer   kcrust_9 ,kcrust_10,kph2o_1
      integer   kph2o_2  ,kph2o_3  ,kph2o_4
      integer   kph2o_5  ,kph2o_6  ,kph2o_7
      integer   kph2o_8  ,kph2o_9  ,kph2o_10
      integer   kpcl_1   ,kpcl_2   ,kpcl_3
      integer   kpcl_4   ,kpcl_5   ,kpcl_6
      integer   kpcl_7   ,kpcl_8   ,kpcl_9
      integer   kpcl_10  ,kna_1    ,kna_2
      integer   kna_3    ,kna_4    ,kna_5
      integer   kna_6    ,kna_7    ,kna_8
      integer   kna_9    ,kna_10   ,kpnh4_1
      integer   kpnh4_2  ,kpnh4_3  ,kpnh4_4
      integer   kpnh4_5  ,kpnh4_6  ,kpnh4_7
      integer   kpnh4_8  ,kpnh4_9  ,kpnh4_10
      integer   kpno3_1  ,kpno3_2  ,kpno3_3
      integer   kpno3_4  ,kpno3_5  ,kpno3_6
      integer   kpno3_7  ,kpno3_8  ,kpno3_9
      integer   kpno3_10 ,kpso4_1  ,kpso4_2
      integer   kpso4_3  ,kpso4_4  ,kpso4_5
      integer   kpso4_6  ,kpso4_7  ,kpso4_8
      integer   kpso4_9  ,kpso4_10 ,kph2o
c
      equivalence (kmap(1), kno  ), (kmap(2), kno2 ), (kmap(3), ko3  ),
     &            (kmap(4), kpan ), (kmap(5), kcres), (kmap(6), kpan2),
     &            (kmap(7), kmpan), (kmap(8), kpbzn), (kmap(9), knphe),
     &            (kmap(10),krno3), (kmap(11),kdcb2), (kmap(12),kdcb3),
     &            (kmap(13),khno4), (kmap(14),kacet), (kmap(15),kald2),
     &            (kmap(16),kalk1), (kmap(17),kalk2), (kmap(18),kalk3),
     &            (kmap(19),kalk4), (kmap(20),kalk5), (kmap(21),karo1),
     &            (kmap(22),karo2), (kmap(23),kbacl), (kmap(24),kbald),
     &            (kmap(25),kbcl1), (kmap(26),kbcl2), (kmap(27),kbuta),
     &            (kmap(28),kccho), (kmap(29),kccrs), 
     &		  (kmap(30),kcpo1), (kmap(31),kcpo2), (kmap(32),kcpo3), 
     &		  (kmap(33),kcpo4), (kmap(34),kcpo5), (kmap(35),kcpo6), 
     &		  (kmap(36),kcpo7), (kmap(37),kcpo8),
     &            (kmap(38),kcoo1), (kmap(39),kcoo2), (kmap(40),kcoo3),
     &            (kmap(41),kcoo4), (kmap(42),kcoo5), (kmap(43),kcoo6),
     &            (kmap(44),kcoo7), (kmap(45),kcoo8),
     &            (kmap(46),kcbs1), (kmap(47),kcbs2), (kmap(48),kcbs3), 
     &		  (kmap(49),kcbs4), (kmap(50),kcbs5),
     &		  (kmap(51),kcas1), (kmap(52),kcas2), (kmap(53),kcas3), 
     &		  (kmap(54),kcas4), (kmap(55),kcas5),
     &            (kmap(56),kcns1), (kmap(57),kcns2), (kmap(58),kcns3),
     &            (kmap(59),kcns4), (kmap(60),kcns5), (kmap(61),kcns6),
     &            (kmap(62),kcns7), (kmap(63),kcns8),
c
     &           (kmap(64),kcl2 ), (kmap(65),kco  ), (kmap(66),kco2h),
     &           (kmap(67),kco3h), (kmap(68),kcooh), (kmap(69),kcprm),
     &           (kmap(70),kdcb1), (kmap(71),keth ), (kmap(72),kethe),
     &           (kmap(73),ketoh), (kmap(74),kfcrs), (kmap(75),kfmcl),
     &           (kmap(76),kform), (kmap(77),kfprm), (kmap(78),kgly ),
     &           (kmap(79),kh2o2), (kmap(80),khc2h), (kmap(81),khcho),
     &           (kmap(82),khcl ), (kmap(83),khono), (kmap(84),khno3),
     &           (kmap(85),kho2h), (kmap(86),khocl), (kmap(87),kicl1),
     &           (kmap(88),kicl2), (kmap(89),kisop), (kmap(90),kispd),
     &           (kmap(91),kmek ), (kmap(92),kmeoh), (kmap(93),kmeth),
     &           (kmap(94),kmgly), (kmap(95),kmvk ), (kmap(96),kna  ),
     &           (kmap(97),knh3 ), (kmap(98),kntr ), (kmap(99),knxoy),
     &           (kmap(100),kole ),(kmap(101),kole1),(kmap(102),kole2),
     &           (kmap(103),kbpin),(kmap(104),klimo),(kmap(105),kmono),
     &           (kmap(106),ksesq),(kmap(107),kopen),(kmap(108),kpar ),
     &           (kmap(109),kpcl ),(kmap(110),kpec ),(kmap(111),kphen),
     &           (kmap(112),kpna ),(kmap(113),kpnh4),(kmap(114),kpno3),
     &           (kmap(115),kpoa ),(kmap(116),kprod),(kmap(117),kpso4),
     &           (kmap(118),krc2h),(kmap(119),krc3h),(kmap(120),krcho),
     &           (kmap(121),krooh),(kmap(122),kso2 ),
     &		 (kmap(123),kapo1),(kmap(124),kapo2),(kmap(125),kapo3),
     &           (kmap(126),kapo4),(kmap(127),kapo5),(kmap(128),kapo6),
     &		 (kmap(129),kapo7),(kmap(130),kapo8),
     &           (kmap(131),kaoo1),(kmap(132),kaoo2),(kmap(133),kaoo3),
     &           (kmap(134),kaoo4),(kmap(135),kaoo5),(kmap(136),kaoo6),
     &           (kmap(137),kaoo7),(kmap(138),kaoo8),
     &           (kmap(139),kabs1),(kmap(140),kabs2),(kmap(141),kabs3), 
     &		 (kmap(142),kabs4),(kmap(143),kabs5),
     &		 (kmap(144),kaas1),(kmap(145),kaas2),(kmap(146),kaas3),
     &		 (kmap(147),kaas4),(kmap(148),kaas5),
     &		 (kmap(149),kans1),(kmap(150),kans2),(kmap(151),kans3),
     &           (kmap(152),kans4),(kmap(153),kans5),(kmap(154),kans6),
     &		 (kmap(155),kans7),(kmap(156),kans8),
     &           (kmap(157),ksulf), (kmap(158),kterp),
     &           (kmap(159),ktol ), (kmap(160),kxn  ),(kmap(161),kxyl ), 
c
     &(kmap(162),kapo1_1 ), (kmap(163),kapo1_2 ), (kmap(164),kapo1_3 ),
     &(kmap(165),kapo1_4 ), (kmap(166),kapo1_5 ), (kmap(167),kapo1_6 ),
     &(kmap(168),kapo1_7 ), (kmap(169),kapo1_8 ), (kmap(170),kapo1_9 ),
     &(kmap(171),kapo1_10 ),(kmap(172),kapo2_1 ), (kmap(173),kapo2_2 ),
     &(kmap(174),kapo2_3 ), (kmap(175),kapo2_4 ), (kmap(176),kapo2_5 ),
     &(kmap(177),kapo2_6 ), (kmap(178),kapo2_7 ), (kmap(179),kapo2_8 ),
     &(kmap(180),kapo2_9 ), (kmap(181),kapo2_10 ),(kmap(182),kapo3_1 ),
     &(kmap(183),kapo3_2 ), (kmap(184),kapo3_3 ), (kmap(185),kapo3_4 ),
     &(kmap(186),kapo3_5 ), (kmap(187),kapo3_6 ), (kmap(188),kapo3_7 ),
     &(kmap(189),kapo3_8 ), (kmap(190),kapo3_9 ), (kmap(191),kapo3_10 ),
     &(kmap(192),kapo4_1 ), (kmap(193),kapo4_2 ), (kmap(194),kapo4_3 ),
     &(kmap(195),kapo4_4 ), (kmap(196),kapo4_5 ), (kmap(197),kapo4_6 ),
     &(kmap(198),kapo4_7 ), (kmap(199),kapo4_8 ), (kmap(200),kapo4_9 ),
     &(kmap(201),kapo4_10 ),(kmap(202),kapo5_1 ), (kmap(203),kapo5_2 ),
     &(kmap(204),kapo5_3 ), (kmap(205),kapo5_4 ), (kmap(206),kapo5_5 ),
     &(kmap(207),kapo5_6 ), (kmap(208),kapo5_7 ), (kmap(209),kapo5_8 ),
     &(kmap(210),kapo5_9 ), (kmap(211),kapo5_10 ),(kmap(212),kapo6_1 ),
     &(kmap(213),kapo6_2 ), (kmap(214),kapo6_3 ), (kmap(215),kapo6_4 ),
     &(kmap(216),kapo6_5 ), (kmap(217),kapo6_6 ), (kmap(218),kapo6_7 ),
     &(kmap(219),kapo6_8 ), (kmap(220),kapo6_9 ), (kmap(221),kapo6_10 ),
     &(kmap(222),kapo7_1 ), (kmap(223),kapo7_2 ), (kmap(224),kapo7_3 ),
     &(kmap(225),kapo7_4 ), (kmap(226),kapo7_5 ), (kmap(227),kapo7_6 ),
     &(kmap(228),kapo7_7 ), (kmap(229),kapo7_8 ), (kmap(230),kapo7_9 ),
     &(kmap(231),kapo7_10 ),(kmap(232),kapo8_1 ), (kmap(233),kapo8_2 ),
     &(kmap(234),kapo8_3 ), (kmap(235),kapo8_4 ), (kmap(236),kapo8_5 ),
     &(kmap(237),kapo8_6 ), (kmap(238),kapo8_7 ), (kmap(239),kapo8_8 ),
     &(kmap(240),kapo8_9 ), (kmap(241),kapo8_10 ),
     &(kmap(242),kaoo1_1 ), (kmap(243),kaoo1_2 ),
     &(kmap(244),kaoo1_3 ), (kmap(245),kaoo1_4 ), (kmap(246),kaoo1_5 ),
     &(kmap(247),kaoo1_6 ), (kmap(248),kaoo1_7 ), (kmap(249),kaoo1_8 ),
     &(kmap(250),kaoo1_9 ), (kmap(251),kaoo1_10 ),(kmap(252),kaoo2_1 ),
     &(kmap(253),kaoo2_2 ), (kmap(254),kaoo2_3 ), (kmap(255),kaoo2_4 ),
     &(kmap(256),kaoo2_5 ), (kmap(257),kaoo2_6 ), (kmap(258),kaoo2_7 ),
     &(kmap(259),kaoo2_8 ), (kmap(260),kaoo2_9 ), (kmap(261),kaoo2_10 ),
     &(kmap(262),kaoo3_1 ), (kmap(263),kaoo3_2 ), (kmap(264),kaoo3_3 ),
     &(kmap(265),kaoo3_4 ), (kmap(266),kaoo3_5 ), (kmap(267),kaoo3_6 ),
     &(kmap(268),kaoo3_7 ), (kmap(269),kaoo3_8 ), (kmap(270),kaoo3_9 ),
     &(kmap(271),kaoo3_10 ),(kmap(272),kaoo4_1 ), (kmap(273),kaoo4_2 ),
     &(kmap(274),kaoo4_3 ), (kmap(275),kaoo4_4 ), (kmap(276),kaoo4_5 ),
     &(kmap(277),kaoo4_6 ), (kmap(278),kaoo4_7 ), (kmap(279),kaoo4_8 ),
     &(kmap(280),kaoo4_9 ), (kmap(281),kaoo4_10 ),(kmap(282),kaoo5_1 ),
     &(kmap(283),kaoo5_2 ), (kmap(284),kaoo5_3 ), (kmap(285),kaoo5_4 ),
     &(kmap(286),kaoo5_5 ), (kmap(287),kaoo5_6 ), (kmap(288),kaoo5_7 ),
     &(kmap(289),kaoo5_8 ), (kmap(290),kaoo5_9 ), (kmap(291),kaoo5_10 ),
     &(kmap(292),kaoo6_1 ), (kmap(293),kaoo6_2 ), (kmap(294),kaoo6_3 ),
     &(kmap(295),kaoo6_4 ), (kmap(296),kaoo6_5 ), (kmap(297),kaoo6_6 ),
     &(kmap(298),kaoo6_7 ), (kmap(299),kaoo6_8 ), (kmap(300),kaoo6_9 ),
     &(kmap(301),kaoo6_10 ),(kmap(302),kaoo7_1 ), (kmap(303),kaoo7_2 ),
     &(kmap(304),kaoo7_3 ), (kmap(305),kaoo7_4 ), (kmap(306),kaoo7_5 ),
     &(kmap(307),kaoo7_6 ), (kmap(308),kaoo7_7 ), (kmap(309),kaoo7_8 ),
     &(kmap(310),kaoo7_9 ), (kmap(311),kaoo7_10 ),(kmap(312),kaoo8_1 ),
     &(kmap(313),kaoo8_2 ), (kmap(314),kaoo8_3 ), (kmap(315),kaoo8_4 ),
     &(kmap(316),kaoo8_5 ), (kmap(317),kaoo8_6 ), (kmap(318),kaoo8_7 ),
     &(kmap(319),kaoo8_8 ), (kmap(320),kaoo8_9 ), (kmap(321),kaoo8_10 ),
     &(kmap(322),kabs1_1 ),
     &(kmap(323),kabs1_2 ), (kmap(324),kabs1_3 ), (kmap(325),kabs1_4 ),
     &(kmap(326),kabs1_5 ), (kmap(327),kabs1_6 ), (kmap(328),kabs1_7 ),
     &(kmap(329),kabs1_8 ), (kmap(330),kabs1_9 ), (kmap(331),kabs1_10 ),
     &(kmap(332),kabs2_1 ), (kmap(333),kabs2_2 ), (kmap(334),kabs2_3 ),
     &(kmap(335),kabs2_4 ), (kmap(336),kabs2_5 ), (kmap(337),kabs2_6 ),
     &(kmap(338),kabs2_7 ), (kmap(339),kabs2_8 ), (kmap(340),kabs2_9 ),
     &(kmap(341),kabs2_10 ),(kmap(342),kabs3_1 ), (kmap(343),kabs3_2 ),
     &(kmap(344),kabs3_3 ), (kmap(345),kabs3_4 ), (kmap(346),kabs3_5 ),
     &(kmap(347),kabs3_6 ), (kmap(348),kabs3_7 ), (kmap(349),kabs3_8 ),
     &(kmap(350),kabs3_9 ), (kmap(351),kabs3_10 ),(kmap(352),kabs4_1 ),
     &(kmap(353),kabs4_2 ), (kmap(354),kabs4_3 ), (kmap(355),kabs4_4 ),
     &(kmap(356),kabs4_5 ), (kmap(357),kabs4_6 ), (kmap(358),kabs4_7 ),
     &(kmap(359),kabs4_8 ), (kmap(360),kabs4_9 ), (kmap(361),kabs4_10 ),
     &(kmap(362),kabs5_1 ),
     &(kmap(363),kabs5_2 ), (kmap(364),kabs5_3 ), (kmap(365),kabs5_4 ),
     &(kmap(366),kabs5_5 ), (kmap(367),kabs5_6 ), (kmap(368),kabs5_7 ),
     &(kmap(369),kabs5_8 ), (kmap(370),kabs5_9 ), (kmap(371),kabs5_10 ),
     &(kmap(372),kaas1_1 ), (kmap(373),kaas1_2 ), (kmap(374),kaas1_3 ),
     &(kmap(375),kaas1_4 ), (kmap(376),kaas1_5 ), (kmap(377),kaas1_6 ),
     &(kmap(378),kaas1_7 ), (kmap(379),kaas1_8 ), (kmap(380),kaas1_9 ),
     &(kmap(381),kaas1_10 ),(kmap(382),kaas2_1 ), (kmap(383),kaas2_2 ),
     &(kmap(384),kaas2_3 ), (kmap(385),kaas2_4 ), (kmap(386),kaas2_5 ),
     &(kmap(387),kaas2_6 ), (kmap(388),kaas2_7 ), (kmap(389),kaas2_8 ),
     &(kmap(390),kaas2_9 ), (kmap(391),kaas2_10 ),(kmap(392),kaas3_1 ),
     &(kmap(393),kaas3_2 ), (kmap(394),kaas3_3 ), (kmap(395),kaas3_4 ),
     &(kmap(396),kaas3_5 ), (kmap(397),kaas3_6 ), (kmap(398),kaas3_7 ),
     &(kmap(399),kaas3_8 ), (kmap(400),kaas3_9 ), (kmap(401),kaas3_10 ),
     &(kmap(402),kaas4_1 ), (kmap(403),kaas4_2 ), (kmap(404),kaas4_3 ),
     &(kmap(405),kaas4_4 ), (kmap(406),kaas4_5 ), (kmap(407),kaas4_6 ),
     &(kmap(408),kaas4_7 ), (kmap(409),kaas4_8 ), (kmap(410),kaas4_9 ),
     &(kmap(411),kaas4_10 ),
     &(kmap(412),kaas5_1 ), (kmap(413),kaas5_2 ), (kmap(414),kaas5_3 ),
     &(kmap(415),kaas5_4 ), (kmap(416),kaas5_5 ), (kmap(417),kaas5_6 ),
     &(kmap(418),kaas5_7 ), (kmap(419),kaas5_8 ), (kmap(420),kaas5_9 ),
     &(kmap(421),kaas5_10 ),
     &(kmap(422),kans1_1 ), (kmap(423),kans1_2 ), (kmap(424),kans1_3 ),
     &(kmap(425),kans1_4 ), (kmap(426),kans1_5 ), (kmap(427),kans1_6 ),
     &(kmap(428),kans1_7 ), (kmap(429),kans1_8 ), (kmap(430),kans1_9 ),
     &(kmap(431),kans1_10 ),(kmap(432),kans2_1 ), (kmap(433),kans2_2 ),
     &(kmap(434),kans2_3 ), (kmap(435),kans2_4 ), (kmap(436),kans2_5 ),
     &(kmap(437),kans2_6 ), (kmap(438),kans2_7 ), (kmap(439),kans2_8 ),
     &(kmap(440),kans2_9 ), (kmap(441),kans2_10 ),(kmap(442),kans3_1 ),
     &(kmap(443),kans3_2 ), (kmap(444),kans3_3 ), (kmap(445),kans3_4 ),
     &(kmap(446),kans3_5 ), (kmap(447),kans3_6 ), (kmap(448),kans3_7 ),
     &(kmap(449),kans3_8 ), (kmap(450),kans3_9 ), (kmap(451),kans3_10 ),
     &(kmap(452),kans4_1 ), (kmap(453),kans4_2 ), (kmap(454),kans4_3 ),
     &(kmap(455),kans4_4 ), (kmap(456),kans4_5 ), (kmap(457),kans4_6 ),
     &(kmap(458),kans4_7 ), (kmap(459),kans4_8 ), (kmap(460),kans4_9 ),
     &(kmap(461),kans4_10 ),(kmap(462),kans5_1 ), (kmap(463),kans5_2 ),
     &(kmap(464),kans5_3 ), (kmap(465),kans5_4 ), (kmap(466),kans5_5 ),
     &(kmap(467),kans5_6 ), (kmap(468),kans5_7 ), (kmap(469),kans5_8 ),
     &(kmap(470),kans5_9 ), (kmap(471),kans5_10 ),(kmap(472),kans6_1 ),
     &(kmap(473),kans6_2 ), (kmap(474),kans6_3 ), (kmap(475),kans6_4 ),
     &(kmap(476),kans6_5 ), (kmap(477),kans6_6 ), (kmap(478),kans6_7 ),
     &(kmap(479),kans6_8 ), (kmap(480),kans6_9 ), (kmap(481),kans6_10 ),
     &(kmap(482),kans7_1 ), (kmap(483),kans7_2 ), (kmap(484),kans7_3 ),
     &(kmap(485),kans7_4 ), (kmap(486),kans7_5 ), (kmap(487),kans7_6 ),
     &(kmap(488),kans7_7 ), (kmap(489),kans7_8 ), (kmap(490),kans7_9 ),
     &(kmap(491),kans7_10 ),(kmap(492),kans8_1 ), (kmap(493),kans8_2 ),
     &(kmap(494),kans8_3 ), (kmap(495),kans8_4 ), (kmap(496),kans8_5 ),
     &(kmap(497),kans8_6 ), (kmap(498),kans8_7 ), (kmap(499),kans8_8 ),
     &(kmap(500),kans8_9 ), (kmap(501),kans8_10 ),

c   scf

     &(kmap(502),kpoc_1  ),
     &(kmap(503),kpoc_2  ),(kmap(504),kpoc_3  ), (kmap(505),kpoc_4  ),
     &(kmap(506),kpoc_5  ),(kmap(507),kpoc_6  ), (kmap(508),kpoc_7  ),
     &(kmap(509),kpoc_8  ),(kmap(510),kpoc_9  ), (kmap(511),kpoc_10 ),
     &(kmap(512),kpec_1  ),(kmap(513),kpec_2  ), (kmap(514),kpec_3  ),
     &(kmap(515),kpec_4  ),(kmap(516),kpec_5  ), (kmap(517),kpec_6  ),
     &(kmap(518),kpec_7  ),(kmap(519),kpec_8  ), (kmap(520),kpec_9  ),
     &(kmap(521),kpec_10 ),(kmap(522),kcrust_1), (kmap(523),kcrust_2),
     &(kmap(524),kcrust_3),(kmap(525),kcrust_4), (kmap(526),kcrust_5),
     &(kmap(527),kcrust_6),(kmap(528),kcrust_7), (kmap(529),kcrust_8),
     &(kmap(530),kcrust_9),(kmap(531),kcrust_10),(kmap(532),kph2o_1),
     &(kmap(533),kph2o_2 ),(kmap(534),kph2o_3 ), (kmap(535),kph2o_4 ),
     &(kmap(536),kph2o_5 ),(kmap(537),kph2o_6 ), (kmap(538),kph2o_7 ),
     &(kmap(539),kph2o_8 ),(kmap(540),kph2o_9 ), (kmap(541),kph2o_10),
     &(kmap(542),kpcl_1  ),(kmap(543),kpcl_2  ), (kmap(544),kpcl_3  ),
     &(kmap(545),kpcl_4  ),(kmap(546),kpcl_5  ), (kmap(547),kpcl_6  ),
     &(kmap(548),kpcl_7  ),(kmap(549),kpcl_8  ), (kmap(550),kpcl_9  ),
     &(kmap(551),kpcl_10 ),(kmap(552),kna_1   ), (kmap(553),kna_2   ),
     &(kmap(554),kna_3   ),(kmap(555),kna_4   ), (kmap(556),kna_5   ),
     &(kmap(557),kna_6   ),(kmap(558),kna_7   ), (kmap(559),kna_8   ),
     &(kmap(560),kna_9   ),(kmap(561),kna_10  ), (kmap(562),kpnh4_1 ),
     &(kmap(563),kpnh4_2 ),(kmap(564),kpnh4_3 ), (kmap(565),kpnh4_4 ),
     &(kmap(566),kpnh4_5 ),(kmap(567),kpnh4_6 ), (kmap(568),kpnh4_7 ),
     &(kmap(569),kpnh4_8 ),(kmap(570),kpnh4_9 ), (kmap(571),kpnh4_10),
     &(kmap(572),kpno3_1 ),(kmap(573),kpno3_2 ), (kmap(574),kpno3_3 ),
     &(kmap(575),kpno3_4 ),(kmap(576),kpno3_5 ), (kmap(577),kpno3_6 ),
     &(kmap(578),kpno3_7 ),(kmap(579),kpno3_8 ), (kmap(580),kpno3_9 ),
     &(kmap(581),kpno3_10),(kmap(582),kpso4_1 ), (kmap(583),kpso4_2 ),
     &(kmap(584),kpso4_3 ),(kmap(585),kpso4_4 ), (kmap(586),kpso4_5 ),
     &(kmap(587),kpso4_6 ),(kmap(588),kpso4_7 ), (kmap(589),kpso4_8 ),
     &(kmap(590),kpso4_9 ),(kmap(591),kpso4_10), (kmap(592),kph2o   )
c
      integer   ko1d  ,ko    ,kclo 
      integer   kcl   ,kn2o5 ,kno3 
      integer   koh   ,kho2  ,kc2o3
      integer   kxo2  ,kxo2n ,kto2 
      integer   kror  ,kcro  ,kro2r
      integer   kr2o2 ,kro2n ,kcco3
      integer   krco3 ,kmco3 ,kbzco
      integer   kcxo2 ,khco3 ,ktbuo
      integer   kbzo  ,kbzno
c
      equivalence (krad(1), ko1d ), (krad(2), ko   ), (krad(3), kclo ),
     &            (krad(4), kcl  ), (krad(5), kn2o5), (krad(6), kno3 ),
     &            (krad(7), koh  ), (krad(8), kho2 ), (krad(9), kc2o3),
     &            (krad(10),kxo2 ), (krad(11),kxo2n), (krad(12),kto2 ),
     &            (krad(13),kror ), (krad(14),kcro ), (krad(15),kro2r),
     &            (krad(16),kr2o2), (krad(17),kro2n), (krad(18),kcco3),
     &            (krad(19),krco3), (krad(20),kmco3), (krad(21),kbzco),
     &            (krad(22),kcxo2), (krad(23),khco3), (krad(24),ktbuo),
     &            (krad(25),kbzo ), (krad(26),kbzno)
c
c-----------------------------------------------------------------------
c     Variables for chemistry lookup tables:
c
c     tempr  -- temperature table
c     presr  -- pressure table
c     rktbl  -- temperature/pressure-dependent rate constant table
c     htint  -- height AGL table
c     zenint -- zenith angle table
c     prkn   -- reaction rate table
c-----------------------------------------------------------------------
c      
      common /tables/ tempr(NTEMPR), presr(NPRESR),
     &                rktbl(MXRXN,NTEMPR,NPRESR),
     &                htint(NHGHT), zenint(NZEN),
     &                prkn(NZEN,MXPHT1,NHGHT,NHAZE,NALB,NOZN)
c
c-----------------------------------------------------------------------
c     Variables to define parameters for each chemical species:
c
c     henry0   -- Henry's Law constant at STP (molar/atm)
c     tfact    -- Temperature dependence of Henry's Law constant (1/K)
c     diffrat  -- Species diffusivity
c     f0       -- Species reactivity parameter
c     rscale   -- Species scaling factor for surface resistance
c     henso20  -- Henry's Law constant at STP for SO2 (molar/atm)
c     tfactso2 -- Temperature dependence of SO2 Henry's Law constant (1/K)
c     nbin     -- Number of aerosol size bins
c     roprt    -- Aerosol density (g/m3)
c     dcut     -- Aerosol size bin cut points (um)
c     cwmin    -- Minimum cloud water threshold (g/m3)
c     tamin    -- Cloud water freezing threshold (K)
c-----------------------------------------------------------------------
c
      real cwmin,tamin
      common /depchm/ henry0(MXSPEC),tfact(MXSPEC),diffrat(MXSPEC),
     &                f0(MXSPEC),rscale(MXSPEC),henso20,tfactso2,cwmin,
     &                tamin
      common /aerochm/ nbin,roprt(MXSPEC),dcut(MXSPEC,2)
c
c-----------------------------------------------------------------------
c     Pointers used to lookup pig chemistry rate constants
c
c     ipigrxn  -- pointers to the nine reactions
c                 (1)   NO2 + O3 -> NO3
c                 (2)         O3 -> O(1D)
c                 (3)      O(1D) -> O(3P)
c                 (4)      O(1D) -> 2 OH
c                 (5)  NO3 + NO2 -> NO + NO2
c                 (6)  NO3 + NO2 -> N2O5
c                 (7) N2O5 + H2O -> 2 HNO3
c                 (8)       N2O5 -> NO3 + NO2
c                 (9)    NO + NO -> 2 NO2
c
      common /pigrxn/ ipigrxn(9)
c
c----------------------------------------------------------------------
c    Variables for controlling calls to aerosol routines 
c    
c     grd_time   -- time for calls to aerosol routines for each grid
c     date_aer   -- Julian date of current grd_time for each grid
c     dtaero     -- time interval between calls to aerosol routines
c     aero_dt    -- time between calls to aerosol routines for each grid
c     dt_aero    -- user input time interval between calls to aerosol routines
c
      real grd_time(MXGRID)
      real aero_dt(MXGRID)
      real dtaero
      integer date_aer(MXGRID)
c
      common /aero_t/ grd_time,aero_dt,date_aer,dtaero,dt_aero

