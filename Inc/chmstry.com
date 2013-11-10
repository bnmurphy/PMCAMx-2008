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
      integer   kcpo9, kcpo10
      integer   kcoo1, kcoo2, kcoo3, kcoo4
      integer   kcoo5, kcoo6, kcoo7, kcoo8
      integer   kcoo9, kcoo10
      integer   kcbs1, kcbs2, kcbs3, kcbs4
      integer   kcas1, kcas2, kcas3, kcas4
c
      integer   kapo1, kapo2, kapo3, kapo4
      integer   kapo5, kapo6, kapo7, kapo8
      integer   kapo9, kapo10
      integer   kaoo1, kaoo2, kaoo3, kaoo4
      integer   kaoo5, kaoo6, kaoo7, kaoo8
      integer   kaoo9, kaoo10
      integer   kabs1, kabs2, kabs3, kabs4
      integer   kaas1, kaas2, kaas3, kaas4
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
      integer   kapo9_1  ,kapo9_2  ,kapo9_3  ,kapo9_4  ,kapo9_5
      integer   kapo9_6  ,kapo9_7  ,kapo9_8  ,kapo9_9  ,kapo9_10
      integer   kapo10_1  ,kapo10_2  ,kapo10_3  ,kapo10_4  ,kapo10_5
      integer   kapo10_6  ,kapo10_7  ,kapo10_8  ,kapo10_9  ,kapo10_10

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
      integer   kaoo9_1  ,kaoo9_2  ,kaoo9_3  ,kaoo9_4  ,kaoo9_5
      integer   kaoo9_6  ,kaoo9_7  ,kaoo9_8  ,kaoo9_9  ,kaoo9_10
      integer   kaoo10_1  ,kaoo10_2  ,kaoo10_3  ,kaoo10_4  ,kaoo10_5
      integer   kaoo10_6  ,kaoo10_7  ,kaoo10_8  ,kaoo10_9  ,kaoo10_10

      integer   kabs1_1  ,kabs1_2  ,kabs1_3  ,kabs1_4  ,kabs1_5
      integer   kabs1_6  ,kabs1_7  ,kabs1_8  ,kabs1_9  ,kabs1_10
      integer   kabs2_1  ,kabs2_2  ,kabs2_3  ,kabs2_4  ,kabs2_5
      integer   kabs2_6  ,kabs2_7  ,kabs2_8  ,kabs2_9  ,kabs2_10
      integer   kabs3_1  ,kabs3_2  ,kabs3_3  ,kabs3_4  ,kabs3_5
      integer   kabs3_6  ,kabs3_7  ,kabs3_8  ,kabs3_9  ,kabs3_10
      integer   kabs4_1  ,kabs4_2  ,kabs4_3  ,kabs4_4  ,kabs4_5
      integer   kabs4_6  ,kabs4_7  ,kabs4_8  ,kabs4_9  ,kabs4_10

      integer   kaas1_1  ,kaas1_2  ,kaas1_3  ,kaas1_4  ,kaas1_5
      integer   kaas1_6  ,kaas1_7  ,kaas1_8  ,kaas1_9  ,kaas1_10
      integer   kaas2_1  ,kaas2_2  ,kaas2_3  ,kaas2_4  ,kaas2_5
      integer   kaas2_6  ,kaas2_7  ,kaas2_8  ,kaas2_9  ,kaas2_10
      integer   kaas3_1  ,kaas3_2  ,kaas3_3  ,kaas3_4  ,kaas3_5
      integer   kaas3_6  ,kaas3_7  ,kaas3_8  ,kaas3_9  ,kaas3_10
      integer   kaas4_1  ,kaas4_2  ,kaas4_3  ,kaas4_4  ,kaas4_5
      integer   kaas4_6  ,kaas4_7  ,kaas4_8  ,kaas4_9  ,kaas4_10

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
     &            (kmap(28),kccho), (kmap(29),kccrs), (kmap(30),kcpo1),
     &            (kmap(31),kcpo2), (kmap(32),kcpo3), (kmap(33),kcpo4),
     &            (kmap(34),kcpo5), (kmap(35),kcpo6), (kmap(36),kcpo7),
     &            (kmap(37),kcpo8), (kmap(38),kcpo9), (kmap(39),kcpo10),
     &            (kmap(40),kcoo1), (kmap(41),kcoo2), (kmap(42),kcoo3),
     &            (kmap(43),kcoo4), (kmap(44),kcoo5), (kmap(45),kcoo6),
     &            (kmap(46),kcoo7), (kmap(47),kcoo8), (kmap(48),kcoo9),
     &            (kmap(49),kcoo10), (kmap(50),kcbs1), (kmap(51),kcbs2),
     &            (kmap(52),kcbs3), (kmap(53),kcbs4), (kmap(54),kcas1),
     &            (kmap(55),kcas2), (kmap(56),kcas3), (kmap(57),kcas4),
c
     &           (kmap(58),kcl2 ), (kmap(59),kco  ), (kmap(60),kco2h),
     &           (kmap(61),kco3h), (kmap(62),kcooh), (kmap(63),kcprm),
     &           (kmap(64),kdcb1), (kmap(65),keth ), (kmap(66),kethe),
     &           (kmap(67),ketoh), (kmap(68),kfcrs), (kmap(69),kfmcl),
     &           (kmap(70),kform), (kmap(71),kfprm), (kmap(72),kgly ),
     &           (kmap(73),kh2o2), (kmap(74),khc2h), (kmap(75),khcho),
     &           (kmap(76),khcl ), (kmap(77),khono), (kmap(78),khno3),
     &           (kmap(79),kho2h), (kmap(80),khocl), (kmap(81),kicl1),
     &           (kmap(82),kicl2), (kmap(83),kisop), (kmap(84),kispd),
     &           (kmap(85),kmek ),(kmap(86),kmeoh),(kmap(87),kmeth),
     &           (kmap(88),kmgly),(kmap(89),kmvk ),(kmap(90),kna  ),
     &           (kmap(91),knh3 ),(kmap(92),kntr ),(kmap(93),knxoy),
     &           (kmap(94),kole ),(kmap(95),kole1),(kmap(96),kole2),
     &           (kmap(97),kbpin),(kmap(98),klimo),(kmap(99),kmono),
     &           (kmap(100),ksesq),(kmap(101),kopen),(kmap(102),kpar ),
     &           (kmap(103),kpcl ),(kmap(104),kpec ),(kmap(105),kphen),
     &           (kmap(106),kpna ),(kmap(107),kpnh4),(kmap(108),kpno3),
     &           (kmap(109),kpoa ),(kmap(110),kprod),(kmap(111),kpso4),
     &           (kmap(112),krc2h),(kmap(113),krc3h),(kmap(114),krcho),
     &           (kmap(115),krooh),(kmap(116),kso2 ),(kmap(117),kapo1),
     &           (kmap(118),kapo2),(kmap(119),kapo3),(kmap(120),kapo4),
     &           (kmap(121),kapo5),(kmap(122),kapo6),(kmap(123),kapo7),
     &           (kmap(124),kapo8),(kmap(125),kapo9),(kmap(126),kapo10),
     &           (kmap(127),kaoo1),(kmap(128),kaoo2),(kmap(129),kaoo3),
     &           (kmap(130),kaoo4),(kmap(131),kaoo5),(kmap(132),kaoo6),
     &           (kmap(133),kaoo7),(kmap(134),kaoo8),(kmap(135),kaoo9),
     &           (kmap(136),kaoo10),(kmap(137),kabs1),(kmap(138),kabs2),
     &           (kmap(139),kabs3),(kmap(140),kabs4),(kmap(141),kaas1),
     &           (kmap(142),kaas2),(kmap(143),kaas3),(kmap(144),kaas4),
     &           (kmap(145),ksulf),(kmap(146),kterp),
     &           (kmap(147),ktol ),(kmap(148),kxn  ),(kmap(149),kxyl ), 
c
     &(kmap(150),kapo1_1 ),(kmap(151),kapo1_2 ),(kmap(152),kapo1_3 ),
     &(kmap(153),kapo1_4 ),(kmap(154),kapo1_5 ),(kmap(155),kapo1_6 ),
     &(kmap(156),kapo1_7 ),(kmap(157),kapo1_8 ),(kmap(158),kapo1_9 ),
     &(kmap(159),kapo1_10 ),(kmap(160),kapo2_1 ),(kmap(161),kapo2_2 ),
     &(kmap(162),kapo2_3 ),(kmap(163),kapo2_4 ),(kmap(164),kapo2_5 ),
     &(kmap(165),kapo2_6 ),(kmap(166),kapo2_7 ),(kmap(167),kapo2_8 ),
     &(kmap(168),kapo2_9 ),(kmap(169),kapo2_10 ),(kmap(170),kapo3_1 ),
     &(kmap(171),kapo3_2 ),(kmap(172),kapo3_3 ),(kmap(173),kapo3_4 ),
     &(kmap(174),kapo3_5 ),(kmap(175),kapo3_6 ),(kmap(176),kapo3_7 ),
     &(kmap(177),kapo3_8 ),(kmap(178),kapo3_9 ),(kmap(179),kapo3_10 ),
     &(kmap(180),kapo4_1 ),(kmap(181),kapo4_2 ),(kmap(182),kapo4_3 ),
     &(kmap(183),kapo4_4 ),(kmap(184),kapo4_5 ),(kmap(185),kapo4_6 ),
     &(kmap(186),kapo4_7 ),(kmap(187),kapo4_8 ),(kmap(188),kapo4_9 ),
     &(kmap(189),kapo4_10 ),(kmap(190),kapo5_1 ),(kmap(191),kapo5_2 ),
     &(kmap(192),kapo5_3 ),(kmap(193),kapo5_4 ),(kmap(194),kapo5_5 ),
     &(kmap(195),kapo5_6 ),(kmap(196),kapo5_7 ),(kmap(197),kapo5_8 ),
     &(kmap(198),kapo5_9 ),(kmap(199),kapo5_10 ),(kmap(200),kapo6_1 ),
     &(kmap(201),kapo6_2 ),(kmap(202),kapo6_3 ),(kmap(203),kapo6_4 ),
     &(kmap(204),kapo6_5 ),(kmap(205),kapo6_6 ),(kmap(206),kapo6_7 ),
     &(kmap(207),kapo6_8 ),(kmap(208),kapo6_9 ),(kmap(209),kapo6_10 ),
     &(kmap(210),kapo7_1 ),(kmap(211),kapo7_2 ),(kmap(212),kapo7_3 ),
     &(kmap(213),kapo7_4 ),(kmap(214),kapo7_5 ),(kmap(215),kapo7_6 ),
     &(kmap(216),kapo7_7 ),(kmap(217),kapo7_8 ),(kmap(218),kapo7_9 ),
     &(kmap(219),kapo7_10 ),(kmap(220),kapo8_1 ),(kmap(221),kapo8_2 ),
     &(kmap(222),kapo8_3 ),(kmap(223),kapo8_4 ),(kmap(224),kapo8_5 ),
     &(kmap(225),kapo8_6 ),(kmap(226),kapo8_7 ),(kmap(227),kapo8_8 ),
     &(kmap(228),kapo8_9 ),(kmap(229),kapo8_10 ),(kmap(230),kapo9_1 ),
     &(kmap(231),kapo9_2 ),(kmap(232),kapo9_3 ),(kmap(233),kapo9_4 ),
     &(kmap(234),kapo9_5 ),(kmap(235),kapo9_6 ),(kmap(236),kapo9_7 ),
     &(kmap(237),kapo9_8 ),(kmap(238),kapo9_9 ),(kmap(239),kapo9_10 ),
     &(kmap(240),kapo10_1 ),(kmap(241),kapo10_2 ),(kmap(242),kapo10_3 ),
     &(kmap(243),kapo10_4 ),(kmap(244),kapo10_5 ),(kmap(245),kapo10_6 ),
     &(kmap(246),kapo10_7 ),(kmap(247),kapo10_8 ),(kmap(248),kapo10_9 ),
     &(kmap(249),kapo10_10 ),(kmap(250),kaoo1_1 ),(kmap(251),kaoo1_2 ),
     &(kmap(252),kaoo1_3 ),(kmap(253),kaoo1_4 ),(kmap(254),kaoo1_5 ),
     &(kmap(255),kaoo1_6 ),(kmap(256),kaoo1_7 ),(kmap(257),kaoo1_8 ),
     &(kmap(258),kaoo1_9 ),(kmap(259),kaoo1_10 ),(kmap(260),kaoo2_1 ),
     &(kmap(261),kaoo2_2 ),(kmap(262),kaoo2_3 ),(kmap(263),kaoo2_4 ),
     &(kmap(264),kaoo2_5 ),(kmap(265),kaoo2_6 ),(kmap(266),kaoo2_7 ),
     &(kmap(267),kaoo2_8 ),(kmap(268),kaoo2_9 ),(kmap(269),kaoo2_10 ),
     &(kmap(270),kaoo3_1 ),(kmap(271),kaoo3_2 ),(kmap(272),kaoo3_3 ),
     &(kmap(273),kaoo3_4 ),(kmap(274),kaoo3_5 ),(kmap(275),kaoo3_6 ),
     &(kmap(276),kaoo3_7 ),(kmap(277),kaoo3_8 ),(kmap(278),kaoo3_9 ),
     &(kmap(279),kaoo3_10 ),(kmap(280),kaoo4_1 ),(kmap(281),kaoo4_2 ),
     &(kmap(282),kaoo4_3 ),(kmap(283),kaoo4_4 ),(kmap(284),kaoo4_5 ),
     &(kmap(285),kaoo4_6 ),(kmap(286),kaoo4_7 ),(kmap(287),kaoo4_8 ),
     &(kmap(288),kaoo4_9 ),(kmap(289),kaoo4_10 ),(kmap(290),kaoo5_1 ),
     &(kmap(291),kaoo5_2 ),(kmap(292),kaoo5_3 ),(kmap(293),kaoo5_4 ),
     &(kmap(294),kaoo5_5 ),(kmap(295),kaoo5_6 ),(kmap(296),kaoo5_7 ),
     &(kmap(297),kaoo5_8 ),(kmap(298),kaoo5_9 ),(kmap(299),kaoo5_10 ),
     &(kmap(300),kaoo6_1 ),(kmap(301),kaoo6_2 ),(kmap(302),kaoo6_3 ),
     &(kmap(303),kaoo6_4 ),(kmap(304),kaoo6_5 ),(kmap(305),kaoo6_6 ),
     &(kmap(306),kaoo6_7 ),(kmap(307),kaoo6_8 ),(kmap(308),kaoo6_9 ),
     &(kmap(309),kaoo6_10 ),(kmap(310),kaoo7_1 ),(kmap(311),kaoo7_2 ),
     &(kmap(312),kaoo7_3 ),(kmap(313),kaoo7_4 ),(kmap(314),kaoo7_5 ),
     &(kmap(315),kaoo7_6 ),(kmap(316),kaoo7_7 ),(kmap(317),kaoo7_8 ),
     &(kmap(318),kaoo7_9 ),(kmap(319),kaoo7_10 ),(kmap(320),kaoo8_1 ),
     &(kmap(321),kaoo8_2 ),(kmap(322),kaoo8_3 ),(kmap(323),kaoo8_4 ),
     &(kmap(324),kaoo8_5 ),(kmap(325),kaoo8_6 ),(kmap(326),kaoo8_7 ),
     &(kmap(327),kaoo8_8 ),(kmap(328),kaoo8_9 ),(kmap(329),kaoo8_10 ),
     &(kmap(330),kaoo9_1 ),(kmap(331),kaoo9_2 ),(kmap(332),kaoo9_3 ),
     &(kmap(333),kaoo9_4 ),(kmap(334),kaoo9_5 ),(kmap(335),kaoo9_6 ),
     &(kmap(336),kaoo9_7 ),(kmap(337),kaoo9_8 ),(kmap(338),kaoo9_9 ),
     &(kmap(339),kaoo9_10 ),(kmap(340),kaoo10_1 ),(kmap(341),kaoo10_2 ),
     &(kmap(342),kaoo10_3 ),(kmap(343),kaoo10_4 ),(kmap(344),kaoo10_5 ),
     &(kmap(345),kaoo10_6 ),(kmap(346),kaoo10_7 ),(kmap(347),kaoo10_8 ),
     &(kmap(348),kaoo10_9 ),(kmap(349),kaoo10_10 ),(kmap(350),kabs1_1 ),
     &(kmap(351),kabs1_2 ),(kmap(352),kabs1_3 ),(kmap(353),kabs1_4 ),
     &(kmap(354),kabs1_5 ),(kmap(355),kabs1_6 ),(kmap(356),kabs1_7 ),
     &(kmap(357),kabs1_8 ),(kmap(358),kabs1_9 ),(kmap(359),kabs1_10 ),
     &(kmap(360),kabs2_1 ),(kmap(361),kabs2_2 ),(kmap(362),kabs2_3 ),
     &(kmap(363),kabs2_4 ),(kmap(364),kabs2_5 ),(kmap(365),kabs2_6 ),
     &(kmap(366),kabs2_7 ),(kmap(367),kabs2_8 ),(kmap(368),kabs2_9 ),
     &(kmap(369),kabs2_10 ),(kmap(370),kabs3_1 ),(kmap(371),kabs3_2 ),
     &(kmap(372),kabs3_3 ),(kmap(373),kabs3_4 ),(kmap(374),kabs3_5 ),
     &(kmap(375),kabs3_6 ),(kmap(376),kabs3_7 ),(kmap(377),kabs3_8 ),
     &(kmap(378),kabs3_9 ),(kmap(379),kabs3_10 ),(kmap(380),kabs4_1 ),
     &(kmap(381),kabs4_2 ),(kmap(382),kabs4_3 ),(kmap(383),kabs4_4 ),
     &(kmap(384),kabs4_5 ),(kmap(385),kabs4_6 ),(kmap(386),kabs4_7 ),
     &(kmap(387),kabs4_8 ),(kmap(388),kabs4_9 ),(kmap(389),kabs4_10 ),
     &(kmap(390),kaas1_1 ),(kmap(391),kaas1_2 ),(kmap(392),kaas1_3 ),
     &(kmap(393),kaas1_4 ),(kmap(394),kaas1_5 ),(kmap(395),kaas1_6 ),
     &(kmap(396),kaas1_7 ),(kmap(397),kaas1_8 ),(kmap(398),kaas1_9 ),
     &(kmap(399),kaas1_10 ),(kmap(400),kaas2_1 ),(kmap(401),kaas2_2 ),
     &(kmap(402),kaas2_3 ),(kmap(403),kaas2_4 ),(kmap(404),kaas2_5 ),
     &(kmap(405),kaas2_6 ),(kmap(406),kaas2_7 ),(kmap(407),kaas2_8 ),
     &(kmap(408),kaas2_9 ),(kmap(409),kaas2_10 ),(kmap(410),kaas3_1 ),
     &(kmap(411),kaas3_2 ),(kmap(412),kaas3_3 ),(kmap(413),kaas3_4 ),
     &(kmap(414),kaas3_5 ),(kmap(415),kaas3_6 ),(kmap(416),kaas3_7 ),
     &(kmap(417),kaas3_8 ),(kmap(418),kaas3_9 ),(kmap(419),kaas3_10 ),
     &(kmap(420),kaas4_1 ),(kmap(421),kaas4_2 ),(kmap(422),kaas4_3 ),
     &(kmap(423),kaas4_4 ),(kmap(424),kaas4_5 ),(kmap(425),kaas4_6 ),
     &(kmap(426),kaas4_7 ),(kmap(427),kaas4_8 ),(kmap(428),kaas4_9 ),
     &(kmap(429),kaas4_10 ),

c   scf

     &(kmap(430),kpoc_1  ),
     &(kmap(431),kpoc_2  ),(kmap(432),kpoc_3  ),(kmap(433),kpoc_4  ),
     &(kmap(434),kpoc_5  ),(kmap(435),kpoc_6  ),(kmap(436),kpoc_7  ),
     &(kmap(437),kpoc_8  ),(kmap(438),kpoc_9  ),(kmap(439),kpoc_10 ),
     &(kmap(440),kpec_1  ),(kmap(441),kpec_2  ),(kmap(442),kpec_3  ),
     &(kmap(443),kpec_4  ),(kmap(444),kpec_5  ),(kmap(445),kpec_6  ),
     &(kmap(446),kpec_7  ),(kmap(447),kpec_8  ),(kmap(448),kpec_9  ),
     &(kmap(449),kpec_10 ),(kmap(450),kcrust_1),(kmap(451),kcrust_2),
     &(kmap(452),kcrust_3),(kmap(453),kcrust_4),(kmap(454),kcrust_5),
     &(kmap(455),kcrust_6),(kmap(456),kcrust_7),(kmap(457),kcrust_8),
     &(kmap(458),kcrust_9),(kmap(459),kcrust_10),(kmap(460),kph2o_1),
     &(kmap(461),kph2o_2 ),(kmap(462),kph2o_3 ),(kmap(463),kph2o_4 ),
     &(kmap(464),kph2o_5 ),(kmap(465),kph2o_6 ),(kmap(466),kph2o_7 ),
     &(kmap(467),kph2o_8 ),(kmap(468),kph2o_9 ),(kmap(469),kph2o_10),
     &(kmap(470),kpcl_1  ),(kmap(471),kpcl_2  ),(kmap(472),kpcl_3  ),
     &(kmap(473),kpcl_4  ),(kmap(474),kpcl_5  ),(kmap(475),kpcl_6  ),
     &(kmap(476),kpcl_7  ),(kmap(477),kpcl_8  ),(kmap(478),kpcl_9  ),
     &(kmap(479),kpcl_10 ),(kmap(480),kna_1   ),(kmap(481),kna_2   ),
     &(kmap(482),kna_3   ),(kmap(483),kna_4   ),(kmap(484),kna_5   ),
     &(kmap(485),kna_6   ),(kmap(486),kna_7   ),(kmap(487),kna_8   ),
     &(kmap(488),kna_9   ),(kmap(489),kna_10  ),(kmap(490),kpnh4_1 ),
     &(kmap(491),kpnh4_2 ),(kmap(492),kpnh4_3 ),(kmap(493),kpnh4_4 ),
     &(kmap(494),kpnh4_5 ),(kmap(495),kpnh4_6 ),(kmap(496),kpnh4_7 ),
     &(kmap(497),kpnh4_8 ),(kmap(498),kpnh4_9 ),(kmap(499),kpnh4_10),
     &(kmap(500),kpno3_1 ),(kmap(501),kpno3_2 ),(kmap(502),kpno3_3 ),
     &(kmap(503),kpno3_4 ),(kmap(504),kpno3_5 ),(kmap(505),kpno3_6 ),
     &(kmap(506),kpno3_7 ),(kmap(507),kpno3_8 ),(kmap(508),kpno3_9 ),
     &(kmap(509),kpno3_10),(kmap(510),kpso4_1 ),(kmap(511),kpso4_2 ),
     &(kmap(512),kpso4_3 ),(kmap(513),kpso4_4 ),(kmap(514),kpso4_5 ),
     &(kmap(515),kpso4_6 ),(kmap(516),kpso4_7 ),(kmap(517),kpso4_8 ),
     &(kmap(518),kpso4_9 ),(kmap(519),kpso4_10),(kmap(520),kph2o   )
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
c
c---------------------------------------------------------------------
cBNM
c    Variables for printing radical concentrations
c
c	nrads -- number of radicals to print output for
c       crad  -- String array of radical names
c	radout -- Mapped variable for radical concentration pointer
	
      integer   nrads
      parameter (nrads = 3)

      character*4  crads(nrads)
      real  bnmradcnc(nrads,MXCOL1,MXROW1,MXLAY1)

      data crads /'OH','NO3','HO2'/  

      common bnmradcnc    
