c-----CAMx v4.02 030709
c
c     DDMCHM.COM sets species pointers for the DDM chemistry
c     These equivalences must be consistent with the internal
c     species lists defined in data statements in READCHM
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
      common /lname/ lmap(NSPNAM), lrad(NRADNM)

C   BNM CHANGED MANY OF THESE ARRAYS TO ADD NTSOA
      integer   lno   ,lno2  ,lo3  
      integer   lpan  ,lcres ,lpan2
      integer   lmpan ,lpbzn ,lnphe
      integer   lrno3 ,ldcb2 ,ldcb3
      integer   lhno4 ,lacet ,lald2
      integer   lalk1 ,lalk2 ,lalk3
      integer   lalk4 ,lalk5 ,laro1
      integer   laro2 ,lbacl ,lbald
      integer   lbcl1 ,lbcl2 ,lbuta
      integer   lccho ,lccrs  
      integer   lcl2  ,lco   ,lco2h
      integer   lco3h ,lcooh ,lcprm
      integer   ldcb1 ,leth  ,lethe
      integer   letoh ,lfcrs ,lfmcl
      integer   lform ,lfprm ,lgly 
      integer   lh2o2 ,lhc2h ,lhcho
      integer   lhcl  ,lhono ,lhno3
      integer   lho2h ,lhocl ,licl1
      integer   licl2 ,lisop ,lispd
      integer   lmek  ,lmeoh ,lmeth
      integer   lmgly ,lmvk  ,lna  
      integer   lnh3  ,lntr  ,lnxoy
      integer   lole  ,lole1 ,lole2
      integer   lbpin ,llimo ,lmono ,lsesq
      integer   lopen ,lpar  ,lpcl 
      integer   lpec  ,lphen ,lpna 
      integer   lpnh4 ,lpno3 ,lpoa 
      integer   lprod ,lpso4 ,lrc2h
      integer   lrc3h ,lrcho ,lrooh
      integer   lso2
      integer   lsulf
      integer   lterp ,ltol  ,lxn  ,lxyl 
      integer   lcpo1, lcpo2, lcpo3, lcpo4
      integer   lcpo5, lcpo6, lcpo7, lcpo8
      integer   lcoo1, lcoo2, lcoo3, lcoo4
      integer   lcoo5, lcoo6, lcoo7, lcoo8
      integer   lcbs1, lcbs2, lcbs3, lcbs4, lcbs5
      integer   lcas1, lcas2, lcas3, lcas4, lcas5
      integer   lcns1, lcns2, lcns3, lcns4
      integer   lcns5, lcns6, lcns7, lcns8
c
      integer   lapo1, lapo2, lapo3, lapo4
      integer   lapo5, lapo6, lapo7, lapo8
      integer   laoo1, laoo2, laoo3, laoo4
      integer   laoo5, laoo6, laoo7, laoo8
      integer   labs1, labs2, labs3, labs4, labs5
      integer   laas1, laas2, laas3, laas4, laas5
      integer   lans1, lans2, lans3, lans4
      integer   lans5, lans6, lans7, lans8
c
      integer   lapo1_1  ,lapo1_2  ,lapo1_3  ,lapo1_4  ,lapo1_5
      integer   lapo1_6  ,lapo1_7  ,lapo1_8  ,lapo1_9  ,lapo1_10
      integer   lapo2_1  ,lapo2_2  ,lapo2_3  ,lapo2_4  ,lapo2_5
      integer   lapo2_6  ,lapo2_7  ,lapo2_8  ,lapo2_9  ,lapo2_10
      integer   lapo3_1  ,lapo3_2  ,lapo3_3  ,lapo3_4  ,lapo3_5
      integer   lapo3_6  ,lapo3_7  ,lapo3_8  ,lapo3_9  ,lapo3_10
      integer   lapo4_1  ,lapo4_2  ,lapo4_3  ,lapo4_4  ,lapo4_5
      integer   lapo4_6  ,lapo4_7  ,lapo4_8  ,lapo4_9  ,lapo4_10
      integer   lapo5_1  ,lapo5_2  ,lapo5_3  ,lapo5_4  ,lapo5_5
      integer   lapo5_6  ,lapo5_7  ,lapo5_8  ,lapo5_9  ,lapo5_10
      integer   lapo6_1  ,lapo6_2  ,lapo6_3  ,lapo6_4  ,lapo6_5
      integer   lapo6_6  ,lapo6_7  ,lapo6_8  ,lapo6_9  ,lapo6_10
      integer   lapo7_1  ,lapo7_2  ,lapo7_3  ,lapo7_4  ,lapo7_5
      integer   lapo7_6  ,lapo7_7  ,lapo7_8  ,lapo7_9  ,lapo7_10
      integer   lapo8_1  ,lapo8_2  ,lapo8_3  ,lapo8_4  ,lapo8_5
      integer   lapo8_6  ,lapo8_7  ,lapo8_8  ,lapo8_9  ,lapo8_10

      integer   laoo1_1  ,laoo1_2  ,laoo1_3  ,laoo1_4  ,laoo1_5
      integer   laoo1_6  ,laoo1_7  ,laoo1_8  ,laoo1_9  ,laoo1_10
      integer   laoo2_1  ,laoo2_2  ,laoo2_3  ,laoo2_4  ,laoo2_5
      integer   laoo2_6  ,laoo2_7  ,laoo2_8  ,laoo2_9  ,laoo2_10
      integer   laoo3_1  ,laoo3_2  ,laoo3_3  ,laoo3_4  ,laoo3_5
      integer   laoo3_6  ,laoo3_7  ,laoo3_8  ,laoo3_9  ,laoo3_10
      integer   laoo4_1  ,laoo4_2  ,laoo4_3  ,laoo4_4  ,laoo4_5
      integer   laoo4_6  ,laoo4_7  ,laoo4_8  ,laoo4_9  ,laoo4_10
      integer   laoo5_1  ,laoo5_2  ,laoo5_3  ,laoo5_4  ,laoo5_5
      integer   laoo5_6  ,laoo5_7  ,laoo5_8  ,laoo5_9  ,laoo5_10
      integer   laoo6_1  ,laoo6_2  ,laoo6_3  ,laoo6_4  ,laoo6_5
      integer   laoo6_6  ,laoo6_7  ,laoo6_8  ,laoo6_9  ,laoo6_10
      integer   laoo7_1  ,laoo7_2  ,laoo7_3  ,laoo7_4  ,laoo7_5
      integer   laoo7_6  ,laoo7_7  ,laoo7_8  ,laoo7_9  ,laoo7_10
      integer   laoo8_1  ,laoo8_2  ,laoo8_3  ,laoo8_4  ,laoo8_5
      integer   laoo8_6  ,laoo8_7  ,laoo8_8  ,laoo8_9  ,laoo8_10

      integer   labs1_1  ,labs1_2  ,labs1_3  ,labs1_4  ,labs1_5
      integer   labs1_6  ,labs1_7  ,labs1_8  ,labs1_9  ,labs1_10
      integer   labs2_1  ,labs2_2  ,labs2_3  ,labs2_4  ,labs2_5
      integer   labs2_6  ,labs2_7  ,labs2_8  ,labs2_9  ,labs2_10
      integer   labs3_1  ,labs3_2  ,labs3_3  ,labs3_4  ,labs3_5
      integer   labs3_6  ,labs3_7  ,labs3_8  ,labs3_9  ,labs3_10
      integer   labs4_1  ,labs4_2  ,labs4_3  ,labs4_4  ,labs4_5
      integer   labs4_6  ,labs4_7  ,labs4_8  ,labs4_9  ,labs4_10
      integer   labs5_1  ,labs5_2  ,labs5_3  ,labs5_4  ,labs5_5
      integer   labs5_6  ,labs5_7  ,labs5_8  ,labs5_9  ,labs5_10

      integer   laas1_1  ,laas1_2  ,laas1_3  ,laas1_4  ,laas1_5
      integer   laas1_6  ,laas1_7  ,laas1_8  ,laas1_9  ,laas1_10
      integer   laas2_1  ,laas2_2  ,laas2_3  ,laas2_4  ,laas2_5
      integer   laas2_6  ,laas2_7  ,laas2_8  ,laas2_9  ,laas2_10
      integer   laas3_1  ,laas3_2  ,laas3_3  ,laas3_4  ,laas3_5
      integer   laas3_6  ,laas3_7  ,laas3_8  ,laas3_9  ,laas3_10
      integer   laas4_1  ,laas4_2  ,laas4_3  ,laas4_4  ,laas4_5
      integer   laas4_6  ,laas4_7  ,laas4_8  ,laas4_9  ,laas4_10
      integer   laas5_1  ,laas5_2  ,laas5_3  ,laas5_4  ,laas5_5
      integer   laas5_6  ,laas5_7  ,laas5_8  ,laas5_9  ,laas5_10

      integer   lans1_1  ,lans1_2  ,lans1_3  ,lans1_4  ,lans1_5
      integer   lans1_6  ,lans1_7  ,lans1_8  ,lans1_9  ,lans1_10
      integer   lans2_1  ,lans2_2  ,lans2_3  ,lans2_4  ,lans2_5
      integer   lans2_6  ,lans2_7  ,lans2_8  ,lans2_9  ,lans2_10
      integer   lans3_1  ,lans3_2  ,lans3_3  ,lans3_4  ,lans3_5
      integer   lans3_6  ,lans3_7  ,lans3_8  ,lans3_9  ,lans3_10
      integer   lans4_1  ,lans4_2  ,lans4_3  ,lans4_4  ,lans4_5
      integer   lans4_6  ,lans4_7  ,lans4_8  ,lans4_9  ,lans4_10
      integer   lans5_1  ,lans5_2  ,lans5_3  ,lans5_4  ,lans5_5
      integer   lans5_6  ,lans5_7  ,lans5_8  ,lans5_9  ,lans5_10
      integer   lans6_1  ,lans6_2  ,lans6_3  ,lans6_4  ,lans6_5
      integer   lans6_6  ,lans6_7  ,lans6_8  ,lans6_9  ,lans6_10
      integer   lans7_1  ,lans7_2  ,lans7_3  ,lans7_4  ,lans7_5
      integer   lans7_6  ,lans7_7  ,lans7_8  ,lans7_9  ,lans7_10
      integer   lans8_1  ,lans8_2  ,lans8_3  ,lans8_4  ,lans8_5
      integer   lans8_6  ,lans8_7  ,lans8_8  ,lans8_9  ,lans8_10
c
      integer   lpoc_1
      integer   lpoc_2   ,lpoc_3   ,lpoc_4
      integer   lpoc_5   ,lpoc_6   ,lpoc_7
      integer   lpoc_8   ,lpoc_9   ,lpoc_10
      integer   lpec_1   ,lpec_2   ,lpec_3
      integer   lpec_4   ,lpec_5   ,lpec_6
      integer   lpec_7   ,lpec_8   ,lpec_9
      integer   lpec_10  ,lcrust_1 ,lcrust_2
      integer   lcrust_3 ,lcrust_4 ,lcrust_5
      integer   lcrust_6 ,lcrust_7 ,lcrust_8
      integer   lcrust_9 ,lcrust_10,lph2o_1
      integer   lph2o_2  ,lph2o_3  ,lph2o_4
      integer   lph2o_5  ,lph2o_6  ,lph2o_7
      integer   lph2o_8  ,lph2o_9  ,lph2o_10
      integer   lpcl_1   ,lpcl_2   ,lpcl_3
      integer   lpcl_4   ,lpcl_5   ,lpcl_6
      integer   lpcl_7   ,lpcl_8   ,lpcl_9
      integer   lpcl_10  ,lna_1    ,lna_2
      integer   lna_3    ,lna_4    ,lna_5
      integer   lna_6    ,lna_7    ,lna_8
      integer   lna_9    ,lna_10   ,lpnh4_1
      integer   lpnh4_2  ,lpnh4_3  ,lpnh4_4
      integer   lpnh4_5  ,lpnh4_6  ,lpnh4_7
      integer   lpnh4_8  ,lpnh4_9  ,lpnh4_10
      integer   lpno3_1  ,lpno3_2  ,lpno3_3
      integer   lpno3_4  ,lpno3_5  ,lpno3_6
      integer   lpno3_7  ,lpno3_8  ,lpno3_9
      integer   lpno3_10 ,lpso4_1  ,lpso4_2
      integer   lpso4_3  ,lpso4_4  ,lpso4_5
      integer   lpso4_6  ,lpso4_7  ,lpso4_8
      integer   lpso4_9  ,lpso4_10 ,lph2o
c
      equivalence (lmap(1), lno  ), (lmap(2), lno2 ), (lmap(3), lo3  ),
     &            (lmap(4), lpan ), (lmap(5), lcres), (lmap(6), lpan2),
     &            (lmap(7), lmpan), (lmap(8), lpbzn), (lmap(9), lnphe),
     &            (lmap(10),lrno3), (lmap(11),ldcb2), (lmap(12),ldcb3),
     &            (lmap(13),lhno4), (lmap(14),lacet), (lmap(15),lald2),
     &            (lmap(16),lalk1), (lmap(17),lalk2), (lmap(18),lalk3),
     &            (lmap(19),lalk4), (lmap(20),lalk5), (lmap(21),laro1),
     &            (lmap(22),laro2), (lmap(23),lbacl), (lmap(24),lbald),
     &            (lmap(25),lbcl1), (lmap(26),lbcl2), (lmap(27),lbuta),
     &            (lmap(28),lccho), (lmap(29),lccrs), 
     &		  (lmap(30),lcpo1), (lmap(31),lcpo2), (lmap(32),lcpo3), 
     &		  (lmap(33),lcpo4), (lmap(34),lcpo5), (lmap(35),lcpo6), 
     &		  (lmap(36),lcpo7), (lmap(37),lcpo8),
     &            (lmap(38),lcoo1), (lmap(39),lcoo2), (lmap(40),lcoo3),
     &            (lmap(41),lcoo4), (lmap(42),lcoo5), (lmap(43),lcoo6),
     &            (lmap(44),lcoo7), (lmap(45),lcoo8),
     &            (lmap(46),lcbs1), (lmap(47),lcbs2), (lmap(48),lcbs3), 
     &		  (lmap(49),lcbs4), (lmap(50),lcbs5),
     &		  (lmap(51),lcas1), (lmap(52),lcas2), (lmap(53),lcas3), 
     &		  (lmap(54),lcas4), (lmap(55),lcas5),
     &            (lmap(56),lcns1), (lmap(57),lcns2), (lmap(58),lcns3),
     &            (lmap(59),lcns4), (lmap(60),lcns5), (lmap(61),lcns6),
     &            (lmap(62),lcns7), (lmap(63),lcns8),
c
     &           (lmap(64),lcl2 ), (lmap(65),lco  ), (lmap(66),lco2h),
     &           (lmap(67),lco3h), (lmap(68),lcooh), (lmap(69),lcprm),
     &           (lmap(70),ldcb1), (lmap(71),leth ), (lmap(72),lethe),
     &           (lmap(73),letoh), (lmap(74),lfcrs), (lmap(75),lfmcl),
     &           (lmap(76),lform), (lmap(77),lfprm), (lmap(78),lgly ),
     &           (lmap(79),lh2o2), (lmap(80),lhc2h), (lmap(81),lhcho),
     &           (lmap(82),lhcl ), (lmap(83),lhono), (lmap(84),lhno3),
     &           (lmap(85),lho2h), (lmap(86),lhocl), (lmap(87),licl1),
     &           (lmap(88),licl2), (lmap(89),lisop), (lmap(90),lispd),
     &           (lmap(91),lmek ), (lmap(92),lmeoh), (lmap(93),lmeth),
     &           (lmap(94),lmgly), (lmap(95),lmvk ), (lmap(96),lna  ),
     &           (lmap(97),lnh3 ), (lmap(98),lntr ), (lmap(99),lnxoy),
     &           (lmap(100),lole ),(lmap(101),lole1),(lmap(102),lole2),
     &           (lmap(103),lbpin),(lmap(104),llimo),(lmap(105),lmono),
     &           (lmap(106),lsesq),(lmap(107),lopen),(lmap(108),lpar ),
     &           (lmap(109),lpcl ),(lmap(110),lpec ),(lmap(111),lphen),
     &           (lmap(112),lpna ),(lmap(113),lpnh4),(lmap(114),lpno3),
     &           (lmap(115),lpoa ),(lmap(116),lprod),(lmap(117),lpso4),
     &           (lmap(118),lrc2h),(lmap(119),lrc3h),(lmap(120),lrcho),
     &           (lmap(121),lrooh),(lmap(122),lso2 ),
     &		 (lmap(123),lapo1),(lmap(124),lapo2),(lmap(125),lapo3),
     &           (lmap(126),lapo4),(lmap(127),lapo5),(lmap(128),lapo6),
     &		 (lmap(129),lapo7),(lmap(130),lapo8),
     &           (lmap(131),laoo1),(lmap(132),laoo2),(lmap(133),laoo3),
     &           (lmap(134),laoo4),(lmap(135),laoo5),(lmap(136),laoo6),
     &           (lmap(137),laoo7),(lmap(138),laoo8),
     &           (lmap(139),labs1),(lmap(140),labs2),(lmap(141),labs3), 
     &		 (lmap(142),labs4),(lmap(143),labs5),
     &		 (lmap(144),laas1),(lmap(145),laas2),(lmap(146),laas3),
     &		 (lmap(147),laas4),(lmap(148),laas5),
     &		 (lmap(149),lans1),(lmap(150),lans2),(lmap(151),lans3),
     &           (lmap(152),lans4),(lmap(153),lans5),(lmap(154),lans6),
     &		 (lmap(155),lans7),(lmap(156),lans8),
     &           (lmap(157),lsulf), (lmap(158),lterp),
     &           (lmap(159),ltol ), (lmap(160),lxn  ),(lmap(161),lxyl ), 
c
     &(lmap(162),lapo1_1 ), (lmap(163),lapo1_2 ), (lmap(164),lapo1_3 ),
     &(lmap(165),lapo1_4 ), (lmap(166),lapo1_5 ), (lmap(167),lapo1_6 ),
     &(lmap(168),lapo1_7 ), (lmap(169),lapo1_8 ), (lmap(170),lapo1_9 ),
     &(lmap(171),lapo1_10 ),(lmap(172),lapo2_1 ), (lmap(173),lapo2_2 ),
     &(lmap(174),lapo2_3 ), (lmap(175),lapo2_4 ), (lmap(176),lapo2_5 ),
     &(lmap(177),lapo2_6 ), (lmap(178),lapo2_7 ), (lmap(179),lapo2_8 ),
     &(lmap(180),lapo2_9 ), (lmap(181),lapo2_10 ),(lmap(182),lapo3_1 ),
     &(lmap(183),lapo3_2 ), (lmap(184),lapo3_3 ), (lmap(185),lapo3_4 ),
     &(lmap(186),lapo3_5 ), (lmap(187),lapo3_6 ), (lmap(188),lapo3_7 ),
     &(lmap(189),lapo3_8 ), (lmap(190),lapo3_9 ), (lmap(191),lapo3_10 ),
     &(lmap(192),lapo4_1 ), (lmap(193),lapo4_2 ), (lmap(194),lapo4_3 ),
     &(lmap(195),lapo4_4 ), (lmap(196),lapo4_5 ), (lmap(197),lapo4_6 ),
     &(lmap(198),lapo4_7 ), (lmap(199),lapo4_8 ), (lmap(200),lapo4_9 ),
     &(lmap(201),lapo4_10 ),(lmap(202),lapo5_1 ), (lmap(203),lapo5_2 ),
     &(lmap(204),lapo5_3 ), (lmap(205),lapo5_4 ), (lmap(206),lapo5_5 ),
     &(lmap(207),lapo5_6 ), (lmap(208),lapo5_7 ), (lmap(209),lapo5_8 ),
     &(lmap(210),lapo5_9 ), (lmap(211),lapo5_10 ),(lmap(212),lapo6_1 ),
     &(lmap(213),lapo6_2 ), (lmap(214),lapo6_3 ), (lmap(215),lapo6_4 ),
     &(lmap(216),lapo6_5 ), (lmap(217),lapo6_6 ), (lmap(218),lapo6_7 ),
     &(lmap(219),lapo6_8 ), (lmap(220),lapo6_9 ), (lmap(221),lapo6_10 ),
     &(lmap(222),lapo7_1 ), (lmap(223),lapo7_2 ), (lmap(224),lapo7_3 ),
     &(lmap(225),lapo7_4 ), (lmap(226),lapo7_5 ), (lmap(227),lapo7_6 ),
     &(lmap(228),lapo7_7 ), (lmap(229),lapo7_8 ), (lmap(230),lapo7_9 ),
     &(lmap(231),lapo7_10 ),(lmap(232),lapo8_1 ), (lmap(233),lapo8_2 ),
     &(lmap(234),lapo8_3 ), (lmap(235),lapo8_4 ), (lmap(236),lapo8_5 ),
     &(lmap(237),lapo8_6 ), (lmap(238),lapo8_7 ), (lmap(239),lapo8_8 ),
     &(lmap(240),lapo8_9 ), (lmap(241),lapo8_10 ),
     &(lmap(242),laoo1_1 ), (lmap(243),laoo1_2 ),
     &(lmap(244),laoo1_3 ), (lmap(245),laoo1_4 ), (lmap(246),laoo1_5 ),
     &(lmap(247),laoo1_6 ), (lmap(248),laoo1_7 ), (lmap(249),laoo1_8 ),
     &(lmap(250),laoo1_9 ), (lmap(251),laoo1_10 ),(lmap(252),laoo2_1 ),
     &(lmap(253),laoo2_2 ), (lmap(254),laoo2_3 ), (lmap(255),laoo2_4 ),
     &(lmap(256),laoo2_5 ), (lmap(257),laoo2_6 ), (lmap(258),laoo2_7 ),
     &(lmap(259),laoo2_8 ), (lmap(260),laoo2_9 ), (lmap(261),laoo2_10 ),
     &(lmap(262),laoo3_1 ), (lmap(263),laoo3_2 ), (lmap(264),laoo3_3 ),
     &(lmap(265),laoo3_4 ), (lmap(266),laoo3_5 ), (lmap(267),laoo3_6 ),
     &(lmap(268),laoo3_7 ), (lmap(269),laoo3_8 ), (lmap(270),laoo3_9 ),
     &(lmap(271),laoo3_10 ),(lmap(272),laoo4_1 ), (lmap(273),laoo4_2 ),
     &(lmap(274),laoo4_3 ), (lmap(275),laoo4_4 ), (lmap(276),laoo4_5 ),
     &(lmap(277),laoo4_6 ), (lmap(278),laoo4_7 ), (lmap(279),laoo4_8 ),
     &(lmap(280),laoo4_9 ), (lmap(281),laoo4_10 ),(lmap(282),laoo5_1 ),
     &(lmap(283),laoo5_2 ), (lmap(284),laoo5_3 ), (lmap(285),laoo5_4 ),
     &(lmap(286),laoo5_5 ), (lmap(287),laoo5_6 ), (lmap(288),laoo5_7 ),
     &(lmap(289),laoo5_8 ), (lmap(290),laoo5_9 ), (lmap(291),laoo5_10 ),
     &(lmap(292),laoo6_1 ), (lmap(293),laoo6_2 ), (lmap(294),laoo6_3 ),
     &(lmap(295),laoo6_4 ), (lmap(296),laoo6_5 ), (lmap(297),laoo6_6 ),
     &(lmap(298),laoo6_7 ), (lmap(299),laoo6_8 ), (lmap(300),laoo6_9 ),
     &(lmap(301),laoo6_10 ),(lmap(302),laoo7_1 ), (lmap(303),laoo7_2 ),
     &(lmap(304),laoo7_3 ), (lmap(305),laoo7_4 ), (lmap(306),laoo7_5 ),
     &(lmap(307),laoo7_6 ), (lmap(308),laoo7_7 ), (lmap(309),laoo7_8 ),
     &(lmap(310),laoo7_9 ), (lmap(311),laoo7_10 ),(lmap(312),laoo8_1 ),
     &(lmap(313),laoo8_2 ), (lmap(314),laoo8_3 ), (lmap(315),laoo8_4 ),
     &(lmap(316),laoo8_5 ), (lmap(317),laoo8_6 ), (lmap(318),laoo8_7 ),
     &(lmap(319),laoo8_8 ), (lmap(320),laoo8_9 ), (lmap(321),laoo8_10 ),
     &(lmap(322),labs1_1 ),
     &(lmap(323),labs1_2 ), (lmap(324),labs1_3 ), (lmap(325),labs1_4 ),
     &(lmap(326),labs1_5 ), (lmap(327),labs1_6 ), (lmap(328),labs1_7 ),
     &(lmap(329),labs1_8 ), (lmap(330),labs1_9 ), (lmap(331),labs1_10 ),
     &(lmap(332),labs2_1 ), (lmap(333),labs2_2 ), (lmap(334),labs2_3 ),
     &(lmap(335),labs2_4 ), (lmap(336),labs2_5 ), (lmap(337),labs2_6 ),
     &(lmap(338),labs2_7 ), (lmap(339),labs2_8 ), (lmap(340),labs2_9 ),
     &(lmap(341),labs2_10 ),(lmap(342),labs3_1 ), (lmap(343),labs3_2 ),
     &(lmap(344),labs3_3 ), (lmap(345),labs3_4 ), (lmap(346),labs3_5 ),
     &(lmap(347),labs3_6 ), (lmap(348),labs3_7 ), (lmap(349),labs3_8 ),
     &(lmap(350),labs3_9 ), (lmap(351),labs3_10 ),(lmap(352),labs4_1 ),
     &(lmap(353),labs4_2 ), (lmap(354),labs4_3 ), (lmap(355),labs4_4 ),
     &(lmap(356),labs4_5 ), (lmap(357),labs4_6 ), (lmap(358),labs4_7 ),
     &(lmap(359),labs4_8 ), (lmap(360),labs4_9 ), (lmap(361),labs4_10 ),
     &(lmap(362),labs5_1 ),
     &(lmap(363),labs5_2 ), (lmap(364),labs5_3 ), (lmap(365),labs5_4 ),
     &(lmap(366),labs5_5 ), (lmap(367),labs5_6 ), (lmap(368),labs5_7 ),
     &(lmap(369),labs5_8 ), (lmap(370),labs5_9 ), (lmap(371),labs5_10 ),
     &(lmap(372),laas1_1 ), (lmap(373),laas1_2 ), (lmap(374),laas1_3 ),
     &(lmap(375),laas1_4 ), (lmap(376),laas1_5 ), (lmap(377),laas1_6 ),
     &(lmap(378),laas1_7 ), (lmap(379),laas1_8 ), (lmap(380),laas1_9 ),
     &(lmap(381),laas1_10 ),(lmap(382),laas2_1 ), (lmap(383),laas2_2 ),
     &(lmap(384),laas2_3 ), (lmap(385),laas2_4 ), (lmap(386),laas2_5 ),
     &(lmap(387),laas2_6 ), (lmap(388),laas2_7 ), (lmap(389),laas2_8 ),
     &(lmap(390),laas2_9 ), (lmap(391),laas2_10 ),(lmap(392),laas3_1 ),
     &(lmap(393),laas3_2 ), (lmap(394),laas3_3 ), (lmap(395),laas3_4 ),
     &(lmap(396),laas3_5 ), (lmap(397),laas3_6 ), (lmap(398),laas3_7 ),
     &(lmap(399),laas3_8 ), (lmap(400),laas3_9 ), (lmap(401),laas3_10 ),
     &(lmap(402),laas4_1 ), (lmap(403),laas4_2 ), (lmap(404),laas4_3 ),
     &(lmap(405),laas4_4 ), (lmap(406),laas4_5 ), (lmap(407),laas4_6 ),
     &(lmap(408),laas4_7 ), (lmap(409),laas4_8 ), (lmap(410),laas4_9 ),
     &(lmap(411),laas4_10 ),
     &(lmap(412),laas5_1 ), (lmap(413),laas5_2 ), (lmap(414),laas5_3 ),
     &(lmap(415),laas5_4 ), (lmap(416),laas5_5 ), (lmap(417),laas5_6 ),
     &(lmap(418),laas5_7 ), (lmap(419),laas5_8 ), (lmap(420),laas5_9 ),
     &(lmap(421),laas5_10 ),
     &(lmap(422),lans1_1 ), (lmap(423),lans1_2 ), (lmap(424),lans1_3 ),
     &(lmap(425),lans1_4 ), (lmap(426),lans1_5 ), (lmap(427),lans1_6 ),
     &(lmap(428),lans1_7 ), (lmap(429),lans1_8 ), (lmap(430),lans1_9 ),
     &(lmap(431),lans1_10 ),(lmap(432),lans2_1 ), (lmap(433),lans2_2 ),
     &(lmap(434),lans2_3 ), (lmap(435),lans2_4 ), (lmap(436),lans2_5 ),
     &(lmap(437),lans2_6 ), (lmap(438),lans2_7 ), (lmap(439),lans2_8 ),
     &(lmap(440),lans2_9 ), (lmap(441),lans2_10 ),(lmap(442),lans3_1 ),
     &(lmap(443),lans3_2 ), (lmap(444),lans3_3 ), (lmap(445),lans3_4 ),
     &(lmap(446),lans3_5 ), (lmap(447),lans3_6 ), (lmap(448),lans3_7 ),
     &(lmap(449),lans3_8 ), (lmap(450),lans3_9 ), (lmap(451),lans3_10 ),
     &(lmap(452),lans4_1 ), (lmap(453),lans4_2 ), (lmap(454),lans4_3 ),
     &(lmap(455),lans4_4 ), (lmap(456),lans4_5 ), (lmap(457),lans4_6 ),
     &(lmap(458),lans4_7 ), (lmap(459),lans4_8 ), (lmap(460),lans4_9 ),
     &(lmap(461),lans4_10 ),(lmap(462),lans5_1 ), (lmap(463),lans5_2 ),
     &(lmap(464),lans5_3 ), (lmap(465),lans5_4 ), (lmap(466),lans5_5 ),
     &(lmap(467),lans5_6 ), (lmap(468),lans5_7 ), (lmap(469),lans5_8 ),
     &(lmap(470),lans5_9 ), (lmap(471),lans5_10 ),(lmap(472),lans6_1 ),
     &(lmap(473),lans6_2 ), (lmap(474),lans6_3 ), (lmap(475),lans6_4 ),
     &(lmap(476),lans6_5 ), (lmap(477),lans6_6 ), (lmap(478),lans6_7 ),
     &(lmap(479),lans6_8 ), (lmap(480),lans6_9 ), (lmap(481),lans6_10 ),
     &(lmap(482),lans7_1 ), (lmap(483),lans7_2 ), (lmap(484),lans7_3 ),
     &(lmap(485),lans7_4 ), (lmap(486),lans7_5 ), (lmap(487),lans7_6 ),
     &(lmap(488),lans7_7 ), (lmap(489),lans7_8 ), (lmap(490),lans7_9 ),
     &(lmap(491),lans7_10 ),(lmap(492),lans8_1 ), (lmap(493),lans8_2 ),
     &(lmap(494),lans8_3 ), (lmap(495),lans8_4 ), (lmap(496),lans8_5 ),
     &(lmap(497),lans8_6 ), (lmap(498),lans8_7 ), (lmap(499),lans8_8 ),
     &(lmap(500),lans8_9 ), (lmap(501),lans8_10 ),

c   scf

     &(lmap(502),lpoc_1  ),
     &(lmap(503),lpoc_2  ),(lmap(504),lpoc_3  ), (lmap(505),lpoc_4  ),
     &(lmap(506),lpoc_5  ),(lmap(507),lpoc_6  ), (lmap(508),lpoc_7  ),
     &(lmap(509),lpoc_8  ),(lmap(510),lpoc_9  ), (lmap(511),lpoc_10 ),
     &(lmap(512),lpec_1  ),(lmap(513),lpec_2  ), (lmap(514),lpec_3  ),
     &(lmap(515),lpec_4  ),(lmap(516),lpec_5  ), (lmap(517),lpec_6  ),
     &(lmap(518),lpec_7  ),(lmap(519),lpec_8  ), (lmap(520),lpec_9  ),
     &(lmap(521),lpec_10 ),(lmap(522),lcrust_1), (lmap(523),lcrust_2),
     &(lmap(524),lcrust_3),(lmap(525),lcrust_4), (lmap(526),lcrust_5),
     &(lmap(527),lcrust_6),(lmap(528),lcrust_7), (lmap(529),lcrust_8),
     &(lmap(530),lcrust_9),(lmap(531),lcrust_10),(lmap(532),lph2o_1),
     &(lmap(533),lph2o_2 ),(lmap(534),lph2o_3 ), (lmap(535),lph2o_4 ),
     &(lmap(536),lph2o_5 ),(lmap(537),lph2o_6 ), (lmap(538),lph2o_7 ),
     &(lmap(539),lph2o_8 ),(lmap(540),lph2o_9 ), (lmap(541),lph2o_10),
     &(lmap(542),lpcl_1  ),(lmap(543),lpcl_2  ), (lmap(544),lpcl_3  ),
     &(lmap(545),lpcl_4  ),(lmap(546),lpcl_5  ), (lmap(547),lpcl_6  ),
     &(lmap(548),lpcl_7  ),(lmap(549),lpcl_8  ), (lmap(550),lpcl_9  ),
     &(lmap(551),lpcl_10 ),(lmap(552),lna_1   ), (lmap(553),lna_2   ),
     &(lmap(554),lna_3   ),(lmap(555),lna_4   ), (lmap(556),lna_5   ),
     &(lmap(557),lna_6   ),(lmap(558),lna_7   ), (lmap(559),lna_8   ),
     &(lmap(560),lna_9   ),(lmap(561),lna_10  ), (lmap(562),lpnh4_1 ),
     &(lmap(563),lpnh4_2 ),(lmap(564),lpnh4_3 ), (lmap(565),lpnh4_4 ),
     &(lmap(566),lpnh4_5 ),(lmap(567),lpnh4_6 ), (lmap(568),lpnh4_7 ),
     &(lmap(569),lpnh4_8 ),(lmap(570),lpnh4_9 ), (lmap(571),lpnh4_10),
     &(lmap(572),lpno3_1 ),(lmap(573),lpno3_2 ), (lmap(574),lpno3_3 ),
     &(lmap(575),lpno3_4 ),(lmap(576),lpno3_5 ), (lmap(577),lpno3_6 ),
     &(lmap(578),lpno3_7 ),(lmap(579),lpno3_8 ), (lmap(580),lpno3_9 ),
     &(lmap(581),lpno3_10),(lmap(582),lpso4_1 ), (lmap(583),lpso4_2 ),
     &(lmap(584),lpso4_3 ),(lmap(585),lpso4_4 ), (lmap(586),lpso4_5 ),
     &(lmap(587),lpso4_6 ),(lmap(588),lpso4_7 ), (lmap(589),lpso4_8 ),
     &(lmap(590),lpso4_9 ),(lmap(591),lpso4_10), (lmap(592),lph2o   )
c
      integer   lo1d  ,lo    ,lclo 
      integer   lcl   ,ln2o5 ,lno3 
      integer   loh   ,lho2  ,lc2o3
      integer   lxo2  ,lxo2n ,lto2 
      integer   lror  ,lcro  ,lro2r
      integer   lr2o2 ,lro2n ,lcco3
      integer   lrco3 ,lmco3 ,lbzco
      integer   lcxo2 ,lhco3 ,ltbuo
      integer   lbzo  ,lbzno
c
      equivalence (lrad(1), lo1d ), (lrad(2), lo   ), (lrad(3), lclo ),
     &            (lrad(4), lcl  ), (lrad(5), ln2o5), (lrad(6), lno3 ),
     &            (lrad(7), loh  ), (lrad(8), lho2 ), (lrad(9), lc2o3),
     &            (lrad(10),lxo2 ), (lrad(11),lxo2n), (lrad(12),lto2 ),
     &            (lrad(13),lror ), (lrad(14),lcro ), (lrad(15),lro2r),
     &            (lrad(16),lr2o2), (lrad(17),lro2n), (lrad(18),lcco3),
     &            (lrad(19),lrco3), (lrad(20),lmco3), (lrad(21),lbzco),
     &            (lrad(22),lcxo2), (lrad(23),lhco3), (lrad(24),ltbuo),
     &            (lrad(25),lbzo ), (lrad(26),lbzno)

