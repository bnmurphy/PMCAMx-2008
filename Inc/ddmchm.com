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
      integer   lcpo9, lcpo10
      integer   lcoo1, lcoo2, lcoo3, lcoo4
      integer   lcoo5, lcoo6, lcoo7, lcoo8
      integer   lcoo9, lcoo10
      integer   lcbs1, lcbs2, lcbs3, lcbs4
      integer   lcas1, lcas2, lcas3, lcas4
c
      integer   lapo1, lapo2, lapo3, lapo4
      integer   lapo5, lapo6, lapo7, lapo8
      integer   lapo9, lapo10
      integer   laoo1, laoo2, laoo3, laoo4
      integer   laoo5, laoo6, laoo7, laoo8
      integer   laoo9, laoo10
      integer   labs1, labs2, labs3, labs4
      integer   laas1, laas2, laas3, laas4
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
      integer   lapo9_1  ,lapo9_2  ,lapo9_3  ,lapo9_4  ,lapo9_5
      integer   lapo9_6  ,lapo9_7  ,lapo9_8  ,lapo9_9  ,lapo9_10
      integer   lapo10_1  ,lapo10_2  ,lapo10_3  ,lapo10_4  ,lapo10_5
      integer   lapo10_6  ,lapo10_7  ,lapo10_8  ,lapo10_9  ,lapo10_10

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
      integer   laoo9_1  ,laoo9_2  ,laoo9_3  ,laoo9_4  ,laoo9_5
      integer   laoo9_6  ,laoo9_7  ,laoo9_8  ,laoo9_9  ,laoo9_10
      integer   laoo10_1  ,laoo10_2  ,laoo10_3  ,laoo10_4  ,laoo10_5
      integer   laoo10_6  ,laoo10_7  ,laoo10_8  ,laoo10_9  ,laoo10_10

      integer   labs1_1  ,labs1_2  ,labs1_3  ,labs1_4  ,labs1_5
      integer   labs1_6  ,labs1_7  ,labs1_8  ,labs1_9  ,labs1_10
      integer   labs2_1  ,labs2_2  ,labs2_3  ,labs2_4  ,labs2_5
      integer   labs2_6  ,labs2_7  ,labs2_8  ,labs2_9  ,labs2_10
      integer   labs3_1  ,labs3_2  ,labs3_3  ,labs3_4  ,labs3_5
      integer   labs3_6  ,labs3_7  ,labs3_8  ,labs3_9  ,labs3_10
      integer   labs4_1  ,labs4_2  ,labs4_3  ,labs4_4  ,labs4_5
      integer   labs4_6  ,labs4_7  ,labs4_8  ,labs4_9  ,labs4_10

      integer   laas1_1  ,laas1_2  ,laas1_3  ,laas1_4  ,laas1_5
      integer   laas1_6  ,laas1_7  ,laas1_8  ,laas1_9  ,laas1_10
      integer   laas2_1  ,laas2_2  ,laas2_3  ,laas2_4  ,laas2_5
      integer   laas2_6  ,laas2_7  ,laas2_8  ,laas2_9  ,laas2_10
      integer   laas3_1  ,laas3_2  ,laas3_3  ,laas3_4  ,laas3_5
      integer   laas3_6  ,laas3_7  ,laas3_8  ,laas3_9  ,laas3_10
      integer   laas4_1  ,laas4_2  ,laas4_3  ,laas4_4  ,laas4_5
      integer   laas4_6  ,laas4_7  ,laas4_8  ,laas4_9  ,laas4_10

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
     &            (lmap(28),lccho), (lmap(29),lccrs), (lmap(30),lcpo1),
     &            (lmap(31),lcpo2), (lmap(32),lcpo3), (lmap(33),lcpo4),
     &            (lmap(34),lcpo5), (lmap(35),lcpo6), (lmap(36),lcpo7),
     &            (lmap(37),lcpo8), (lmap(38),lcpo9), (lmap(39),lcpo10),
     &            (lmap(40),lcoo1), (lmap(41),lcoo2), (lmap(42),lcoo3),
     &            (lmap(43),lcoo4), (lmap(44),lcoo5), (lmap(45),lcoo6),
     &            (lmap(46),lcoo7), (lmap(47),lcoo8), (lmap(48),lcoo9),
     &            (lmap(49),lcoo10), (lmap(50),lcbs1), (lmap(51),lcbs2),
     &            (lmap(52),lcbs3), (lmap(53),lcbs4), (lmap(54),lcas1),
     &            (lmap(55),lcas2), (lmap(56),lcas3), (lmap(57),lcas4),
c
     &           (lmap(58),lcl2 ), (lmap(59),lco  ), (lmap(60),lco2h),
     &           (lmap(61),lco3h), (lmap(62),lcooh), (lmap(63),lcprm),
     &           (lmap(64),ldcb1), (lmap(65),leth ), (lmap(66),lethe),
     &           (lmap(67),letoh), (lmap(68),lfcrs), (lmap(69),lfmcl),
     &           (lmap(70),lform), (lmap(71),lfprm), (lmap(72),lgly ),
     &           (lmap(73),lh2o2), (lmap(74),lhc2h), (lmap(75),lhcho),
     &           (lmap(76),lhcl ), (lmap(77),lhono), (lmap(78),lhno3),
     &           (lmap(79),lho2h), (lmap(80),lhocl), (lmap(81),licl1),
     &           (lmap(82),licl2), (lmap(83),lisop), (lmap(84),lispd),
     &           (lmap(85),lmek ),(lmap(86),lmeoh),(lmap(87),lmeth),
     &           (lmap(88),lmgly),(lmap(89),lmvk ),(lmap(90),lna  ),
     &           (lmap(91),lnh3 ),(lmap(92),lntr ),(lmap(93),lnxoy),
     &           (lmap(94),lole ),(lmap(95),lole1),(lmap(96),lole2),
     &           (lmap(97),lbpin),(lmap(98),llimo),(lmap(99),lmono),
     &           (lmap(100),lsesq),(lmap(101),lopen),(lmap(102),lpar ),
     &           (lmap(103),lpcl ),(lmap(104),lpec ),(lmap(105),lphen),
     &           (lmap(106),lpna ),(lmap(107),lpnh4),(lmap(108),lpno3),
     &           (lmap(109),lpoa ),(lmap(110),lprod),(lmap(111),lpso4),
     &           (lmap(112),lrc2h),(lmap(113),lrc3h),(lmap(114),lrcho),
     &           (lmap(115),lrooh),(lmap(116),lso2 ),(lmap(117),lapo1),
     &           (lmap(118),lapo2),(lmap(119),lapo3),(lmap(120),lapo4),
     &           (lmap(121),lapo5),(lmap(122),lapo6),(lmap(123),lapo7),
     &           (lmap(124),lapo8),(lmap(125),lapo9),(lmap(126),lapo10),
     &           (lmap(127),laoo1),(lmap(128),laoo2),(lmap(129),laoo3),
     &           (lmap(130),laoo4),(lmap(131),laoo5),(lmap(132),laoo6),
     &           (lmap(133),laoo7),(lmap(134),laoo8),(lmap(135),laoo9),
     &           (lmap(136),laoo10),(lmap(137),labs1),(lmap(138),labs2),
     &           (lmap(139),labs3),(lmap(140),labs4),(lmap(141),laas1),
     &           (lmap(142),laas2),(lmap(143),laas3),(lmap(144),laas4),
     &           (lmap(145),lsulf),(lmap(146),lterp),
     &           (lmap(147),ltol ),(lmap(148),lxn  ),(lmap(149),lxyl ), 
c
     &(lmap(150),lapo1_1 ),(lmap(151),lapo1_2 ),(lmap(152),lapo1_3 ),
     &(lmap(153),lapo1_4 ),(lmap(154),lapo1_5 ),(lmap(155),lapo1_6 ),
     &(lmap(156),lapo1_7 ),(lmap(157),lapo1_8 ),(lmap(158),lapo1_9 ),
     &(lmap(159),lapo1_10 ),(lmap(160),lapo2_1 ),(lmap(161),lapo2_2 ),
     &(lmap(162),lapo2_3 ),(lmap(163),lapo2_4 ),(lmap(164),lapo2_5 ),
     &(lmap(165),lapo2_6 ),(lmap(166),lapo2_7 ),(lmap(167),lapo2_8 ),
     &(lmap(168),lapo2_9 ),(lmap(169),lapo2_10 ),(lmap(170),lapo3_1 ),
     &(lmap(171),lapo3_2 ),(lmap(172),lapo3_3 ),(lmap(173),lapo3_4 ),
     &(lmap(174),lapo3_5 ),(lmap(175),lapo3_6 ),(lmap(176),lapo3_7 ),
     &(lmap(177),lapo3_8 ),(lmap(178),lapo3_9 ),(lmap(179),lapo3_10 ),
     &(lmap(180),lapo4_1 ),(lmap(181),lapo4_2 ),(lmap(182),lapo4_3 ),
     &(lmap(183),lapo4_4 ),(lmap(184),lapo4_5 ),(lmap(185),lapo4_6 ),
     &(lmap(186),lapo4_7 ),(lmap(187),lapo4_8 ),(lmap(188),lapo4_9 ),
     &(lmap(189),lapo4_10 ),(lmap(190),lapo5_1 ),(lmap(191),lapo5_2 ),
     &(lmap(192),lapo5_3 ),(lmap(193),lapo5_4 ),(lmap(194),lapo5_5 ),
     &(lmap(195),lapo5_6 ),(lmap(196),lapo5_7 ),(lmap(197),lapo5_8 ),
     &(lmap(198),lapo5_9 ),(lmap(199),lapo5_10 ),(lmap(200),lapo6_1 ),
     &(lmap(201),lapo6_2 ),(lmap(202),lapo6_3 ),(lmap(203),lapo6_4 ),
     &(lmap(204),lapo6_5 ),(lmap(205),lapo6_6 ),(lmap(206),lapo6_7 ),
     &(lmap(207),lapo6_8 ),(lmap(208),lapo6_9 ),(lmap(209),lapo6_10 ),
     &(lmap(210),lapo7_1 ),(lmap(211),lapo7_2 ),(lmap(212),lapo7_3 ),
     &(lmap(213),lapo7_4 ),(lmap(214),lapo7_5 ),(lmap(215),lapo7_6 ),
     &(lmap(216),lapo7_7 ),(lmap(217),lapo7_8 ),(lmap(218),lapo7_9 ),
     &(lmap(219),lapo7_10 ),(lmap(220),lapo8_1 ),(lmap(221),lapo8_2 ),
     &(lmap(222),lapo8_3 ),(lmap(223),lapo8_4 ),(lmap(224),lapo8_5 ),
     &(lmap(225),lapo8_6 ),(lmap(226),lapo8_7 ),(lmap(227),lapo8_8 ),
     &(lmap(228),lapo8_9 ),(lmap(229),lapo8_10 ),(lmap(230),lapo9_1 ),
     &(lmap(231),lapo9_2 ),(lmap(232),lapo9_3 ),(lmap(233),lapo9_4 ),
     &(lmap(234),lapo9_5 ),(lmap(235),lapo9_6 ),(lmap(236),lapo9_7 ),
     &(lmap(237),lapo9_8 ),(lmap(238),lapo9_9 ),(lmap(239),lapo9_10 ),
     &(lmap(240),lapo10_1 ),(lmap(241),lapo10_2 ),(lmap(242),lapo10_3 ),
     &(lmap(243),lapo10_4 ),(lmap(244),lapo10_5 ),(lmap(245),lapo10_6 ),
     &(lmap(246),lapo10_7 ),(lmap(247),lapo10_8 ),(lmap(248),lapo10_9 ),
     &(lmap(249),lapo10_10 ),(lmap(250),laoo1_1 ),(lmap(251),laoo1_2 ),
     &(lmap(252),laoo1_3 ),(lmap(253),laoo1_4 ),(lmap(254),laoo1_5 ),
     &(lmap(255),laoo1_6 ),(lmap(256),laoo1_7 ),(lmap(257),laoo1_8 ),
     &(lmap(258),laoo1_9 ),(lmap(259),laoo1_10 ),(lmap(260),laoo2_1 ),
     &(lmap(261),laoo2_2 ),(lmap(262),laoo2_3 ),(lmap(263),laoo2_4 ),
     &(lmap(264),laoo2_5 ),(lmap(265),laoo2_6 ),(lmap(266),laoo2_7 ),
     &(lmap(267),laoo2_8 ),(lmap(268),laoo2_9 ),(lmap(269),laoo2_10 ),
     &(lmap(270),laoo3_1 ),(lmap(271),laoo3_2 ),(lmap(272),laoo3_3 ),
     &(lmap(273),laoo3_4 ),(lmap(274),laoo3_5 ),(lmap(275),laoo3_6 ),
     &(lmap(276),laoo3_7 ),(lmap(277),laoo3_8 ),(lmap(278),laoo3_9 ),
     &(lmap(279),laoo3_10 ),(lmap(280),laoo4_1 ),(lmap(281),laoo4_2 ),
     &(lmap(282),laoo4_3 ),(lmap(283),laoo4_4 ),(lmap(284),laoo4_5 ),
     &(lmap(285),laoo4_6 ),(lmap(286),laoo4_7 ),(lmap(287),laoo4_8 ),
     &(lmap(288),laoo4_9 ),(lmap(289),laoo4_10 ),(lmap(290),laoo5_1 ),
     &(lmap(291),laoo5_2 ),(lmap(292),laoo5_3 ),(lmap(293),laoo5_4 ),
     &(lmap(294),laoo5_5 ),(lmap(295),laoo5_6 ),(lmap(296),laoo5_7 ),
     &(lmap(297),laoo5_8 ),(lmap(298),laoo5_9 ),(lmap(299),laoo5_10 ),
     &(lmap(300),laoo6_1 ),(lmap(301),laoo6_2 ),(lmap(302),laoo6_3 ),
     &(lmap(303),laoo6_4 ),(lmap(304),laoo6_5 ),(lmap(305),laoo6_6 ),
     &(lmap(306),laoo6_7 ),(lmap(307),laoo6_8 ),(lmap(308),laoo6_9 ),
     &(lmap(309),laoo6_10 ),(lmap(310),laoo7_1 ),(lmap(311),laoo7_2 ),
     &(lmap(312),laoo7_3 ),(lmap(313),laoo7_4 ),(lmap(314),laoo7_5 ),
     &(lmap(315),laoo7_6 ),(lmap(316),laoo7_7 ),(lmap(317),laoo7_8 ),
     &(lmap(318),laoo7_9 ),(lmap(319),laoo7_10 ),(lmap(320),laoo8_1 ),
     &(lmap(321),laoo8_2 ),(lmap(322),laoo8_3 ),(lmap(323),laoo8_4 ),
     &(lmap(324),laoo8_5 ),(lmap(325),laoo8_6 ),(lmap(326),laoo8_7 ),
     &(lmap(327),laoo8_8 ),(lmap(328),laoo8_9 ),(lmap(329),laoo8_10 ),
     &(lmap(330),laoo9_1 ),(lmap(331),laoo9_2 ),(lmap(332),laoo9_3 ),
     &(lmap(333),laoo9_4 ),(lmap(334),laoo9_5 ),(lmap(335),laoo9_6 ),
     &(lmap(336),laoo9_7 ),(lmap(337),laoo9_8 ),(lmap(338),laoo9_9 ),
     &(lmap(339),laoo9_10 ),(lmap(340),laoo10_1 ),(lmap(341),laoo10_2 ),
     &(lmap(342),laoo10_3 ),(lmap(343),laoo10_4 ),(lmap(344),laoo10_5 ),
     &(lmap(345),laoo10_6 ),(lmap(346),laoo10_7 ),(lmap(347),laoo10_8 ),
     &(lmap(348),laoo10_9 ),(lmap(349),laoo10_10 ),(lmap(350),labs1_1 ),
     &(lmap(351),labs1_2 ),(lmap(352),labs1_3 ),(lmap(353),labs1_4 ),
     &(lmap(354),labs1_5 ),(lmap(355),labs1_6 ),(lmap(356),labs1_7 ),
     &(lmap(357),labs1_8 ),(lmap(358),labs1_9 ),(lmap(359),labs1_10 ),
     &(lmap(360),labs2_1 ),(lmap(361),labs2_2 ),(lmap(362),labs2_3 ),
     &(lmap(363),labs2_4 ),(lmap(364),labs2_5 ),(lmap(365),labs2_6 ),
     &(lmap(366),labs2_7 ),(lmap(367),labs2_8 ),(lmap(368),labs2_9 ),
     &(lmap(369),labs2_10 ),(lmap(370),labs3_1 ),(lmap(371),labs3_2 ),
     &(lmap(372),labs3_3 ),(lmap(373),labs3_4 ),(lmap(374),labs3_5 ),
     &(lmap(375),labs3_6 ),(lmap(376),labs3_7 ),(lmap(377),labs3_8 ),
     &(lmap(378),labs3_9 ),(lmap(379),labs3_10 ),(lmap(380),labs4_1 ),
     &(lmap(381),labs4_2 ),(lmap(382),labs4_3 ),(lmap(383),labs4_4 ),
     &(lmap(384),labs4_5 ),(lmap(385),labs4_6 ),(lmap(386),labs4_7 ),
     &(lmap(387),labs4_8 ),(lmap(388),labs4_9 ),(lmap(389),labs4_10 ),
     &(lmap(390),laas1_1 ),(lmap(391),laas1_2 ),(lmap(392),laas1_3 ),
     &(lmap(393),laas1_4 ),(lmap(394),laas1_5 ),(lmap(395),laas1_6 ),
     &(lmap(396),laas1_7 ),(lmap(397),laas1_8 ),(lmap(398),laas1_9 ),
     &(lmap(399),laas1_10 ),(lmap(400),laas2_1 ),(lmap(401),laas2_2 ),
     &(lmap(402),laas2_3 ),(lmap(403),laas2_4 ),(lmap(404),laas2_5 ),
     &(lmap(405),laas2_6 ),(lmap(406),laas2_7 ),(lmap(407),laas2_8 ),
     &(lmap(408),laas2_9 ),(lmap(409),laas2_10 ),(lmap(410),laas3_1 ),
     &(lmap(411),laas3_2 ),(lmap(412),laas3_3 ),(lmap(413),laas3_4 ),
     &(lmap(414),laas3_5 ),(lmap(415),laas3_6 ),(lmap(416),laas3_7 ),
     &(lmap(417),laas3_8 ),(lmap(418),laas3_9 ),(lmap(419),laas3_10 ),
     &(lmap(420),laas4_1 ),(lmap(421),laas4_2 ),(lmap(422),laas4_3 ),
     &(lmap(423),laas4_4 ),(lmap(424),laas4_5 ),(lmap(425),laas4_6 ),
     &(lmap(426),laas4_7 ),(lmap(427),laas4_8 ),(lmap(428),laas4_9 ),
     &(lmap(429),laas4_10 ),

c   scf

     &(lmap(430),lpoc_1  ),
     &(lmap(431),lpoc_2  ),(lmap(432),lpoc_3  ),(lmap(433),lpoc_4  ),
     &(lmap(434),lpoc_5  ),(lmap(435),lpoc_6  ),(lmap(436),lpoc_7  ),
     &(lmap(437),lpoc_8  ),(lmap(438),lpoc_9  ),(lmap(439),lpoc_10 ),
     &(lmap(440),lpec_1  ),(lmap(441),lpec_2  ),(lmap(442),lpec_3  ),
     &(lmap(443),lpec_4  ),(lmap(444),lpec_5  ),(lmap(445),lpec_6  ),
     &(lmap(446),lpec_7  ),(lmap(447),lpec_8  ),(lmap(448),lpec_9  ),
     &(lmap(449),lpec_10 ),(lmap(450),lcrust_1),(lmap(451),lcrust_2),
     &(lmap(452),lcrust_3),(lmap(453),lcrust_4),(lmap(454),lcrust_5),
     &(lmap(455),lcrust_6),(lmap(456),lcrust_7),(lmap(457),lcrust_8),
     &(lmap(458),lcrust_9),(lmap(459),lcrust_10),(lmap(460),lph2o_1),
     &(lmap(461),lph2o_2 ),(lmap(462),lph2o_3 ),(lmap(463),lph2o_4 ),
     &(lmap(464),lph2o_5 ),(lmap(465),lph2o_6 ),(lmap(466),lph2o_7 ),
     &(lmap(467),lph2o_8 ),(lmap(468),lph2o_9 ),(lmap(469),lph2o_10),
     &(lmap(470),lpcl_1  ),(lmap(471),lpcl_2  ),(lmap(472),lpcl_3  ),
     &(lmap(473),lpcl_4  ),(lmap(474),lpcl_5  ),(lmap(475),lpcl_6  ),
     &(lmap(476),lpcl_7  ),(lmap(477),lpcl_8  ),(lmap(478),lpcl_9  ),
     &(lmap(479),lpcl_10 ),(lmap(480),lna_1   ),(lmap(481),lna_2   ),
     &(lmap(482),lna_3   ),(lmap(483),lna_4   ),(lmap(484),lna_5   ),
     &(lmap(485),lna_6   ),(lmap(486),lna_7   ),(lmap(487),lna_8   ),
     &(lmap(488),lna_9   ),(lmap(489),lna_10  ),(lmap(490),lpnh4_1 ),
     &(lmap(491),lpnh4_2 ),(lmap(492),lpnh4_3 ),(lmap(493),lpnh4_4 ),
     &(lmap(494),lpnh4_5 ),(lmap(495),lpnh4_6 ),(lmap(496),lpnh4_7 ),
     &(lmap(497),lpnh4_8 ),(lmap(498),lpnh4_9 ),(lmap(499),lpnh4_10),
     &(lmap(500),lpno3_1 ),(lmap(501),lpno3_2 ),(lmap(502),lpno3_3 ),
     &(lmap(503),lpno3_4 ),(lmap(504),lpno3_5 ),(lmap(505),lpno3_6 ),
     &(lmap(506),lpno3_7 ),(lmap(507),lpno3_8 ),(lmap(508),lpno3_9 ),
     &(lmap(509),lpno3_10),(lmap(510),lpso4_1 ),(lmap(511),lpso4_2 ),
     &(lmap(512),lpso4_3 ),(lmap(513),lpso4_4 ),(lmap(514),lpso4_5 ),
     &(lmap(515),lpso4_6 ),(lmap(516),lpso4_7 ),(lmap(517),lpso4_8 ),
     &(lmap(518),lpso4_9 ),(lmap(519),lpso4_10),(lmap(520),lph2o   )
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
