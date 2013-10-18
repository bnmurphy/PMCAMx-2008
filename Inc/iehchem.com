c-----CAMx v4.02 030709
c
c     IEHCHEM.COM contains all chemistry variables
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
      common /iname/ imap(NSPNAM), irad(NRADNM)
c
      integer   ino   ,ino2  ,io3  
      integer   ipan  ,icres ,ipan2
      integer   impan ,ipbzn ,inphe
      integer   irno3 ,idcb2 ,idcb3
      integer   ihno4 ,iacet ,iald2
      integer   ialk1 ,ialk2 ,ialk3
      integer   ialk4 ,ialk5 ,iaro1
      integer   iaro2 ,ibacl ,ibald
      integer   ibcl1 ,ibcl2 ,ibuta
      integer   iccho ,iccrs  
      integer   icl2  ,ico   ,ico2h
      integer   ico3h ,icooh ,icprm
      integer   idcb1 ,ieth  ,iethe
      integer   ietoh ,ifcrs ,ifmcl
      integer   iform ,ifprm ,igly 
      integer   ih2o2 ,ihc2h ,ihcho
      integer   ihcl  ,ihono ,ihno3
      integer   iho2h ,ihocl ,iicl1
      integer   iicl2 ,iisop ,iispd
      integer   imek  ,imeoh ,imeth
      integer   imgly ,imvk  ,ina  
      integer   inh3  ,intr  ,inxoy
      integer   iole  ,iole1 ,iole2
      integer   ibpin ,ilimo ,imono ,isesq
      integer   iopen ,ipar  ,ipcl 
      integer   ipec  ,iphen ,ipna 
      integer   ipnh4 ,ipno3 ,ipoa 
      integer   iprod ,ipso4 ,irc2h
      integer   irc3h ,ircho ,irooh
      integer   iso2
      integer   isulf
      integer   iterp ,itol  ,ixn  ,ixyl 
      integer   icpo1, icpo2, icpo3, icpo4
      integer   icpo5, icpo6, icpo7, icpo8
      integer   icpo9, icpo10
      integer   icoo1, icoo2, icoo3, icoo4
      integer   icoo5, icoo6, icoo7, icoo8
      integer   icoo9, icoo10
      integer   icbs1, icbs2, icbs3, icbs4
      integer   icas1, icas2, icas3, icas4
c
      integer   iapo1, iapo2, iapo3, iapo4
      integer   iapo5, iapo6, iapo7, iapo8
      integer   iapo9, iapo10
      integer   iaoo1, iaoo2, iaoo3, iaoo4
      integer   iaoo5, iaoo6, iaoo7, iaoo8
      integer   iaoo9, iaoo10
      integer   iabs1, iabs2, iabs3, iabs4
      integer   iaas1, iaas2, iaas3, iaas4
c
      integer   iapo1_1  ,iapo1_2  ,iapo1_3  ,iapo1_4  ,iapo1_5
      integer   iapo1_6  ,iapo1_7  ,iapo1_8  ,iapo1_9  ,iapo1_10
      integer   iapo2_1  ,iapo2_2  ,iapo2_3  ,iapo2_4  ,iapo2_5
      integer   iapo2_6  ,iapo2_7  ,iapo2_8  ,iapo2_9  ,iapo2_10
      integer   iapo3_1  ,iapo3_2  ,iapo3_3  ,iapo3_4  ,iapo3_5
      integer   iapo3_6  ,iapo3_7  ,iapo3_8  ,iapo3_9  ,iapo3_10
      integer   iapo4_1  ,iapo4_2  ,iapo4_3  ,iapo4_4  ,iapo4_5
      integer   iapo4_6  ,iapo4_7  ,iapo4_8  ,iapo4_9  ,iapo4_10
      integer   iapo5_1  ,iapo5_2  ,iapo5_3  ,iapo5_4  ,iapo5_5
      integer   iapo5_6  ,iapo5_7  ,iapo5_8  ,iapo5_9  ,iapo5_10
      integer   iapo6_1  ,iapo6_2  ,iapo6_3  ,iapo6_4  ,iapo6_5
      integer   iapo6_6  ,iapo6_7  ,iapo6_8  ,iapo6_9  ,iapo6_10
      integer   iapo7_1  ,iapo7_2  ,iapo7_3  ,iapo7_4  ,iapo7_5
      integer   iapo7_6  ,iapo7_7  ,iapo7_8  ,iapo7_9  ,iapo7_10
      integer   iapo8_1  ,iapo8_2  ,iapo8_3  ,iapo8_4  ,iapo8_5
      integer   iapo8_6  ,iapo8_7  ,iapo8_8  ,iapo8_9  ,iapo8_10
      integer   iapo9_1  ,iapo9_2  ,iapo9_3  ,iapo9_4  ,iapo9_5
      integer   iapo9_6  ,iapo9_7  ,iapo9_8  ,iapo9_9  ,iapo9_10
      integer   iapo10_1  ,iapo10_2  ,iapo10_3  ,iapo10_4  ,iapo10_5
      integer   iapo10_6  ,iapo10_7  ,iapo10_8  ,iapo10_9  ,iapo10_10

      integer   iaoo1_1  ,iaoo1_2  ,iaoo1_3  ,iaoo1_4  ,iaoo1_5
      integer   iaoo1_6  ,iaoo1_7  ,iaoo1_8  ,iaoo1_9  ,iaoo1_10
      integer   iaoo2_1  ,iaoo2_2  ,iaoo2_3  ,iaoo2_4  ,iaoo2_5
      integer   iaoo2_6  ,iaoo2_7  ,iaoo2_8  ,iaoo2_9  ,iaoo2_10
      integer   iaoo3_1  ,iaoo3_2  ,iaoo3_3  ,iaoo3_4  ,iaoo3_5
      integer   iaoo3_6  ,iaoo3_7  ,iaoo3_8  ,iaoo3_9  ,iaoo3_10
      integer   iaoo4_1  ,iaoo4_2  ,iaoo4_3  ,iaoo4_4  ,iaoo4_5
      integer   iaoo4_6  ,iaoo4_7  ,iaoo4_8  ,iaoo4_9  ,iaoo4_10
      integer   iaoo5_1  ,iaoo5_2  ,iaoo5_3  ,iaoo5_4  ,iaoo5_5
      integer   iaoo5_6  ,iaoo5_7  ,iaoo5_8  ,iaoo5_9  ,iaoo5_10
      integer   iaoo6_1  ,iaoo6_2  ,iaoo6_3  ,iaoo6_4  ,iaoo6_5
      integer   iaoo6_6  ,iaoo6_7  ,iaoo6_8  ,iaoo6_9  ,iaoo6_10
      integer   iaoo7_1  ,iaoo7_2  ,iaoo7_3  ,iaoo7_4  ,iaoo7_5
      integer   iaoo7_6  ,iaoo7_7  ,iaoo7_8  ,iaoo7_9  ,iaoo7_10
      integer   iaoo8_1  ,iaoo8_2  ,iaoo8_3  ,iaoo8_4  ,iaoo8_5
      integer   iaoo8_6  ,iaoo8_7  ,iaoo8_8  ,iaoo8_9  ,iaoo8_10
      integer   iaoo9_1  ,iaoo9_2  ,iaoo9_3  ,iaoo9_4  ,iaoo9_5
      integer   iaoo9_6  ,iaoo9_7  ,iaoo9_8  ,iaoo9_9  ,iaoo9_10
      integer   iaoo10_1  ,iaoo10_2  ,iaoo10_3  ,iaoo10_4  ,iaoo10_5
      integer   iaoo10_6  ,iaoo10_7  ,iaoo10_8  ,iaoo10_9  ,iaoo10_10

      integer   iabs1_1  ,iabs1_2  ,iabs1_3  ,iabs1_4  ,iabs1_5
      integer   iabs1_6  ,iabs1_7  ,iabs1_8  ,iabs1_9  ,iabs1_10
      integer   iabs2_1  ,iabs2_2  ,iabs2_3  ,iabs2_4  ,iabs2_5
      integer   iabs2_6  ,iabs2_7  ,iabs2_8  ,iabs2_9  ,iabs2_10
      integer   iabs3_1  ,iabs3_2  ,iabs3_3  ,iabs3_4  ,iabs3_5
      integer   iabs3_6  ,iabs3_7  ,iabs3_8  ,iabs3_9  ,iabs3_10
      integer   iabs4_1  ,iabs4_2  ,iabs4_3  ,iabs4_4  ,iabs4_5
      integer   iabs4_6  ,iabs4_7  ,iabs4_8  ,iabs4_9  ,iabs4_10

      integer   iaas1_1  ,iaas1_2  ,iaas1_3  ,iaas1_4  ,iaas1_5
      integer   iaas1_6  ,iaas1_7  ,iaas1_8  ,iaas1_9  ,iaas1_10
      integer   iaas2_1  ,iaas2_2  ,iaas2_3  ,iaas2_4  ,iaas2_5
      integer   iaas2_6  ,iaas2_7  ,iaas2_8  ,iaas2_9  ,iaas2_10
      integer   iaas3_1  ,iaas3_2  ,iaas3_3  ,iaas3_4  ,iaas3_5
      integer   iaas3_6  ,iaas3_7  ,iaas3_8  ,iaas3_9  ,iaas3_10
      integer   iaas4_1  ,iaas4_2  ,iaas4_3  ,iaas4_4  ,iaas4_5
      integer   iaas4_6  ,iaas4_7  ,iaas4_8  ,iaas4_9  ,iaas4_10

c
      integer   ipoc_1
      integer   ipoc_2   ,ipoc_3   ,ipoc_4
      integer   ipoc_5   ,ipoc_6   ,ipoc_7
      integer   ipoc_8   ,ipoc_9   ,ipoc_10
      integer   ipec_1   ,ipec_2   ,ipec_3
      integer   ipec_4   ,ipec_5   ,ipec_6
      integer   ipec_7   ,ipec_8   ,ipec_9
      integer   ipec_10  ,icrust_1 ,icrust_2
      integer   icrust_3 ,icrust_4 ,icrust_5
      integer   icrust_6 ,icrust_7 ,icrust_8
      integer   icrust_9 ,icrust_10,iph2o_1
      integer   iph2o_2  ,iph2o_3  ,iph2o_4
      integer   iph2o_5  ,iph2o_6  ,iph2o_7
      integer   iph2o_8  ,iph2o_9  ,iph2o_10
      integer   ipcl_1   ,ipcl_2   ,ipcl_3
      integer   ipcl_4   ,ipcl_5   ,ipcl_6
      integer   ipcl_7   ,ipcl_8   ,ipcl_9
      integer   ipcl_10  ,ina_1    ,ina_2
      integer   ina_3    ,ina_4    ,ina_5
      integer   ina_6    ,ina_7    ,ina_8
      integer   ina_9    ,ina_10   ,ipnh4_1
      integer   ipnh4_2  ,ipnh4_3  ,ipnh4_4
      integer   ipnh4_5  ,ipnh4_6  ,ipnh4_7
      integer   ipnh4_8  ,ipnh4_9  ,ipnh4_10
      integer   ipno3_1  ,ipno3_2  ,ipno3_3
      integer   ipno3_4  ,ipno3_5  ,ipno3_6
      integer   ipno3_7  ,ipno3_8  ,ipno3_9
      integer   ipno3_10 ,ipso4_1  ,ipso4_2
      integer   ipso4_3  ,ipso4_4  ,ipso4_5
      integer   ipso4_6  ,ipso4_7  ,ipso4_8
      integer   ipso4_9  ,ipso4_10 ,iph2o
c
      equivalence (imap(1), ino  ), (imap(2), ino2 ), (imap(3), io3  ),
     &            (imap(4), ipan ), (imap(5), icres), (imap(6), ipan2),
     &            (imap(7), impan), (imap(8), ipbzn), (imap(9), inphe),
     &            (imap(10),irno3), (imap(11),idcb2), (imap(12),idcb3),
     &            (imap(13),ihno4), (imap(14),iacet), (imap(15),iald2),
     &            (imap(16),ialk1), (imap(17),ialk2), (imap(18),ialk3),
     &            (imap(19),ialk4), (imap(20),ialk5), (imap(21),iaro1),
     &            (imap(22),iaro2), (imap(23),ibacl), (imap(24),ibald),
     &            (imap(25),ibcl1), (imap(26),ibcl2), (imap(27),ibuta),
     &            (imap(28),iccho), (imap(29),iccrs), (imap(30),icpo1),
     &            (imap(31),icpo2), (imap(32),icpo3), (imap(33),icpo4),
     &            (imap(34),icpo5), (imap(35),icpo6), (imap(36),icpo7),
     &            (imap(37),icpo8), (imap(38),icpo9), (imap(39),icpo10),
     &            (imap(40),icoo1), (imap(41),icoo2), (imap(42),icoo3),
     &            (imap(43),icoo4), (imap(44),icoo5), (imap(45),icoo6),
     &            (imap(46),icoo7), (imap(47),icoo8), (imap(48),icoo9),
     &            (imap(49),icoo10), (imap(50),icbs1), (imap(51),icbs2),
     &            (imap(52),icbs3), (imap(53),icbs4), (imap(54),icas1),
     &            (imap(55),icas2), (imap(56),icas3), (imap(57),icas4),
c
     &           (imap(58),icl2 ), (imap(59),ico  ), (imap(60),ico2h),
     &           (imap(61),ico3h), (imap(62),icooh), (imap(63),icprm),
     &           (imap(64),idcb1), (imap(65),ieth ), (imap(66),iethe),
     &           (imap(67),ietoh), (imap(68),ifcrs), (imap(69),ifmcl),
     &           (imap(70),iform), (imap(71),ifprm), (imap(72),igly ),
     &           (imap(73),ih2o2), (imap(74),ihc2h), (imap(75),ihcho),
     &           (imap(76),ihcl ), (imap(77),ihono), (imap(78),ihno3),
     &           (imap(79),iho2h), (imap(80),ihocl), (imap(81),iicl1),
     &           (imap(82),iicl2), (imap(83),iisop), (imap(84),iispd),
     &           (imap(85),imek ),(imap(86),imeoh),(imap(87),imeth),
     &           (imap(88),imgly),(imap(89),imvk ),(imap(90),ina  ),
     &           (imap(91),inh3 ),(imap(92),intr ),(imap(93),inxoy),
     &           (imap(94),iole ),(imap(95),iole1),(imap(96),iole2),
     &           (imap(97),ibpin),(imap(98),ilimo),(imap(99),imono),
     &           (imap(100),isesq),(imap(101),iopen),(imap(102),ipar ),
     &           (imap(103),ipcl ),(imap(104),ipec ),(imap(105),iphen),
     &           (imap(106),ipna ),(imap(107),ipnh4),(imap(108),ipno3),
     &           (imap(109),ipoa ),(imap(110),iprod),(imap(111),ipso4),
     &           (imap(112),irc2h),(imap(113),irc3h),(imap(114),ircho),
     &           (imap(115),irooh),(imap(116),iso2 ),(imap(117),iapo1),
     &           (imap(118),iapo2),(imap(119),iapo3),(imap(120),iapo4),
     &           (imap(121),iapo5),(imap(122),iapo6),(imap(123),iapo7),
     &           (imap(124),iapo8),(imap(125),iapo9),(imap(126),iapo10),
     &           (imap(127),iaoo1),(imap(128),iaoo2),(imap(129),iaoo3),
     &           (imap(130),iaoo4),(imap(131),iaoo5),(imap(132),iaoo6),
     &           (imap(133),iaoo7),(imap(134),iaoo8),(imap(135),iaoo9),
     &           (imap(136),iaoo10),(imap(137),iabs1),(imap(138),iabs2),
     &           (imap(139),iabs3),(imap(140),iabs4),(imap(141),iaas1),
     &           (imap(142),iaas2),(imap(143),iaas3),(imap(144),iaas4),
     &           (imap(145),isulf),(imap(146),iterp),
     &           (imap(147),itol ),(imap(148),ixn  ),(imap(149),ixyl ), 
c
     &(imap(150),iapo1_1 ),(imap(151),iapo1_2 ),(imap(152),iapo1_3 ),
     &(imap(153),iapo1_4 ),(imap(154),iapo1_5 ),(imap(155),iapo1_6 ),
     &(imap(156),iapo1_7 ),(imap(157),iapo1_8 ),(imap(158),iapo1_9 ),
     &(imap(159),iapo1_10 ),(imap(160),iapo2_1 ),(imap(161),iapo2_2 ),
     &(imap(162),iapo2_3 ),(imap(163),iapo2_4 ),(imap(164),iapo2_5 ),
     &(imap(165),iapo2_6 ),(imap(166),iapo2_7 ),(imap(167),iapo2_8 ),
     &(imap(168),iapo2_9 ),(imap(169),iapo2_10 ),(imap(170),iapo3_1 ),
     &(imap(171),iapo3_2 ),(imap(172),iapo3_3 ),(imap(173),iapo3_4 ),
     &(imap(174),iapo3_5 ),(imap(175),iapo3_6 ),(imap(176),iapo3_7 ),
     &(imap(177),iapo3_8 ),(imap(178),iapo3_9 ),(imap(179),iapo3_10 ),
     &(imap(180),iapo4_1 ),(imap(181),iapo4_2 ),(imap(182),iapo4_3 ),
     &(imap(183),iapo4_4 ),(imap(184),iapo4_5 ),(imap(185),iapo4_6 ),
     &(imap(186),iapo4_7 ),(imap(187),iapo4_8 ),(imap(188),iapo4_9 ),
     &(imap(189),iapo4_10 ),(imap(190),iapo5_1 ),(imap(191),iapo5_2 ),
     &(imap(192),iapo5_3 ),(imap(193),iapo5_4 ),(imap(194),iapo5_5 ),
     &(imap(195),iapo5_6 ),(imap(196),iapo5_7 ),(imap(197),iapo5_8 ),
     &(imap(198),iapo5_9 ),(imap(199),iapo5_10 ),(imap(200),iapo6_1 ),
     &(imap(201),iapo6_2 ),(imap(202),iapo6_3 ),(imap(203),iapo6_4 ),
     &(imap(204),iapo6_5 ),(imap(205),iapo6_6 ),(imap(206),iapo6_7 ),
     &(imap(207),iapo6_8 ),(imap(208),iapo6_9 ),(imap(209),iapo6_10 ),
     &(imap(210),iapo7_1 ),(imap(211),iapo7_2 ),(imap(212),iapo7_3 ),
     &(imap(213),iapo7_4 ),(imap(214),iapo7_5 ),(imap(215),iapo7_6 ),
     &(imap(216),iapo7_7 ),(imap(217),iapo7_8 ),(imap(218),iapo7_9 ),
     &(imap(219),iapo7_10 ),(imap(220),iapo8_1 ),(imap(221),iapo8_2 ),
     &(imap(222),iapo8_3 ),(imap(223),iapo8_4 ),(imap(224),iapo8_5 ),
     &(imap(225),iapo8_6 ),(imap(226),iapo8_7 ),(imap(227),iapo8_8 ),
     &(imap(228),iapo8_9 ),(imap(229),iapo8_10 ),(imap(230),iapo9_1 ),
     &(imap(231),iapo9_2 ),(imap(232),iapo9_3 ),(imap(233),iapo9_4 ),
     &(imap(234),iapo9_5 ),(imap(235),iapo9_6 ),(imap(236),iapo9_7 ),
     &(imap(237),iapo9_8 ),(imap(238),iapo9_9 ),(imap(239),iapo9_10 ),
     &(imap(240),iapo10_1 ),(imap(241),iapo10_2 ),(imap(242),iapo10_3 ),
     &(imap(243),iapo10_4 ),(imap(244),iapo10_5 ),(imap(245),iapo10_6 ),
     &(imap(246),iapo10_7 ),(imap(247),iapo10_8 ),(imap(248),iapo10_9 ),
     &(imap(249),iapo10_10 ),(imap(250),iaoo1_1 ),(imap(251),iaoo1_2 ),
     &(imap(252),iaoo1_3 ),(imap(253),iaoo1_4 ),(imap(254),iaoo1_5 ),
     &(imap(255),iaoo1_6 ),(imap(256),iaoo1_7 ),(imap(257),iaoo1_8 ),
     &(imap(258),iaoo1_9 ),(imap(259),iaoo1_10 ),(imap(260),iaoo2_1 ),
     &(imap(261),iaoo2_2 ),(imap(262),iaoo2_3 ),(imap(263),iaoo2_4 ),
     &(imap(264),iaoo2_5 ),(imap(265),iaoo2_6 ),(imap(266),iaoo2_7 ),
     &(imap(267),iaoo2_8 ),(imap(268),iaoo2_9 ),(imap(269),iaoo2_10 ),
     &(imap(270),iaoo3_1 ),(imap(271),iaoo3_2 ),(imap(272),iaoo3_3 ),
     &(imap(273),iaoo3_4 ),(imap(274),iaoo3_5 ),(imap(275),iaoo3_6 ),
     &(imap(276),iaoo3_7 ),(imap(277),iaoo3_8 ),(imap(278),iaoo3_9 ),
     &(imap(279),iaoo3_10 ),(imap(280),iaoo4_1 ),(imap(281),iaoo4_2 ),
     &(imap(282),iaoo4_3 ),(imap(283),iaoo4_4 ),(imap(284),iaoo4_5 ),
     &(imap(285),iaoo4_6 ),(imap(286),iaoo4_7 ),(imap(287),iaoo4_8 ),
     &(imap(288),iaoo4_9 ),(imap(289),iaoo4_10 ),(imap(290),iaoo5_1 ),
     &(imap(291),iaoo5_2 ),(imap(292),iaoo5_3 ),(imap(293),iaoo5_4 ),
     &(imap(294),iaoo5_5 ),(imap(295),iaoo5_6 ),(imap(296),iaoo5_7 ),
     &(imap(297),iaoo5_8 ),(imap(298),iaoo5_9 ),(imap(299),iaoo5_10 ),
     &(imap(300),iaoo6_1 ),(imap(301),iaoo6_2 ),(imap(302),iaoo6_3 ),
     &(imap(303),iaoo6_4 ),(imap(304),iaoo6_5 ),(imap(305),iaoo6_6 ),
     &(imap(306),iaoo6_7 ),(imap(307),iaoo6_8 ),(imap(308),iaoo6_9 ),
     &(imap(309),iaoo6_10 ),(imap(310),iaoo7_1 ),(imap(311),iaoo7_2 ),
     &(imap(312),iaoo7_3 ),(imap(313),iaoo7_4 ),(imap(314),iaoo7_5 ),
     &(imap(315),iaoo7_6 ),(imap(316),iaoo7_7 ),(imap(317),iaoo7_8 ),
     &(imap(318),iaoo7_9 ),(imap(319),iaoo7_10 ),(imap(320),iaoo8_1 ),
     &(imap(321),iaoo8_2 ),(imap(322),iaoo8_3 ),(imap(323),iaoo8_4 ),
     &(imap(324),iaoo8_5 ),(imap(325),iaoo8_6 ),(imap(326),iaoo8_7 ),
     &(imap(327),iaoo8_8 ),(imap(328),iaoo8_9 ),(imap(329),iaoo8_10 ),
     &(imap(330),iaoo9_1 ),(imap(331),iaoo9_2 ),(imap(332),iaoo9_3 ),
     &(imap(333),iaoo9_4 ),(imap(334),iaoo9_5 ),(imap(335),iaoo9_6 ),
     &(imap(336),iaoo9_7 ),(imap(337),iaoo9_8 ),(imap(338),iaoo9_9 ),
     &(imap(339),iaoo9_10 ),(imap(340),iaoo10_1 ),(imap(341),iaoo10_2 ),
     &(imap(342),iaoo10_3 ),(imap(343),iaoo10_4 ),(imap(344),iaoo10_5 ),
     &(imap(345),iaoo10_6 ),(imap(346),iaoo10_7 ),(imap(347),iaoo10_8 ),
     &(imap(348),iaoo10_9 ),(imap(349),iaoo10_10 ),(imap(350),iabs1_1 ),
     &(imap(351),iabs1_2 ),(imap(352),iabs1_3 ),(imap(353),iabs1_4 ),
     &(imap(354),iabs1_5 ),(imap(355),iabs1_6 ),(imap(356),iabs1_7 ),
     &(imap(357),iabs1_8 ),(imap(358),iabs1_9 ),(imap(359),iabs1_10 ),
     &(imap(360),iabs2_1 ),(imap(361),iabs2_2 ),(imap(362),iabs2_3 ),
     &(imap(363),iabs2_4 ),(imap(364),iabs2_5 ),(imap(365),iabs2_6 ),
     &(imap(366),iabs2_7 ),(imap(367),iabs2_8 ),(imap(368),iabs2_9 ),
     &(imap(369),iabs2_10 ),(imap(370),iabs3_1 ),(imap(371),iabs3_2 ),
     &(imap(372),iabs3_3 ),(imap(373),iabs3_4 ),(imap(374),iabs3_5 ),
     &(imap(375),iabs3_6 ),(imap(376),iabs3_7 ),(imap(377),iabs3_8 ),
     &(imap(378),iabs3_9 ),(imap(379),iabs3_10 ),(imap(380),iabs4_1 ),
     &(imap(381),iabs4_2 ),(imap(382),iabs4_3 ),(imap(383),iabs4_4 ),
     &(imap(384),iabs4_5 ),(imap(385),iabs4_6 ),(imap(386),iabs4_7 ),
     &(imap(387),iabs4_8 ),(imap(388),iabs4_9 ),(imap(389),iabs4_10 ),
     &(imap(390),iaas1_1 ),(imap(391),iaas1_2 ),(imap(392),iaas1_3 ),
     &(imap(393),iaas1_4 ),(imap(394),iaas1_5 ),(imap(395),iaas1_6 ),
     &(imap(396),iaas1_7 ),(imap(397),iaas1_8 ),(imap(398),iaas1_9 ),
     &(imap(399),iaas1_10 ),(imap(400),iaas2_1 ),(imap(401),iaas2_2 ),
     &(imap(402),iaas2_3 ),(imap(403),iaas2_4 ),(imap(404),iaas2_5 ),
     &(imap(405),iaas2_6 ),(imap(406),iaas2_7 ),(imap(407),iaas2_8 ),
     &(imap(408),iaas2_9 ),(imap(409),iaas2_10 ),(imap(410),iaas3_1 ),
     &(imap(411),iaas3_2 ),(imap(412),iaas3_3 ),(imap(413),iaas3_4 ),
     &(imap(414),iaas3_5 ),(imap(415),iaas3_6 ),(imap(416),iaas3_7 ),
     &(imap(417),iaas3_8 ),(imap(418),iaas3_9 ),(imap(419),iaas3_10 ),
     &(imap(420),iaas4_1 ),(imap(421),iaas4_2 ),(imap(422),iaas4_3 ),
     &(imap(423),iaas4_4 ),(imap(424),iaas4_5 ),(imap(425),iaas4_6 ),
     &(imap(426),iaas4_7 ),(imap(427),iaas4_8 ),(imap(428),iaas4_9 ),
     &(imap(429),iaas4_10 ),

c   scf

     &(imap(430),ipoc_1  ),
     &(imap(431),ipoc_2  ),(imap(432),ipoc_3  ),(imap(433),ipoc_4  ),
     &(imap(434),ipoc_5  ),(imap(435),ipoc_6  ),(imap(436),ipoc_7  ),
     &(imap(437),ipoc_8  ),(imap(438),ipoc_9  ),(imap(439),ipoc_10 ),
     &(imap(440),ipec_1  ),(imap(441),ipec_2  ),(imap(442),ipec_3  ),
     &(imap(443),ipec_4  ),(imap(444),ipec_5  ),(imap(445),ipec_6  ),
     &(imap(446),ipec_7  ),(imap(447),ipec_8  ),(imap(448),ipec_9  ),
     &(imap(449),ipec_10 ),(imap(450),icrust_1),(imap(451),icrust_2),
     &(imap(452),icrust_3),(imap(453),icrust_4),(imap(454),icrust_5),
     &(imap(455),icrust_6),(imap(456),icrust_7),(imap(457),icrust_8),
     &(imap(458),icrust_9),(imap(459),icrust_10),(imap(460),iph2o_1),
     &(imap(461),iph2o_2 ),(imap(462),iph2o_3 ),(imap(463),iph2o_4 ),
     &(imap(464),iph2o_5 ),(imap(465),iph2o_6 ),(imap(466),iph2o_7 ),
     &(imap(467),iph2o_8 ),(imap(468),iph2o_9 ),(imap(469),iph2o_10),
     &(imap(470),ipcl_1  ),(imap(471),ipcl_2  ),(imap(472),ipcl_3  ),
     &(imap(473),ipcl_4  ),(imap(474),ipcl_5  ),(imap(475),ipcl_6  ),
     &(imap(476),ipcl_7  ),(imap(477),ipcl_8  ),(imap(478),ipcl_9  ),
     &(imap(479),ipcl_10 ),(imap(480),ina_1   ),(imap(481),ina_2   ),
     &(imap(482),ina_3   ),(imap(483),ina_4   ),(imap(484),ina_5   ),
     &(imap(485),ina_6   ),(imap(486),ina_7   ),(imap(487),ina_8   ),
     &(imap(488),ina_9   ),(imap(489),ina_10  ),(imap(490),ipnh4_1 ),
     &(imap(491),ipnh4_2 ),(imap(492),ipnh4_3 ),(imap(493),ipnh4_4 ),
     &(imap(494),ipnh4_5 ),(imap(495),ipnh4_6 ),(imap(496),ipnh4_7 ),
     &(imap(497),ipnh4_8 ),(imap(498),ipnh4_9 ),(imap(499),ipnh4_10),
     &(imap(500),ipno3_1 ),(imap(501),ipno3_2 ),(imap(502),ipno3_3 ),
     &(imap(503),ipno3_4 ),(imap(504),ipno3_5 ),(imap(505),ipno3_6 ),
     &(imap(506),ipno3_7 ),(imap(507),ipno3_8 ),(imap(508),ipno3_9 ),
     &(imap(509),ipno3_10),(imap(510),ipso4_1 ),(imap(511),ipso4_2 ),
     &(imap(512),ipso4_3 ),(imap(513),ipso4_4 ),(imap(514),ipso4_5 ),
     &(imap(515),ipso4_6 ),(imap(516),ipso4_7 ),(imap(517),ipso4_8 ),
     &(imap(518),ipso4_9 ),(imap(519),ipso4_10),(imap(520),iph2o   )
c
      integer   io1d  ,io    ,iclo 
      integer   icl   ,in2o5 ,ino3 
      integer   ioh   ,iho2  ,ic2o3
      integer   ixo2  ,ixo2n ,ito2 
      integer   iror  ,icro  ,iro2r
      integer   ir2o2 ,iro2n ,icco3
      integer   irco3 ,imco3 ,ibzco
      integer   icxo2 ,ihco3 ,itbuo
      integer   ibzo  ,ibzno
c
      equivalence (irad(1), io1d ), (irad(2), io   ), (irad(3), iclo ),
     &            (irad(4), icl  ), (irad(5), in2o5), (irad(6), ino3 ),
     &            (irad(7), ioh  ), (irad(8), iho2 ), (irad(9), ic2o3),
     &            (irad(10),ixo2 ), (irad(11),ixo2n), (irad(12),ito2 ),
     &            (irad(13),iror ), (irad(14),icro ), (irad(15),iro2r),
     &            (irad(16),ir2o2), (irad(17),iro2n), (irad(18),icco3),
     &            (irad(19),irco3), (irad(20),imco3), (irad(21),ibzco),
     &            (irad(22),icxo2), (irad(23),ihco3), (irad(24),itbuo),
     &            (irad(25),ibzo ), (irad(26),ibzno)
