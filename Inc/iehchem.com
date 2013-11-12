c-----CAMx v4.02 030709
c
c     IEHCHEM.COM contains all chemistry variables
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
      common /iname/ imap(NSPNAM), irad(NRADNM)
c
c   BNM  CHANGED MANY OF THESE ARRAYS TO ADD NTSOA
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
      integer   icoo1, icoo2, icoo3, icoo4
      integer   icoo5, icoo6, icoo7, icoo8
      integer   icbs1, icbs2, icbs3, icbs4, icbs5
      integer   icas1, icas2, icas3, icas4, icas5
      integer   icns1, icns2, icns3, icns4
      integer   icns5, icns6, icns7, icns8
c
      integer   iapo1, iapo2, iapo3, iapo4
      integer   iapo5, iapo6, iapo7, iapo8
      integer   iaoo1, iaoo2, iaoo3, iaoo4
      integer   iaoo5, iaoo6, iaoo7, iaoo8
      integer   iabs1, iabs2, iabs3, iabs4, iabs5
      integer   iaas1, iaas2, iaas3, iaas4, iaas5
      integer   ians1, ians2, ians3, ians4
      integer   ians5, ians6, ians7, ians8
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

      integer   iabs1_1  ,iabs1_2  ,iabs1_3  ,iabs1_4  ,iabs1_5
      integer   iabs1_6  ,iabs1_7  ,iabs1_8  ,iabs1_9  ,iabs1_10
      integer   iabs2_1  ,iabs2_2  ,iabs2_3  ,iabs2_4  ,iabs2_5
      integer   iabs2_6  ,iabs2_7  ,iabs2_8  ,iabs2_9  ,iabs2_10
      integer   iabs3_1  ,iabs3_2  ,iabs3_3  ,iabs3_4  ,iabs3_5
      integer   iabs3_6  ,iabs3_7  ,iabs3_8  ,iabs3_9  ,iabs3_10
      integer   iabs4_1  ,iabs4_2  ,iabs4_3  ,iabs4_4  ,iabs4_5
      integer   iabs4_6  ,iabs4_7  ,iabs4_8  ,iabs4_9  ,iabs4_10
      integer   iabs5_1  ,iabs5_2  ,iabs5_3  ,iabs5_4  ,iabs5_5
      integer   iabs5_6  ,iabs5_7  ,iabs5_8  ,iabs5_9  ,iabs5_10

      integer   iaas1_1  ,iaas1_2  ,iaas1_3  ,iaas1_4  ,iaas1_5
      integer   iaas1_6  ,iaas1_7  ,iaas1_8  ,iaas1_9  ,iaas1_10
      integer   iaas2_1  ,iaas2_2  ,iaas2_3  ,iaas2_4  ,iaas2_5
      integer   iaas2_6  ,iaas2_7  ,iaas2_8  ,iaas2_9  ,iaas2_10
      integer   iaas3_1  ,iaas3_2  ,iaas3_3  ,iaas3_4  ,iaas3_5
      integer   iaas3_6  ,iaas3_7  ,iaas3_8  ,iaas3_9  ,iaas3_10
      integer   iaas4_1  ,iaas4_2  ,iaas4_3  ,iaas4_4  ,iaas4_5
      integer   iaas4_6  ,iaas4_7  ,iaas4_8  ,iaas4_9  ,iaas4_10
      integer   iaas5_1  ,iaas5_2  ,iaas5_3  ,iaas5_4  ,iaas5_5
      integer   iaas5_6  ,iaas5_7  ,iaas5_8  ,iaas5_9  ,iaas5_10

      integer   ians1_1  ,ians1_2  ,ians1_3  ,ians1_4  ,ians1_5
      integer   ians1_6  ,ians1_7  ,ians1_8  ,ians1_9  ,ians1_10
      integer   ians2_1  ,ians2_2  ,ians2_3  ,ians2_4  ,ians2_5
      integer   ians2_6  ,ians2_7  ,ians2_8  ,ians2_9  ,ians2_10
      integer   ians3_1  ,ians3_2  ,ians3_3  ,ians3_4  ,ians3_5
      integer   ians3_6  ,ians3_7  ,ians3_8  ,ians3_9  ,ians3_10
      integer   ians4_1  ,ians4_2  ,ians4_3  ,ians4_4  ,ians4_5
      integer   ians4_6  ,ians4_7  ,ians4_8  ,ians4_9  ,ians4_10
      integer   ians5_1  ,ians5_2  ,ians5_3  ,ians5_4  ,ians5_5
      integer   ians5_6  ,ians5_7  ,ians5_8  ,ians5_9  ,ians5_10
      integer   ians6_1  ,ians6_2  ,ians6_3  ,ians6_4  ,ians6_5
      integer   ians6_6  ,ians6_7  ,ians6_8  ,ians6_9  ,ians6_10
      integer   ians7_1  ,ians7_2  ,ians7_3  ,ians7_4  ,ians7_5
      integer   ians7_6  ,ians7_7  ,ians7_8  ,ians7_9  ,ians7_10
      integer   ians8_1  ,ians8_2  ,ians8_3  ,ians8_4  ,ians8_5
      integer   ians8_6  ,ians8_7  ,ians8_8  ,ians8_9  ,ians8_10
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
     &            (imap(28),iccho), (imap(29),iccrs), 
     &		  (imap(30),icpo1), (imap(31),icpo2), (imap(32),icpo3), 
     &		  (imap(33),icpo4), (imap(34),icpo5), (imap(35),icpo6), 
     &		  (imap(36),icpo7), (imap(37),icpo8),
     &            (imap(38),icoo1), (imap(39),icoo2), (imap(40),icoo3),
     &            (imap(41),icoo4), (imap(42),icoo5), (imap(43),icoo6),
     &            (imap(44),icoo7), (imap(45),icoo8),
     &            (imap(46),icbs1), (imap(47),icbs2), (imap(48),icbs3), 
     &		  (imap(49),icbs4), (imap(50),icbs5),
     &		  (imap(51),icas1), (imap(52),icas2), (imap(53),icas3), 
     &		  (imap(54),icas4), (imap(55),icas5),
     &            (imap(56),icns1), (imap(57),icns2), (imap(58),icns3),
     &            (imap(59),icns4), (imap(60),icns5), (imap(61),icns6),
     &            (imap(62),icns7), (imap(63),icns8),
c
     &           (imap(64),icl2 ), (imap(65),ico  ), (imap(66),ico2h),
     &           (imap(67),ico3h), (imap(68),icooh), (imap(69),icprm),
     &           (imap(70),idcb1), (imap(71),ieth ), (imap(72),iethe),
     &           (imap(73),ietoh), (imap(74),ifcrs), (imap(75),ifmcl),
     &           (imap(76),iform), (imap(77),ifprm), (imap(78),igly ),
     &           (imap(79),ih2o2), (imap(80),ihc2h), (imap(81),ihcho),
     &           (imap(82),ihcl ), (imap(83),ihono), (imap(84),ihno3),
     &           (imap(85),iho2h), (imap(86),ihocl), (imap(87),iicl1),
     &           (imap(88),iicl2), (imap(89),iisop), (imap(90),iispd),
     &           (imap(91),imek ), (imap(92),imeoh), (imap(93),imeth),
     &           (imap(94),imgly), (imap(95),imvk ), (imap(96),ina  ),
     &           (imap(97),inh3 ), (imap(98),intr ), (imap(99),inxoy),
     &           (imap(100),iole ),(imap(101),iole1),(imap(102),iole2),
     &           (imap(103),ibpin),(imap(104),ilimo),(imap(105),imono),
     &           (imap(106),isesq),(imap(107),iopen),(imap(108),ipar ),
     &           (imap(109),ipcl ),(imap(110),ipec ),(imap(111),iphen),
     &           (imap(112),ipna ),(imap(113),ipnh4),(imap(114),ipno3),
     &           (imap(115),ipoa ),(imap(116),iprod),(imap(117),ipso4),
     &           (imap(118),irc2h),(imap(119),irc3h),(imap(120),ircho),
     &           (imap(121),irooh),(imap(122),iso2 ),
     &		 (imap(123),iapo1),(imap(124),iapo2),(imap(125),iapo3),
     &           (imap(126),iapo4),(imap(127),iapo5),(imap(128),iapo6),
     &		 (imap(129),iapo7),(imap(130),iapo8),
     &           (imap(131),iaoo1),(imap(132),iaoo2),(imap(133),iaoo3),
     &           (imap(134),iaoo4),(imap(135),iaoo5),(imap(136),iaoo6),
     &           (imap(137),iaoo7),(imap(138),iaoo8),
     &           (imap(139),iabs1),(imap(140),iabs2),(imap(141),iabs3), 
     &		 (imap(142),iabs4),(imap(143),iabs5),
     &		 (imap(144),iaas1),(imap(145),iaas2),(imap(146),iaas3),
     &		 (imap(147),iaas4),(imap(148),iaas5),
     &		 (imap(149),ians1),(imap(150),ians2),(imap(151),ians3),
     &           (imap(152),ians4),(imap(153),ians5),(imap(154),ians6),
     &		 (imap(155),ians7),(imap(156),ians8),
     &           (imap(157),isulf), (imap(158),iterp),
     &           (imap(159),itol ), (imap(160),ixn  ),(imap(161),ixyl ), 
c
     &(imap(162),iapo1_1 ), (imap(163),iapo1_2 ), (imap(164),iapo1_3 ),
     &(imap(165),iapo1_4 ), (imap(166),iapo1_5 ), (imap(167),iapo1_6 ),
     &(imap(168),iapo1_7 ), (imap(169),iapo1_8 ), (imap(170),iapo1_9 ),
     &(imap(171),iapo1_10 ),(imap(172),iapo2_1 ), (imap(173),iapo2_2 ),
     &(imap(174),iapo2_3 ), (imap(175),iapo2_4 ), (imap(176),iapo2_5 ),
     &(imap(177),iapo2_6 ), (imap(178),iapo2_7 ), (imap(179),iapo2_8 ),
     &(imap(180),iapo2_9 ), (imap(181),iapo2_10 ),(imap(182),iapo3_1 ),
     &(imap(183),iapo3_2 ), (imap(184),iapo3_3 ), (imap(185),iapo3_4 ),
     &(imap(186),iapo3_5 ), (imap(187),iapo3_6 ), (imap(188),iapo3_7 ),
     &(imap(189),iapo3_8 ), (imap(190),iapo3_9 ), (imap(191),iapo3_10 ),
     &(imap(192),iapo4_1 ), (imap(193),iapo4_2 ), (imap(194),iapo4_3 ),
     &(imap(195),iapo4_4 ), (imap(196),iapo4_5 ), (imap(197),iapo4_6 ),
     &(imap(198),iapo4_7 ), (imap(199),iapo4_8 ), (imap(200),iapo4_9 ),
     &(imap(201),iapo4_10 ),(imap(202),iapo5_1 ), (imap(203),iapo5_2 ),
     &(imap(204),iapo5_3 ), (imap(205),iapo5_4 ), (imap(206),iapo5_5 ),
     &(imap(207),iapo5_6 ), (imap(208),iapo5_7 ), (imap(209),iapo5_8 ),
     &(imap(210),iapo5_9 ), (imap(211),iapo5_10 ),(imap(212),iapo6_1 ),
     &(imap(213),iapo6_2 ), (imap(214),iapo6_3 ), (imap(215),iapo6_4 ),
     &(imap(216),iapo6_5 ), (imap(217),iapo6_6 ), (imap(218),iapo6_7 ),
     &(imap(219),iapo6_8 ), (imap(220),iapo6_9 ), (imap(221),iapo6_10 ),
     &(imap(222),iapo7_1 ), (imap(223),iapo7_2 ), (imap(224),iapo7_3 ),
     &(imap(225),iapo7_4 ), (imap(226),iapo7_5 ), (imap(227),iapo7_6 ),
     &(imap(228),iapo7_7 ), (imap(229),iapo7_8 ), (imap(230),iapo7_9 ),
     &(imap(231),iapo7_10 ),(imap(232),iapo8_1 ), (imap(233),iapo8_2 ),
     &(imap(234),iapo8_3 ), (imap(235),iapo8_4 ), (imap(236),iapo8_5 ),
     &(imap(237),iapo8_6 ), (imap(238),iapo8_7 ), (imap(239),iapo8_8 ),
     &(imap(240),iapo8_9 ), (imap(241),iapo8_10 ),
     &(imap(242),iaoo1_1 ), (imap(243),iaoo1_2 ),
     &(imap(244),iaoo1_3 ), (imap(245),iaoo1_4 ), (imap(246),iaoo1_5 ),
     &(imap(247),iaoo1_6 ), (imap(248),iaoo1_7 ), (imap(249),iaoo1_8 ),
     &(imap(250),iaoo1_9 ), (imap(251),iaoo1_10 ),(imap(252),iaoo2_1 ),
     &(imap(253),iaoo2_2 ), (imap(254),iaoo2_3 ), (imap(255),iaoo2_4 ),
     &(imap(256),iaoo2_5 ), (imap(257),iaoo2_6 ), (imap(258),iaoo2_7 ),
     &(imap(259),iaoo2_8 ), (imap(260),iaoo2_9 ), (imap(261),iaoo2_10 ),
     &(imap(262),iaoo3_1 ), (imap(263),iaoo3_2 ), (imap(264),iaoo3_3 ),
     &(imap(265),iaoo3_4 ), (imap(266),iaoo3_5 ), (imap(267),iaoo3_6 ),
     &(imap(268),iaoo3_7 ), (imap(269),iaoo3_8 ), (imap(270),iaoo3_9 ),
     &(imap(271),iaoo3_10 ),(imap(272),iaoo4_1 ), (imap(273),iaoo4_2 ),
     &(imap(274),iaoo4_3 ), (imap(275),iaoo4_4 ), (imap(276),iaoo4_5 ),
     &(imap(277),iaoo4_6 ), (imap(278),iaoo4_7 ), (imap(279),iaoo4_8 ),
     &(imap(280),iaoo4_9 ), (imap(281),iaoo4_10 ),(imap(282),iaoo5_1 ),
     &(imap(283),iaoo5_2 ), (imap(284),iaoo5_3 ), (imap(285),iaoo5_4 ),
     &(imap(286),iaoo5_5 ), (imap(287),iaoo5_6 ), (imap(288),iaoo5_7 ),
     &(imap(289),iaoo5_8 ), (imap(290),iaoo5_9 ), (imap(291),iaoo5_10 ),
     &(imap(292),iaoo6_1 ), (imap(293),iaoo6_2 ), (imap(294),iaoo6_3 ),
     &(imap(295),iaoo6_4 ), (imap(296),iaoo6_5 ), (imap(297),iaoo6_6 ),
     &(imap(298),iaoo6_7 ), (imap(299),iaoo6_8 ), (imap(300),iaoo6_9 ),
     &(imap(301),iaoo6_10 ),(imap(302),iaoo7_1 ), (imap(303),iaoo7_2 ),
     &(imap(304),iaoo7_3 ), (imap(305),iaoo7_4 ), (imap(306),iaoo7_5 ),
     &(imap(307),iaoo7_6 ), (imap(308),iaoo7_7 ), (imap(309),iaoo7_8 ),
     &(imap(310),iaoo7_9 ), (imap(311),iaoo7_10 ),(imap(312),iaoo8_1 ),
     &(imap(313),iaoo8_2 ), (imap(314),iaoo8_3 ), (imap(315),iaoo8_4 ),
     &(imap(316),iaoo8_5 ), (imap(317),iaoo8_6 ), (imap(318),iaoo8_7 ),
     &(imap(319),iaoo8_8 ), (imap(320),iaoo8_9 ), (imap(321),iaoo8_10 ),
     &(imap(322),iabs1_1 ),
     &(imap(323),iabs1_2 ), (imap(324),iabs1_3 ), (imap(325),iabs1_4 ),
     &(imap(326),iabs1_5 ), (imap(327),iabs1_6 ), (imap(328),iabs1_7 ),
     &(imap(329),iabs1_8 ), (imap(330),iabs1_9 ), (imap(331),iabs1_10 ),
     &(imap(332),iabs2_1 ), (imap(333),iabs2_2 ), (imap(334),iabs2_3 ),
     &(imap(335),iabs2_4 ), (imap(336),iabs2_5 ), (imap(337),iabs2_6 ),
     &(imap(338),iabs2_7 ), (imap(339),iabs2_8 ), (imap(340),iabs2_9 ),
     &(imap(341),iabs2_10 ),(imap(342),iabs3_1 ), (imap(343),iabs3_2 ),
     &(imap(344),iabs3_3 ), (imap(345),iabs3_4 ), (imap(346),iabs3_5 ),
     &(imap(347),iabs3_6 ), (imap(348),iabs3_7 ), (imap(349),iabs3_8 ),
     &(imap(350),iabs3_9 ), (imap(351),iabs3_10 ),(imap(352),iabs4_1 ),
     &(imap(353),iabs4_2 ), (imap(354),iabs4_3 ), (imap(355),iabs4_4 ),
     &(imap(356),iabs4_5 ), (imap(357),iabs4_6 ), (imap(358),iabs4_7 ),
     &(imap(359),iabs4_8 ), (imap(360),iabs4_9 ), (imap(361),iabs4_10 ),
     &(imap(362),iabs5_1 ),
     &(imap(363),iabs5_2 ), (imap(364),iabs5_3 ), (imap(365),iabs5_4 ),
     &(imap(366),iabs5_5 ), (imap(367),iabs5_6 ), (imap(368),iabs5_7 ),
     &(imap(369),iabs5_8 ), (imap(370),iabs5_9 ), (imap(371),iabs5_10 ),
     &(imap(372),iaas1_1 ), (imap(373),iaas1_2 ), (imap(374),iaas1_3 ),
     &(imap(375),iaas1_4 ), (imap(376),iaas1_5 ), (imap(377),iaas1_6 ),
     &(imap(378),iaas1_7 ), (imap(379),iaas1_8 ), (imap(380),iaas1_9 ),
     &(imap(381),iaas1_10 ),(imap(382),iaas2_1 ), (imap(383),iaas2_2 ),
     &(imap(384),iaas2_3 ), (imap(385),iaas2_4 ), (imap(386),iaas2_5 ),
     &(imap(387),iaas2_6 ), (imap(388),iaas2_7 ), (imap(389),iaas2_8 ),
     &(imap(390),iaas2_9 ), (imap(391),iaas2_10 ),(imap(392),iaas3_1 ),
     &(imap(393),iaas3_2 ), (imap(394),iaas3_3 ), (imap(395),iaas3_4 ),
     &(imap(396),iaas3_5 ), (imap(397),iaas3_6 ), (imap(398),iaas3_7 ),
     &(imap(399),iaas3_8 ), (imap(400),iaas3_9 ), (imap(401),iaas3_10 ),
     &(imap(402),iaas4_1 ), (imap(403),iaas4_2 ), (imap(404),iaas4_3 ),
     &(imap(405),iaas4_4 ), (imap(406),iaas4_5 ), (imap(407),iaas4_6 ),
     &(imap(408),iaas4_7 ), (imap(409),iaas4_8 ), (imap(410),iaas4_9 ),
     &(imap(411),iaas4_10 ),
     &(imap(412),iaas5_1 ), (imap(413),iaas5_2 ), (imap(414),iaas5_3 ),
     &(imap(415),iaas5_4 ), (imap(416),iaas5_5 ), (imap(417),iaas5_6 ),
     &(imap(418),iaas5_7 ), (imap(419),iaas5_8 ), (imap(420),iaas5_9 ),
     &(imap(421),iaas5_10 ),
     &(imap(422),ians1_1 ), (imap(423),ians1_2 ), (imap(424),ians1_3 ),
     &(imap(425),ians1_4 ), (imap(426),ians1_5 ), (imap(427),ians1_6 ),
     &(imap(428),ians1_7 ), (imap(429),ians1_8 ), (imap(430),ians1_9 ),
     &(imap(431),ians1_10 ),(imap(432),ians2_1 ), (imap(433),ians2_2 ),
     &(imap(434),ians2_3 ), (imap(435),ians2_4 ), (imap(436),ians2_5 ),
     &(imap(437),ians2_6 ), (imap(438),ians2_7 ), (imap(439),ians2_8 ),
     &(imap(440),ians2_9 ), (imap(441),ians2_10 ),(imap(442),ians3_1 ),
     &(imap(443),ians3_2 ), (imap(444),ians3_3 ), (imap(445),ians3_4 ),
     &(imap(446),ians3_5 ), (imap(447),ians3_6 ), (imap(448),ians3_7 ),
     &(imap(449),ians3_8 ), (imap(450),ians3_9 ), (imap(451),ians3_10 ),
     &(imap(452),ians4_1 ), (imap(453),ians4_2 ), (imap(454),ians4_3 ),
     &(imap(455),ians4_4 ), (imap(456),ians4_5 ), (imap(457),ians4_6 ),
     &(imap(458),ians4_7 ), (imap(459),ians4_8 ), (imap(460),ians4_9 ),
     &(imap(461),ians4_10 ),(imap(462),ians5_1 ), (imap(463),ians5_2 ),
     &(imap(464),ians5_3 ), (imap(465),ians5_4 ), (imap(466),ians5_5 ),
     &(imap(467),ians5_6 ), (imap(468),ians5_7 ), (imap(469),ians5_8 ),
     &(imap(470),ians5_9 ), (imap(471),ians5_10 ),(imap(472),ians6_1 ),
     &(imap(473),ians6_2 ), (imap(474),ians6_3 ), (imap(475),ians6_4 ),
     &(imap(476),ians6_5 ), (imap(477),ians6_6 ), (imap(478),ians6_7 ),
     &(imap(479),ians6_8 ), (imap(480),ians6_9 ), (imap(481),ians6_10 ),
     &(imap(482),ians7_1 ), (imap(483),ians7_2 ), (imap(484),ians7_3 ),
     &(imap(485),ians7_4 ), (imap(486),ians7_5 ), (imap(487),ians7_6 ),
     &(imap(488),ians7_7 ), (imap(489),ians7_8 ), (imap(490),ians7_9 ),
     &(imap(491),ians7_10 ),(imap(492),ians8_1 ), (imap(493),ians8_2 ),
     &(imap(494),ians8_3 ), (imap(495),ians8_4 ), (imap(496),ians8_5 ),
     &(imap(497),ians8_6 ), (imap(498),ians8_7 ), (imap(499),ians8_8 ),
     &(imap(500),ians8_9 ), (imap(501),ians8_10 ),

c   scf

     &(imap(502),ipoc_1  ),
     &(imap(503),ipoc_2  ),(imap(504),ipoc_3  ), (imap(505),ipoc_4  ),
     &(imap(506),ipoc_5  ),(imap(507),ipoc_6  ), (imap(508),ipoc_7  ),
     &(imap(509),ipoc_8  ),(imap(510),ipoc_9  ), (imap(511),ipoc_10 ),
     &(imap(512),ipec_1  ),(imap(513),ipec_2  ), (imap(514),ipec_3  ),
     &(imap(515),ipec_4  ),(imap(516),ipec_5  ), (imap(517),ipec_6  ),
     &(imap(518),ipec_7  ),(imap(519),ipec_8  ), (imap(520),ipec_9  ),
     &(imap(521),ipec_10 ),(imap(522),icrust_1), (imap(523),icrust_2),
     &(imap(524),icrust_3),(imap(525),icrust_4), (imap(526),icrust_5),
     &(imap(527),icrust_6),(imap(528),icrust_7), (imap(529),icrust_8),
     &(imap(530),icrust_9),(imap(531),icrust_10),(imap(532),iph2o_1),
     &(imap(533),iph2o_2 ),(imap(534),iph2o_3 ), (imap(535),iph2o_4 ),
     &(imap(536),iph2o_5 ),(imap(537),iph2o_6 ), (imap(538),iph2o_7 ),
     &(imap(539),iph2o_8 ),(imap(540),iph2o_9 ), (imap(541),iph2o_10),
     &(imap(542),ipcl_1  ),(imap(543),ipcl_2  ), (imap(544),ipcl_3  ),
     &(imap(545),ipcl_4  ),(imap(546),ipcl_5  ), (imap(547),ipcl_6  ),
     &(imap(548),ipcl_7  ),(imap(549),ipcl_8  ), (imap(550),ipcl_9  ),
     &(imap(551),ipcl_10 ),(imap(552),ina_1   ), (imap(553),ina_2   ),
     &(imap(554),ina_3   ),(imap(555),ina_4   ), (imap(556),ina_5   ),
     &(imap(557),ina_6   ),(imap(558),ina_7   ), (imap(559),ina_8   ),
     &(imap(560),ina_9   ),(imap(561),ina_10  ), (imap(562),ipnh4_1 ),
     &(imap(563),ipnh4_2 ),(imap(564),ipnh4_3 ), (imap(565),ipnh4_4 ),
     &(imap(566),ipnh4_5 ),(imap(567),ipnh4_6 ), (imap(568),ipnh4_7 ),
     &(imap(569),ipnh4_8 ),(imap(570),ipnh4_9 ), (imap(571),ipnh4_10),
     &(imap(572),ipno3_1 ),(imap(573),ipno3_2 ), (imap(574),ipno3_3 ),
     &(imap(575),ipno3_4 ),(imap(576),ipno3_5 ), (imap(577),ipno3_6 ),
     &(imap(578),ipno3_7 ),(imap(579),ipno3_8 ), (imap(580),ipno3_9 ),
     &(imap(581),ipno3_10),(imap(582),ipso4_1 ), (imap(583),ipso4_2 ),
     &(imap(584),ipso4_3 ),(imap(585),ipso4_4 ), (imap(586),ipso4_5 ),
     &(imap(587),ipso4_6 ),(imap(588),ipso4_7 ), (imap(589),ipso4_8 ),
     &(imap(590),ipso4_9 ),(imap(591),ipso4_10), (imap(592),iph2o   )
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

