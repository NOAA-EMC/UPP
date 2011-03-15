    module RQSTFLD_mod
!--------------------------------------------------------------------
      implicit none
!
      INTEGER, PARAMETER :: MXFLD=510,MXLVL=70
      CHARACTER*20 AVBL(MXFLD),FIELD(MXFLD)
      CHARACTER*6 DATSET      
!
      LOGICAL RITEHD,RITE2
!
      integer :: KGTYPE,IOUTYP,SVALUE,NFLD,IGET(MXFLD),           &
                 IQ(MXFLD),IS(MXFLD),ISMSTG(MXFLD),               &
                 ISMFUL(MXFLD),ISMOUT(MXFLD),LVLS(MXLVL,MXFLD),   &
                 IDENT(MXFLD),IFILV(MXFLD),ID(25),IGDS(18)
      real    :: DEC(MXFLD)

!initialization
!
!     THIS FILE CONTAINS ALL THE UNIQUE FIELDS THE
!     ETA POST PROCESSOR CAN CURRENTLY GENERATE.  
!
!	IFILV IS FLAG FOR IDENTIFYING MASS OR VELOCITY POINT
!	   =0 DATA IS VELOCITY POINT
!	   =1 DATA IS MASS POINT
!	AVBL IS CHARACTER STRING IDENTIFYING THE FIELD.
!	IQ  IS THE GRIB PDS OCTET 9 - PARAMETER (TABLE 2)
!	IS  IS THE GRIB PDS OCTET 10 - LEVEL TYPE (TABLE 3 & 3a)
!
!     WANT MORE/DIFFERENT FIELDS? 
!	(1) ADD CODE TO CALCULATE FIELD(S) IN APPROPRIATE ROUTINE(S),
!       (2) ADD FIELD(S) TO THIS LIST WITH A UNIQUE ITAG TAG,
!       (3) EDIT INPUT (CONTROL) FILE ACCORDINGLY,
!       (3) INCREASE PARAMETER MXFLD IN COMMON BLOCK RQSTFLD.comm.
!
!     CURRENT NUMBER OF FIELDS LISTED:  180
!
!0       1         2         3         4         5         6         7
!234567890123456789012345678901234567890123456789012345678901234567890
!
      DATA IFILV(001),AVBL(001),IQ(001),IS(001)       &
     &                      /1,'PRESS ON MDL SFCS   ',001,109/
      DATA IFILV(077),AVBL(077),IQ(077),IS(077)       &
     &                      /1,'HEIGHT ON MDL SFCS  ',007,109/
      DATA IFILV(002),AVBL(002),IQ(002),IS(002)       &
     &                      /1,'TEMP ON MDL SFCS    ',011,109/
      DATA IFILV(003),AVBL(003),IQ(003),IS(003)       &
     &                      /1,'POT TEMP ON MDL SFCS',013,109/
      DATA IFILV(004),AVBL(004),IQ(004),IS(004)       &
     &                      /1,'DWPT TEMP ON MDL SFC',017,109/
      DATA IFILV(005),AVBL(005),IQ(005),IS(005)       &
     &                      /1,'SPEC HUM ON MDL SFCS',051,109/
      DATA IFILV(006),AVBL(006),IQ(006),IS(006)       &
     &                      /1,'REL HUM ON MDL SFCS ',052,109/
      DATA IFILV(083),AVBL(083),IQ(083),IS(083)       &
     &                      /1,'MST CNVG ON MDL SFCS',135,109/
      DATA IFILV(007),AVBL(007),IQ(007),IS(007)       &
     &                      /0,'U WIND ON MDL SFCS  ',033,109/
      DATA IFILV(008),AVBL(008),IQ(008),IS(008)       &
     &                      /0,'V WIND ON MDL SFCS  ',034,109/
      DATA IFILV(009),AVBL(009),IQ(009),IS(009)       &
     &                      /1,'OMEGA ON MDL SFCS   ',039,109/
      DATA IFILV(010),AVBL(010),IQ(010),IS(010)       &
     &                      /1,'ABS VORT ON MDL SFCS',041,109/
      DATA IFILV(084),AVBL(084),IQ(084),IS(084)       &
     &                      /1,'STRMFUNC ON MDL SFCS',035,109/
      DATA IFILV(011),AVBL(011),IQ(011),IS(011)       &
     &                      /1,'TRBLNT KE ON MDL SFC',158,109/
      DATA IFILV(111),AVBL(111),IQ(111),IS(111)       &
     &                      /1,'RCHDSN NO ON MDL SFC',254,109/
      DATA IFILV(146),AVBL(146),IQ(146),IS(146)       &
     &                      /1,'MASTER LENGTH SCALE ',226,109/
      DATA IFILV(147),AVBL(147),IQ(147),IS(147)       &
     &                      /1,'ASYMPT MSTR LEN SCL ',227,109/
      DATA IFILV(012),AVBL(012),IQ(012),IS(012)       &
     &                      /1,'HEIGHT OF PRESS SFCS',007,100/
      DATA IFILV(013),AVBL(013),IQ(013),IS(013)       &
     &                      /1,'TEMP ON PRESS SFCS  ',011,100/
      DATA IFILV(014),AVBL(014),IQ(014),IS(014)       &
     &                      /1,'POT TEMP ON P SFCS  ',013,100/
      DATA IFILV(015),AVBL(015),IQ(015),IS(015)       &
     &                      /1,'DWPT TEMP ON P SFCS ',017,100/
      DATA IFILV(016),AVBL(016),IQ(016),IS(016)       &
     &                      /1,'SPEC HUM ON P SFCS  ',051,100/
      DATA IFILV(017),AVBL(017),IQ(017),IS(017)       &
     &                      /1,'REL HUMID ON P SFCS ',052,100/
      DATA IFILV(085),AVBL(085),IQ(085),IS(085)       &
     &                      /1,'MST CNVG ON P SFCS  ',135,100/
      DATA IFILV(018),AVBL(018),IQ(018),IS(018)       &
     &                      /1,'U WIND ON PRESS SFCS',033,100/
      DATA IFILV(019),AVBL(019),IQ(019),IS(019)       &
     &                      /1,'V WIND ON PRESS SFCS',034,100/
      DATA IFILV(020),AVBL(020),IQ(020),IS(020)       &
     &                      /1,'OMEGA ON PRESS SFCS ',039,100/
      DATA IFILV(021),AVBL(021),IQ(021),IS(021)       &
     &                      /1,'ABS VORT ON P SFCS  ',041,100/
      DATA IFILV(086),AVBL(086),IQ(086),IS(086)       &
     &                      /1,'STRMFUNC ON P SFCS  ',035,100/
      DATA IFILV(022),AVBL(022),IQ(022),IS(022)       &
     &                      /1,'TRBLNT KE ON P SFCS ',158,100/
      DATA IFILV(153),AVBL(153),IQ(153),IS(153)       &
     &                      /1,'CLOUD WATR ON P SFCS',153,100/
      DATA IFILV(166),AVBL(166),IQ(166),IS(166)       &
     &                      /1,'CLOUD ICE ON P SFCS ',058,100/
      DATA IFILV(023),AVBL(023),IQ(023),IS(023)       &
     &                      /1,'MESINGER MEAN SLP   ',130,102/
      DATA IFILV(105),AVBL(105),IQ(105),IS(105)       &
     &                      /1,'SHUELL MEAN SLP     ',002,102/
      DATA IFILV(445),AVBL(445),IQ(445),IS(445)       &
     &                      /1,'MAPS SLP            ',002,102/
      DATA IFILV(138),AVBL(138),IQ(138),IS(138)       &
     &                      /1,'SHELTER PRESSURE    ',001,105/
      DATA IFILV(106),AVBL(106),IQ(106),IS(106)       &
     &                      /1,'SHELTER TEMPERATURE ',011,105/
      DATA IFILV(112),AVBL(112),IQ(112),IS(112)       &
     &                      /1,'SHELTER SPEC HUMID  ',051,105/
      DATA IFILV(414),AVBL(414),IQ(414),IS(414)       &
     &                      /1,'SHELTER MIX RATIO   ',053,105/
      DATA IFILV(113),AVBL(113),IQ(113),IS(113)       &
     &                      /1,'SHELTER DEWPOINT    ',017,105/
      DATA IFILV(114),AVBL(114),IQ(114),IS(114)       &
     &                      /1,'SHELTER REL HUMID   ',052,105/
      DATA IFILV(064),AVBL(064),IQ(064),IS(064)       &
     &                      /1,'U WIND AT ANEMOM HT ',033,105/
      DATA IFILV(065),AVBL(065),IQ(065),IS(065)       &
     &                      /1,'V WIND AT ANEMOM HT ',034,105/
      DATA IFILV(158),AVBL(158),IQ(158),IS(158)       &
     &                      /1,'POT TEMP AT 10 M    ',013,105/
      DATA IFILV(159),AVBL(159),IQ(159),IS(159)       &
     &                      /1,'SPEC HUM AT 10 M    ',051,105/
      DATA IFILV(024),AVBL(024),IQ(024),IS(024)       &
     &                      /1,'SURFACE PRESSURE    ',001,001/
      DATA IFILV(025),AVBL(025),IQ(025),IS(025)       &
     &                      /1,'SURFACE HEIGHT      ',007,001/
      DATA IFILV(027),AVBL(027),IQ(027),IS(027)       &
     &                      /1,'SURFACE POT TEMP    ',013,001/
      DATA IFILV(028),AVBL(028),IQ(028),IS(028)       &
     &                      /1,'SURFACE SPEC HUMID  ',051,001/
      DATA IFILV(029),AVBL(029),IQ(029),IS(029)       &
     &                      /1,'SURFACE DEWPOINT    ',017,001/
      DATA IFILV(076),AVBL(076),IQ(076),IS(076)       &
     &                      /1,'SURFACE REL HUMID   ',052,001/
      DATA IFILV(026),AVBL(026),IQ(026),IS(026)       &
     &                      /1,'SFC (SKIN) TEMPRATUR',011,001/
      DATA IFILV(115),AVBL(115),IQ(115),IS(115)       &
     &                      /1,'BOTTOM SOIL TEMP    ',085,111/
      DATA IFILV(116),AVBL(116),IQ(116),IS(116)       &
     &                      /1,'SOIL TEMPERATURE    ',085,112/
      DATA IFILV(117),AVBL(117),IQ(117),IS(117)       &
     &                      /1,'SOIL MOISTURE       ',144,112/
      DATA IFILV(036),AVBL(036),IQ(036),IS(036)       &
     &                      /1,'TOTAL SOIL MOISTURE ',086,112/
      DATA IFILV(118),AVBL(118),IQ(118),IS(118)       &
     &                      /1,'PLANT CANOPY SFC WTR',223,001/
      DATA IFILV(119),AVBL(119),IQ(119),IS(119)       &
     &                      /1,'SNOW WATER EQUIVALNT',065,001/
      DATA IFILV(120),AVBL(120),IQ(120),IS(120)       &
     &                      /1,'PERCENT SNOW COVER  ',238,001/
      DATA IFILV(169),AVBL(169),IQ(169),IS(169)       &
     &                      /1,'SFC EXCHANGE COEF   ',208,001/
      DATA IFILV(170),AVBL(170),IQ(170),IS(170)       &
     &                      /1,'GREEN VEG COVER     ',087,001/
      DATA IFILV(171),AVBL(171),IQ(171),IS(171)       &
     &                      /1,'SOIL MOISTURE AVAIL ',207,112/
      DATA IFILV(152),AVBL(152),IQ(152),IS(152)       &
     &                      /1,'INST GROUND HEAT FLX',155,001/
      DATA IFILV(030),AVBL(030),IQ(030),IS(030)       &
     &                      /1,'LIFTED INDEX--SURFCE',131,101/
      DATA IFILV(031),AVBL(031),IQ(031),IS(031)       &
     &                      /1,'LIFTED INDEX--BEST  ',132,116/
      DATA IFILV(075),AVBL(075),IQ(075),IS(075)       &
     &                      /1,'LIFTED INDEX--BNDLYR',024,116/
      DATA IFILV(032),AVBL(032),IQ(032),IS(032)       &
     &                      /1,'CNVCT AVBL POT ENRGY',157,001/
      DATA IFILV(107),AVBL(107),IQ(107),IS(107)       &
     &                      /1,'CNVCT INHIBITION    ',156,001/
      DATA IFILV(080),AVBL(080),IQ(080),IS(080)       &
     &                      /1,'PRECIPITABLE WATER  ',054,200/
      DATA IFILV(162),AVBL(162),IQ(162),IS(162)       &
     &                      /1,'STORM REL HELICITY  ',190,106/
      DATA IFILV(163),AVBL(163),IQ(163),IS(163)       &
     &                      /1,'U COMP STORM MOTION ',196,106/
      DATA IFILV(164),AVBL(164),IQ(164),IS(164)       &
     &                      /1,'V COMP STORM MOTION ',197,106/
      DATA IFILV(087),AVBL(087),IQ(087),IS(087)       &
     &                      /1,'ACM TOTAL PRECIP    ',061,001/
      DATA IFILV(033),AVBL(033),IQ(033),IS(033)       &
     &                      /1,'ACM CONVCTIVE PRECIP',063,001/
      DATA IFILV(034),AVBL(034),IQ(034),IS(034)       &
     &                      /1,'ACM GRD SCALE PRECIP',062,001/
      DATA IFILV(035),AVBL(035),IQ(035),IS(035)       &
     &                      /1,'ACM SNOWFALL        ',065,001/
      DATA IFILV(244),AVBL(244),IQ(244),IS(244)       &
     &                      /1,'ACM GRD SCALE SW ICE',079,001/
      DATA IFILV(121),AVBL(121),IQ(121),IS(121)       &
     &                      /1,'ACM SNOW TOTAL/MELT ',099,001/
      DATA IFILV(122),AVBL(122),IQ(122),IS(122)       &
     &                      /1,'ACM STORM SFC RNOFF ',235,001/
      DATA IFILV(123),AVBL(123),IQ(123),IS(123)       &
     &                      /1,'ACM BSFL-GDWR RNOFF ',234,001/
      DATA IFILV(160),AVBL(160),IQ(160),IS(160)       &
     &                      /1,'INSTANT PRECIP TYPE ',140,001/
      DATA IFILV(407),AVBL(407),IQ(407),IS(407)       &
     &                      /1,'GSD PRECIP TYPE     ',140,001/
      DATA IFILV(167),AVBL(167),IQ(167),IS(167)       &
     &                      /1,'INSTANT PRECIP RATE ',059,001/
      DATA IFILV(172),AVBL(172),IQ(172),IS(172)       &
     &                      /1,'FROZEN FRAC CLD SCHM',194,001/
      DATA IFILV(124),AVBL(124),IQ(124),IS(124)       &
     &                      /1,'CLD WTR ON MDL SFCS ',153,109/
      DATA IFILV(125),AVBL(125),IQ(125),IS(125)       &
     &                      /1,'CLD ICE ON MDL SFCS ',058,109/
      DATA IFILV(145),AVBL(145),IQ(145),IS(145)       &
     &                      /1,'CLD FRAC ON MDL SFCS',071,109/
      DATA IFILV(037),AVBL(037),IQ(037),IS(037)       &
     &                      /1,'LOW CLOUD FRACTION  ',073,214/
      DATA IFILV(038),AVBL(038),IQ(038),IS(038)       &
     &                      /1,'MID CLOUD FRACTION  ',074,224/
      DATA IFILV(039),AVBL(039),IQ(039),IS(039)       &
     &                      /1,'HIGH CLOUD FRACTION ',075,234/
      DATA IFILV(161),AVBL(161),IQ(161),IS(161)       &
     &                      /1,'TOTAL CLD FRACTION  ',071,200/
      DATA IFILV(144),AVBL(144),IQ(144),IS(144)       &
     &                      /1,'AVG TOTAL CLD FRAC  ',071,200/
      DATA IFILV(139),AVBL(139),IQ(139),IS(139)       &
     &                      /1,'AVG STRAT CLD FRAC  ',213,200/
      DATA IFILV(143),AVBL(143),IQ(143),IS(143)       &
     &                      /1,'AVG CNVCT CLD FRAC  ',072,200/
      DATA IFILV(148),AVBL(148),IQ(148),IS(148)       &
     &                      /1,'CLOUD BOT PRESSURE  ',001,002/
      DATA IFILV(487),AVBL(487),IQ(487),IS(487)       &
     &                      /1,'GSD CLD BOT PRESSURE',001,002/
      DATA IFILV(149),AVBL(149),IQ(149),IS(149)       &
     &                      /1,'CLOUD TOP PRESSURE  ',001,003/
      DATA IFILV(406),AVBL(406),IQ(406),IS(406)       &
     &                      /1,'GSD CLD TOP PRESSURE',001,003/
      DATA IFILV(109),AVBL(109),IQ(109),IS(109)       &
     &                      /1,'LCL AGL HEIGHT      ',007,005/
      DATA IFILV(110),AVBL(110),IQ(110),IS(110)       &
     &                      /1,'LCL PRESSURE        ',001,005/
      DATA IFILV(078),AVBL(078),IQ(078),IS(078)       &
     &                      /1,'AVE GRDSCL RN TMPTDY',241,109/
      DATA IFILV(079),AVBL(079),IQ(079),IS(079)       &
     &                      /1,'AVE CNVCT RN TMPTDY ',242,109/
      DATA IFILV(168),AVBL(168),IQ(168),IS(168)     &
     &                      /1,'CLOUD TOP TEMPS     ',011,003/
      DATA IFILV(140),AVBL(140),IQ(140),IS(140)     &
     &                      /1,'RADFLX CNVG TMP TNDY',216,109/
      DATA IFILV(040),AVBL(040),IQ(040),IS(040)     &
     &                      /1,'SW RAD TEMP TNDY    ',250,109/
      DATA IFILV(041),AVBL(041),IQ(041),IS(041)     &
     &                      /1,'LW RAD TEMP TNDY    ',251,109/
      DATA IFILV(141),AVBL(141),IQ(141),IS(141)     &
     &                      /1,'INSTN OUT SFC SW RAD',211,001/
      DATA IFILV(142),AVBL(142),IQ(142),IS(142)     &
     &                      /1,'INSTN OUT SFC LW RAD',212,001/
      DATA IFILV(126),AVBL(126),IQ(126),IS(126)     &
     &                      /1,'AVE INCMG SFC SW RAD',204,001/
      DATA IFILV(127),AVBL(127),IQ(127),IS(127)     &
     &                      /1,'AVE INCMG SFC LW RAD',205,001/
      DATA IFILV(128),AVBL(128),IQ(128),IS(128)     &
     &                      /1,'AVE OUTGO SFC SW RAD',211,001/
      DATA IFILV(129),AVBL(129),IQ(129),IS(129)     &
     &                      /1,'AVE OUTGO SFC LW RAD',212,001/
      DATA IFILV(130),AVBL(130),IQ(130),IS(130)     &
     &                      /1,'AVE OUTGO TOA SW RAD',211,008/
      DATA IFILV(131),AVBL(131),IQ(131),IS(131)     &
     &                      /1,'AVE OUTGO TOA LW RAD',212,008/
      DATA IFILV(156),AVBL(156),IQ(156),IS(156)     &
     &                      /1,'INSTN INC SFC SW RAD',204,001/
      DATA IFILV(157),AVBL(157),IQ(157),IS(157)     &
     &                      /1,'INSTN INC SFC LW RAD',205,001/
      DATA IFILV(044),AVBL(044),IQ(044),IS(044)     &
     &                      /1,'ROUGHNESS LENGTH    ',083,001/
      DATA IFILV(045),AVBL(045),IQ(045),IS(045)     &
     &                      /1,'FRICTION VELOCITY   ',253,001/
      DATA IFILV(132),AVBL(132),IQ(132),IS(132)     &
     &                      /1,'SFC DRAG COEFFICIENT',252,001/
      DATA IFILV(133),AVBL(133),IQ(133),IS(133)     &
     &                      /1,'SFC U WIND STRESS   ',124,001/
      DATA IFILV(134),AVBL(134),IQ(134),IS(134)     &
     &                      /1,'SFC V WIND STRESS   ',125,001/
      DATA IFILV(043),AVBL(043),IQ(043),IS(043)     &
     &                      /1,'AVE SFC SENHEAT FX  ',122,001/
      DATA IFILV(135),AVBL(135),IQ(135),IS(135)     &
     &                      /1,'AVE GROUND HEAT FX  ',155,001/
      DATA IFILV(136),AVBL(136),IQ(136),IS(136)     &
     &                      /1,'AVE SNO PHSCNG HT FX',229,001/
      DATA IFILV(042),AVBL(042),IQ(042),IS(042)     &
     &                      /1,'AVE SFC LATHEAT FX  ',121,001/
      DATA IFILV(046),AVBL(046),IQ(046),IS(046)     &
     &                      /1,'AVE SFC MOMENTUM FX ',172,001/
      DATA IFILV(047),AVBL(047),IQ(047),IS(047)     &
     &                      /1,'ACC SFC EVAPORATION ',057,001/
      DATA IFILV(137),AVBL(137),IQ(137),IS(137)     &
     &                      /1,'ACC POT EVAPORATION ',228,001/
      DATA IFILV(154),AVBL(154),IQ(154),IS(154)     &
     &                      /1,'INST SFC SENHEAT FX ',122,001/
      DATA IFILV(155),AVBL(155),IQ(155),IS(155)     &
     &                      /1,'INST SFC LATHEAT FX ',121,001/
      DATA IFILV(048),AVBL(048),IQ(048),IS(048)     &
     &                      /1,'LATITUDE            ',176,001/
      DATA IFILV(049),AVBL(049),IQ(049),IS(049)     &
     &                      /1,'LONGITUDE           ',177,001/
      DATA IFILV(050),AVBL(050),IQ(050),IS(050)     &
     &                      /1,'LAND/SEA MASK       ',081,001/
      DATA IFILV(051),AVBL(051),IQ(051),IS(051)     &
     &                      /1,'SEA ICE MASK        ',091,001/
      DATA IFILV(052),AVBL(052),IQ(052),IS(052)     &
     &                      /1,'MASS POINT MDL SFC  ',173,001/
      DATA IFILV(053),AVBL(053),IQ(053),IS(053)     &
     &                      /1,'VEL POINT MDL SFC   ',174,001/
      DATA IFILV(150),AVBL(150),IQ(150),IS(150)     &
     &                      /1,'SFC MIDDAY ALBEDO   ',084,001/
      DATA IFILV(151),AVBL(151),IQ(151),IS(151)     &
     &                      /1,'SEA SFC TEMPERATURE ',080,001/
      DATA IFILV(054),AVBL(054),IQ(054),IS(054)     &
     &                      /1,'PRESS AT TROPOPAUSE ',001,007/
      DATA IFILV(055),AVBL(055),IQ(055),IS(055)     &
     &                      /1,'TEMP AT TROPOPAUSE  ',011,007/
      DATA IFILV(108),AVBL(108),IQ(108),IS(108)     &
     &                      /1,'POTENTL TEMP AT TROP',013,007/
      DATA IFILV(056),AVBL(056),IQ(056),IS(056)     &
     &                      /1,'U WIND AT TROPOPAUSE',033,007/
      DATA IFILV(057),AVBL(057),IQ(057),IS(057)     &
     &                      /1,'V WIND AT TROPOPAUSE',034,007/
      DATA IFILV(058),AVBL(058),IQ(058),IS(058)     &
     &                      /1,'SHEAR AT TROPOPAUSE ',136,007/
      DATA IFILV(059),AVBL(059),IQ(059),IS(059)     &
     &                      /1,'TEMP AT FD HEIGHTS  ',011,103/
      DATA IFILV(060),AVBL(060),IQ(060),IS(060)     &
     &                      /1,'U WIND AT FD HEIGHTS',033,103/
      DATA IFILV(061),AVBL(061),IQ(061),IS(061)     &
     &                      /1,'V WIND AT FD HEIGHTS',034,103/
      DATA IFILV(062),AVBL(062),IQ(062),IS(062)     &
     &                      /1,'HEIGHT OF FRZ LVL   ',007,004/
      DATA IFILV(063),AVBL(063),IQ(063),IS(063)     &
     &                      /1,'REL HUMID AT FRZ LVL',052,004/
      DATA IFILV(165),AVBL(165),IQ(165),IS(165)     &
     &                      /1,'HIGHEST FREEZE LVL  ',007,204/
      DATA IFILV(350),AVBL(350),IQ(350),IS(350)     &
     &                      /1,'HIGHEST FRZ LVL RH  ',052,204/
      DATA IFILV(067),AVBL(067),IQ(067),IS(067)     &
     &                      /1,'PRESS IN BNDRY LYR  ',001,116/
      DATA IFILV(068),AVBL(068),IQ(068),IS(068)     &
     &                      /1,'TEMP IN BNDRY LYR   ',011,116/
      DATA IFILV(069),AVBL(069),IQ(069),IS(069)     &
     &                      /1,'POT TMP IN BNDRY LYR',013,116/
      DATA IFILV(070),AVBL(070),IQ(070),IS(070)     &
     &                      /1,'DWPT IN BNDRY LYR   ',017,116/
      DATA IFILV(071),AVBL(071),IQ(071),IS(071)     &
     &                      /1,'SPC HUM IN BNDRY LYR',051,116/
      DATA IFILV(072),AVBL(072),IQ(072),IS(072)     &
     &                      /1,'REL HUM IN BNDRY LYR',052,116/
      DATA IFILV(088),AVBL(088),IQ(088),IS(088)     &
     &                      /1,'MST CNV IN BNDRY LYR',135,116/
      DATA IFILV(089),AVBL(089),IQ(089),IS(089)     &
     &                      /1,'P WATER IN BNDRY LYR',054,116/
      DATA IFILV(073),AVBL(073),IQ(073),IS(073)     &
     &                      /1,'U WIND IN BNDRY LYR ',033,116/
      DATA IFILV(074),AVBL(074),IQ(074),IS(074)     &
     &                      /1,'V WIND IN BNDRY LYR ',034,116/
      DATA IFILV(090),AVBL(090),IQ(090),IS(090)     &
     &                      /1,'OMEGA IN BNDRY LYR  ',039,116/
      DATA IFILV(066),AVBL(066),IQ(066),IS(066)     &
     &                      /1,'LFM 0.33-1.00 RELHUM',052,108/
      DATA IFILV(081),AVBL(081),IQ(081),IS(081)     &
     &                      /1,'LFM 0.66-1.00 RELHUM',052,108/
      DATA IFILV(082),AVBL(082),IQ(082),IS(082)     &
     &                      /1,'LFM 0.33-0.66 RELHUM',052,108/
      DATA IFILV(104),AVBL(104),IQ(104),IS(104)     &
     &                      /1,'LFM 0.33-1.00 PWAT  ',054,108/
      DATA IFILV(091),AVBL(091),IQ(091),IS(091)     &
     &                      /1,'NGM 0.98230 PRESSURE',001,107/
      DATA IFILV(092),AVBL(092),IQ(092),IS(092)     &
     &                      /1,'NGM 0.98230 TMPRATUR',011,107/
      DATA IFILV(093),AVBL(093),IQ(093),IS(093)     &
     &                      /1,'NGM 0.98230 SPC HUM ',051,107/
      DATA IFILV(094),AVBL(094),IQ(094),IS(094)     &
     &                      /1,'NGM 0.98230 REL HUM ',052,107/
      DATA IFILV(095),AVBL(095),IQ(095),IS(095)     &
     &                      /1,'NGM 0.98230 U WIND  ',033,107/
      DATA IFILV(096),AVBL(096),IQ(096),IS(096)     &
     &                      /1,'NGM 0.98230 V WIND  ',034,107/
      DATA IFILV(097),AVBL(097),IQ(097),IS(097)     &
     &                      /1,'NGM 0.89671 TMPRATUR',011,107/
      DATA IFILV(098),AVBL(098),IQ(098),IS(098)     &
     &                      /1,'NGM 0.78483 TMPRATUR',011,107/
      DATA IFILV(099),AVBL(099),IQ(099),IS(099)     &
     &                      /1,'NGM 0.47-1.00 RELHUM',052,108/
      DATA IFILV(100),AVBL(100),IQ(100),IS(100)     &
     &                      /1,'NGM 0.47-0.96 RELHUM',052,108/
      DATA IFILV(101),AVBL(101),IQ(101),IS(101)     &
     &                      /1,'NGM 0.18-0.47 RELHUM',052,108/
      DATA IFILV(102),AVBL(102),IQ(102),IS(102)     &
     &                      /1,'NGM 0.84-0.98 RELHUM',052,108/
      DATA IFILV(103),AVBL(103),IQ(103),IS(103)     &
     &                      /1,'NGM 0.85-1.00 QCONVG',135,108/
      DATA IFILV(173),AVBL(173),IQ(173),IS(173)     &
     &                      /1,'MAX WIND PRESS LEVEL',001,006/
      DATA IFILV(174),AVBL(174),IQ(174),IS(174)     &
     &                      /1,'MAX WIND HGHT LEVEL ',007,006/
      DATA IFILV(175),AVBL(175),IQ(175),IS(175)     &
     &                      /1,'U COMP MAX WIND     ',033,006/
      DATA IFILV(176),AVBL(176),IQ(176),IS(176)     &
     &                      /1,'V COMP MAX WIND     ',034,006/
      DATA IFILV(177),AVBL(177),IQ(177),IS(177)     &
     &                      /1,'HEIGHT AT TROPOPAUSE',007,007/
      DATA IFILV(178),AVBL(178),IQ(178),IS(178)     &
     &                      /1,'CLOUD BOTTOM HEIGHT ',007,002/
      DATA IFILV(408),AVBL(408),IQ(408),IS(408)     &
     &                      /1,'GSD CLD BOT HEIGHT  ',007,002/
      DATA IFILV(179),AVBL(179),IQ(179),IS(179)     &
     &                      /1,'CLOUD TOP HEIGHT    ',007,003/
      DATA IFILV(409),AVBL(409),IQ(409),IS(409)     &
     &                      /1,'GSD CLD TOP HEIGHT  ',007,003/
      DATA IFILV(180),AVBL(180),IQ(180),IS(180)     &
     &                      /1,'VISIBILITY          ',020,001/
      DATA IFILV(410),AVBL(410),IQ(410),IS(410)     &
     &                      /1,'GSD VISIBILITY      ',020,001/
! CRA
      DATA IFILV(411),AVBL(411),IQ(411),IS(411)     &
     &                      /1,'INSTN WIND POWER AGL',126,105/
      DATA IFILV(412),AVBL(412),IQ(412),IS(412)     &
     &                      /1,'U WIND AT 80M AGL   ',049,105/
      DATA IFILV(413),AVBL(413),IQ(413),IS(413)     &
     &                      /1,'V WIND AT 80M AGL   ',050,105/
! CRA
      DATA IFILV(181),AVBL(181),IQ(181),IS(181)     &
     &                      /1,'RAIN ON MDL SFCS    ',170,109/
      DATA IFILV(182),AVBL(182),IQ(182),IS(182)     &
     &                      /1,'SNOW ON MDL SFCS    ',171,109/
      DATA IFILV(183),AVBL(183),IQ(183),IS(183)     &
     &                      /1,'RAIN ON P SFCS      ',170,100/
      DATA IFILV(184),AVBL(184),IQ(184),IS(184)     &
     &                      /1,'SNOW ON P SFCS      ',171,100/
      DATA IFILV(415),AVBL(415),IQ(415),IS(415)     &
     &                      /1,'GRAUPEL ON MDL SFCS ',179,109/
      DATA IFILV(416),AVBL(416),IQ(416),IS(416)     &
     &                      /1,'GRAUPEL ON P SFCS   ',179,100/

! SRD
      DATA IFILV(420),AVBL(420),IQ(420),IS(420)     &
     &                      /1,'MAX UPDRAFT HELICITY',215,106/
      DATA IFILV(421),AVBL(421),IQ(421),IS(421)     &
     &                      /1,'MAX 1km REFLECTIVITY',217,105/
      DATA IFILV(422),AVBL(422),IQ(422),IS(422)     &
     &                      /1,'MAX 10m WIND SPEED  ',229,105/
      DATA IFILV(423),AVBL(423),IQ(423),IS(423)     &
     &                      /1,'MAX UPDRAFT VERT VEL',220,106/
      DATA IFILV(424),AVBL(424),IQ(424),IS(424)     &
     &                      /1,'MAX DNDRAFT VERT VEL',223,106/
      DATA IFILV(425),AVBL(425),IQ(425),IS(425)     &
     &                      /1,'MEAN VERT VEL       ',040,108/
      DATA IFILV(426),AVBL(426),IQ(426),IS(426)     &
     &                      /1,'ECHO TOPS IN KFT    ',007,105/
      DATA IFILV(427),AVBL(427),IQ(427),IS(427)     &
     &                      /1,'UPDRAFT HELICITY PRM',227,106/
      DATA IFILV(428),AVBL(428),IQ(428),IS(428)     &
     &                      /1,'VERT INTEG GRAUP    ',179,200/
      DATA IFILV(429),AVBL(429),IQ(429),IS(429)     & 
     &                      /1,'MAX VERT INTEG GRAUP',228,200/
! SRD
! CRA
      DATA IFILV(430),AVBL(430),IQ(430),IS(430)     &
     &                      /1,'U COMP 0-1 KM SHEAR ',230,106/
      DATA IFILV(431),AVBL(431),IQ(431),IS(431)     &
     &                      /1,'V COMP 0-1 KM SHEAR ',238,106/
      DATA IFILV(432),AVBL(432),IQ(432),IS(432)     &
     &                      /1,'U COMP 0-6 KM SHEAR ',239,106/
      DATA IFILV(433),AVBL(433),IQ(433),IS(433)     &
     &                      /1,'V COMP 0-6 KM SHEAR ',241,106/
! CRA

! Add precipitation buckets between outputs
      DATA IFILV(434),AVBL(434),IQ(434),IS(434)       &
     &                      /1,'BUCKET TOTAL PRECIP ',061,001/
      DATA IFILV(435),AVBL(435),IQ(435),IS(435)       &
     &                      /1,'BUCKET CONV PRECIP  ',063,001/
      DATA IFILV(436),AVBL(436),IQ(436),IS(436)       &
     &                      /1,'BUCKET GRDSCALE PRCP',062,001/
      DATA IFILV(437),AVBL(437),IQ(437),IS(437)       &
     &                      /1,'BUCKET SNOW  PRECIP ',065,001/

!
!--- Added new cloud microphysics fields & displaying more
!    convective cloud properties  (Jin, '01;  Ferrier, Feb '02)     
!
!
!--- The following fields have been added to the post under
!    PDS Octet 4 = 129.  All other fields above are with PDS Octet
!    4 = 2.  Most of the fields below, except for the cloud top
!    and cloud base pressures, have PDS Octet 4 = 129.  These new
!    grib parameters are listed in Table 129 of the GRIB documentation.
!    See Table 2 in Office Note 388 (ON388) for more details.
!
!--- F_rain, F_ice, F_RimeF => PDS Octet 4 = 129
!
      DATA IFILV(185),AVBL(185),IQ(185),IS(185)     &
     &                      /1,'F_rain on MDL SFCS  ',131,109/
      DATA IFILV(186),AVBL(186),IQ(186),IS(186)     &
     &                      /1,'F_ice ON MDL SFCS   ',132,109/
      DATA IFILV(187),AVBL(187),IQ(187),IS(187)     &
     &                      /1,'F_RimeF ON MDL SFCS ',133,109/
!
!--- The following cloud pressure fields have PDS Octet 4 = 2
!
      DATA IFILV(188),AVBL(188),IQ(188),IS(188)     &
     &                      /1,'CONV CLOUD BOT PRESS',001,242/
      DATA IFILV(189),AVBL(189),IQ(189),IS(189)     &
     &                      /1,'CONV CLOUD TOP PRESS',001,243/
      DATA IFILV(190),AVBL(190),IQ(190),IS(190)     &
     &                      /1,'SHAL CU CLD BOT PRES',001,248/
      DATA IFILV(191),AVBL(191),IQ(191),IS(191)     &
     &                      /1,'SHAL CU CLD TOP PRES',001,249/
      DATA IFILV(192),AVBL(192),IQ(192),IS(192)     &
     &                      /1,'DEEP CU CLD BOT PRES',001,251/
      DATA IFILV(193),AVBL(193),IQ(193),IS(193)     &
     &                      /1,'DEEP CU CLD TOP PRES',001,252/
      DATA IFILV(194),AVBL(194),IQ(194),IS(194)     &
     &                      /1,'GRID CLOUD BOT PRESS',001,206/
      DATA IFILV(195),AVBL(195),IQ(195),IS(195)     &
     &                      /1,'GRID CLOUD TOP PRESS',001,207/
      DATA IFILV(196),AVBL(196),IQ(196),IS(196)     &
     &                      /1,'CONV CLOUD FRACTION ',072,200/
!
!--- These remaining fields have PDS Octet 4 = 129 (Table 129, ON388)      
!
      DATA IFILV(197),AVBL(197),IQ(197),IS(197)     &
     &                      /1,'CU CLOUD EFFICIENCY ',134,200/
      DATA IFILV(198),AVBL(198),IQ(198),IS(198)     &
     &                      /1,'CONDENSATE ON P SFCS',135,100/
      DATA IFILV(199),AVBL(199),IQ(199),IS(199)     &
     &                      /1,'CONDENSATE MDL SFCS ',135,109/
      DATA IFILV(200),AVBL(200),IQ(200),IS(200)     &
     &                      /1,'TOTAL COLUMN CLD WTR',136,200/
      DATA IFILV(201),AVBL(201),IQ(201),IS(201)     &
     &                      /1,'TOTAL COLUMN CLD ICE',137,200/
      DATA IFILV(202),AVBL(202),IQ(202),IS(202)     &
     &                      /1,'TOTAL COLUMN RAIN   ',138,200/
      DATA IFILV(203),AVBL(203),IQ(203),IS(203)     &
     &                      /1,'TOTAL COLUMN SNOW   ',139,200/
      DATA IFILV(204),AVBL(204),IQ(204),IS(204)     &
     &                      /1,'TOTAL COL CONDENSATE',140,200/
! See below for total supercooled liquid & melting ice ... IFILV(285)     
! H CHUANG--ADD INTERPOLATED FIELDS ON SIGMA LEVELS
      DATA IFILV(205),AVBL(205),IQ(205),IS(205)     &
     &                      /1,'HEIGHT OF SIGMA SFCS',007,107/
      DATA IFILV(206),AVBL(206),IQ(206),IS(206)     &
     &                      /1,'TEMP ON SIGMA SFCS  ',011,107/
      DATA IFILV(207),AVBL(207),IQ(207),IS(207)     &
     &                      /1,'SPEC HUM ON S SFCS  ',051,107/
      DATA IFILV(208),AVBL(208),IQ(208),IS(208)     &
     &                      /0,'U WIND ON SIGMA SFCS',033,107/
      DATA IFILV(209),AVBL(209),IQ(209),IS(209)     &
     &                      /0,'V WIND ON SIGMA SFCS',034,107/
      DATA IFILV(210),AVBL(210),IQ(210),IS(210)     &
     &                      /1,'OMEGA ON SIGMA SFCS ',039,107/
      DATA IFILV(211),AVBL(211),IQ(211),IS(211)     &
     &                      /1,'CLOUD WATR ON S SFCS',153,107/
      DATA IFILV(212),AVBL(212),IQ(212),IS(212)     &
     &                      /1,'CLOUD ICE ON S SFCS ',058,107/
      DATA IFILV(213),AVBL(213),IQ(213),IS(213)     &
     &                      /1,'RAIN ON S SFCS      ',170,107/
      DATA IFILV(214),AVBL(214),IQ(214),IS(214)     &
     &                      /1,'SNOW ON S SFCS      ',171,107/
      DATA IFILV(215),AVBL(215),IQ(215),IS(215)     &
     &                      /1,'CONDENSATE ON S SFCS',135,107/
      DATA IFILV(216),AVBL(216),IQ(216),IS(216)     &
     &                      /1,'PRESS ON SIG SFCS   ',001,107/
      DATA IFILV(217),AVBL(217),IQ(217),IS(217)     &
     &                      /1,'TRBLNT KE ON S SFCS ',158,107/
      DATA IFILV(222),AVBL(222),IQ(222),IS(222)     &
     &                      /1,'CLD FRAC ON SIG SFCS',071,107/
      DATA IFILV(255),AVBL(255),IQ(255),IS(255)     &
     &                      /1,'GRAUPEL ON S SFCS   ',179,107/
! H CHUANG--ADD FIXED AND LSM FIELDS
      DATA IFILV(218),AVBL(218),IQ(218),IS(218)     &
     &                      /1,'VEGETATION TYPE     ',225,001/
      DATA IFILV(219),AVBL(219),IQ(219),IS(219)     &
     &                      /1,'SOIL TYPE           ',224,001/
      DATA IFILV(220),AVBL(220),IQ(220),IS(220)     &
     &                      /1,'CANOPY CONDUCTANCE  ',181,001/
      DATA IFILV(221),AVBL(221),IQ(221),IS(221)     &
     &                      /1,'PBL HEIGHT          ',221,001/
      DATA IFILV(223),AVBL(223),IQ(223),IS(223)     &
     &                      /1,'SLOPE TYPE          ',222,001/
      DATA IFILV(224),AVBL(224),IQ(224),IS(224)     &
     &                      /1,'SNOW DEPTH          ',066,001/
      DATA IFILV(225),AVBL(225),IQ(225),IS(225)     &
     &                      /1,'LIQUID SOIL MOISTURE',160,112/
      DATA IFILV(226),AVBL(226),IQ(226),IS(226)     &
     &                      /1,'SNOW FREE ALBEDO    ',170,001/
      DATA IFILV(227),AVBL(227),IQ(227),IS(227)     &
     &                      /1,'MAXIMUM SNOW ALBEDO ',159,001/
      DATA IFILV(228),AVBL(228),IQ(228),IS(228)     &
     &                      /1,'CANOPY WATER EVAP   ',200,001/
      DATA IFILV(229),AVBL(229),IQ(229),IS(229)     &
     &                      /1,'DIRECT SOIL EVAP    ',199,001/
      DATA IFILV(230),AVBL(230),IQ(230),IS(230)     &
     &                      /1,'PLANT TRANSPIRATION ',210,001/
      DATA IFILV(231),AVBL(231),IQ(231),IS(231)     &
     &                      /1,'SNOW SUBLIMATION    ',198,001/
      DATA IFILV(232),AVBL(232),IQ(232),IS(232)     &
     &                      /1,'AIR DRY SOIL MOIST  ',231,001/
      DATA IFILV(233),AVBL(233),IQ(233),IS(233)     &
     &                      /1,'SOIL MOIST POROSITY ',240,001/
      DATA IFILV(234),AVBL(234),IQ(234),IS(234)     &
     &                      /1,'MIN STOMATAL RESIST ',203,001/
      DATA IFILV(235),AVBL(235),IQ(235),IS(235)     &
     &                      /1,'NO OF ROOT LAYERS   ',171,001/
      DATA IFILV(236),AVBL(236),IQ(236),IS(236)     &
     &                      /1,'SOIL MOIST WILT PT  ',219,001/
      DATA IFILV(237),AVBL(237),IQ(237),IS(237)     &
     &                      /1,'SOIL MOIST REFERENCE',230,001/
      DATA IFILV(238),AVBL(238),IQ(238),IS(238)     &
     &                      /1,'CANOPY COND SOLAR   ',246,001/
      DATA IFILV(239),AVBL(239),IQ(239),IS(239)     &
     &                      /1,'CANOPY COND TEMP    ',247,001/
      DATA IFILV(240),AVBL(240),IQ(240),IS(240)     &
     &                      /1,'CANOPY COND HUMID   ',248,001/
      DATA IFILV(241),AVBL(241),IQ(241),IS(241)     &
     &                      /1,'CANOPY COND SOILM   ',249,001/
      DATA IFILV(242),AVBL(242),IQ(242),IS(242)     &
     &                      /1,'POTENTIAL EVAP      ',145,001/
      DATA IFILV(243),AVBL(243),IQ(243),IS(243)     &
     &                      /1,'DIFFUSION H RATE S S',182,107/ 
!
      DATA IFILV(245),AVBL(245),IQ(245),IS(245)     &
     &                      /1,'SFC WIND GUST       ',180,001/
      DATA IFILV(246),AVBL(246),IQ(246),IS(246)     &
     &                      /1,'LIFT PCL LVL PRESS  ',141,116/
      DATA IFILV(247),AVBL(247),IQ(247),IS(247)     &
     &                      /1,'LOW WET BULB ZERO HT',007,245/

      DATA IFILV(248),AVBL(248),IQ(248),IS(248)     &
     &                      /1,'EMISSIVITY          ',193,001/

      DATA IFILV(249),AVBL(249),IQ(249),IS(249)     &
     &                      /1,'CONV PRECIP RATE    ',214,001/
!--- USING Table 129
!
      DATA IFILV(250),AVBL(250),IQ(250),IS(250)     &
     &                      /1,'RADAR REFL MDL SFCS ',211,109/
      DATA IFILV(251),AVBL(251),IQ(251),IS(251)     &
     &                      /1,'RADAR REFL ON P SFCS',211,100/
      DATA IFILV(252),AVBL(252),IQ(252),IS(252)     &
     &                      /1,'COMPOSITE RADAR REFL',212,200/
      DATA IFILV(253),AVBL(253),IQ(253),IS(253)     &
     &                      /1,'RADAR REFL AGL      ',211,105/
      DATA IFILV(254),AVBL(254),IQ(254),IS(254)     &
     &                      /1,'LEAF AREA INDEX     ',182,001/
!
      DATA IFILV(256),AVBL(256),IQ(256),IS(256)     &
     &                      /1,'ACM LSM PRECIP      ',154,001/
!
!--- FOLLOWINGS ARE AVIATION-RELATED FIELDS: ADDED BY BINBIN ZHOU
!
      DATA IFILV(257),AVBL(257),IQ(257),IS(257)     &
     &                      /1,'IN-FLIGHT ICING     ',186,100/
      DATA IFILV(258),AVBL(258),IQ(258),IS(258)     &
     &                      /1,'CLEAR AIR TURBULENCE',185,100/
      DATA IFILV(259),AVBL(259),IQ(259),IS(259)     &
     &                      /1,'0-2000FT LLWS       ',136,106/
      DATA IFILV(260),AVBL(260),IQ(260),IS(260)     &
     &                      /1,'CEILING             ',007,215/
      DATA IFILV(261),AVBL(261),IQ(261),IS(261)     &
     &                      /1,'FLIGHT RESTRICTION  ',020,002/
!
      DATA IFILV(262),AVBL(262),IQ(262),IS(262)     &
     &                      /1,'INSTN CLR INC SFC SW',161,001/
      DATA IFILV(263),AVBL(263),IQ(263),IS(263)     &
     &                      /1,'F_RimeF ON P SFCS   ',133,100/
      DATA IFILV(264),AVBL(264),IQ(264),IS(264)     &
     &                      /1,'W WIND ON MDL SFCS  ',040,109/
!      DATA IFILV(265),AVBL(265),IQ(265),IS(265)     &
!     &                      /1,'BRIGHTNESS TEMP     ',118,008/
      DATA IFILV(265),AVBL(265),IQ(265),IS(265)     &
     &                      /1,'BRIGHTNESS TEMP     ',213,008/     
! H CHUANG--ADD GFS products
      DATA IFILV(266),AVBL(266),IQ(266),IS(266)     &
     &                      /1,'AVE ALBEDO          ',084,001/
      DATA IFILV(267),AVBL(267),IQ(267),IS(267)     &
     &                      /1,'OZONE ON MDL SFCS   ',154,109/
      DATA IFILV(268),AVBL(268),IQ(268),IS(268)     &
     &                      /1,'OZONE ON P SFCS     ',154,100/
      DATA IFILV(269),AVBL(269),IQ(269),IS(269)     &
     &                      /1,'SFC ZONAL MOMEN FX  ',124,001/
      DATA IFILV(270),AVBL(270),IQ(270),IS(270)     &
     &                      /1,'SFC MERID MOMEN FX  ',125,001/
      DATA IFILV(271),AVBL(271),IQ(271),IS(271)     &
     &                      /1,'AVE PRECIP RATE     ',059,001/
      DATA IFILV(272),AVBL(272),IQ(272),IS(272)     &
     &                      /1,'AVE CONV PRECIP RATE',214,001/
! CMAQ requested fields     
      DATA IFILV(273),AVBL(273),IQ(273),IS(273)     &
     &                      /1,'HYBRID SIGMA DP     ',001,110/
      DATA IFILV(274),AVBL(274),IQ(274),IS(274)     &
     &                      /1,'INSTN OUT TOA LW RAD',212,008/
!      DATA IFILV(275),AVBL(275),IQ(275),IS(275)     &
!     &                      /1,'BRIGHTNESS TEMP NCAR',213,008/
      DATA IFILV(275),AVBL(275),IQ(275),IS(275)     &
     &                      /1,'BRIGHTNESS TEMP NCAR',118,008/
      DATA IFILV(282),AVBL(282),IQ(282),IS(282)     &
     &                      /1,'MODEL TOP PRESSURE  ',001,008/
      DATA IFILV(283),AVBL(283),IQ(283),IS(283)     &
     &                      /1,'HYBRID PRESSURE DP  ',001,110/
!
!--- USING Table 129
!
      DATA IFILV(276),AVBL(276),IQ(276),IS(276)     &
     &                      /1,'COMPOSITE RAIN REFL ',165,200/
      DATA IFILV(277),AVBL(277),IQ(277),IS(277)     &
     &                      /1,'COMPOSITE ICE REFL  ',166,200/
      DATA IFILV(278),AVBL(278),IQ(278),IS(278)     &
     &                      /1,'COMPOSITE CONV REFL ',167,200/
      DATA IFILV(279),AVBL(279),IQ(279),IS(279)     &
     &                      /1,'RAIN RADAR REFL AGL ',165,105/
      DATA IFILV(280),AVBL(280),IQ(280),IS(280)     &
     &                      /1,'ICE RADAR REFL AGL  ',166,105/
      DATA IFILV(281),AVBL(281),IQ(281),IS(281)     &
     &                      /1,'CONV RADAR REFL AGL ',167,105/
!
!--- USING Table 2
!
      DATA IFILV(284),AVBL(284),IQ(284),IS(284)     &
     &                      /1,'W WIND ON P SFCS    ',040,100/
!
!--- USING Table 129
!
      DATA IFILV(285),AVBL(285),IQ(285),IS(285)     &
     &                      /1,'TOTAL COLD LIQUID   ',168,200/
      DATA IFILV(286),AVBL(286),IQ(286),IS(286)     &
     &                      /1,'TOTAL MELTING ICE   ',169,200/
!
!--- USING Table 2
!
      DATA IFILV(287),AVBL(287),IQ(287),IS(287)     &
     &                      /1,'COLD LIQ BOT HEIGHT ',007,253/
      DATA IFILV(288),AVBL(288),IQ(288),IS(288)     &
     &                      /1,'COLD LIQ TOP HEIGHT ',007,254/     
      DATA IFILV(289),AVBL(289),IQ(289),IS(289)     &
     &                      /1,'RICH NO PBL HEIGHT  ',007,220/
!
!---- New Column-integrated fields
      DATA IFILV(290),AVBL(290),IQ(290),IS(290)     &
     &                      /1,'TOT COL SW T TNDY   ',250,200/
      DATA IFILV(291),AVBL(291),IQ(291),IS(291)     &
     &                      /1,'TOT COL LW T TNDY   ',251,200/
      DATA IFILV(292),AVBL(292),IQ(292),IS(292)     &
     &                      /1,'TOT COL GRD T TNDY  ',241,200/
      DATA IFILV(293),AVBL(293),IQ(293),IS(293)     &
     &                      /1,'TOT COL CNVCT T TNDY',242,200/
      DATA IFILV(294),AVBL(294),IQ(294),IS(294)     &
     &                      /1,'RADFLX TMP TNDY ON P',216,100/
      DATA IFILV(295),AVBL(295),IQ(295),IS(295)     &
     &                      /1,'TOT COL MST CNVG    ',135,200/
      DATA IFILV(296),AVBL(296),IQ(296),IS(296)     &
     &                      /1,'HPC T ON SIGMA SFCS ',011,107/
! H CHUANG--ADD GFS products
      DATA IFILV(297),AVBL(297),IQ(297),IS(297)     &
     &                      /1,'AVE CLR INC UV-B SW ',201,001/
      DATA IFILV(298),AVBL(298),IQ(298),IS(298)     &
     &                      /1,'AVE INC UV-B SW     ',200,001/
      DATA IFILV(299),AVBL(299),IQ(299),IS(299)     &
     &                      /1,'TOT COL OZONE       ',010,200/
      DATA IFILV(300),AVBL(300),IQ(300),IS(300)     &
     &                      /1,'AVE LOW CLOUD FRAC  ',071,214/
      DATA IFILV(301),AVBL(301),IQ(301),IS(301)     &
     &                      /1,'AVE MID CLOUD FRAC  ',071,224/
      DATA IFILV(302),AVBL(302),IQ(302),IS(302)     &
     &                      /1,'AVE HIGH CLOUD FRAC ',071,234/
      DATA IFILV(303),AVBL(303),IQ(303),IS(303)     &
     &                      /1,'AVE LOW CLOUD BOT P ',001,212/
      DATA IFILV(304),AVBL(304),IQ(304),IS(304)     &
     &                      /1,'AVE LOW CLOUD TOP P ',001,213/
      DATA IFILV(305),AVBL(305),IQ(305),IS(305)     &
     &                      /1,'AVE LOW CLOUD TOP T ',011,213/
      DATA IFILV(306),AVBL(306),IQ(306),IS(306)     &
     &                      /1,'AVE MID CLOUD BOT P ',001,222/
      DATA IFILV(307),AVBL(307),IQ(307),IS(307)     &
     &                      /1,'AVE MID CLOUD TOP P ',001,223/
      DATA IFILV(308),AVBL(308),IQ(308),IS(308)     &
     &                      /1,'AVE MID CLOUD TOP T ',011,223/
      DATA IFILV(309),AVBL(309),IQ(309),IS(309)     &
     &                      /1,'AVE HIGH CLOUD BOT P',001,232/
      DATA IFILV(310),AVBL(310),IQ(310),IS(310)     &
     &                      /1,'AVE HIGH CLOUD TOP P',001,233/
      DATA IFILV(311),AVBL(311),IQ(311),IS(311)     &
     &                      /1,'AVE HIGH CLOUD TOP T',011,233/
      DATA IFILV(312),AVBL(312),IQ(312),IS(312)     &
     &                      /1,'TOT COL REL HUM     ',052,200/
      DATA IFILV(313),AVBL(313),IQ(313),IS(313)     &
     &                      /1,'CLOUD WORK FUNCTION ',146,200/
      DATA IFILV(314),AVBL(314),IQ(314),IS(314)     &
     &                      /1,'MAX WIND TEMPERATURE',011,006/
      DATA IFILV(315),AVBL(315),IQ(315),IS(315)     &
     &                      /1,'AVE Z GRAVITY STRESS',147,001/
      DATA IFILV(316),AVBL(316),IQ(316),IS(316)     &
     &                      /1,'AVE M GRAVITY STRESS',148,001/
      DATA IFILV(317),AVBL(317),IQ(317),IS(317)     &
     &                      /1,'AVE PRECIP TYPE     ',140,001/
      DATA IFILV(318),AVBL(318),IQ(318),IS(318)     &
     &                      /1,'LFM 0.44-1.00 RELHUM',052,108/
      DATA IFILV(319),AVBL(319),IQ(319),IS(319)     &
     &                      /1,'LFM 0.72-0.94 RELHUM',052,108/
      DATA IFILV(320),AVBL(320),IQ(320),IS(320)     &
     &                      /1,'LFM 0.44-0.72 RELHUM',052,108/
      DATA IFILV(321),AVBL(321),IQ(321),IS(321)     &
     &                      /1,'NGM 0.9950 TEMP     ',011,107/
      DATA IFILV(322),AVBL(322),IQ(322),IS(322)     &
     &                      /1,'NGM 0.9950 POT TEMP ',013,107/
      DATA IFILV(323),AVBL(323),IQ(323),IS(323)     &
     &                      /1,'NGM 0.9950 REL HUM  ',052,107/
      DATA IFILV(324),AVBL(324),IQ(324),IS(324)     &
     &                      /1,'NGM 0.9950 U WIND   ',033,107/
      DATA IFILV(325),AVBL(325),IQ(325),IS(325)     &
     &                      /1,'NGM 0.9950 V WIND   ',034,107/ 
      DATA IFILV(326),AVBL(326),IQ(326),IS(326)     &
     &                      /1,'NGM 0.9950 OMEGA    ',039,107/
      DATA IFILV(327),AVBL(327),IQ(327),IS(327)     &
     &                      /1,'GOES TB - CH 2      ',213,008/ !table 129
      DATA IFILV(328),AVBL(328),IQ(328),IS(328)     &
     &                      /1,'GOES TB - CH 3      ',214,008/ !table 129
      DATA IFILV(329),AVBL(329),IQ(329),IS(329)     &
     &                      /1,'GOES TB - CH 4      ',215,008/ !table 129
      DATA IFILV(330),AVBL(330),IQ(330),IS(330)     &
     &                      /1,'GOES TB - CH 5      ',216,008/ !table 129     
      DATA IFILV(331),AVBL(331),IQ(331),IS(331)     &
     &                      /1,'CLD FRAC ON P SFCS  ',071,100/
      DATA IFILV(332),AVBL(332),IQ(332),IS(332)     &
     &                      /1,'U WIND ON THETA SFCS',033,113/
      DATA IFILV(333),AVBL(333),IQ(333),IS(333)     &
     &                      /1,'V WIND ON THETA SFCS',034,113/
      DATA IFILV(334),AVBL(334),IQ(334),IS(334)     &
     &                      /1,'TEMP ON THETA SFCS  ',011,113/
      DATA IFILV(335),AVBL(335),IQ(335),IS(335)     &
     &                      /1,'PV ON THETA SFCS    ',004,113/
      DATA IFILV(353),AVBL(353),IQ(353),IS(353)     &
     &                      /1,'M STRMFUNC ON THETA ',037,113/
      DATA IFILV(351),AVBL(351),IQ(351),IS(351)     &
     &                      /1,'S STAB ON THETA SFCS',019,113/
      DATA IFILV(352),AVBL(352),IQ(352),IS(352)     &
     &                      /1,'RH ON THETA SFCS    ',052,113/
      DATA IFILV(336),AVBL(336),IQ(336),IS(336)     &
     &                      /1,'U WIND ON PV SFCS   ',033,117/
      DATA IFILV(337),AVBL(337),IQ(337),IS(337)     &
     &                      /1,'V WIND ON PV SFCS   ',034,117/
      DATA IFILV(338),AVBL(338),IQ(338),IS(338)     &
     &                      /1,'TEMP ON PV SFCS     ',011,117/
      DATA IFILV(339),AVBL(339),IQ(339),IS(339)     &
     &                      /1,'HEIGHT ON PV SFCS   ',007,117/
      DATA IFILV(340),AVBL(340),IQ(340),IS(340)     &
     &                      /1,'PRESSURE ON PV SFCS ',001,117/
      DATA IFILV(341),AVBL(341),IQ(341),IS(341)     &
     &                      /1,'SHEAR ON PV SFCS    ',136,117/
      DATA IFILV(342),AVBL(342),IQ(342),IS(342)     &
     &                      /1,'PBL CLD FRACTION    ',071,211/
      DATA IFILV(343),AVBL(343),IQ(343),IS(343)     &
     &                      /1,'AVE WATER RUNOFF    ',090,001/
      DATA IFILV(344),AVBL(344),IQ(344),IS(344)     &
     &                      /1,'PBL REGIME          ',220,001/
      DATA IFILV(345),AVBL(345),IQ(345),IS(345)     &
     &                      /1,'MAX SHELTER TEMP    ',015,105/
      DATA IFILV(346),AVBL(346),IQ(346),IS(346)     &
     &                      /1,'MIN SHELTER TEMP    ',016,105/
      DATA IFILV(347),AVBL(347),IQ(347),IS(347)     &
     &                      /1,'MAX SHELTER RH      ',218,105/ !table129
      DATA IFILV(348),AVBL(348),IQ(348),IS(348)     &
     &                      /1,'MIN SHELTER RH      ',217,105/ !table129
      DATA IFILV(349),AVBL(349),IQ(349),IS(349)     &
     &                      /1,'ICE THICKNESS       ',092,001/
      DATA IFILV(354),AVBL(354),IQ(354),IS(354)     &
     &                      /1,'SW TNDY ON P SFCS   ',250,100/
      DATA IFILV(355),AVBL(355),IQ(355),IS(355)     &
     &                      /1,'LW TNDY ON P SFCS   ',251,100/
      DATA IFILV(356),AVBL(356),IQ(356),IS(356)     &
     &                      /1,'VDIFF TNDY ON P SFCS',246,100/
      DATA IFILV(357),AVBL(357),IQ(357),IS(357)     &
     &                      /1,'D CNVCT TNDY ON P SF',242,100/
      DATA IFILV(358),AVBL(358),IQ(358),IS(358)     &
     &                      /1,'S CNVCT TNDY ON P SF',244,100/
      DATA IFILV(359),AVBL(359),IQ(359),IS(359)     &
     &                      /1,'GRDSCL TNDY ON P SFC',241,100/
      DATA IFILV(360),AVBL(360),IQ(360),IS(360)     &
     &                      /1,'VDIFF MOIS ON P SFCS',249,100/
      DATA IFILV(361),AVBL(361),IQ(361),IS(361)     &
     &                      /1,'D CNVCT MOIS ON P SF',243,100/
      DATA IFILV(362),AVBL(362),IQ(362),IS(362)     &
     &                      /1,'S CNVCT MOIS ON P SF',245,100/
      DATA IFILV(363),AVBL(363),IQ(363),IS(363)     &
     &                      /1,'N RAD TNDY ON P SFCS',173,100/
      DATA IFILV(364),AVBL(364),IQ(364),IS(364)     &
     &                      /1,'OZONE VDIFF ON P SFC',174,100/
      DATA IFILV(365),AVBL(365),IQ(365),IS(365)     &
     &                      /1,'OZONE PROD ON P SFCS',175,100/
      DATA IFILV(366),AVBL(366),IQ(366),IS(366)     &
     &                      /1,'OZONE TNDY ON P SFCS',188,100/
      DATA IFILV(367),AVBL(367),IQ(367),IS(367)     &
     &                      /1,'MASS WEIGHTED PV    ',139,100/
      DATA IFILV(368),AVBL(368),IQ(368),IS(368)     &
     &                      /1,'UNKNOWN D3D ON P SFC',239,100/
      DATA IFILV(369),AVBL(369),IQ(369),IS(369)     &
     &                      /1,'VDIFF Z ACCE ON P SF',247,100/
      DATA IFILV(370),AVBL(370),IQ(370),IS(370)     &
     &                      /1,'G DRAG Z ACCE ON P S',181,100/
      DATA IFILV(371),AVBL(371),IQ(371),IS(371)     &
     &                      /1,'CNVCT U M MIX ON P S',183,100/
      DATA IFILV(372),AVBL(372),IQ(372),IS(372)     &
     &                      /1,'VDIFF M ACCE ON P SF',248,100/
      DATA IFILV(373),AVBL(373),IQ(373),IS(373)     &
     &                      /1,'G DRAG M ACCE ON P S',182,100/
      DATA IFILV(374),AVBL(374),IQ(374),IS(374)     &
     &                      /1,'CNVCT V M MIX ON P S',184,100/
      DATA IFILV(375),AVBL(375),IQ(375),IS(375)     &
     &                      /1,'N CNVCT CLD FRA ON P',213,100/
      DATA IFILV(376),AVBL(376),IQ(376),IS(376)     &
     &                      /1,'GOES BRIGHTNESS-CH 3',221,008/ !table 129
      DATA IFILV(377),AVBL(377),IQ(377),IS(377)     &
     &                      /1,'GOES BRIGHTNESS-CH 4',222,008/ !table 129 
      DATA IFILV(378),AVBL(378),IQ(378),IS(378)     &
     &                      /1,'OMEGA ON THETA SFCS ',039,113/
! D3D fields
      DATA IFILV(379),AVBL(379),IQ(379),IS(379)     &
     &                      /1,'T DIAB TNDY ON P SFC',215,100/
      DATA IFILV(391),AVBL(391),IQ(391),IS(391)     &
     &                      /1,'CNVCT U M FLX ON P S',202,100/
      DATA IFILV(392),AVBL(392),IQ(392),IS(392)     &
     &                      /1,'CNVCT D M FLX ON P S',209,100/
      DATA IFILV(393),AVBL(393),IQ(393),IS(393)     &
     &                      /1,'CNVCT DET M FLX ON P',219,100/
      DATA IFILV(394),AVBL(394),IQ(394),IS(394)     &
     &                      /1,'CNVCT Z G DRAG ON P ',196,100/
      DATA IFILV(395),AVBL(395),IQ(395),IS(395)     &
     &                      /1,'CNVCT M G DRAG ON P ',197,100/
!---- Using table 129   !aqm PLee 1/07
      DATA IFILV(380),AVBL(380),IQ(380),IS(380)     &
     &                      /1,'DIFFUSION H RATE MDL',182,109/
!---- Using table 2     !aqm PLee 3/07
      DATA IFILV(381),AVBL(381),IQ(381),IS(381)     &
     &                      /1,'MIXHT HEIGHT        ',067,001/
! NEW GFS FLUX FILE FIELDS
      DATA IFILV(382),AVBL(382),IQ(382),IS(382)     &
     &                      /1,'AVE CLR INC SFC LW  ',163,001/
      DATA IFILV(383),AVBL(383),IQ(383),IS(383)     &
     &                      /1,'AVE CLR INC SFC SW  ',161,001/
      DATA IFILV(384),AVBL(384),IQ(384),IS(384)     &
     &                      /1,'AVE CLR OUT SFC LW  ',162,001/
      DATA IFILV(385),AVBL(385),IQ(385),IS(385)     &
     &                      /1,'AVE CLR OUT TOA LW  ',162,008/
      DATA IFILV(386),AVBL(386),IQ(386),IS(386)     &
     &                      /1,'AVE CLR OUT SFC SW  ',160,001/
      DATA IFILV(387),AVBL(387),IQ(387),IS(387)     &
     &                      /1,'AVE CLR OUT TOA SW  ',160,008/
      DATA IFILV(388),AVBL(388),IQ(388),IS(388)     &
     &                      /1,'AVE INC TOA SW      ',204,008/
      DATA IFILV(389),AVBL(389),IQ(389),IS(389)     &  
     &                      /1,'TRANSPORT U WIND    ',033,220/
      DATA IFILV(390),AVBL(390),IQ(390),IS(390)     &
     &                      /1,'TRANSPORT V WIND    ',034,220/
! Add TIGGE FIELDS
      DATA IFILV(396),AVBL(396),IQ(396),IS(396)     & 
     &                      /1,'SUNSHINE DURATION   ',191,001/ !table 133
      DATA IFILV(397),AVBL(397),IQ(397),IS(397)     &
     &                      /1,'FIELD CAPACITY      ',220,001/ !table 130
! Add ICAO FIELDS
      DATA IFILV(398),AVBL(398),IQ(398),IS(398)     & 
     &                      /1,'ICAO HGHT MAX WIND  ',005,006/ 
      DATA IFILV(399),AVBL(399),IQ(399),IS(399)     &
     &                      /1,'ICAO HGHT AT TROP   ',005,007/
      DATA IFILV(400),AVBL(400),IQ(400),IS(400)     &
     &                      /1,'RADAR ECHO TOP      ',240,200/
! Add MORE CFSRR FIELDS
! surface Visible beam downward solar flux
      DATA IFILV(401),AVBL(401),IQ(401),IS(401) & 
     &                      /1,'AVE IN SFC VIS SW BE',166,001/
!surface Visible diffuse downward solar flux
      DATA IFILV(402),AVBL(402),IQ(402),IS(402) & 
     &                      /1,'AVE IN SFC VIS SW DF',167,001/
!surface Near IR beam downward solar flux
      DATA IFILV(403),AVBL(403),IQ(403),IS(403) & 
     &                      /1,'AVE IN SFC IR SW BE ',168,001/
!surface Near IR diffuse downward solar flux
      DATA IFILV(404),AVBL(404),IQ(404),IS(404) & 
     &                      /1,'AVE IN SFC IR SW DF ',169,001/
! SNOWFALL RATE
      DATA IFILV(405),AVBL(405),IQ(405),IS(405) &
     &                      /1,'AVE SNOWFALL RATE   ',064,001/
! ADD DUST FIELDS
      DATA IFILV(438),AVBL(438),IQ(438),IS(438) &
     &                      /1,'DUST 1 ON P SFCS    ',240,100/
      DATA IFILV(439),AVBL(439),IQ(439),IS(439) &
     &                      /1,'DUST 2 ON P SFCS    ',241,100/
      DATA IFILV(440),AVBL(440),IQ(440),IS(440) &
     &                      /1,'DUST 3 ON P SFCS    ',242,100/
      DATA IFILV(441),AVBL(441),IQ(441),IS(441) &
     &                      /1,'DUST 4 ON P SFCS    ',243,100/
      DATA IFILV(442),AVBL(442),IQ(442),IS(442) &
     &                      /1,'DUST 5 ON P SFCS    ',244,100/
!
      DATA IFILV(443),AVBL(443),IQ(443),IS(443) &
     &                      /1,'EQUIL LEVEL HEIGHT  ',007,247/
      DATA IFILV(444),AVBL(444),IQ(444),IS(444) &
     &                      /1,'LIGHTNING           ',187,001/
! GOES WEST
      DATA IFILV(446),AVBL(446),IQ(446),IS(446)     &
     &                      /1,'GOES W TB - CH 2    ',241,008/ !Table 130
      DATA IFILV(447),AVBL(447),IQ(447),IS(447)     &
     &                      /1,'GOES W TB - CH 3    ',242,008/ !Table 130
      DATA IFILV(448),AVBL(448),IQ(448),IS(448)     &
     &                      /1,'GOES W TB - CH 4    ',243,008/ !Table 130
      DATA IFILV(449),AVBL(449),IQ(449),IS(449)     &
     &                      /1,'GOES W TB - CH 5    ',244,008/ !Table 130
! NCAR GFIP
      DATA IFILV(450),AVBL(450),IQ(450),IS(450)     &
     &                      /1,'NCAR IN-FLIGHT ICING',168,100/
! Flight level Q
      DATA IFILV(451),AVBL(451),IQ(451),IS(451)     &
     &                      /1,'SPE HUM AT FD HEIGHT',051,103/
! Virtual T based CAPE
      DATA IFILV(452),AVBL(452),IQ(452),IS(452)       &
     &                      /1,'TV CNVCT AVBL POT EN',202,001/
      DATA IFILV(453),AVBL(453),IQ(453),IS(453)       &
     &                      /1,'TV CNVCT INHIBITION ',201,001/
      DATA IFILV(454),AVBL(454),IQ(454),IS(454)     &
     &                      /1,'VENTILATION RATE    ',241,220/
      DATA IFILV(455),AVBL(455),IQ(455),IS(455)     &
     &                      /1,'HAINES INDEX        ',250,001/
      DATA IFILV(456),AVBL(456),IQ(456),IS(456)     &
     &                      /1,'GOESE TB-2 NON NADIR',213,008/ !table 129
      DATA IFILV(457),AVBL(457),IQ(457),IS(457)     &
     &                      /1,'GOESE TB-3 NON NADIR',214,008/ !table 129
      DATA IFILV(458),AVBL(458),IQ(458),IS(458)     &
     &                      /1,'GOESE TB-4 NON NADIR',215,008/ !table 129
      DATA IFILV(459),AVBL(459),IQ(459),IS(459)     &
     &                      /1,'GOESE TB-5 NON NADIR',216,008/ !table 129 
      DATA IFILV(460),AVBL(460),IQ(460),IS(460)     &
     &                      /1,'GOESW TB-2 NON NADIR',241,008/ !table 130
      DATA IFILV(461),AVBL(461),IQ(461),IS(461)     &
     &                      /1,'GOESW TB-3 NON NADIR',242,008/ !table 130
      DATA IFILV(462),AVBL(462),IQ(462),IS(462)     &
     &                      /1,'GOESW TB-4 NON NADIR',243,008/ !table 130
      DATA IFILV(463),AVBL(463),IQ(463),IS(463)     &
     &                      /1,'GOESW TB-5 NON NADIR',244,008/ !table 130   
      DATA IFILV(482),AVBL(482),IQ(482),IS(482)     &
     &                      /1,'PRESS AT FD HEIGHTS ',001,103/
      DATA IFILV(483),AVBL(483),IQ(483),IS(483)     &
     &                      /1,'AMSRE TB - CH 9     ',176,008/ !table 133   
      DATA IFILV(484),AVBL(484),IQ(484),IS(484)     &
     &                      /1,'AMSRE TB - CH 10    ',177,008/ !table 133
      DATA IFILV(485),AVBL(485),IQ(485),IS(485)     &
     &                      /1,'AMSRE TB - CH 11    ',178,008/ !table 133
      DATA IFILV(486),AVBL(486),IQ(486),IS(486)     &
     &                      /1,'AMSRE TB - CH 12    ',179,008/ !table 133
!end initialization
!
   end module RQSTFLD_mod
