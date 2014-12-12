CDECK  ID>, PHOCDE. 
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOtos CDE's
C.
C.    Purpose:  Keep definitions  for PHOTOS QED correction Monte Carlo.
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  29/11/89
C.                                                Last Update: 10/01/92
C.
C.----------------------------------------------------------------------
CDECK  ID>, PHOINI. 
      SUBROUTINE PHOINI
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays INItialisation
C.
C.    Purpose:  Initialisation  routine  for  the  PHOTOS  QED radiation
C.              package.  Should be called  at least once  before a call
C.              to the steering program 'PHOTOS' is made.
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
C.                                                Last Update: 12/04/90
C.
C.----------------------------------------------------------------------
      INTEGER INIT
      SAVE INIT
      DATA INIT/ 0/
C--
C--   Return if already initialized...
      IF (INIT.NE.0) RETURN
      INIT=1
C--
C--   Preset parameters in PHOTOS commons
      CALL PHOCIN
C--
C--   Print info
      CALL PHOINF
C--
C--   Initialize Marsaglia and Zaman random number generator
      CALL PHORIN
      RETURN
      END
CDECK  ID>, PHOCIN. 
      SUBROUTINE PHOCIN
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton Common INitialisation
C.
C.    Purpose:  Initialisation of parameters in common blocks.
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  Commons /PHOLUN/, /PHOPHO/, /PHOCOP/, /PHPICO/
C.                                and /PHSEED/.
C.
C.    Author(s):  B. van Eijk                     Created at:  26/11/89
C.                                                Last Update: 10/01/92
C.
C.----------------------------------------------------------------------
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      LOGICAL QEDRAD
      COMMON/PHOQED/QEDRAD(NMXHEP)
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHCUT,XPHMAX,XPHOTO
      COMMON/PHOPHD/COSTHG,SINTHG
      COMMON/PHOPHS/XPHCUT,XPHMAX,XPHOTO
      REAL ALPHA
      COMMON/PHOCOP/ALPHA
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      INTEGER ISEED,I97,J97
      REAL URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      INTEGER INIT,I
      SAVE INIT
      DATA INIT/ 0/
C--
C--   Return if already initialized...
      IF (INIT.NE.0) RETURN
      INIT=1
C--
C--   Preset switch  for  photon emission to 'TRUE' for each particle in
C--   /HEPEVT/, this interface is needed for KORALB and KORALZ...
      DO 10 I=1,NMXHEP
   10 QEDRAD(I)=.TRUE.
C--
C--   Logical output unit for printing of PHOTOS error messages
      PHLUN=6
C--
C--   Set cut parameter for photon radiation
      XPHCUT=0.001
C--
C--   Define some constants
      ALPHA=0.00729735039
      PI=3.14159265358979324
      TWOPI=6.28318530717958648
C--
C--   Default seeds Marsaglia and Zaman random number generator
      ISEED(1)=1802
      ISEED(2)=9373
C--
C--   Initialise status counter for warning messages
      DO 20 I=1,PHOMES
   20 STATUS(I)=0
      RETURN
      END
CDECK  ID>, PHOINF. 
      SUBROUTINE PHOINF
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays general INFo
C.
C.    Purpose:  Print PHOTOS info
C.
C.    Input Parameters:   PHOLUN
C.
C.    Output Parameters:  PHOVN1, PHOVN2
C.
C.    Author(s):  B. van Eijk                     Created at:  12/04/90
C.                                                Last Update: 10/01/92
C.
C.----------------------------------------------------------------------
      INTEGER IV1,IV2,IV3
      INTEGER PHOVN1,PHOVN2
      COMMON/PHOVER/PHOVN1,PHOVN2
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
C--
C--   PHOTOS version number and release date
      PHOVN1=107
      PHOVN2=260693
C--
C--   Print info
      WRITE(PHLUN,9000)
      WRITE(PHLUN,9020)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9030)
      IV1=PHOVN1/100
      IV2=PHOVN1-IV1*100
      WRITE(PHLUN,9040) IV1,IV2
      IV1=PHOVN2/10000
      IV2=(PHOVN2-IV1*10000)/100
      IV3=PHOVN2-IV1*10000-IV2*100
      WRITE(PHLUN,9050) IV1,IV2,IV3
      WRITE(PHLUN,9030)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9060)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9070)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9020)
      RETURN
 9000 FORMAT(1H1)
 9010 FORMAT(1H ,'*',T81,'*')
 9020 FORMAT(1H ,80('*'))
 9030 FORMAT(1H ,'*',26X,26('='),T81,'*')
 9040 FORMAT(1H ,'*',28X,'PHOTOS, Version: ',I2,'.',I2,T81,'*')
 9050 FORMAT(1H ,'*',28X,'Released at:  ',I2,'/',I2,'/',I2,T81,'*')
 9060 FORMAT(1H ,'*',18X,'PHOTOS QED Corrections in Particle Decays',
     &T81,'*')
 9070 FORMAT(1H ,'*',9X,'Monte Carlo Program - by E. Barberio, B. van Ei
     &jk and Z. Was',T81,'*')
      END
CDECK  ID>, PHOTOS. 
      SUBROUTINE PHOTOS(IPARR)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   Photon radiation in decays
C.
C.    Purpose:  Order (alpha) radiative corrections  are  calculated in
C.              the decay of the IPPAR-th particle in the standard  HEP
C.              common /HEPEVT/.  Photon radiation takes place from one
C.              of the charged daughters of the decaying particle IPPAR
C.              Beware that the kinematics  of both parent and daughter
C.              should preferably be defined in the resframe of the pa-
C.              rent.  If not, boosts etc. will  considerably slow down
C.              the program.  Casacade decays  will  be  treated  until
C.              the end of the cascade,  photons may be  added  at each
C.              stage in the decay chain.
C.
C.    Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
C.                                /HEPEVT/ and the common itself,
C.
C.    Output Parameters:  Common  /HEPEVT/, either  with  or  without a
C.                                photon(s) added.
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
C.                                                Last Update: 26/03/93
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION MINMAS,BET(3),GAM,PB,MPASQR,MCHREN
      DOUBLE PRECISION BETA,EPS,DEL1,DEL2,DATA
      REAL PCHARG(5),PHOCHA,PHOSPI,PHORAN,PHOCOR,PHOTON(5),MASSUM
      INTEGER IP,IPARR,IPPAR,I,J,K,L,ME,NCHARG,NEUPOI,NLAST,THEDUM
      INTEGER IDABS,MOTHER,POSPHO,WTDUM
      LOGICAL CASCAD
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL QEDRAD
      COMMON/PHOQED/QEDRAD(NMXHEP)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(NMXHEP)
      INTEGER PDMMAX
      PARAMETER (PDMMAX=100)
      INTEGER CHAPOI
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5),CHAPOI(PDMMAX)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHCUT,XPHMAX,XPHOTO
      COMMON/PHOPHD/COSTHG,SINTHG
      COMMON/PHOPHS/XPHCUT,XPHMAX,XPHOTO
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
C--
      IPPAR=ABS(IPARR)
      CALL PHOCHK(IPPAR)
C--   Store pointers for cascade treatement...
      IP=IPPAR
      NLAST=NHEP
      CASCAD=.FALSE.
C--
C--   Check decay multiplicity..
      IF (JDAHEP(1,IP).EQ.0) RETURN
C--
C--   Loop over daughters, determine charge multiplicity
   10 NCHARG=0
      IREP=0
      MINMAS=0.
      MASSUM=0.
      DO 20 I=JDAHEP(1,IP),JDAHEP(2,IP)
C--
C--   Cascade decay ?
        IF (JDAHEP(1,I).NE.0) CASCAD=.TRUE.
        IF (IPARR.LE.0) CASCAD=.FALSE.
C--
C--   Exclude marked particles, quarks and gluons etc...
        IDABS=ABS(IDHEP(I))
        IF (QEDRAD(I).AND.CHKIF(I)) THEN
          IF (PHOCHA(IDHEP(I)).NE.0) THEN
            NCHARG=NCHARG+1
            IF (NCHARG.GT.PDMMAX) THEN
              DATA=NCHARG
              CALL PHOERR(1,'PHOTOS',DATA)
            ENDIF
            CHAPOI(NCHARG)=I
          ENDIF
          MINMAS=MINMAS+PHEP(5,I)**2
        ENDIF
        MASSUM=MASSUM+PHEP(5,I)
   20 CONTINUE
      IF (NCHARG.NE.0) THEN
C--
C--   Check that sum of daughter masses does not exceed parent mass
        IF ((PHEP(5,IP)-MASSUM)/PHEP(5,IP).GT.2.*XPHCUT) THEN
C--
C--   Order  charged  particles  according  to decreasing mass, this  to
C--   increase efficiency (smallest mass is treated first).
          IF (NCHARG.GT.1) CALL PHOOMA(1,NCHARG,CHAPOI)
C--
C--   Check whether parent is in its rest frame...
   30     IF (ABS(PHEP(4,IP)-PHEP(5,IP))/PHEP(5,IP).GT.1.E-8) THEN
C--
C--   Boost charged daughter to rest frame of parent...
            DO 40 J=1,3
   40       BET(J)=-PHEP(J,IP)/PHEP(5,IP)
            GAM=PHEP(4,IP)/PHEP(5,IP)
            PB=BET(1)*PHEP(1,CHAPOI(NCHARG))+BET(2)*PHEP(2,
     &      CHAPOI(NCHARG))+BET(3)*PHEP(3,CHAPOI(NCHARG))
            DO 50 J=1,3
   50       PCHARG(J)=PHEP(J,CHAPOI(NCHARG))+BET(J)*(PHEP(4,CHAPOI
     &      (NCHARG))+PB/(GAM+1.))
            PCHARG(4)=GAM*PHEP(4,CHAPOI(NCHARG))+PB
C--
C--   ...and reconstruct boosted resultant neutral system
            DO 60 J=1,3
   60       PNEUTR(J)=-PCHARG(J)
            PNEUTR(4)=PHEP(5,IP)-PCHARG(4)
          ELSE
            DO 70 J=1,3
   70       PNEUTR(J)=-PHEP(J,CHAPOI(NCHARG))
            PNEUTR(4)=PHEP(5,IP)-PHEP(4,CHAPOI(NCHARG))
          ENDIF
C--
C--   Calculate  invariant  mass of 'neutral' etc. systems
          MPASQR=PHEP(5,IP)**2
          MCHSQR=PHEP(5,CHAPOI(NCHARG))**2
          IF ((JDAHEP(2,IP)-JDAHEP(1,IP)).EQ.1) THEN
            NEUPOI=JDAHEP(1,IP)
            IF (NEUPOI.EQ.CHAPOI(NCHARG)) NEUPOI=JDAHEP(2,IP)
            MNESQR=PHEP(5,NEUPOI)**2
            PNEUTR(5)=PHEP(5,NEUPOI)
          ELSE
            MNESQR=PNEUTR(4)**2-PNEUTR(1)**2-PNEUTR(2)**2-PNEUTR(3)**2
            MNESQR=MAX(MNESQR,MINMAS-MCHSQR)
            PNEUTR(5)=SQRT(MNESQR)
          ENDIF
C--
C--   Determine kinematical limit...
          XPHMAX=(MPASQR-(PNEUTR(5)+PHEP(5,CHAPOI(NCHARG)))**2)/MPASQR
C--
C--   Photon energy fraction...
          CALL PHOENE(MPASQR,MCHREN,BETA,IDHEP(CHAPOI(NCHARG)))
C--
C--   Energy fraction not too large (very seldom) ? Define angle.
          IF ((XPHOTO.LT.XPHCUT).OR.(XPHOTO.GT.XPHMAX)) THEN
C--
C--   No radiation was accepted, check  for more daughters  that may ra-
C--   diate and correct radiation probability...
            NCHARG=NCHARG-1
            IF (NCHARG.GT.0) THEN
              IREP=IREP+1
              GOTO 30
            ENDIF
          ELSE
C--
C--   Angle is generated  in  the  frame defined  by  charged vector and
C--   PNEUTR, distribution is taken in the infrared limit...
            EPS=MCHREN/(1.+BETA)
C--
C--   Calculate sin(theta) and cos(theta) from interval variables
            DEL1=(2.-EPS)*(EPS/(2.-EPS))**PHORAN(THEDUM)
            DEL2=2.-DEL1
            COSTHG=(1.-DEL1)/BETA
            SINTHG=SQRT(DEL1*DEL2-MCHREN)/BETA
C--
C--   Determine spin of  particle and construct code  for matrix element
            ME=2.*PHOSPI(IDHEP(CHAPOI(NCHARG)))+1.
C--
C--   Weighting procedure with 'exact' matrix element, reconstruct kine-
C--   matics for photon, neutral and charged system and update /HEPEVT/.
            IF (PHORAN(WTDUM).LE.PHOCOR(MPASQR,MCHREN,ME))
     &      CALL PHOKIN(IP,NCHARG)
          ENDIF
        ELSE
          DATA=PHEP(5,IP)-MASSUM
          CALL PHOERR(10,'PHOTOS',DATA)
        ENDIF
      ENDIF
C--
C--   Check for cascade decays...
      IF (CASCAD) THEN
        DO 80 I=IP+1,NLAST
          IF (JDAHEP(1,I).NE.0) THEN
            IP=I
            GOTO 10
          ENDIF
   80   CONTINUE
      ENDIF
C--
C--   rearrange  /HEPEVT/  to get correct order...
        IF (NHEP.GT.NLAST) THEN
          DO 160 I=NLAST+1,NHEP
C--
C--   Photon mother and position...
            MOTHER=JMOHEP(1,I)
            POSPHO=JDAHEP(2,MOTHER)+1
C--
C--   Exclude photon in sequence !
            IF (POSPHO.NE.NHEP) THEN
C--
C--   Intermediate save of photon energy/momentum
              DO 90 J=1,5
   90         PHOTON(J)=PHEP(J,I)
C--
C--   Order /HEPEVT/
              DO 120 K=I,POSPHO+1,-1
                ISTHEP(K)=ISTHEP(K-1)
                QEDRAD(K)=QEDRAD(K-1)
                IDHEP(K)=IDHEP(K-1)
                DO 100 L=1,2
                  JMOHEP(L,K)=JMOHEP(L,K-1)
  100           JDAHEP(L,K)=JDAHEP(L,K-1)
                DO 110 L=1,5
  110           PHEP(L,K)=PHEP(L,K-1)
                DO 120 L=1,4
  120         VHEP(L,K)=VHEP(L,K-1)
C--
C--   Correct pointers assuming most dirty /HEPEVT/...
              DO 130 K=1,NHEP
                DO 130 L=1,2
                  IF ((JMOHEP(L,K).NE.0).AND.(JMOHEP(L,K).GE.
     &            POSPHO)) JMOHEP(L,K)=JMOHEP(L,K)+1
                  IF ((JDAHEP(L,K).NE.0).AND.(JDAHEP(L,K).GE.
     &            POSPHO)) JDAHEP(L,K)=JDAHEP(L,K)+1
  130         CONTINUE
C--
C--   Store photon energy/momentum
              DO 140 J=1,5
  140         PHEP(J,POSPHO)=PHOTON(J)
            ENDIF
C--
C--   Store pointers for the photon...
            JDAHEP(2,MOTHER)=POSPHO
            ISTHEP(POSPHO)=1
            IDHEP(POSPHO)=22
            JMOHEP(1,POSPHO)=MOTHER
            JMOHEP(2,POSPHO)=0
            JDAHEP(1,POSPHO)=0
            JDAHEP(2,POSPHO)=0
C--
C--   Get photon production vertex position
            DO 150 J=1,4
  150       VHEP(J,POSPHO)=VHEP(J,POSPHO-1)
  160     CONTINUE
        ENDIF
      RETURN
      END
CDECK  ID>, PHZODE. 
      LOGICAL FUNCTION PHZODE(IPPAR,I)
C.----------------------------------------------------------------------
C.
C.    PHZODE:   checkin family relations
C.
C.    Purpose:  checks whether particle I originates from particle
C.              IPPAR
C.
C.    Input Parameters:  IPPAR  :  position of particle starting cascade
C.                       I      :  position to be checked
C.
C.    Output Parameter:  .true. :  particles are related
C.                       .false.: particles are not related
C.
C.    Author(s):  Z. Was                          Created at:  12/10/92
C.                                                Last Update:
C.
C.----------------------------------------------------------------------
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER IPPAR,I,POTOM
      PHZODE=.FALSE.
      POTOM=I
 10   CONTINUE
        POTOM=JMOHEP(1,POTOM)
        PHZODE=POTOM.EQ.IPPAR
      IF (POTOM.GT.IPPAR) GOTO 10
      END
CDECK  ID>, PHOCHK. 
      SUBROUTINE PHOCHK(IPPAR)
C.----------------------------------------------------------------------
C.
C.    PHOCHK:   checkin family correctness
C.
C.    Purpose:  checks whether particles in the common block /HEPEVT/
C.              can be served by PHOTOS
C.
C.    Author(s):  Z. Was                           Created at: 22/10/92
C.                                                Last Update: 27/10/92
C.
C.----------------------------------------------------------------------
C     ********************
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(NMXHEP)
  
      LOGICAL F,PHZODE
      INTEGER IDABS,NLAST,I,IP,IJ,IPPAR
C these are OK .... if you do not like somebody else, add here.
      F(IDABS)=
     &     ( ((IDABS.GT.9).AND.(IDABS.LE.40)) .OR. (IDABS.GT.100) )
     & .AND.(IDABS.NE.21)
     $ .AND.(IDABS.NE.2101).AND.(IDABS.NE.3101).AND.(IDABS.NE.3201)
     & .AND.(IDABS.NE.1103).AND.(IDABS.NE.2103).AND.(IDABS.NE.2203)
     & .AND.(IDABS.NE.3103).AND.(IDABS.NE.3203).AND.(IDABS.NE.3303)
C
      NLAST = NHEP
C
C checking for good particles
      DO 10 I=IPPAR,NLAST
      IDABS    = ABS(IDHEP(I))
      CHKIF(I)= F(IDABS).AND.PHZODE(IPPAR,I)
 10   CONTINUE
c  checking daughters of wrong mothers
      DO 20 IP=IPPAR,NLAST
       IDABS=ABS(IDHEP(IP))
       IF (.NOT.F(IDABS)) THEN
         DO 30 IJ=JDAHEP(1,IP),JDAHEP(2,IP)
 30      CHKIF(IJ)=.FALSE.
        ENDIF
 20   CONTINUE
      END
CDECK  ID>, PHOENE. 
      SUBROUTINE PHOENE(MPASQR,MCHREN,BETA,IDENT)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays calculation  of photon ENErgy
C.              fraction
C.
C.    Purpose:  Subroutine  returns  photon  energy fraction (in (parent
C.              mass)/2 units) for the decay bremsstrahlung.
C.
C.    Input Parameters:  MPASQR:  Mass of decaying system squared,
C.                       XPHCUT:  Minimum energy fraction of photon,
C.                       XPHMAX:  Maximum energy fraction of photon.
C.
C.    Output Parameter:  MCHREN:  Renormalised mass squared,
C.                       BETA:    Beta factor due to renormalisation,
C.                       XPHOTO:  Photon energy fraction,
C.                       XF:      Correction factor for PHOFAC.
C.
C.    Author(s):  S. Jadach, Z. Was               Created at:  01/01/89
C.                B. van Eijk                     Last Update: 26/03/93
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION MPASQR,MCHREN,BIGLOG,BETA,DATA
      INTEGER IWT1,IRN,IWT2
      REAL PRSOFT,PRHARD,PHORAN,PHOFAC
      INTEGER PDMMAX
      PARAMETER (PDMMAX=100)
      INTEGER CHAPOI
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      INTEGER IDENT
      REAL PHOCHA
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5),CHAPOI(PDMMAX)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHCUT,XPHMAX,XPHOTO
      COMMON/PHOPHD/COSTHG,SINTHG
      COMMON/PHOPHS/XPHCUT,XPHMAX,XPHOTO
      REAL ALPHA
      COMMON/PHOCOP/ALPHA
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
C--
      IF (XPHMAX.LE.XPHCUT) THEN
        XPHOTO=0.0
        RETURN
      ENDIF
C--   Probabilities for hard and soft bremstrahlung...
      MCHREN=4.*MCHSQR/MPASQR/(1.+MCHSQR/MPASQR)**2
      BETA=SQRT(1.-MCHREN)
      BIGLOG=LOG(MPASQR/MCHSQR*(1.+BETA)**2/4.*(1.+MCHSQR/MPASQR)**2)
      PRHARD=ALPHA/PI/BETA*BIGLOG*(LOG(XPHMAX/XPHCUT)-.75+XPHCUT/
     &XPHMAX-.25*XPHCUT**2/XPHMAX**2)
      PRHARD=PRHARD*PHOCHA(IDENT)**2
      IF (IREP.EQ.0) PROBH=0.
      PRHARD=PRHARD*PHOFAC(0)
      PROBH=PRHARD
      PRSOFT=1.-PRHARD
C--
C--   Check on kinematical bounds
      IF (PRSOFT.LT.0.1) THEN
        DATA=PRSOFT
        CALL PHOERR(2,'PHOENE',DATA)
      ENDIF
      IF (PHORAN(IWT1).LT.PRSOFT) THEN
C--
C--   No photon... (ie. photon too soft)
        XPHOTO=0.
      ELSE
C--
C--   Hard  photon... (ie.  photon  hard enough).
C--   Calculate  Altarelli-Parisi Kernel
   10   XPHOTO=EXP(PHORAN(IRN)*LOG(XPHCUT/XPHMAX))
        XPHOTO=XPHOTO*XPHMAX
        IF (PHORAN(IWT2).GT.((1.+(1.-XPHOTO/XPHMAX)**2)/2.)) GOTO 10
      ENDIF
C--
C--   Calculate parameter for PHOFAC function
      XF=4.*MCHSQR*MPASQR/(MPASQR+MCHSQR-MNESQR)**2
      RETURN
      END
CDECK  ID>, PHOCOR. 
      FUNCTION PHOCOR(MPASQR,MCHREN,ME)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays CORrection weight from
C.              matrix elements
C.
C.    Purpose:  Calculate  photon  angle.  The reshaping functions  will
C.              have  to  depend  on the spin S of the charged particle.
C.              We define:  ME = 2 * S + 1 !
C.
C.    Input Parameters:  MPASQR:  Parent mass squared,
C.                       MCHREN:  Renormalised mass of charged system,
C.                       ME:      2 * spin + 1 determines matrix element
C.
C.    Output Parameter:  Function value.
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
C.                                                Last Update: 21/03/93
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION MPASQR,MCHREN,BETA,XX,YY,DATA
      INTEGER ME
      REAL PHOCOR,PHOFAC,WT1,WT2,WT3
      INTEGER PDMMAX
      PARAMETER (PDMMAX=100)
      INTEGER CHAPOI
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5),CHAPOI(PDMMAX)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHCUT,XPHMAX,XPHOTO
      COMMON/PHOPHD/COSTHG,SINTHG
      COMMON/PHOPHS/XPHCUT,XPHMAX,XPHOTO
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
C--
C--   Shaping (modified by ZW)...
      XX=4.*MCHSQR/MPASQR*(1.-XPHOTO)/(1.-XPHOTO+(MCHSQR-MNESQR)/
     &MPASQR)**2
      IF (ME.EQ.1) THEN
        YY=1.
        WT3=(1.-XPHOTO/XPHMAX)/((1.+(1.-XPHOTO/XPHMAX)**2)/2.)
      ELSEIF (ME.EQ.2) THEN
        YY=0.5*(1.-XPHOTO/XPHMAX+1./(1.-XPHOTO/XPHMAX))
        WT3=1.
      ELSEIF ((ME.EQ.3).OR.(ME.EQ.4).OR.(ME.EQ.5)) THEN
        YY=1.
        WT3=(1.+(1.-XPHOTO/XPHMAX)**2-(XPHOTO/XPHMAX)**3)/(1.+(1.
     &  -XPHOTO/XPHMAX)** 2)
      ELSE
        DATA=(ME-1.)/2.
        CALL PHOERR(6,'PHOCOR',DATA)
        YY=1.
        WT3=1.
      ENDIF
      BETA=SQRT(1.-XX)
      WT1=(1.-COSTHG*SQRT(1.-MCHREN))/(1.-COSTHG*BETA)
      WT2=(1.-XX/YY/(1.-BETA**2*COSTHG**2))*(1.+COSTHG*BETA)/2.
      WT2=WT2*PHOFAC(1)
      PHOCOR=WT1*WT2*WT3
      CORWT=PHOCOR
      IF (PHOCOR.GT.1.) THEN
        DATA=PHOCOR
        CALL PHOERR(3,'PHOCOR',DATA)
      ENDIF
      RETURN
      END
CDECK  ID>, PHOFAC. 
      FUNCTION PHOFAC(MODE)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays control FACtor
C.
C.    Purpose:  This is the control function for the photon spectrum and
C.              final weighting.  It is  called  from PHOENE for genera-
C.              ting the raw photon energy spectrum (MODE=0) and in PHO-
C.              COR to scale the final weight (MODE=1).  The factor con-
C.              sists of 3 terms.  Addition of  the factor FF which mul-
C.              tiplies PHOFAC for MODE=0 and divides PHOFAC for MODE=1,
C.              does not affect  the results for  the MC generation.  An
C.              appropriate choice  for FF can speed up the calculation.
C.              Note that a too small value of FF may cause weight over-
C.              flow in PHOCOR  and will generate a warning, halting the
C.              execution.  PRX  should  be  included for repeated calls
C.              for  the  same event, allowing more particles to radiate
C.              photons.  At  the  first  call IREP=0, for  more  than 1
C.              charged  decay  products, IREP >= 1.  Thus,  PRSOFT  (no
C.              photon radiation  probability  in  the  previous  calls)
C.              appropriately scales the strength of the bremsstrahlung.
C.
C.    Input Parameters:  MODE, PROBH, XF
C.
C.    Output Parameter:  Function value
C.
C.    Author(s):  S. Jadach, Z. Was               Created at:  01/01/89
C.                B. van Eijk                     Last Update: 13/02/90
C.
C.----------------------------------------------------------------------
      REAL PHOFAC,FF,PRX
      INTEGER MODE
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
      SAVE PRX,FF
      DATA PRX,FF/ 0., 0./
      IF (MODE.EQ.0) THEN
        IF (IREP.EQ.0) PRX=1.
        PRX=PRX/(1.-PROBH)
        FF=1.
C--
C--   Following options are not considered for the time being...
C--   (1) Good choice, but does not save very much time:
C--       FF=(1.0-SQRT(XF)/2.0)/(1.0+SQRT(XF)/2.0)
C--   (2) Taken from the blue, but works without weight overflows...
C--       FF=(1.-XF/(1-(1-SQRT(XF))**2))*(1+(1-SQRT(XF))/SQRT(1-XF))/2
        PHOFAC=FF*PRX
      ELSE
        PHOFAC=1./FF
      ENDIF
      END
CDECK  ID>, PHOKIN. 
      SUBROUTINE PHOKIN(IP,NCHARG)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in  decays reconstruction of KINematics
C.
C.    Purpose:  Starting  from   the  charged  particle energy/momentum,
C.              PNEUTR, photon  energy  fraction and photon  angle  with
C.              respect  to  the axis formed by charged particle energy/
C.              momentum  vector  and PNEUTR, scale the energy/momentum,
C.              keeping the original direction of the neutral system  in
C.              the lab. frame untouched.
C.
C.    Input Parameters:   IP:      Pointer  to   decaying  particle   in
C.                                 /HEPEVT/  and   the   common   itself
C.                        NCHARG:  Element  of the array CHAPOI  contai-
C.                                 ning pointer to the charged radiating
C.                                 daughter in /HEPEVT/.
C.
C.    Output Parameters:  Common /HEPEVT/, with photon added.
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
C.                                                Last Update: 21/12/92
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION PHOAN1,PHOAN2,ANGLE,FI1,FI3,FI4,FI5,TH1,TH3,TH4
      DOUBLE PRECISION BET(3),PB,GAM,PARNE,QNEW,QOLD,DATA
      INTEGER IP,FI3DUM,NCHARG,I,II,J,K,NEUDAU,FIRST,LAST
      REAL EPHOTO,PMAVIR,PHOTRI
      REAL GNEUT,PHORAN,CCOSTH,SSINTH,PVEC(4),PBOOST(5)
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER PDMMAX
      PARAMETER (PDMMAX=100)
      INTEGER CHAPOI
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5),CHAPOI(PDMMAX)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHCUT,XPHMAX,XPHOTO
      COMMON/PHOPHD/COSTHG,SINTHG
      COMMON/PHOPHS/XPHCUT,XPHMAX,XPHOTO
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      LOGICAL BOOST
      BOOST=.FALSE.
C--
C--   Check whether parent is in its rest frame...
      IF (ABS(PHEP(4,IP)-PHEP(5,IP))/PHEP(5,IP).GT.1.E-8) THEN
        BOOST=.TRUE.
C--
C--   Boost daughter particles to rest frame of parent...
C--   Resultant neutral system already calculated in rest frame !
        DO 10 J=1,3
   10   BET(J)=-PHEP(J,IP)/PHEP(5,IP)
        GAM=PHEP(4,IP)/PHEP(5,IP)
        DO 30 I=JDAHEP(1,IP),JDAHEP(2,IP)
          PB=BET(1)*PHEP(1,I)+BET(2)*PHEP(2,I)+BET(3)*PHEP(3,I)
          DO 20 J=1,3
   20     PHEP(J,I)=PHEP(J,I)+BET(J)*(PHEP(4,I)+PB/(GAM+1.))
   30   PHEP(4,I)=GAM*PHEP(4,I)+PB
      ENDIF
      EPHOTO=XPHOTO*PHEP(5,IP)/2.
      PMAVIR=SQRT(PHEP(5,IP)*(PHEP(5,IP)-2.*EPHOTO))
C--
C--   Find pointer to the first component of 'neutral' system
      DO 40 I=JDAHEP(1,IP),JDAHEP(2,IP)
        IF (I.NE.CHAPOI(NCHARG)) THEN
          NEUDAU=I
          GOTO 50
        ENDIF
   40 CONTINUE
C--
C--   Pointer not found...
      DATA=NCHARG
      CALL PHOERR(5,'PHOKIN',DATA)
C--
C--   Reconstruct  kinematics  of  charged particle  and  neutral system
   50 FI1=PHOAN1(PNEUTR(1),PNEUTR(2))
C--
C--   Choose axis along  z of  PNEUTR, calculate  angle  between x and y
C--   components  and z  and x-y plane and  perform Lorentz transform...
      TH1=PHOAN2(PNEUTR(3),SQRT(PNEUTR(1)**2+PNEUTR(2)**2))
      CALL PHORO3(-FI1,PNEUTR(1))
      CALL PHORO2(-TH1,PNEUTR(1))
C--
C--   Take  away  photon energy from charged particle and PNEUTR !  Thus
C--   the onshell charged particle  decays into virtual charged particle
C--   and photon.  The virtual charged  particle mass becomes:
C--   SQRT(PHEP(5,IP)*(PHEP(5,IP)-2*EPHOTO)).  Construct  new PNEUTR mo-
C--   mentum in the rest frame of the parent:
C--   1) Scaling parameters...
      QNEW=PHOTRI(PMAVIR,PNEUTR(5),PHEP(5,CHAPOI(NCHARG)))
      QOLD=PNEUTR(3)
      GNEUT=(QNEW**2+QOLD**2+MNESQR)/(QNEW*QOLD+SQRT((QNEW**2+MNESQR)*
     &(QOLD**2+MNESQR)))
      IF (GNEUT.LT.1.) THEN
        DATA=0.
        CALL PHOERR(4,'PHOKIN',DATA)
      ENDIF
      PARNE=GNEUT-SQRT(MAX(GNEUT**2-1.0,0.))
C--
C--   2) ...reductive boost...
      CALL PHOBO3(PARNE,PNEUTR)
C--
C--   ...calculate photon energy in the reduced system...
      NHEP=NHEP+1
      PHEP(4,NHEP)=EPHOTO*PHEP(5,IP)/PMAVIR
C--
C--   ...and photon momenta
      CCOSTH=-COSTHG
      SSINTH=SINTHG
      TH3=PHOAN2(CCOSTH,SSINTH)
      FI3=TWOPI*PHORAN(FI3DUM)
      PHEP(1,NHEP)=PHEP(4,NHEP)*SINTHG*COS(FI3)
      PHEP(2,NHEP)=PHEP(4,NHEP)*SINTHG*SIN(FI3)
C--
C--   Minus sign because axis opposite direction of charged particle !
      PHEP(3,NHEP)=-PHEP(4,NHEP)*COSTHG
      PHEP(5,NHEP)=0.
C--
C--   Rotate in order to get photon along z-axis
      CALL PHORO3(-FI3,PNEUTR(1))
      CALL PHORO3(-FI3,PHEP(1,NHEP))
      CALL PHORO2(-TH3,PNEUTR(1))
      CALL PHORO2(-TH3,PHEP(1,NHEP))
      ANGLE=EPHOTO/PHEP(4,NHEP)
C--
C--   Boost to the rest frame of decaying particle
      CALL PHOBO3(ANGLE,PNEUTR(1))
      CALL PHOBO3(ANGLE,PHEP(1,NHEP))
C--
C--   Back in the parent rest frame but PNEUTR not yet oriented !
      FI4=PHOAN1(PNEUTR(1),PNEUTR(2))
      TH4=PHOAN2(PNEUTR(3),SQRT(PNEUTR(1)**2+PNEUTR(2)**2))
      CALL PHORO3(FI4,PNEUTR(1))
      CALL PHORO3(FI4,PHEP(1,NHEP))
C--
        DO 60 I=2,4
   60   PVEC(I)=0.
        PVEC(1)=1.
        CALL PHORO3(-FI3,PVEC)
        CALL PHORO2(-TH3,PVEC)
        CALL PHOBO3(ANGLE,PVEC)
        CALL PHORO3(FI4,PVEC)
        CALL PHORO2(-TH4,PNEUTR)
        CALL PHORO2(-TH4,PHEP(1,NHEP))
        CALL PHORO2(-TH4,PVEC)
        FI5=PHOAN1(PVEC(1),PVEC(2))
C--
C--   Charged particle restores original direction
        CALL PHORO3(-FI5,PNEUTR)
        CALL PHORO3(-FI5,PHEP(1,NHEP))
        CALL PHORO2(TH1,PNEUTR(1))
        CALL PHORO2(TH1,PHEP(1,NHEP))
        CALL PHORO3(FI1,PNEUTR)
        CALL PHORO3(FI1,PHEP(1,NHEP))
C--   See whether neutral system has multiplicity larger than 1...
      IF ((JDAHEP(2,IP)-JDAHEP(1,IP)).GT.1) THEN
C--   Find pointers to components of 'neutral' system
C--
        FIRST=NEUDAU
        LAST=JDAHEP(2,IP)
        DO 70 I=FIRST,LAST
          IF (I.NE.CHAPOI(NCHARG).AND.(JMOHEP(1,I).EQ.IP)) THEN
C--
C--   Reconstruct kinematics...
            CALL PHORO3(-FI1,PHEP(1,I))
            CALL PHORO2(-TH1,PHEP(1,I))
C--
C--   ...reductive boost
            CALL PHOBO3(PARNE,PHEP(1,I))
C--
C--   Rotate in order to get photon along z-axis
            CALL PHORO3(-FI3,PHEP(1,I))
            CALL PHORO2(-TH3,PHEP(1,I))
C--
C--   Boost to the rest frame of decaying particle
            CALL PHOBO3(ANGLE,PHEP(1,I))
C--
C--   Back in the parent rest-frame but PNEUTR not yet oriented.
            CALL PHORO3(FI4,PHEP(1,I))
            CALL PHORO2(-TH4,PHEP(1,I))
C--
C--   Charged particle restores original direction
            CALL PHORO3(-FI5,PHEP(1,I))
            CALL PHORO2(TH1,PHEP(1,I))
            CALL PHORO3(FI1,PHEP(1,I))
          ENDIF
   70   CONTINUE
      ELSE
C--
C--   ...only one 'neutral' particle in addition to photon!
C        CALL PHORO2(TH1-TH4,PNEUTR(1))
C        CALL PHORO2(TH1-TH4,PHEP(1,NHEP))
C        CALL PHORO3(FI1,PNEUTR(1))
C        CALL PHORO3(FI1,PHEP(1,NHEP))
        DO 80 J=1,4
   80   PHEP(J,NEUDAU)=PNEUTR(J)
      ENDIF
C--
C--   All 'neutrals' treated, fill /HEPEVT/ for charged particle...
      DO 90 J=1,3
   90 PHEP(J,CHAPOI(NCHARG))=-(PHEP(J,NHEP)+PNEUTR(J))
      PHEP(4,CHAPOI(NCHARG))=PHEP(5,IP)-(PHEP(4,NHEP)+PNEUTR(4))
C--
C--   When parent was not in its rest-frame, boost back...
      IF (BOOST) THEN
        DO 110 J=JDAHEP(1,IP),JDAHEP(2,IP)
          PB=-BET(1)*PHEP(1,J)-BET(2)*PHEP(2,J)-BET(3)*PHEP(3,J)
          DO 100 K=1,3
  100     PHEP(K,J)=PHEP(K,J)-BET(K)*(PHEP(4,J)+PB/(GAM+1.))
  110   PHEP(4,J)=GAM*PHEP(4,J)+PB
C--
C--   ...boost photon
        PB=-BET(1)*PHEP(1,NHEP)-BET(2)*PHEP(2,NHEP)-BET(3)*PHEP(3,NHEP)
        DO 120 K=1,3
  120   PHEP(K,NHEP)=PHEP(K,NHEP)-BET(K)*(PHEP(4,NHEP)+PB/(GAM+1.))
        PHEP(4,NHEP)=GAM*PHEP(4,NHEP)+PB
        BOOST=.FALSE.
      ENDIF
C--
C--   Photon mother and daughter pointers !
      JMOHEP(1,NHEP)=IP
      JMOHEP(2,NHEP)=0
      JDAHEP(1,NHEP)=0
      JDAHEP(2,NHEP)=0
C--
C--   Modify momentum/energy for cascade decays
      DO 150 I=JDAHEP(1,IP),JDAHEP(2,IP)
        IF (JDAHEP(1,I).NE.0) THEN
C--
C--   Reconstruct original energy/momentum of daughter
          DO 130 J=1,5
  130     PBOOST(J)=0.
          FIRST=JDAHEP(1,I)
          LAST=JDAHEP(2,I)
          DO 140 K=FIRST,LAST
            DO 140 J=1,4
  140     PBOOST(J)=PBOOST(J)+PHEP(J,K)
          PBOOST(5)=PHEP(5,I)
C--
C--   Correct energy/momentum of cascade daughters
          IF (JMOHEP(1,I).EQ.IP) THEN
            II=I
            CALL PHOBOS(II,PBOOST,PHEP(1,I),FIRST,LAST)
          ENDIF
        ENDIF
  150 CONTINUE
      END
CDECK  ID>, PHOBOS. 
      SUBROUTINE PHOBOS(IP,PBOOS1,PBOOS2,FIRST,LAST)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays BOoSt routine
C.
C.    Purpose:  Boost particles  in  cascade decay  to parent rest frame
C.              and boost back with modified boost vector.
C.
C.    Input Parameters:       IP:  pointer of particle starting chain
C.                                 to be boosted
C.                        PBOOS1:  Boost vector to rest frame,
C.                        PBOOS2:  Boost vector to modified frame,
C.                        FIRST:   Pointer to first particle to be boos-
C.                                 ted (/HEPEVT/),
C.                        LAST:    Pointer to last  particle to be boos-
C.                                 ted (/HEPEVT/).
C.
C.    Output Parameters:  Common /HEPEVT/.
C.
C.    Author(s):  B. van Eijk                     Created at:  13/02/90
C.                Z. Was                          Last Update: 21/12/92
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION BET1(3),BET2(3),GAM1,GAM2,PB,DATA
      INTEGER I,J,FIRST,LAST,MAXSTA,NSTACK,IP
      PARAMETER (MAXSTA=100)
      INTEGER STACK(MAXSTA)
      REAL PBOOS1(5),PBOOS2(5)
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      NSTACK=0
      DO 10 J=1,3
        BET1(J)=-PBOOS1(J)/PBOOS1(5)
   10 BET2(J)=PBOOS2(J)/PBOOS2(5)
      GAM1=PBOOS1(4)/PBOOS1(5)
      GAM2=PBOOS2(4)/PBOOS2(5)
C--
C--   Boost vector to parent rest frame...
   20 DO 50 I=FIRST,LAST
        PB=BET1(1)*PHEP(1,I)+BET1(2)*PHEP(2,I)+BET1(3)*PHEP(3,I)
        IF (JMOHEP(1,I).EQ.IP) THEN
         DO 30 J=1,3
   30    PHEP(J,I)=PHEP(J,I)+BET1(J)*(PHEP(4,I)+PB/(GAM1+1.))
         PHEP(4,I)=GAM1*PHEP(4,I)+PB
C--
C--    ...and boost back to modified parent frame.
         PB=BET2(1)*PHEP(1,I)+BET2(2)*PHEP(2,I)+BET2(3)*PHEP(3,I)
         DO 40 J=1,3
   40    PHEP(J,I)=PHEP(J,I)+BET2(J)*(PHEP(4,I)+PB/(GAM2+1.))
         PHEP(4,I)=GAM2*PHEP(4,I)+PB
         IF (JDAHEP(1,I).NE.0) THEN
           NSTACK=NSTACK+1
C--
C--    Check on stack length...
           IF (NSTACK.GT.MAXSTA) THEN
             DATA=NSTACK
             CALL PHOERR(7,'PHOBOS',DATA)
           ENDIF
           STACK(NSTACK)=I
         ENDIF
        ENDIF
   50 CONTINUE
      IF (NSTACK.NE.0) THEN
C--
C--   Now go one step further in the decay tree...
        FIRST=JDAHEP(1,STACK(NSTACK))
        LAST=JDAHEP(2,STACK(NSTACK))
        IP=STACK(NSTACK)
        NSTACK=NSTACK-1
        GOTO 20
      ENDIF
      RETURN
      END
CDECK  ID>, PHOTRI. 
      FUNCTION PHOTRI(A,B,C)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays calculation of TRIangle fie
C.
C.    Purpose:  Calculation of triangle function for phase space.
C.
C.    Input Parameters:  A, B, C (Virtual) particle masses.
C.
C.    Output Parameter:  Function value =
C.                       SQRT(LAMBDA(A**2,B**2,C**2))/(2*A)
C.
C.    Author(s):  B. van Eijk                     Created at:  15/11/89
C.                                                Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION DA,DB,DC,DAPB,DAMB,DTRIAN
      REAL A,B,C,PHOTRI
      DA=A
      DB=B
      DC=C
      DAPB=DA+DB
      DAMB=DA-DB
      DTRIAN=SQRT((DAMB-DC)*(DAPB+DC)*(DAMB+DC)*(DAPB-DC))
      PHOTRI=DTRIAN/(DA+DA)
      RETURN
      END
CDECK  ID>, PHOAN1. 
      FUNCTION PHOAN1(X,Y)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays calculation of ANgle '1'
C.
C.    Purpose:  Calculate angle from X and Y
C.
C.    Input Parameters:  X, Y
C.
C.    Output Parameter:  Function value
C.
C.    Author(s):  S. Jadach                       Created at:  01/01/89
C.                B. van Eijk                     Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION PHOAN1
      REAL X,Y
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      IF (ABS(Y).LT.ABS(X)) THEN
        PHOAN1=ATAN(ABS(Y/X))
        IF (X.LE.0.) PHOAN1=PI-PHOAN1
      ELSE
        PHOAN1=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      IF (Y.LT.0.) PHOAN1=TWOPI-PHOAN1
      RETURN
      END
CDECK  ID>, PHOAN2. 
      FUNCTION PHOAN2(X,Y)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays calculation of ANgle '2'
C.
C.    Purpose:  Calculate angle from X and Y
C.
C.    Input Parameters:  X, Y
C.
C.    Output Parameter:  Function value
C.
C.    Author(s):  S. Jadach                       Created at:  01/01/89
C.                B. van Eijk                     Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION PHOAN2
      REAL X,Y
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      IF (ABS(Y).LT.ABS(X)) THEN
        PHOAN2=ATAN(ABS(Y/X))
        IF (X.LE.0.) PHOAN2=PI-PHOAN2
      ELSE
        PHOAN2=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      RETURN
      END
CDECK  ID>, PHOBO3. 
      SUBROUTINE PHOBO3(ANGLE,PVEC)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays BOost routine '3'
C.
C.    Purpose:  Boost  vector PVEC  along z-axis where ANGLE = EXP(ETA),
C.              ETA is the hyperbolic velocity.
C.
C.    Input Parameters:  ANGLE, PVEC
C.
C.    Output Parameter:  PVEC
C.
C.    Author(s):  S. Jadach                       Created at:  01/01/89
C.                B. van Eijk                     Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION QPL,QMI,ANGLE
      REAL PVEC(4)
      QPL=(PVEC(4)+PVEC(3))*ANGLE
      QMI=(PVEC(4)-PVEC(3))/ANGLE
      PVEC(3)=(QPL-QMI)/2.
      PVEC(4)=(QPL+QMI)/2.
      RETURN
      END
CDECK  ID>, PHORO2. 
      SUBROUTINE PHORO2(ANGLE,PVEC)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays ROtation routine '2'
C.
C.    Purpose:  Rotate  x and z components  of vector PVEC  around angle
C.              'ANGLE'.
C.
C.    Input Parameters:  ANGLE, PVEC
C.
C.    Output Parameter:  PVEC
C.
C.    Author(s):  S. Jadach                       Created at:  01/01/89
C.                B. van Eijk                     Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION CS,SN,ANGLE
      REAL PVEC(4)
      CS=COS(ANGLE)*PVEC(1)+SIN(ANGLE)*PVEC(3)
      SN=-SIN(ANGLE)*PVEC(1)+COS(ANGLE)*PVEC(3)
      PVEC(1)=CS
      PVEC(3)=SN
      RETURN
      END
CDECK  ID>, PHORO3. 
      SUBROUTINE PHORO3(ANGLE,PVEC)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays ROtation routine '3'
C.
C.    Purpose:  Rotate  x and y components  of vector PVEC  around angle
C.              'ANGLE'.
C.
C.    Input Parameters:  ANGLE, PVEC
C.
C.    Output Parameter:  PVEC
C.
C.    Author(s):  S. Jadach                       Created at:  01/01/89
C.                B. van Eijk                     Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION CS,SN,ANGLE
      REAL PVEC(4)
      CS=COS(ANGLE)*PVEC(1)-SIN(ANGLE)*PVEC(2)
      SN=SIN(ANGLE)*PVEC(1)+COS(ANGLE)*PVEC(2)
      PVEC(1)=CS
      PVEC(2)=SN
      RETURN
      END
CDECK  ID>, PHORIN. 
      SUBROUTINE PHORIN
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation  in decays RANdom number generator init
C.
C.    Purpose:  Initialse PHORAN  with  the user  specified seeds in the
C.              array ISEED.  For details  see also:  F. James  CERN DD-
C.              Report November 1988.
C.
C.    Input Parameters:   ISEED(*)
C.
C.    Output Parameters:  URAN, CRAN, CDRAN, CMRAN, I97, J97
C.
C.    Author(s):  B. van Eijk and F. James        Created at:  27/09/89
C.                                                Last Update: 22/02/90
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION DATA
      REAL S,T
      INTEGER I,IS1,IS2,IS3,IS4,IS5,J
      INTEGER ISEED,I97,J97
      REAL URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
C--
C--   Check value range of seeds
      IF ((ISEED(1).LT.0).OR.(ISEED(1).GE.31328)) THEN
        DATA=ISEED(1)
        CALL PHOERR(8,'PHORIN',DATA)
      ENDIF
      IF ((ISEED(2).LT.0).OR.(ISEED(2).GE.30081)) THEN
        DATA=ISEED(2)
        CALL PHOERR(9,'PHORIN',DATA)
      ENDIF
C--
C--   Calculate Marsaglia and Zaman seeds (by F. James)
      IS1=MOD(ISEED(1)/177,177)+2
      IS2=MOD(ISEED(1),177)+2
      IS3=MOD(ISEED(2)/169,178)+1
      IS4=MOD(ISEED(2),169)
      DO 20 I=1,97
        S=0.
        T=0.5
        DO 10 J=1,24
          IS5=MOD (MOD(IS1*IS2,179)*IS3,179)
          IS1=IS2
          IS2=IS3
          IS3=IS5
          IS4=MOD(53*IS4+1,169)
          IF (MOD(IS4*IS5,64).GE.32) S=S+T
   10   T=0.5*T
   20 URAN(I)=S
      CRAN=362436./16777216.
      CDRAN=7654321./16777216.
      CMRAN=16777213./16777216.
      I97=97
      J97=33
      RETURN
      END
CDECK  ID>, PHORAN. 
      FUNCTION PHORAN(IDUM)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays RANdom number generator based
C.              on Marsaglia Algorithm
C.
C.    Purpose:  Generate  uniformly  distributed  random numbers between
C.              0 and 1.  Super long period:  2**144.  See also:
C.              G. Marsaglia and A. Zaman,  FSU-SCR-87-50,  for seed mo-
C.              difications  to  this version  see:  F. James DD-Report,
C.              November 1988.  The generator  has  to be initialized by
C.              a call to PHORIN.
C.
C.    Input Parameters:   IDUM (integer dummy)
C.
C.    Output Parameters:  Function value
C.
C.    Author(s):  B. van Eijk, G. Marsaglia and   Created at:  27/09/89
C.                A. Zaman                        Last Update: 27/09/89
C.
C.----------------------------------------------------------------------
      REAL PHORAN
      INTEGER IDUM
      INTEGER ISEED,I97,J97
      REAL URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
   10 PHORAN=URAN(I97)-URAN(J97)
      IF (PHORAN.LT.0.) PHORAN=PHORAN+1.
      URAN(I97)=PHORAN
      I97=I97-1
      IF (I97.EQ.0) I97=97
      J97=J97-1
      IF (J97.EQ.0) J97=97
      CRAN=CRAN-CDRAN
      IF (CRAN.LT.0.) CRAN=CRAN+CMRAN
      PHORAN=PHORAN-CRAN
      IF (PHORAN.LT.0.) PHORAN=PHORAN+1.
      IF (PHORAN.LE.0.) GOTO 10
      RETURN
      END
CDECK  ID>, PHOCHA. 
      FUNCTION PHOCHA(IDHEP)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays CHArge determination
C.
C.    Purpose:  Calculate the charge  of particle  with code IDHEP.  The
C.              code  of the  particle  is  defined by the Particle Data
C.              Group in Phys. Lett. B204 (1988) 1.
C.
C.    Input Parameter:   IDHEP
C.
C.    Output Parameter:  Funtion value = charge  of  particle  with code
C.                       IDHEP
C.
C.    Author(s):  E. Barberio and B. van Eijk     Created at:  29/11/89
C.                                                Last update: 02/01/90
C.
C.----------------------------------------------------------------------
      REAL PHOCHA
      INTEGER IDHEP,IDABS,Q1,Q2,Q3
C--
C--   Array 'CHARGE' contains the charge  of the first 101 particles ac-
C--   cording  to  the PDG particle code... (0 is added for convenience)
      REAL CHARGE(0:100)
      DATA CHARGE/ 0.,
     &-0.3333333333,  0.6666666667, -0.3333333333, 0.6666666667,
     &-0.3333333333,  0.6666666667, -0.3333333333, 0.6666666667,
     & 2*0., -1., 0., -1., 0., -1., 0., -1., 6*0., 1., 12*0., 1., 63*0./
      IDABS=ABS(IDHEP)
      IF (IDABS.LE.100) THEN
C--
C--   Charge of quark, lepton, boson etc....
        PHOCHA = CHARGE(IDABS)
      ELSE
C--
C--   Check on particle build out of quarks, unpack its code...
        Q3=MOD(IDABS/1000,10)
        Q2=MOD(IDABS/100,10)
        Q1=MOD(IDABS/10,10)
        IF (Q3.EQ.0) THEN
C--
C--   ...meson...
          IF(MOD(Q2,2).EQ.0) THEN
            PHOCHA=CHARGE(Q2)-CHARGE(Q1)
          ELSE
            PHOCHA=CHARGE(Q1)-CHARGE(Q2)
          ENDIF
        ELSE
C--
C--   ...diquarks or baryon.
          PHOCHA=CHARGE(Q1)+CHARGE(Q2)+CHARGE(Q3)
        ENDIF
      ENDIF
C--
C--   Find the sign of the charge...
      IF (IDHEP.LT.0.) PHOCHA=-PHOCHA
      RETURN
      END
CDECK  ID>, PHOSPI. 
      FUNCTION PHOSPI(IDHEP)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation  in decays function for SPIn determina-
C.              tion
C.
C.    Purpose:  Calculate  the spin  of particle  with  code IDHEP.  The
C.              code  of the particle  is  defined  by the Particle Data
C.              Group in Phys. Lett. B204 (1988) 1.
C.
C.    Input Parameter:   IDHEP
C.
C.    Output Parameter:  Funtion  value = spin  of  particle  with  code
C.                       IDHEP
C.
C.    Author(s):  E. Barberio and B. van Eijk     Created at:  29/11/89
C.                                                Last update: 02/01/90
C.
C.----------------------------------------------------------------------
      REAL PHOSPI
      INTEGER IDHEP,IDABS
C--
C--   Array 'SPIN' contains the spin  of  the first 100 particles accor-
C--   ding to the PDG particle code...
      REAL SPIN(100)
      DATA SPIN/ 8*.5, 1., 0., 8*.5, 2*0., 4*1., 76*0./
      IDABS=ABS(IDHEP)
C--
C--   Spin of quark, lepton, boson etc....
      IF (IDABS.LE.100) THEN
        PHOSPI=SPIN(IDABS)
      ELSE
C--
C--   ...other particles, however...
        PHOSPI=(MOD(IDABS,10)-1.)/2.
C--
C--   ...K_short and K_long are special !!
        PHOSPI=MAX(PHOSPI,0.)
      ENDIF
      RETURN
      END
CDECK  ID>, PHOOMA. 
      SUBROUTINE PHOOMA(IFIRST,ILAST,POINTR)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays Order MAss vector
C.
C.    Purpose:  Order  the  contents  of array 'POINTR' according to the
C.              decreasing value in the array 'MASS'.
C.
C.    Input Parameters:  IFIRST, ILAST:  Pointers  to  the  vector loca-
C.                                       tion be sorted,
C.                       POINTR:         Unsorted array with pointers to
C.                                       /HEPEVT/.
C.
C.    Output Parameter:  POINTR:         Sorted arrays  with  respect to
C.                                       particle mass 'PHEP(5,*)'.
C.
C.    Author(s):  B. van Eijk                     Created at:  28/11/89
C.                                                Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER PDMMAX
      PARAMETER (PDMMAX=100)
      INTEGER IFIRST,ILAST,I,J,BUFPOI,POINTR(PDMMAX)
      REAL BUFMAS,MASS(PDMMAX)
      IF (IFIRST.EQ.ILAST) RETURN
C--
C--   Copy particle masses
      DO 10 I=IFIRST,ILAST
   10 MASS(I)=PHEP(5,POINTR(I))
C--
C--   Order the masses in a decreasing series
      DO 30 I=IFIRST,ILAST-1
        DO 20 J=I+1,ILAST
          IF (MASS(J).LE.MASS(I)) GOTO 20
          BUFPOI=POINTR(J)
          POINTR(J)=POINTR(I)
          POINTR(I)=BUFPOI
          BUFMAS=MASS(J)
          MASS(J)=MASS(I)
          MASS(I)=BUFMAS
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
CDECK  ID>, PHOERR. 
      SUBROUTINE PHOERR(IMES,TEXT,DATA)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays ERRror handling
C.
C.    Purpose:  Inform user  about (fatal) errors and warnings generated
C.              by either the user or the program.
C.
C.    Input Parameters:   IMES, TEXT, DATA
C.
C.    Output Parameters:  None
C.
C.    Author(s):  B. van Eijk                     Created at:  29/11/89
C.                                                Last Update: 10/01/92
C.
C.----------------------------------------------------------------------
      DOUBLE PRECISION DATA
      INTEGER DATA1,DATA2,IMES,IERROR
      REAL SDATA
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      CHARACTER TEXT*(*)
      SAVE IERROR
      DATA IERROR/ 0/
      IF (IMES.LE.PHOMES) STATUS(IMES)=STATUS(IMES)+1
C--
C--   Count number of non-fatal errors...
      IF ((IMES.EQ. 6).AND.(STATUS(IMES).GE.2)) RETURN
      IF ((IMES.EQ.10).AND.(STATUS(IMES).GE.2)) RETURN
      SDATA=DATA
      WRITE(PHLUN,9000)
      WRITE(PHLUN,9120)
      GOTO (10,20,30,40,50,60,70,80,90,100),IMES
      WRITE(PHLUN,9130) IMES
      GOTO 120
   10 WRITE(PHLUN,9010) TEXT,INT(SDATA)
      GOTO 110
   20 WRITE(PHLUN,9020) TEXT,SDATA
      GOTO 110
   30 WRITE(PHLUN,9030) TEXT,SDATA
      GOTO 110
   40 WRITE(PHLUN,9040) TEXT
      GOTO 110
   50 WRITE(PHLUN,9050) TEXT,INT(SDATA)
      GOTO 110
   60 WRITE(PHLUN,9060) TEXT,SDATA
      GOTO 130
   70 WRITE(PHLUN,9070) TEXT,INT(SDATA)
      GOTO 110
   80 WRITE(PHLUN,9080) TEXT,INT(SDATA)
      GOTO 110
   90 WRITE(PHLUN,9090) TEXT,INT(SDATA)
      GOTO 110
  100 WRITE(PHLUN,9100) TEXT,SDATA
      GOTO 130
  110 CONTINUE
      WRITE(PHLUN,9140)
      WRITE(PHLUN,9120)
      WRITE(PHLUN,9000)
      STOP
  120 IERROR=IERROR+1
      IF (IERROR.GE.10) THEN
        WRITE(PHLUN,9150)
        WRITE(PHLUN,9120)
        WRITE(PHLUN,9000)
        STOP
      ENDIF
  130 WRITE(PHLUN,9120)
      WRITE(PHLUN,9000)
      RETURN
 9000 FORMAT(1H ,80('*'))
 9010 FORMAT(1H ,'* ',A,': Too many charged Particles, NCHARG =',I6,T81,
     &'*')
 9020 FORMAT(1H ,'* ',A,': Too much Bremsstrahlung required, PRSOFT = ',
     &F15.6,T81,'*')
 9030 FORMAT(1H ,'* ',A,': Combined Weight is exceeding 1., Weight = ',
     &F15.6,T81,'*')
 9040 FORMAT(1H ,'* ',A,
     &': Error in Rescaling charged and neutral Vectors',T81,'*')
 9050 FORMAT(1H ,'* ',A,
     &': Non matching charged Particle Pointer, NCHARG = ',I5,T81,'*')
 9060 FORMAT(1H ,'* ',A,
     &': Do you really work with a Particle of Spin: ',F4.1,' ?',T81,
     &'*')
 9070 FORMAT(1H ,'* ',A, ': Stack Length exceeded, NSTACK = ',I5 ,T81,
     &'*')
 9080 FORMAT(1H ,'* ',A,
     &': Random Number Generator Seed(1) out of Range: ',I8,T81,'*')
 9090 FORMAT(1H ,'* ',A,
     &': Random Number Generator Seed(2) out of Range: ',I8,T81,'*')
 9100 FORMAT(1H ,'* ',A,
     &': Available Phase Space below Cut-off: ',F15.6,' GeV/c^2',T81,
     &'*')
 9120 FORMAT(1H ,'*',T81,'*')
 9130 FORMAT(1H ,'* Funny Error Message: ',I4,' ! What to do ?',T81,'*')
 9140 FORMAT(1H ,'* Fatal Error Message, I stop this Run !',T81,'*')
 9150 FORMAT(1H ,'* 10 Error Messages generated, I stop this Run !',T81,
     &'*')
      END
CDECK  ID>, PHOREP. 
      SUBROUTINE PHOREP
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays run summary REPort
C.
C.    Purpose:  Inform user about success and/or restrictions of PHOTOS
C.              encountered during execution.
C.
C.    Input Parameters:   Common /PHOSTA/
C.
C.    Output Parameters:  None
C.
C.    Author(s):  B. van Eijk                     Created at:  10/01/92
C.                                                Last Update: 10/01/92
C.
C.----------------------------------------------------------------------
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      INTEGER I
      LOGICAL ERROR
      ERROR=.FALSE.
      WRITE(PHLUN,9000)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9020)
      WRITE(PHLUN,9030)
      WRITE(PHLUN,9040)
      WRITE(PHLUN,9030)
      WRITE(PHLUN,9020)
      DO 10 I=1,PHOMES
        IF (STATUS(I).EQ.0) GOTO 10
        IF ((I.EQ.6).OR.(I.EQ.10)) THEN
          WRITE(PHLUN,9050) I,STATUS(I)
        ELSE
          ERROR=.TRUE.
          WRITE(PHLUN,9060) I,STATUS(I)
        ENDIF
   10 CONTINUE
      IF (.NOT.ERROR) WRITE(PHLUN,9070)
      WRITE(PHLUN,9020)
      WRITE(PHLUN,9010)
      RETURN
 9000 FORMAT(1H1)
 9010 FORMAT(1H ,80('*'))
 9020 FORMAT(1H ,'*',T81,'*')
 9030 FORMAT(1H ,'*',26X,25('='),T81,'*')
 9040 FORMAT(1H ,'*',30X,'PHOTOS Run Summary',T81,'*')
 9050 FORMAT(1H ,'*',22X,'Warning #',I2,' occured',I6,' times',T81,'*')
 9060 FORMAT(1H ,'*',23X,'Error #',I2,' occured',I6,' times',T81,'*')
 9070 FORMAT(1H ,'*',16X,'PHOTOS Execution has successfully terminated',
     &T81,'*')
      END
