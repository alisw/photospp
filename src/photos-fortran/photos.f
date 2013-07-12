      SUBROUTINE PHTYPE(ID)
C.----------------------------------------------------------------------
C.
C.    PHTYPE:   Central manadgement routine.              
C.
C.    Purpose:   defines what kind of the 
C.              actions will be performed at point ID. 
C.
C.    Input Parameters:       ID:  pointer of particle starting branch
C.                                 in /PH_HEPEVT/ to be treated.
C.
C.    Output Parameters:  Common /PH_HEPEVT/.
C.
C.    Author(s):  Z. Was                          Created at:  24/05/93
C.                P. Golonka                      Last Update: 27/06/04
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      REAL*8 FINT,FSEC,EXPEPS
      COMMON /PHOKEY/ FSEC,FINT,EXPEPS,INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      LOGICAL EXPINI
      INTEGER NX,K,NCHAN
      PARAMETER (NX=10)
      REAL*8 PRO,PRSUM,ESU
      COMMON /PHOEXP/ PRO(NX),NCHAN,EXPINI

      INTEGER ID,NHEP0
      LOGICAL IPAIR
      REAL*8 RN,PHORANC,SUM
      INTEGER WTDUM
      LOGICAL IFOUR
C--
      IFOUR=(.TRUE.).AND.(ITRE) ! we can make internal choice whether 
                                ! we want 3 or four photons at most.
      IPAIR=.TRUE.
C--   Check decay multiplicity..
      IF (JDAHEP(1,ID).EQ.0) RETURN
C      IF (JDAHEP(1,ID).EQ.JDAHEP(2,ID)) RETURN
C--
      NHEP0=NHEP
C--
      IF    (IEXP)  THEN
         EXPINI=.TRUE.   ! Initialization/cleaning
         DO NCHAN=1,NX
           PRO(NCHAN)=0.D0
         ENDDO
         NCHAN=0
         
         FSEC=1.0D0
         CALL PHOMAK(ID,NHEP0)! Initialization/crude formfactors into 
                                                   ! PRO(NCHAN)
         EXPINI=.FALSE.
         RN=PHORANC(WTDUM)
         PRSUM=0
         DO K=1,NX
          PRSUM=PRSUM+PRO(K)
         ENDDO
         ESU=EXP(-PRSUM) ! exponent for crude Poissonian multiplicity 
                         ! distribution, will be later overwritten 
                         ! to give probability for k
         SUM=ESU         ! distribuant for the crude Poissonian 
                         ! at first for k=0
         DO K=1,100      ! hard coded max (photon) multiplicity is 100
           IF(RN.LT.SUM) GOTO 100
           ESU=ESU*PRSUM/K  ! we get at K ESU=EXP(-PRSUM)*PRSUM**K/K!
           SUM=SUM+ESU      ! thus we get distribuant at K.
           NCHAN=0
           CALL PHOMAK(ID,NHEP0) ! LOOPING
           IF(SUM.GT.1D0-EXPEPS) GOTO 100
         ENDDO
 100     CONTINUE
      ELSEIF(IFOUR) THEN
C-- quatro photon emission
        FSEC=1.0D0
        RN=PHORANC(WTDUM)
        IF (RN.GE.23.D0/24D0) THEN
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
        ELSEIF (RN.GE.17.D0/24D0) THEN
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
        ELSEIF (RN.GE.9.D0/24D0) THEN
          CALL PHOMAK(ID,NHEP0)
        ENDIF
      ELSEIF(ITRE) THEN
C-- triple photon emission
        FSEC=1.0D0
        RN=PHORANC(WTDUM)
        IF (RN.GE.5.D0/6D0) THEN
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
        ELSEIF (RN.GE.2.D0/6D0) THEN
          CALL PHOMAK(ID,NHEP0)
        ENDIF
      ELSEIF(ISEC) THEN
C-- double photon emission
        FSEC=1.0D0
        RN=PHORANC(WTDUM)
        IF (RN.GE.0.5D0) THEN
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
        ENDIF
      ELSE
C-- single photon emission
        FSEC=1.0D0
        CALL PHOMAK(ID,NHEP0)
      ENDIF
C--
C-- electron positron pair (coomented out for a while
C      IF (IPAIR) CALL PHOPAR(ID,NHEP0)
      END  
      SUBROUTINE PHOMAK(IPPAR,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOMAK:   PHOtos MAKe
C.
C.    Purpose:  Single or double bremstrahlung radiative corrections  
C.              are generated in  the decay of the IPPAR-th particle in 
C.              the  HEP common /PH_HEPEVT/. Example of the use of 
C.              general tools.
C.
C.    Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
C.                                /PH_HEPEVT/ and the common itself
C.
C.    Output Parameters:  Common  /PH_HEPEVT/, either  with  or  without
C.                                particles added.
C.
C.    Author(s):  Z. Was,                         Created at:  26/05/93
C.                                                Last Update: 29/01/05
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION DATA
      REAL*8 PHORANC
      INTEGER IP,IPPAR,NCHARG,IDME
      INTEGER WTDUM,IDUM,NHEP0
      INTEGER NCHARB,NEUDAU
      REAL*8 RN,WT,PHINT,XDUMM,PHwtNLO
      LOGICAL BOOST
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      REAL*8 FINT,FSEC,EXPEPS
      COMMON /PHOKEY/ FSEC,FINT,EXPEPS,INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
C--
      IP=IPPAR
      IDUM=1
      NCHARG=0
C--
        CALL PHOIN(IP,BOOST,NHEP0)
        CALL PHOCHK(JDAHEP(1,IP))
        WT=0.0D0
        CALL PHOPRE(1,WT,NEUDAU,NCHARB)

        IF (WT.EQ.0.0D0) RETURN
        RN=PHORANC(WTDUM)
C PHODO is caling PHORANC, thus change of series if it is moved before if
        CALL PHODO(1,NCHARB,NEUDAU)
C we eliminate /FINT in variant B.
C get ID of channel dependent ME, ID=0 means no 

        CALL ME_CHANNEL(IDME)                 ! corrections for matrix elements
                                              ! controlled by IDME
                                              ! write(*,*) 'KANALIK IDME=',IDME

        IF(     IDME.EQ.0) THEN               ! default 

          IF (INTERF) WT=WT*PHINT(IDUM)/FINT  ! FINT must be in variant A
          IF (IFW) CALL PHOBW   (WT)          ! extra weight for leptonic W decay 

        ELSEIF (IDME.EQ.2) THEN               ! ME weight for leptonic W decay

          CALL PHOBWnlo(WT)
          WT=WT*2D0/FINT

        ELSEIF (IDME.EQ.1) THEN               ! ME weight for leptonic Z decay

         xdumm=0.5D0
         WT=WT*PHwtnlo(xdumm)/FINT

        ELSE
         write(*,*) 'problem with ME_CHANNEL  IDME=',IDME
         stop
        ENDIF


        DATA=WT 
        IF (WT.GT.1.0D0) CALL PHOERR(3,'WT_INT',DATA)
C weighting
      IF (RN.LE.WT) THEN 
        CALL PHOOUT(IP,BOOST,NHEP0)
      ENDIF
      RETURN
      END
