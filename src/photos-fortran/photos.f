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
