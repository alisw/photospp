C********************************************************
C Example of use of Photos C++ interface.
C D+ D- -> tau+ tau- HEPEVT events are constructed.
C Taus are subsequently decayed via Photos++.
C
C @author Tomasz Przedzinski
C @date 21 August 2013
C********************************************************
      PROGRAM PH_HEPVT_TEST

C     INITIALIZE PHOTOS++
      CALL PHOTOS_INIT

C     PREPARE SIMPLE EVENT in HEPEVT
      CALL SIMPLE_EVENT
      WRITE(*,*) "##############"
      WRITE(*,*) "PHODMP: BEFORE"
      WRITE(*,*) "##############"
      CALL PHODMP

      CALL PHOTOS_PROCESS
C      CALL PHOTOS_PROCESS_PARTICLE(4)
C      CALL PHOTOS_PROCESS_BRANCH(4)

      WRITE(*,*) "##############"
      WRITE(*,*) "PHODMP: AFTER"
      WRITE(*,*) "##############"
      CALL PHODMP
      END

      SUBROUTINE SIMPLE_EVENT
      REAL*8 AMELL

      AMELL = 0.0005111;

C     CREATE SIMPLE EVENT: e+ e- -> Z -> tau+ tau- -> pi nu_tau,  pi nu_tau
C     ---------------------------------------------------------------------
      CALL ADD_PARTICLE(  11, 6,  1.7763568394002505D-15, -3.5565894425761324D-15,  4.5521681043409913D+01, 4.5521681043409934D+01,
     $                    AMELL,                   0,  0,  3,  3)
      CALL ADD_PARTICLE( -11, 6, -1.7763568394002505D-15,  3.5488352204797800D-15, -4.5584999071936601D+01, 4.5584999071936622D+01,
     $                    AMELL,                   0,  0,  3,  3)
      CALL ADD_PARTICLE(  23, 5,  0D0,                     0D0,                    -6.3318028526687442D-02, 9.1106680115346506D+01,
     $                    9.1106658112716090D+01,  1,  2,  4,  5)
      CALL ADD_PARTICLE(  15, 2, -2.3191595992562256D+01, -2.6310500920665142D+01, -2.9046412466624929D+01, 4.5573504956498098D+01,
     $                    1.7769900000002097D+00,  3,  0,  6,  7)
      CALL ADD_PARTICLE( -15, 2,  2.3191595992562256D+01,  2.6310500920665142D+01,  2.8983094438098242D+01, 4.5533175158848429D+01,
     $                    1.7769900000000818D+00,  3,  0,  8,  9)
      CALL ADD_PARTICLE(  16, 1, -1.2566536214715378D+00, -1.7970251138317268D+00, -1.3801323581022720D+00, 2.5910119010468553D+00,
     $                    9.9872238934040070D-03,  4,  0,  0, 0)
      CALL ADD_PARTICLE(-211, 1, -2.1935073012334062D+01, -2.4513624017269400D+01, -2.7666443730700312D+01, 4.2982749776866747D+01,
     $                    1.3956783711910248D-01,  4,  0,  0, 0)
      CALL ADD_PARTICLE( -16, 1,  8.4364531743909055D+00,  8.3202830831667836D+00,  9.6202800273055971D+00, 1.5262723881157640D+01,
     $                    9.9829332903027534D-03,  5,  0,  0, 0)
      CALL ADD_PARTICLE( 211, 1,  1.4755273459419701D+01,  1.7990366047940022D+01,  1.9362977676297948D+01, 3.0270707771933196D+01,
     $                    1.3956753909587860D-01,  5,  0,  0, 0)
      END

      SUBROUTINE ADD_PARTICLE(ID,STATUS,PX,PY,PZ,E,M,MOTHER1,MOTHER2,DAUGHTER1,DAUGHTER2)
      INTEGER ID,STATUS,MOTHER1,MOTHER2,DAUGHTER1,DAUGHTER2
      REAL*8  PX,PY,PZ,E,M

      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      REAL*8  phep,  vhep ! to be real*4/ *8  depending on host
      INTEGER nevhep,nhep,isthep,idhep,jmohep,
     $        jdahep
      COMMON /hepevt/
     $      nevhep,               ! serial number
     $      nhep,                 ! number of particles
     $      isthep(nmxhep),   ! status code
     $      idhep(nmxhep),    ! particle ident KF
     $      jmohep(2,nmxhep), ! parent particles
     $      jdahep(2,nmxhep), ! childreen particles
     $      phep(5,nmxhep),   ! four-momentum, mass [GeV]
     $      vhep(4,nmxhep)    ! vertex [mm]
      SAVE hepevt
      NHEP=NHEP+1
      IDHEP(NHEP) =ID
      ISTHEP(NHEP)=STATUS
      PHEP(1,NHEP)=PX
      PHEP(2,NHEP)=PY
      PHEP(3,NHEP)=PZ
      PHEP(4,NHEP)=E
      PHEP(5,NHEP)=M
      JMOHEP(1,NHEP)=MOTHER1
      JMOHEP(2,NHEP)=MOTHER2
      JDAHEP(1,NHEP)=DAUGHTER1
      JDAHEP(2,NHEP)=DAUGHTER2
      END

      SUBROUTINE PHODMP
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays event DuMP routine
C.
C.    Purpose:  Print event record.
C.
C.    Input Parameters:   Common /HEPEVT/
C.
C.    Output Parameters:  None
C.
C.    Author(s):  B. van Eijk                     Created at:  05/06/90
C.                                                Last Update: 05/06/90
C.
C.----------------------------------------------------------------------
C      IMPLICIT NONE
      DOUBLE PRECISION  SUMVEC(5)
      INTEGER I,J
C this is the hepevt class in old style. No d_h_ class prD-name
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      REAL*8  phep,  vhep ! to be real*4/ *8  depending on host
      INTEGER nevhep,nhep,isthep,idhep,jmohep,
     $        jdahep
      COMMON /hepevt/
     $      nevhep,               ! serial number
     $      nhep,                 ! number of particles
     $      isthep(nmxhep),   ! status code
     $      idhep(nmxhep),    ! particle ident KF
     $      jmohep(2,nmxhep), ! parent particles
     $      jdahep(2,nmxhep), ! childreen particles
     $      phep(5,nmxhep),   ! four-momentum, mass [GeV]
     $      vhep(4,nmxhep)    ! vertex [mm]

      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      DO 10 I=1,5
   10 SUMVEC(I)=0.
C--
C--   Print event number...
      WRITE(PHLUN,9000)
      WRITE(PHLUN,9010) NEVHEP
      WRITE(PHLUN,9080)
      WRITE(PHLUN,9020)
      DO 30 I=1,NHEP
C--
C--   For 'stable particle' calculate vector momentum sum
        IF (JDAHEP(1,I).EQ.0) THEN
          DO 20 J=1,4
   20     SUMVEC(J)=SUMVEC(J)+PHEP(J,I)
          IF (JMOHEP(2,I).EQ.0) THEN
            WRITE(PHLUN,9030) I,IDHEP(I),JMOHEP(1,I),(PHEP(J,I),J=1,5)
          ELSE
            WRITE(PHLUN,9040) I,IDHEP(I),JMOHEP(1,I),JMOHEP(2,I),(PHEP
     &      (J,I),J=1,5)
          ENDIF
        ELSE
          IF (JMOHEP(2,I).EQ.0) THEN
            WRITE(PHLUN,9050) I,IDHEP(I),JMOHEP(1,I),JDAHEP(1,I),
     &      JDAHEP(2,I),(PHEP(J,I),J=1,5)
          ELSE
            WRITE(PHLUN,9060) I,IDHEP(I),JMOHEP(1,I),JMOHEP(2,I),
     &      JDAHEP(1,I),JDAHEP(2,I),(PHEP(J,I),J=1,5)
          ENDIF
        ENDIF
   30 CONTINUE
      SUMVEC(5)=SQRT(SUMVEC(4)**2-SUMVEC(1)**2-SUMVEC(2)**2-
     &SUMVEC(3)**2)
      WRITE(PHLUN,9070) (SUMVEC(J),J=1,5)
      RETURN
 9000 FORMAT(1H0,80('='))
 9010 FORMAT(1H ,29X,'Event No.:',I10)
 9020 FORMAT(1H0,1X,'Nr',3X,'Type',3X,'Parent(s)',2X,'Daughter(s)',6X,
     &'Px',7X,'Py',7X,'Pz',7X,'E',4X,'Inv. M.')
 9030 FORMAT(1H ,I4,I7,3X,I4,9X,'Stable',2X,5F9.2)
 9040 FORMAT(1H ,I4,I7,I4,' - ',I4,5X,'Stable',2X,5F9.2)
 9050 FORMAT(1H ,I4,I7,3X,I4,6X,I4,' - ',I4,5F9.2)
 9060 FORMAT(1H ,I4,I7,I4,' - ',I4,2X,I4,' - ',I4,5F9.2)
 9070 FORMAT(1H0,23X,'Vector Sum: ', 5F9.2)
 9080 FORMAT(1H0,6X,'Particle Parameters')
      END
