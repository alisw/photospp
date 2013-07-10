
*///////////////////////////////////////////////////////////////////////
*//                                                                     
*//  !!!!!!! WARNING!!!!!   This source may be  agressive !!!!          
*//                                                                     
*//  Due to short common block names it may owerwrite variables in other
*//  parts of the code.                                                 
*//                                                                     
*//  One should add suffix c_Photos_ to names of all commons as soon as 
*//  possible!!  
*///////////////////////////////////////////////////////////////////////

C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOtos CDE-s
C.
C.    Purpose:  Keep definitions  for PHOTOS QED correction Monte Carlo.
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  29/11/89
C.                                                Last Update: 10/08/93
C.
C. =========================================================
C.    General Structure Information:                       =
C. =========================================================
C:   ROUTINES:
C.             1) INITIALIZATION (all in C++ now)
C.             2) GENERAL INTERFACE:
C.                                      PHOBOS
C.                                      PHOIN
C.                                      PHOTWO (specific interface
C.                                      PHOOUT
C.                                      PHOCHK
C.                                      PHTYPE (specific interface
C.                                      PHOMAK (specific interface
C.             3) QED PHOTON GENERATION:
C.                                      PHINT
C.                                      PHOBW
C.                                      PHOPRE
C.                                      PHOOMA
C.                                      PHOENE
C.                                      PHOCOR
C.                                      PHOFAC
C.                                      PHODO
C.             4) UTILITIES:
C.                                      PHOTRI
C.                                      PHOAN1
C.                                      PHOAN2
C.                                      PHOBO3
C.                                      PHORO2
C.                                      PHORO3
C.                                      PHOCHA
C.                                      PHOSPI
C.                                      PHOERR
C.                                      PHOREP
C.                                      PHLUPA
C.                                      PHCORK
C.                                      IPHQRK
C.                                      IPHEKL
C.   COMMONS:
C.   NAME     USED IN SECT. # OF OCC.     Comment
C.   PHOQED   1) 2)            3      Flags whether emisson to be gen. 
C.   PHOLUN   1) 4)            6      Output device number
C.   PHOCOP   1) 3)            4      photon coupling & min energy
C.   PHPICO   1) 3) 4)         5      PI & 2*PI
C.   PHOSTA   1) 4)            3      Status information
C.   PHOKEY   1) 2) 3)         7      Keys for nonstandard application
C.   PHOVER   1)               1      Version info for outside
C.   HEPEVT   2)               2      PDG common
C.   PH_HEPEVT2)               8      PDG common internal
C.   PHOEVT   2) 3)           10      PDG branch
C.   PHOIF    2) 3)            2      emission flags for PDG branch 
C.   PHOMOM   3)               5      param of char-neutr system
C.   PHOPHS   3)               5      photon momentum parameters
C.   PHOPRO   3)               4      var. for photon rep. (in branch)
C.   PHOCMS   2)               3      parameters of boost to branch CMS
C.   PHNUM    4)               1      event number from outside         
C.----------------------------------------------------------------------
      SUBROUTINE PHOIN(IP,BOOST,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOIN:   PHOtos INput
C.
C.    Purpose:  copies IP branch of the common /PH_HEPEVT/ into /PHOEVT/
C.              moves branch into its CMS system.
C.
C.    Input Parameters:       IP:  pointer of particle starting branch
C.                                 to be copied
C.                        BOOST:   Flag whether boost to CMS was or was 
C     .                            not performed.
C.
C.    Output Parameters:  Commons: /PHOEVT/, /PHOCMS/
C.
C.    Author(s):  Z. Was                          Created at:  24/05/93
C.                                                Last Update: 16/11/93
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER IP,IP2,I,FIRST,LAST,LL,NA
      LOGICAL BOOST
      INTEGER J,NHEP0
      DOUBLE PRECISION BET(3),GAM,PB
      COMMON /PHOCMS/ BET,GAM
      LOGICAL INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      REAL*8 FINT,FSEC,EXPEPS
      COMMON /PHOKEY/ FSEC,FINT,EXPEPS,INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
C--
C let-s calculate size of the little common entry
        FIRST=JDAHEP(1,IP)
        LAST =JDAHEP(2,IP)
        NPHO=3+LAST-FIRST+NHEP-NHEP0
        NEVPHO=NPHO
C let-s take in decaying particle
           IDPHO(1)=IDHEP(IP)
           JDAPHO(1,1)=3
           JDAPHO(2,1)=3+LAST-FIRST
           DO I=1,5
             PPHO(I,1)=PHEP(I,IP)
           ENDDO
C let-s take in eventual second mother
         IP2=JMOHEP(2,JDAHEP(1,IP))
         IF((IP2.NE.0).AND.(IP2.NE.IP)) THEN 
           IDPHO(2)=IDHEP(IP2)
           JDAPHO(1,2)=3
           JDAPHO(2,2)=3+LAST-FIRST
           DO I=1,5
             PPHO(I,2)=PHEP(I,IP2)
           ENDDO
         ELSE
           IDPHO(2)=0
           DO I=1,5
             PPHO(I,2)=0.0D0
           ENDDO
         ENDIF
C let-s take in daughters
        DO LL=0,LAST-FIRST
           IDPHO(3+LL)=IDHEP(FIRST+LL)
           JMOPHO(1,3+LL)=JMOHEP(1,FIRST+LL)
           IF (JMOHEP(1,FIRST+LL).EQ.IP) JMOPHO(1,3+LL)=1
           DO I=1,5
             PPHO(I,3+LL)=PHEP(I,FIRST+LL)
           ENDDO
        ENDDO
        IF (NHEP.GT.NHEP0) THEN
C let-s take in illegitimate daughters
        NA=3+LAST-FIRST 
        DO LL=1,NHEP-NHEP0
           IDPHO(NA+LL)=IDHEP(NHEP0+LL)
           JMOPHO(1,NA+LL)=JMOHEP(1,NHEP0+LL)
           IF (JMOHEP(1,NHEP0+LL).EQ.IP) JMOPHO(1,NA+LL)=1
           DO I=1,5
             PPHO(I,NA+LL)=PHEP(I,NHEP0+LL)
           ENDDO
        ENDDO
C--        there is NHEP-NHEP0 daugters more.
           JDAPHO(2,1)=3+LAST-FIRST+NHEP-NHEP0
        ENDIF
        IF(IDPHO(NPHO).EQ.22)CALL PHLUPA(100001)
!        IF(IDPHO(NPHO).EQ.22) stop
        CALL PHCORK(0)
        IF(IDPHO(NPHO).EQ.22)CALL PHLUPA(100002)
C special case of t tbar production process
        IF(IFTOP) CALL PHOTWO(0)
        BOOST=.FALSE.
C--   Check whether parent is in its rest frame...
C ZBW ND  27.07.2009:
C bug reported by Vladimir Savinov localized and fixed.
C protection against rounding error was back-firing if soft
C momentum of mother was physical. Consequence was that PHCORK was
C messing up masses of final state particles in vertex of the decay.
C Only configurations with previously generated photons of energy fraction
C smaller than 0.0001 were affected. Effect was numerically insignificant. 

C      IF (     (ABS(PPHO(4,1)-PPHO(5,1)).GT.PPHO(5,1)*1.D-8)
C     $    .AND.(PPHO(5,1).NE.0))                            THEN

      IF ((ABS(PPHO(1,1)+ABS(PPHO(2,1))+ABS(PPHO(3,1))).GT.
     $    PPHO(5,1)*1.D-8).AND.(PPHO(5,1).NE.0))    THEN

        BOOST=.TRUE.

C--
C--   Boost daughter particles to rest frame of parent...
C--   Resultant neutral system already calculated in rest frame !
        DO 10 J=1,3
   10   BET(J)=-PPHO(J,1)/PPHO(5,1)
        GAM=PPHO(4,1)/PPHO(5,1)
        DO 30 I=JDAPHO(1,1),JDAPHO(2,1)
          PB=BET(1)*PPHO(1,I)+BET(2)*PPHO(2,I)+BET(3)*PPHO(3,I)
          DO 20 J=1,3
   20     PPHO(J,I)=PPHO(J,I)+BET(J)*(PPHO(4,I)+PB/(GAM+1.D0))
   30   PPHO(4,I)=GAM*PPHO(4,I)+PB
C--    Finally boost mother as well
          I=1   
          PB=BET(1)*PPHO(1,I)+BET(2)*PPHO(2,I)+BET(3)*PPHO(3,I)
          DO J=1,3
            PPHO(J,I)=PPHO(J,I)+BET(J)*(PPHO(4,I)+PB/(GAM+1.D0))
          ENDDO
          PPHO(4,I)=GAM*PPHO(4,I)+PB
      ENDIF
C special case of t tbar production process
        IF(IFTOP) CALL PHOTWO(1)
      CALL PHLUPA(2)
        IF(IDPHO(NPHO).EQ.22) CALL PHLUPA(10000)
!        IF(IDPHO(NPHO-1).EQ.22) stop
      END 
      SUBROUTINE PHOOUT(IP,BOOST,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOOUT:   PHOtos OUTput
C.
C.    Purpose:  copies back IP branch of the common /PH_HEPEVT/ from 
C.              /PHOEVT/ moves branch back from its CMS system.
C.
C.    Input Parameters:       IP:  pointer of particle starting branch
C.                                 to be given back.
C.                        BOOST:   Flag whether boost to CMS was or was 
C     .                            not performed.
C.
C.    Output Parameters:  Common /PHOEVT/, 
C.
C.    Author(s):  Z. Was                          Created at:  24/05/93
C.                                                Last Update:
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER IP,LL,FIRST,LAST,I
      LOGICAL BOOST
      INTEGER NN,J,K,NHEP0,NA
      DOUBLE PRECISION BET(3),GAM,PB
      COMMON /PHOCMS/ BET,GAM
      IF(NPHO.EQ.NEVPHO) RETURN
C--   When parent was not in its rest-frame, boost back...
      CALL PHLUPA(10)
      IF (BOOST) THEN
        DO 110 J=JDAPHO(1,1),JDAPHO(2,1)
          PB=-BET(1)*PPHO(1,J)-BET(2)*PPHO(2,J)-BET(3)*PPHO(3,J)
          DO 100 K=1,3
  100     PPHO(K,J)=PPHO(K,J)-BET(K)*(PPHO(4,J)+PB/(GAM+1.D0))
  110   PPHO(4,J)=GAM*PPHO(4,J)+PB
C--   ...boost photon, or whatever else has shown up
        DO NN=NEVPHO+1,NPHO
          PB=-BET(1)*PPHO(1,NN)-BET(2)*PPHO(2,NN)-BET(3)*PPHO(3,NN)
          DO 120 K=1,3
  120     PPHO(K,NN)=PPHO(K,NN)-BET(K)*(PPHO(4,NN)+PB/(GAM+1.D0))
          PPHO(4,NN)=GAM*PPHO(4,NN)+PB
        ENDDO
      ENDIF
        CALL PHCORK(0)   ! we have to use it because it clears input 
                         ! for grandaughters modified in C++
        FIRST=JDAHEP(1,IP)
        LAST =JDAHEP(2,IP)
C let-s take in original daughters
        DO LL=0,LAST-FIRST
         IDHEP(FIRST+LL) = IDPHO(3+LL)
           DO I=1,5
             PHEP(I,FIRST+LL) = PPHO(I,3+LL)
           ENDDO
        ENDDO
C let-s take newcomers to the end of HEPEVT.
        NA=3+LAST-FIRST
        DO LL=1,NPHO-NA
         IDHEP(NHEP0+LL) = IDPHO(NA+LL)
         ISTHEP(NHEP0+LL)=ISTPHO(NA+LL)
         JMOHEP(1,NHEP0+LL)=IP
         JMOHEP(2,NHEP0+LL)=JMOHEP(2,JDAHEP(1,IP))
         JDAHEP(1,NHEP0+LL)=0
         JDAHEP(2,NHEP0+LL)=0
           DO I=1,5
             PHEP(I,NHEP0+LL) = PPHO(I,NA+LL)
           ENDDO
        ENDDO
        NHEP=NHEP+NPHO-NEVPHO
        CALL PHLUPA(20)
      END 
      SUBROUTINE PHOCHK(JFIRST)
C.----------------------------------------------------------------------
C.
C.    PHOCHK:   checking branch.
C.
C.    Purpose:  checks whether particles in the common block /PHOEVT/
C.              can be served by PHOMAK. 
C.              JFIRST is the position in /PH_HEPEVT/ (!) of the first 
C.              daughter of sub-branch under action.
C.
C.
C.    Author(s):  Z. Was                           Created at: 22/10/92
C.                                                Last Update: 11/12/00
C.
C.----------------------------------------------------------------------
C     ********************
      IMPLICIT NONE
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(NMXPHO)
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      LOGICAL QEDRAD
      COMMON/PH_PHOQED/QEDRAD(NMXHEP)
      INTEGER JFIRST
      LOGICAL F
      INTEGER IDABS,NLAST,I,IPPAR
      LOGICAL INTERF,ISEC,ITRE,IEXP,IFTOP,IFW,IFNPI0,IFKL
      REAL*8 FINT,FSEC,EXPEPS
      COMMON /PHOKEY/ FSEC,FINT,EXPEPS,INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      LOGICAL IFRAD
      INTEGER IDENT,K,IQRK,IPHQRK,IEKL,IPHEKL
      NLAST = NPHO
C
      IPPAR=1

      DO 10 I=IPPAR,NLAST
      IDABS    = ABS(IDPHO(I))
C possibly call on PHZODE is a dead (to be omitted) code. 
      CHKIF(I)= F(IDABS)       .AND.F(ABS(IDPHO(1)))
     &  .AND.   (IDPHO(2).EQ.0)
      IF(I.GT.2) CHKIF(I)=CHKIF(I).AND.QEDRAD(JFIRST+I-IPPAR-2)
 10   CONTINUE
C--
C now we go to special cases, where CHKIF(I) will be overwritten
C--
      IF(IFTOP) THEN
C special case of top pair production
        DO  K=JDAPHO(2,1),JDAPHO(1,1),-1
           IF(IDPHO(K).NE.22) THEN
             IDENT=K
             GOTO 15
           ENDIF
        ENDDO
 15     CONTINUE
        IFRAD=((IDPHO(1).EQ.21).AND.(IDPHO(2).EQ.21))
     &  .OR. ((ABS(IDPHO(1)).LE.6).AND.((IDPHO(2)).EQ.(-IDPHO(1))))
        IFRAD=IFRAD
     &        .AND.(ABS(IDPHO(3)).EQ.6).AND.((IDPHO(4)).EQ.(-IDPHO(3)))
     &        .AND.(IDENT.EQ.4)   
        IF(IFRAD) THEN    
           DO 20 I=IPPAR,NLAST
           CHKIF(I)= .TRUE.
           IF(I.GT.2) CHKIF(I)=CHKIF(I).AND.QEDRAD(JFIRST+I-IPPAR-2)
 20        CONTINUE
        ENDIF
      ENDIF
C--
C--
      IF(IFTOP) THEN
C special case of top decay
        DO  K=JDAPHO(2,1),JDAPHO(1,1),-1
           IF(IDPHO(K).NE.22) THEN
             IDENT=K
             GOTO 25
           ENDIF
        ENDDO
 25     CONTINUE
        IFRAD=((ABS(IDPHO(1)).EQ.6).AND.(IDPHO(2).EQ.0))
        IFRAD=IFRAD
     &        .AND.((ABS(IDPHO(3)).EQ.24).AND.(ABS(IDPHO(4)).EQ.5)
     &        .OR.(ABS(IDPHO(3)).EQ.5).AND.(ABS(IDPHO(4)).EQ.24))
     &        .AND.(IDENT.EQ.4)   
        IF(IFRAD) THEN    
           DO 30 I=IPPAR,NLAST
           CHKIF(I)= .TRUE.
           IF(I.GT.2) CHKIF(I)=CHKIF(I).AND.QEDRAD(JFIRST+I-IPPAR-2)
 30        CONTINUE
        ENDIF
      ENDIF
C--
C--
      END
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

      SUBROUTINE PHOPRE(IPARR,WT,NEUDAU,NCHARB)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   Photon radiation in decays
C.
C.    Purpose:  Order (alpha) radiative corrections  are  generated  in
C.              the decay of the IPPAR-th particle in the HEP-like
C.              common /PHOEVT/.  Photon radiation takes place from one
C.              of the charged daughters of the decaying particle IPPAR
C.              WT is calculated, eventual rejection will be performed
C.              later after inclusion of interference weight.
C.
C.    Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
C.                                /PHOEVT/ and the common itself,
C.
C.    Output Parameters:  Common  /PHOEVT/, either  with  or  without a
C.                                photon(s) added.
C.                        WT      weight of the configuration 
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
C.                                                Last Update: 29/01/05
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION MINMAS,MPASQR,MCHREN
      DOUBLE PRECISION BETA,EPS,DEL1,DEL2,DATA,BIGLOG
      REAL*8 PHOCHA,PHOSPI,PHORANC,PHOCOR,PHOCORN,MASSUM
      INTEGER IP,IPARR,IPPAR,I,J,ME,NCHARG,NEUPOI,NLAST,THEDUM
      INTEGER IDABS,IDUM
      INTEGER NCHARB,NEUDAU
      REAL*8 WT,WGT
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(NMXPHO)
      INTEGER CHAPOI(NMXPHO)
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL*8 PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL*8 XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL*8 ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      INTEGER IREP,IDME
      REAL*8 PROBH,CORWT,XF
      COMMON/PHOPRO/PROBH,CORWT,XF,IREP
C may be it is not the best place, but ...
      LOGICAL INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      REAL*8 FINT,FSEC,EXPEPS
      COMMON /PHOKEY/ FSEC,FINT,EXPEPS,INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      REAL*8             WT1,WT2,WT3
      COMMON /PHWT/ BETA,WT1,WT2,WT3
      DOUBLE PRECISION phocorWT3,phocorWT2,phocorWT1
      common/phocorWT/phocorWT3,phocorWT2,phocorWT1
      real*8 a,b
C--
      IPPAR=IPARR
C--   Store pointers for cascade treatement...
      IP=IPPAR
      NLAST=NPHO
      IDUM=1
C--
C--   Check decay multiplicity..
      IF (JDAPHO(1,IP).EQ.0) RETURN
C--
C--   Loop over daughters, determine charge multiplicity
   10 NCHARG=0
      IREP=0
      MINMAS=0.D0
      MASSUM=0.D0
      DO 20 I=JDAPHO(1,IP),JDAPHO(2,IP)
C--
C--
C--   Exclude marked particles, quarks and gluons etc...
        IDABS=ABS(IDPHO(I))
        IF (CHKIF(I-JDAPHO(1,IP)+3)) THEN
          IF (PHOCHA(IDPHO(I)).NE.0) THEN
            NCHARG=NCHARG+1
            IF (NCHARG.GT.NMXPHO) THEN
              DATA=NCHARG
              CALL PHOERR(1,'PHOTOS',DATA)
            ENDIF
            CHAPOI(NCHARG)=I
          ENDIF
          MINMAS=MINMAS+PPHO(5,I)**2
        ENDIF
        MASSUM=MASSUM+PPHO(5,I)
   20 CONTINUE
      IF (NCHARG.NE.0) THEN
C--
C--   Check that sum of daughter masses does not exceed parent mass
        IF ((PPHO(5,IP)-MASSUM)/PPHO(5,IP).GT.2.D0*XPHCUT) THEN
C--
   30       CONTINUE
            DO 70 J=1,3
   70       PNEUTR(J)=-PPHO(J,CHAPOI(NCHARG))
            PNEUTR(4)=PPHO(5,IP)-PPHO(4,CHAPOI(NCHARG))
C--
C--   Calculate  invariant  mass of 'neutral' etc. systems
          MPASQR=PPHO(5,IP)**2
          MCHSQR=PPHO(5,CHAPOI(NCHARG))**2
          IF ((JDAPHO(2,IP)-JDAPHO(1,IP)).EQ.1) THEN
            NEUPOI=JDAPHO(1,IP)
            IF (NEUPOI.EQ.CHAPOI(NCHARG)) NEUPOI=JDAPHO(2,IP)
            MNESQR=PPHO(5,NEUPOI)**2
            PNEUTR(5)=PPHO(5,NEUPOI)
          ELSE
            MNESQR=PNEUTR(4)**2-PNEUTR(1)**2-PNEUTR(2)**2-PNEUTR(3)**2
            MNESQR=MAX(MNESQR,MINMAS-MCHSQR)
            PNEUTR(5)=SQRT(MNESQR)
          ENDIF
C--
C--   Determine kinematical limit...
          XPHMAX=(MPASQR-(PNEUTR(5)+PPHO(5,CHAPOI(NCHARG)))**2)/MPASQR
C--
C--   Photon energy fraction...
          CALL PHOENE(MPASQR,MCHREN,BETA,BIGLOG,IDPHO(CHAPOI(NCHARG)))
C--
         IF (XPHOTO.LT.-4D0) THEN
            NCHARG=0  ! we really stop trials
            XPHOTO=0d0! in this case !!
C--   Energy fraction not too large (very seldom) ? Define angle.
          ELSEIF ((XPHOTO.LT.XPHCUT).OR.(XPHOTO.GT.XPHMAX)) THEN
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
            EPS=MCHREN/(1.D0+BETA)
C--
C--   Calculate sin(theta) and cos(theta) from interval variables
            DEL1=(2.D0-EPS)*(EPS/(2.D0-EPS))**PHORANC(THEDUM)
            DEL2=2.D0-DEL1

C ----------- VARIANT B ------------------
CC corrections for more efiicient interference correction,
CC instead of doubling crude distribution, we add flat parallel channel
C           IF (PHORANC(THEDUM).LT.BIGLOG/BETA/(BIGLOG/BETA+2*FINT)) THEN
C              COSTHG=(1.D0-DEL1)/BETA
C              SINTHG=SQRT(DEL1*DEL2-MCHREN)/BETA
C           ELSE
C             COSTHG=-1D0+2*PHORANC(THEDUM)
C             SINTHG= SQRT(1D0-COSTHG**2)
C           ENDIF
C
C           IF (FINT.GT.1.0D0) THEN
C
C              WGT=1D0/(1D0-BETA*COSTHG)
C              WGT=WGT/(WGT+FINT)
C       !       WGT=1D0   ! ??
C
C           ELSE
C              WGT=1D0
C           ENDIF
C
C ----------- END OF VARIANT B ------------------

C ----------- VARIANT A ------------------
              COSTHG=(1.D0-DEL1)/BETA
              SINTHG=SQRT(DEL1*DEL2-MCHREN)/BETA
              WGT=1D0
C ----------- END OF VARIANT A ------------------

C--
C--   Determine spin of  particle and construct code  for matrix element
            ME=2.D0*PHOSPI(IDPHO(CHAPOI(NCHARG)))+1.D0
C--
C--   Weighting procedure with 'exact' matrix element, reconstruct kine-
C--   matics for photon, neutral and charged system and update /PHOEVT/.
C--   Find pointer to the first component of 'neutral' system
      DO  I=JDAPHO(1,IP),JDAPHO(2,IP)
        IF (I.NE.CHAPOI(NCHARG)) THEN
          NEUDAU=I
          GOTO 51
        ENDIF
      ENDDO
C--
C--   Pointer not found...
      DATA=NCHARG
      CALL PHOERR(5,'PHOKIN',DATA)
 51   CONTINUE
      NCHARB=CHAPOI(NCHARG)
      NCHARB=NCHARB-JDAPHO(1,IP)+3
      NEUDAU=NEUDAU-JDAPHO(1,IP)+3

      CALL ME_CHANNEL(IDME)  !  two options introduced temporarily. 
                             !  In future always PHOCOR-->PHOCORN
                             !  Tests and adjustment of wts for Znlo needed.
                             !  otherwise simple change. PHOCORN implements
                             !  exact ME for scalar to 2 scalar decays.
      IF(IDME.EQ.2) THEN
        b=PHOCORN(MPASQR,MCHREN,ME)
        WT=b*WGT
        WT=WT/(1-xphoto/xphmax+0.5*(xphoto/xphmax)**2)*(1-xphoto/xphmax)/2 ! factor to go to WnloWT
      ELSEIF(IDME.EQ.1) THEN

        a=PHOCOR(MPASQR,MCHREN,ME)
        b=PHOCORN(MPASQR,MCHREN,ME)
        WT=b*WGT 
        WT=WT*wt1*wt2*wt3/phocorwt1/phocorwt2/phocorwt3 ! factor to go to ZnloWT
!        write(*,*) ' -----------'
!        write(*,*)   wt1,' ',wt2,' ',wt3
!        write(*,*)   phocorwt1,' ',phocorwt2,' ',phocorwt3
      ELSE
        a=PHOCOR(MPASQR,MCHREN,ME)
        WT=a*WGT
!        WT=b*WGT!/(1-xphoto/xphmax+0.5*(xphoto/xphmax)**2)*(1-xphoto/xphmax)/2
      ENDIF


          ENDIF
        ELSE
          DATA=PPHO(5,IP)-MASSUM
          CALL PHOERR(10,'PHOTOS',DATA)
        ENDIF
      ENDIF
C--
      RETURN
      END
      SUBROUTINE PHOENE(MPASQR,MCHREN,BETA,BIGLOG,IDENT)
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
C.                B. van Eijk, P.Golonka          Last Update: 29/01/05
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION MPASQR,MCHREN,BIGLOG,BETA,DATA
      INTEGER IWT1,IRN,IWT2
      REAL*8 PRSOFT,PRHARD,PHORANC,PHOFAC
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL*8 PNEUTR
      INTEGER IDENT
      REAL*8 PHOCHA,PRKILL,RRR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL*8 XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL*8 ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      REAL*8 PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      INTEGER IREP
      REAL*8 PROBH,CORWT,XF
      COMMON/PHOPRO/PROBH,CORWT,XF,IREP
      LOGICAL INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      REAL*8 FINT,FSEC,EXPEPS
      COMMON /PHOKEY/ FSEC,FINT,EXPEPS,INTERF,ISEC,ITRE,IEXP,IFTOP,IFW
      INTEGER NX,NCHAN,K,IDME
      PARAMETER (NX=10)
      LOGICAL EXPINI
      REAL*8 PRO,PRSUM
      COMMON /PHOEXP/ PRO(NX),NCHAN,EXPINI
C--
      IF (XPHMAX.LE.XPHCUT) THEN
        BETA=PHOFAC(-1)  ! to zero counter, here beta is dummy
        XPHOTO=0.0D0
        RETURN
      ENDIF
C--   Probabilities for hard and soft bremstrahlung...
      MCHREN=4.D0*MCHSQR/MPASQR/(1.D0+MCHSQR/MPASQR)**2
      BETA=SQRT(1.D0-MCHREN)

C ----------- VARIANT B ------------------
CC we replace 1D0/BETA*BIGLOG with (1D0/BETA*BIGLOG+2*FINT) 
CC for integral of new crude
C      BIGLOG=LOG(MPASQR/MCHSQR*(1.D0+BETA)**2/4.D0*
C     &          (1.D0+MCHSQR/MPASQR)**2)
C      PRHARD=ALPHA/PI*(1D0/BETA*BIGLOG+2*FINT)*(LOG(XPHMAX/XPHCUT)
C     &-.75D0+XPHCUT/XPHMAX-.25D0*XPHCUT**2/XPHMAX**2)
C      PRHARD=PRHARD*PHOCHA(IDENT)**2*FSEC
C ----------- END OF VARIANT B ------------------

C ----------- VARIANT A ------------------
      BIGLOG=LOG(MPASQR/MCHSQR*(1.D0+BETA)**2/4.D0*
     &          (1.D0+MCHSQR/MPASQR)**2)
      PRHARD=ALPHA/PI*(1D0/BETA*BIGLOG)*
     &(LOG(XPHMAX/XPHCUT)-.75D0+XPHCUT/XPHMAX-.25D0*XPHCUT**2/XPHMAX**2)
      PRHARD=PRHARD*PHOCHA(IDENT)**2*FSEC*FINT
        CALL ME_CHANNEL(IDME)
!        write(*,*) 'KANALIK IDME=',IDME
        IF (IDME.EQ.0) THEN  ! default
           continue
        ELSEIF (IDME.EQ.1) THEN
           PRHARD=PRHARD/(1d0+0.75*ALPHA/PI) !  NLO
        ELSEIF (IDME.EQ.2) THEN
           continue ! work on virtual crrections in W decay to be done.
        ELSE
         write(*,*) 'problem with ME_CHANNEL  IDME=',IDME
         stop
        ENDIF

C ----------- END OF VARIANT A ------------------
      IF (IREP.EQ.0) PROBH=0.D0
      PRKILL=0d0
      IF (IEXP) THEN           ! IEXP
       NCHAN=NCHAN+1
       IF (EXPINI) THEN     ! EXPINI
          PRO(NCHAN)=PRHARD+0.05*(1.0+FINT) ! we store hard photon emission prob 
                                      !for leg NCHAN
          PRHARD=0D0         ! to kill emission at initialization call
          PROBH=PRHARD
       ELSE                 ! EXPINI
        PRSUM=0
        DO K=NCHAN,NX
         PRSUM=PRSUM+PRO(K)
        ENDDO
        PRHARD=PRHARD/PRSUM ! note that PRHARD may be smaller than 
                            !PRO(NCHAN) because it is calculated
                            ! for kinematical configuartion as is 
                            ! (with effects of previous photons)
        PRKILL=PRO(NCHAN)/PRSUM-PRHARD !

       ENDIF                ! EXPINI
        PRSOFT=1.D0-PRHARD
      ELSE                     ! IEXP
       PRHARD=PRHARD*PHOFAC(0) ! PHOFAC is used to control eikonal 
                               ! formfactors for non exp version only
                               ! here PHOFAC(0)=1 at least now.
       PROBH=PRHARD
      ENDIF                    ! IEXP
      PRSOFT=1.D0-PRHARD
C--
C--   Check on kinematical bounds
      IF (IEXP) THEN
       IF (PRSOFT.LT.-5.0D-8) THEN
         DATA=PRSOFT
         CALL PHOERR(2,'PHOENE',DATA)
       ENDIF
      ELSE
       IF (PRSOFT.LT.0.1D0) THEN
         DATA=PRSOFT
         CALL PHOERR(2,'PHOENE',DATA)
       ENDIF
      ENDIF

      RRR=PHORANC(IWT1)
      IF (RRR.LT.PRSOFT) THEN
C--
C--   No photon... (ie. photon too soft)
        XPHOTO=0.D0
        IF (RRR.LT.PRKILL) XPHOTO=-5d0 ! No photon...no further trials
      ELSE
C--
C--   Hard  photon... (ie.  photon  hard enough).
C--   Calculate  Altarelli-Parisi Kernel
   10   XPHOTO=EXP(PHORANC(IRN)*LOG(XPHCUT/XPHMAX))
        XPHOTO=XPHOTO*XPHMAX
        IF (PHORANC(IWT2).GT.((1.D0+(1.D0-XPHOTO/XPHMAX)**2)/2.D0)) 
     &                            GOTO 10
      ENDIF
C--
C--   Calculate parameter for PHOFAC function
      XF=4.D0*MCHSQR*MPASQR/(MPASQR+MCHSQR-MNESQR)**2
      RETURN
      END

      FUNCTION PHOCORN(MPASQR,MCHREN,ME)
c ###### 
c  replace with, 
c ######

C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays CORrection weight from
C.              matrix elements This version for spin 1/2 is verified for
C.              W decay only
C.    Purpose:  Calculate  photon  angle.  The reshaping functions  will
C.              have  to  depend  on the spin S of the charged particle.
C.              We define:  ME = 2 * S + 1 !
C.              THIS IS POSSIBLY ALWAYS BETTER THAN PHOCOR()
C.
C.    Input Parameters:  MPASQR:  Parent mass squared,
C.                       MCHREN:  Renormalised mass of charged system,
C.                       ME:      2 * spin + 1 determines matrix element
C.
C.    Output Parameter:  Function value.
C.
C.    Author(s):  Z. Was, B. van Eijk, G. Nanava  Created at:  26/11/89
C.                                                Last Update: 01/11/12
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION MPASQR,MCHREN,BETA,BETA0,BETA1,XX,YY,DATA
      INTEGER ME
      REAL*8 PHOCOR,PHOFAC,WT1,WT2,WT3,PHOTRI,S1,PHOCORN
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL*8 PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG,phocorWT3,phocorWT2,phocorWT1
      REAL*8 XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      common/phocorWT/phocorWT3,phocorWT2,phocorWT1
      INTEGER IREP
      REAL*8 PROBH,CORWT,XF
      COMMON/PHOPRO/PROBH,CORWT,XF,IREP

C--
C--   Shaping (modified by ZW)...
      XX=4.D0*MCHSQR/MPASQR*(1.D0-XPHOTO)/(1.D0-XPHOTO+(MCHSQR-MNESQR)/
     &MPASQR)**2
      IF (ME.EQ.1) THEN
        S1=MPASQR  * (1.D0-XPHOTO)
        BETA0=2*PHOTRI(1D0,dsqrt(MCHSQR/MPASQR),dsqrt(MNESQR/MPASQR))
        BETA1=2*PHOTRI(1D0,dsqrt(MCHSQR/S1),dsqrt(MNESQR/S1))
        WT1= (1.D0-COSTHG*SQRT(1.D0-MCHREN))
     $      /((1D0+(1D0-XPHOTO/XPHMAX)**2)/2.D0)*XPHOTO          ! de-presampler
     $     
        WT2= beta1/beta0*XPHOTO                                  !phase space jacobians
        WT3=  beta1**2* (1D0-COSTHG**2) * (1D0-XPHOTO)/XPHOTO**2 ! matrix element
     $    /(1D0 +MCHSQR/S1-MNESQR/S1-BETA1*COSTHG)**2/2D0 
      ELSEIF (ME.EQ.2) THEN
        S1=MPASQR  * (1.D0-XPHOTO)
        BETA0=2*PHOTRI(1D0,dsqrt(MCHSQR/MPASQR),dsqrt(MNESQR/MPASQR))
        BETA1=2*PHOTRI(1D0,dsqrt(MCHSQR/S1),dsqrt(MNESQR/S1))
        WT1= (1.D0-COSTHG*SQRT(1.D0-MCHREN))
     $      /((1D0+(1D0-XPHOTO/XPHMAX)**2)/2.D0)*XPHOTO          ! de-presampler
         
        WT2= beta1/beta0*XPHOTO                                  ! phase space jacobians

        WT3= beta1**2* (1D0-COSTHG**2) * (1D0-XPHOTO)/XPHOTO**2  ! matrix element
     $       /(1D0 +MCHSQR/S1-MNESQR/S1-BETA1*COSTHG)**2/2D0 
        WT3=WT3*(1-xphoto/xphmax+0.5*(xphoto/xphmax)**2)/(1-xphoto/xphmax)
c       print*,"WT3=",wt3
        phocorWT3=WT3
        phocorWT2=WT2
        phocorWT1=WT1

c       YY=0.5D0*(1.D0-XPHOTO/XPHMAX+1.D0/(1.D0-XPHOTO/XPHMAX))
c       BETA=SQRT(1.D0-XX)
c       WT1=(1.D0-COSTHG*SQRT(1.D0-MCHREN))/(1.D0-COSTHG*BETA)
c       WT2=(1.D0-XX/YY/(1.D0-BETA**2*COSTHG**2))*(1.D0+COSTHG*BETA)/2.D0
c       WT3=1.D0
      ELSEIF ((ME.EQ.3).OR.(ME.EQ.4).OR.(ME.EQ.5)) THEN
        YY=1.D0
        BETA=SQRT(1.D0-XX)
        WT1=(1.D0-COSTHG*SQRT(1.D0-MCHREN))/(1.D0-COSTHG*BETA)
        WT2=(1.D0-XX/YY/(1.D0-BETA**2*COSTHG**2))*(1.D0+COSTHG*BETA)/2.D0
        WT3=(1.D0+(1.D0-XPHOTO/XPHMAX)**2-(XPHOTO/XPHMAX)**3)/
     &  (1.D0+(1.D0-XPHOTO/XPHMAX)** 2)
      ELSE
        DATA=(ME-1.D0)/2.D0
        CALL PHOERR(6,'PHOCOR',DATA)
        YY=1.D0
        BETA=SQRT(1.D0-XX)
        WT1=(1.D0-COSTHG*SQRT(1.D0-MCHREN))/(1.D0-COSTHG*BETA)
        WT2=(1.D0-XX/YY/(1.D0-BETA**2*COSTHG**2))*(1.D0+COSTHG*BETA)/2.D0
        WT3=1.D0
      ENDIF
      WT2=WT2*PHOFAC(1)
      PHOCOR=WT1*WT2*WT3
      PHOCORN=PHOCOR
      CORWT=PHOCOR
      IF (PHOCOR.GT.1.D0) THEN
        DATA=PHOCOR
        CALL PHOERR(3,'PHOCOR',DATA)
      ENDIF
      RETURN
      END

