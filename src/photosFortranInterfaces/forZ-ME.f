      function phwtnlo(xdumm) 
C.----------------------------------------------------------------------
C.
C.    PHWTNLO:   PHotosWTatNLO
C.
C.    Purpose:  calculates instead of interference weight
C.              complete NLO weight for vector boson decays
C.              of pure vector (or pseudovector) couplings
C.              Proper orientation of beams required.
C.              Uses Zphwtnlo encapsulating actual matrix element
C.              At NLO more tuning on kinematical conf.
C.              than in standard is needed.
C.              Works with KORALZ and KKMC. 
C.              Note some commented out commons from MUSTAAL, KORALZ
C.
C.    Input Parameters:   Common /PHOEVT/ /PHOPS/ /PHOREST/ /PHOPRO/
C.
C.    Output Parameters:  Function value
C.
C.    Author(s):  Z. Was                          Created at:  08/12/05
C.                                                Last Update: 
C.
C.----------------------------------------------------------------------

      IMPLICIT NONE 
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      COMMON/PHOEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      real*8 phep,vhep
      INTEGER NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      SAVE  /PHOEVT/
C      COMMON / UTIL3 / xFIG,xFI,IT,IBREX
      DOUBLE PRECISION COSTHG,SINTHG
      REAL*8 XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
C fi3 orientation of photon, fi1,th1 orientation of neutral
      REAL*8 FI3,fi1,th1
      COMMON /PHOREST/ FI3,fi1,th1
      REAL*8 BETA,WT1,WT2,WT3
      COMMON /PHWT/ BETA,WT1,WT2,WT3
      REAL*8 PROBH,CORWT,XF
      COMMON/PHOPRO/PROBH,CORWT,XF,IREP
      REAL*8 PI
      DATA PI /3.141592653589793238462643D0/
      REAL*8  QP(4),QM(4),PH(4),QQ(4),PP(4),PM(4),QQS(4)
      REAL*8 s,c,svar,xkaM,xkaP,xk,phwtnlo,xdumm,PHINT
      REAL*8 ENE,a,t,u,t1,u1,wagan2,waga,PHASYZ,BT,BU,ENEB
      INTEGER IBREM,K,L,IREP,IDUM,IDHEP3
      integer icont,ide,idf
      REAL*8 delta
      REAL*8 PNEUTR,MCHSQR,MNESQR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION Zphwtnlo
*/////////////////////
!        call phlupa(299 500)


*/////////////////////
!        call phlupa(299 500)

        XK=2*PHEP(4,nhep)/PHEP(4,1)

!        XK=2*PHEP(4,nhep)/PHEP(4,1)/xphmax  ! it is not used becuse here
                                             ! order of emissions is meaningless
      IF(NHEP.LE.4) XK=0.0D0
C the mother must be Z or gamma*  !!!!
      
      IF (XK.GT.1d-10.AND.(IDHEP(1).EQ.22.OR.IDHEP(1).EQ.23)) THEN

!        write(*,*) 'nhep=',nhep
!      DO K=1,3 ENDDO
!      IF (K.EQ.1) IBREM= 1
!      IF (K.EQ.2) IBREM=-1
!      ICONT=ICONT+1
!      IBREM=IBREX        ! that will be input parameter.
!      IBREM=IBREY        ! that IS now   input parameter.

C We initialize twice 4-vectors, here and again later after boost 
C must be the same way. Important is how the reduction procedure will work.
C It seems at present that the beams must be translated to be back to back.
C this may be done after initialising, thus on 4-vectors.

        DO K=1,4
         PP(K)=PHEP(K,1)
         PM(K)=PHEP(K,2)
         QP(K)=PHEP(K,3)
         QM(K)=PHEP(K,4)
         PH(K)=PHEP(K,nhep)
         QQ(K)=0.0
         QQS(K)=QP(K)+QM(K)
        ENDDO


        PP(4)=(PHEP(4,1)+PHEP(4,2))/2
        PM(4)=(PHEP(4,1)+PHEP(4,2))/2
        PP(3)= PP(4)
        PM(3)=-PP(4)
        
         DO L=5,NHEP-1
         DO K=1,4
          QQ(K)=QQ(K)+ PHEP(K,L)
         QQS(K)=QQS(K)+ PHEP(K,L)
         ENDDO
        ENDDO        

C go to the restframe of 3        
        CALL PHOB(1,QQS,qp)
        CALL PHOB(1,QQS,qm)
        CALL PHOB(1,QQS,QQ)
        ENE=(qp(4)+qm(4)+QQ(4))/2

C preserve direction of emitting particle and wipeout QQ 
        IF (IREP.EQ.1) THEN
         a=sqrt(ENE**2-PHEP(5,3)**2)/sqrt(qm(4)**2-PHEP(5,3)**2)
         qm(1)= qm(1)*a
         qm(2)= qm(2)*a
         qm(3)= qm(3)*a
         qp(1)=-qm(1)
         qp(2)=-qm(2)
         qp(3)=-qm(3)
        ELSE
         a=sqrt(ENE**2-PHEP(5,3)**2)/sqrt(qp(4)**2-PHEP(5,3)**2)
         qP(1)= qP(1)*a
         qP(2)= qP(2)*a
         qP(3)= qP(3)*a
         qm(1)=-qP(1)
         qm(2)=-qP(2)
         qm(3)=-qP(3)
        ENDIF
         qp(4)=ENE
         qm(4)=ENE
C go back to reaction frame (QQ eliminated) 
        CALL PHOB(-1,QQS,qp)
        CALL PHOB(-1,QQS,qm)
        CALL PHOB(-1,QQS,QQ)

        svar=PHEP(4,1)**2

        IDHEP3=IDHEP(3)
        phwtnlo=Zphwtnlo
     $        (svar,xk,IDHEP3,IREP,IBREM,qp,qm,ph,pp,pm,COSTHG,BETA,th1)

      ELSE
C in other cases we just use default setups.
        phwtnlo= PHINT(IDUM)
      ENDIF
      end
