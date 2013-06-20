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



      DOUBLE PRECISION FUNCTION Zphwtnlo
     $       (svar,xk,IDHEP3,IREP,IBREM,qp,qm,ph,pp,pm,COSTHG,BETA,th1)
      IMPLICIT NONE
      REAL*8  QP(4),QM(4),PH(4),PP(4),PM(4)
      INTEGER IDHEP3,IREP,IBREM
      INTEGER IDE,IDF
      DOUBLE PRECISION svar,xk
      DOUBLE PRECISION COSTHG,BETA,th1 ! variables of some common blocks
                                       ! passed though function parameters
      DOUBLE PRECISION C,s,xkaM,xkaP,t,u,t1,u1,BT,BU
      DOUBLE PRECISION waga,wagan2
      DOUBLE PRECISION PHASYZ


C IBREM is spurious but it numbers branches of MUSTRAAL

        IBREM=1
       IF (IREP.EQ.1)  IBREM=-1

C we calculate C and S, note that TH1 exists in MUSTRAAL as well. 

        C=COS(TH1) ! this parameter is calculated outside of the class

C from off line application we had:
        IF(IBREM.EQ.-1) C=-C
C ... we need to check it. 
        s=sqrt(1D0-C**2)

        IF (IBREM.EQ.1) THEN
         xkaM=(qp(4)*PH(4)-qp(3)*PH(3)-qp(2)*PH(2)-qp(1)*PH(1))/XK
         xkaP=(qM(4)*PH(4)-qM(3)*PH(3)-qM(2)*PH(2)-qM(1)*PH(1))/XK
        ELSE
         xkap=(qp(4)*PH(4)-qp(3)*PH(3)-qp(2)*PH(2)-qp(1)*PH(1))/XK
         xkaM=(qM(4)*PH(4)-qM(3)*PH(3)-qM(2)*PH(2)-qM(1)*PH(1))/XK
        ENDIF    

!        XK=2*PHEP(4,nhep)/PHEP(4,1)/xphmax   ! it is not used becuse here
                                              ! order of emissions is meaningless
!
!        DELTA=2*PHEP(5,4)**2/svar/(1+(1-XK)**2)*(xKAP/xKAM+xKAM/xKAP)
!        waga=SVAR/4./xkap
!        waga=waga*(1.D0-COSTHG*BETA) ! sprawdzone 1= svar/xKAp/4   * (1.D0-COSTHG*BETA)
!        waga=waga*(1-delta) /wt2 ! sprawdzone ze to jest =2/(1.D0+COSTHG*BETA)
!                                 ! czyli ubija de-interferencje
 

C this is true only for intermediate resonances with afb=0!
        t =2*(qp(4)*PP(4)-qp(3)*PP(3)-qp(2)*PP(2)-qp(1)*PP(1))
        u =2*(qM(4)*PP(4)-qM(3)*PP(3)-qM(2)*PP(2)-qM(1)*PP(1))
        u1=2*(qp(4)*Pm(4)-qp(3)*Pm(3)-qp(2)*Pm(2)-qp(1)*Pm(1))
        t1=2*(qM(4)*Pm(4)-qM(3)*Pm(3)-qM(2)*Pm(2)-qM(1)*Pm(1))

C basically irrelevant lines  ...
        t =t - (qp(4)**2-qp(3)**2-qp(2)**2-qp(1)**2)
        u =u - (qm(4)**2-qm(3)**2-qm(2)**2-qm(1)**2)
        u1=u1- (qp(4)**2-qp(3)**2-qp(2)**2-qp(1)**2)
        t1=t1- (qm(4)**2-qm(3)**2-qm(2)**2-qm(1)**2)




      call GETIDEIDF(IDE,IDF)   ! we adjust to what is f-st,s-nd beam flavour 
       IF (IDE*IDHEP3.GT.0) THEN
        BT=1+PHASYZ(SVAR)
        BU=1-PHASYZ(SVAR)
       ELSE
        BT=1-PHASYZ(SVAR)
        BU=1+PHASYZ(SVAR)
       ENDIF  
        wagan2=2*(BT*t**2+BU*u**2+BT*t1**2+BU*u1**2)
     $        /(1+(1-xk)**2)* 2.0/(BT*(1-c)**2+BU*(1+c)**2)/svar**2

!        waga=waga*wagan2
!        waga=waga*(1-delta) /wt2 ! sprawdzone ze to jest =2/(1.D0+COSTHG*BETA)
        waga=2/(1.D0+COSTHG*BETA)*wagan2  
!     %       * SVAR/4./xkap*(1.D0-COSTHG*BETA)*sqrt(1.0-xk)

        Zphwtnlo=waga
        IF(wagan2.gt.3.8) THEN
!         write(*,*) 'phwtnlo= ',phwtnlo
!         write(*,*) 'idhepy= ',IDHEP(1),IDHEP(2),IDHEP(3),IDHEP(4),IDHEP(5)
         write(*,*) 'IDE=',IDE,'  IDF=',IDF
         write(*,*) 'bt,bu,bt+bu= ',bt,bu,bt+bu
         call PHODMP
         write(*,*) ' '
         write (*,*) IREP,IBREM, '<-- IREP,IBREM '
 !        write(*,*) 'pneutr= ',pneutr
         write(*,*) 'qp    = ',qp
         write(*,*) 'qm    = ',qm
         write(*,*) ' '
         write(*,*) 'ph    = ',ph
 !        write(*,*) 'p1= ',PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1)
 !        write(*,*) 'p2= ',PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2)
 !        write(*,*) 'p3= ',PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3)
 !        write(*,*) 'p4= ',PHEP(1,4),PHEP(2,4),PHEP(3,4),PHEP(4,4)
 !        write(*,*) 'p5= ',PHEP(1,5),PHEP(2,5),PHEP(3,5),PHEP(4,5)

         write (*,*) ' c= ',c,' theta= ',th1
!         write(*,*)  'photos waga daje ... IBREM=',IBREM,' waga=',waga
!         write(*,*) 'xk,COSTHG,c',xk,COSTHG,c
!         write(*,*) SVAR/4./xkap*(1.D0-COSTHG*BETA), 
!     $   (1-delta) /wt2 *(1.D0+COSTHG*BETA)/2, wagan2
!         write(*,*) ' delta, wt2,beta',  delta, wt2,beta
         write(*,*) '  -  '
         write(*,*) 't,u       = ',t,u
         write(*,*) 't1,u1     = ',t1,u1
         write(*,*) 'sredniaki = ',svar*(1-c)/2
     $                            ,svar*(1+c)/2
!         write(*,*) 'xk= ',xk,' c= ',c ,' COSTHG= ' ,COSTHG
         write(*,*) 'PHASYZ(SVAR)=',PHASYZ(SVAR),' svar= ',svar,' waga=',waga
         write(*,*) '  -  '
         write(*,*) 'BT-part= ',2*(BT*t**2+BT*t1**2)
     $        /(1+(1-xk)**2)* 2.0/(BT*(1-c)**2)/svar**2,
     $ ' BU-part= ',2*(BU*u**2+BU*u1**2)
     $        /(1+(1-xk)**2)* 2.0/(BU*(1+c)**2)/svar**2
         write(*,*) 'BT-part*BU-part= ',2*(BT*t**2+BT*t1**2)
     $        /(1+(1-xk)**2)* 2.0/(BT*(1-c)**2)/svar**2
     $         *2*(BU*u**2+BU*u1**2)
     $        /(1+(1-xk)**2)* 2.0/(BU*(1+c)**2)/svar**2, 'wagan2=',wagan2

         write(*,*) 'wagan2=',wagan2
         write(*,*) ' ###################  '
         wagan2=3.8  ! overwrite 
         waga=2/(1.D0+COSTHG*BETA)*wagan2  
!     %       * SVAR/4./xkap*(1.D0-COSTHG*BETA)*sqrt(1.0-xk)

        Zphwtnlo=waga

        ENDIF


      RETURN
      END
