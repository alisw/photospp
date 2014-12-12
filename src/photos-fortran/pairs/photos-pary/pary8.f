
C====================================================================
C Shelved project -- frozem. Basically code of the year 93, but
C some comments added in Nov 1998.
C implementation of electron positron pairs distributions are not OK
C that means never checked
C====================================================================
      SUBROUTINE PHOPAR(IPARR,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOPAR:   PHOtos PAiR radiation in decays
C.
C.    Purpose:  e+e- pairs  are  generated  in
C.              the decay of the IPPAR-th particle in the HEP-like
C.              common /PHOEVT/.  Radiation from one after another
C.              of the charged daughters of the decaying particle IPPAR
C.           
C.           
C.
C.    Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
C.                                /PHOEVT/ and the common itself,
C.                                NHEP0 length of the /HEPEVT/ entry
C.                                before starting any activity on this
C.                                IPPAR decay.
C.    Output Parameters:  Common  /PHOEVT/, either  with  or  without a
C.                                e+e-(s) added.
C.                    
C.
C.    Author(s):  Z. Was,                         Created at:  01/06/93
C.                                                Last Update: 
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION MINMAS
      DOUBLE PRECISION STRENG,PCHAR(4),PNEU(4),PELE(4),PPOZ(4),BUF(4)
      DOUBLE PRECISION DATA
      DOUBLE PRECISION PHOCHA,MASSUM
      INTEGER IP,IPARR,IPPAR,I,J,NCHARG,NLAST,K
      INTEGER IDABS,IDUM
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      DOUBLE PRECISION PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER PDMMAX
      PARAMETER (PDMMAX=2000)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(PDMMAX)
      INTEGER CHAPOI(PDMMAX)
      INTEGER IREP
      INTEGER NHEP0
      LOGICAL BOOST,JESLI
      INTEGER KEYRND
c-- type of RN generator
      COMMON / PHRNDP / KEYRND
      save   / PHRNDP /
C--
      KEYRND=1
      IPPAR=IPARR
C--   Store pointers for cascade treatement...
      IP=1
      NLAST=NPHO
      IDUM=1
C--
C--   Check decay multiplicity..
           CALL PHOIN(IPPAR,BOOST,NHEP0)
           CALL PHLUPA(100)
      IF (JDAPHO(1,IP).EQ.0) RETURN
      IF (JDAPHO(1,IP).EQ.JDAPHO(2,IP)) RETURN
C--
C--   Loop over daughters, determine charge multiplicity
   10 NCHARG=0
      IREP=0
      MINMAS=0.
      MASSUM=0.
      DO 20 I=JDAPHO(1,IP),JDAPHO(2,IP)
C--
C--
C--   Exclude marked particles, quarks and gluons etc...
        IDABS=ABS(IDPHO(I))
        IF (CHKIF(I-JDAPHO(1,IP)+3)) THEN
          IF (PHOCHA(IDPHO(I)).NE.0) THEN
            NCHARG=NCHARG+1
            IF (NCHARG.GT.PDMMAX) THEN
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
C--
C--   Order  charged  particles  according  to decreasing mass, this  to
C--   increase efficiency (smallest mass is treated first).
          IF (NCHARG.GT.1) CALL PHOOMA(1,NCHARG,CHAPOI)
C--
   30       CONTINUE
            DO 70 J=1,3
            PCHAR (J)= PPHO(J,CHAPOI(NCHARG))
   70       PNEU(J)=-PPHO(J,CHAPOI(NCHARG))
            PNEU(4)=PPHO(5,IP)-PPHO(4,CHAPOI(NCHARG))
            PCHAR (4)=PPHO(4,CHAPOI(NCHARG))
            STRENG=1D0
            CALL PHPAR(JESLI,STRENG,PCHAR,PNEU,PELE,PPOZ)
          IF (JESLI) THEN
           DO I=JDAPHO(1,IP),JDAPHO(2,IP)
            DO K=1,4
              BUF(K)=PPHO(K,I)
            ENDDO
            IF (I.EQ.CHAPOI(NCHARG)) THEN
               CALL PHPRTR( 1,BUF)
            ELSE
               CALL PHPRTR(-1,BUF)
            ENDIF
            DO K=1,4
              PPHO(K,I)=BUF(K)
            ENDDO 
           ENDDO 
c--  electron
           CALL PHLUPA(1011)
           NPHO=NPHO+1
           ISTPHO(NPHO)=1
           IDPHO(NPHO) =11
           JMOPHO(1,NPHO)=IP
           JMOPHO(2,NPHO)=0
           JDAPHO(1,NPHO)=0
           JDAPHO(2,NPHO)=0
           DO K=1,4
            PPHO(K,NPHO)=PELE(K)
           ENDDO
           PPHO(5,NPHO)=0.000 511
c--  positron
           NPHO=NPHO+1
           ISTPHO(NPHO)=1
           IDPHO(NPHO) =-11
           JMOPHO(1,NPHO)=IP
           JMOPHO(2,NPHO)=0
           JDAPHO(1,NPHO)=0
           JDAPHO(2,NPHO)=0
           DO K=1,4
            PPHO(K,NPHO)=PPOZ(K)
           ENDDO
           PPHO(5,NPHO)=0.000 511
C-- write in
           CALL PHLUPA(1012)
           CALL PHOOUT(IPPAR,BOOST,NHEP0) 
           CALL PHOIN(IPPAR,BOOST,NHEP0)
           CALL PHLUPA(1013)
          ENDIF
C--
C--   check  for more daughters  that may ra-
C--   diate and correct radiation probability...
C--   probability of pair emission is assumed to be negligible !!!
            NCHARG=NCHARG-1
            IF (NCHARG.GT.0) THEN
              IREP=IREP+1
              GOTO 30
            ENDIF
          ELSE
      ENDIF
C--
      RETURN
      END
      SUBROUTINE PHSPAJ(KUDA,PA,PB,PP,PE)     
C.----------------------------------------------------------------------
C.
C.    PHSPAJ:   PHotos SPying on four momenta  PA,PB,PP,PE
C.
C.    Purpose:  PHotos SPying on four momenta  PA,PB,PP,PE
C.              of the elementary branch with the pair creation
C.              Normally it is desactivated KLUCZ=0.
C.              For KLUCZ=1 printouts will show up.
C.
C.    Input Parameter:    KUDA:   position in program
C.                                PA,PB, charged neutral in generation
C.                                PP,PE pair (to be) generated
C.                    
C.
C.    Author(s):  Z. Was,                         Created at:  01/06/93
C.                                                Last Update: 
C.
C.----------------------------------------------------------------------
*     **********************     
      IMPLICIT NONE   
      REAL*8 SUM(4),PA(4),PB(4),PP(4),PE(4)
      INTEGER KLUCZ,K,NOUT,KUDA
      DATA KLUCZ /0/
      IF (KLUCZ.EQ.0) RETURN
      NOUT=56     
      WRITE(NOUT,*) KUDA,'===================PHSPAJ====================' 
      WRITE(NOUT,3100) ' P2',(PA(K),K=1,4)   
      WRITE(NOUT,3100) ' Q2',(PB(K),K=1,4)   
      WRITE(NOUT,3100) ' PE',(PE(K),K=1,4)   
      WRITE(NOUT,3100) ' PP',(PP(K),K=1,4)   
      DO 200 K=1,4      
  200 SUM(K)=PA(K)+PB(K)+PE(K)+PP(K)         
      WRITE(NOUT,3100) 'SUM',(SUM(K),K=1,4)           
      NOUT=16     
      WRITE(NOUT,*) KUDA,'===================PHSPAJ====================' 
      WRITE(NOUT,3100) ' P2',(PA(K),K=1,4)   
      WRITE(NOUT,3100) ' Q2',(PB(K),K=1,4)   
      WRITE(NOUT,3100) ' PE',(PE(K),K=1,4)   
      WRITE(NOUT,3100) ' PP',(PP(K),K=1,4)   
      WRITE(NOUT,3100) 'SUM',(SUM(K),K=1,4)           
 3100 FORMAT(1X,A3,1X,5F18.13)   
      END   
   
      SUBROUTINE PHPAR(JESLI,STRENG,PA,PB,PE,PP)       
C.----------------------------------------------------------------------
C.
C.    PHPAR:   PHotos PAiR radiation in decays
C.
C.    Purpose:  actual generation of e+e- pair
C.
C.    Input Parameter:    PA,PB  4-momenta of chosen as charged,neutral
C.                               objects before generation
C.    Output Parameters:  JESLI   flags telling whether pair is arriving
C.                        PE,PP   4-momenta of e+/e- generated
C.                        /PHKINP/ generated/calculated kinem var.
C.
C.    Author(s):  Z. Was,                         Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
!      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT NONE 
      REAL*8                                      
     $ FI0,FI1,FI2,FI3,FI4,FI5,TH0,TH1,TH3,TH4
     $,PARNEU,PARCH,BPAR,BSTA,BSTB
      COMMON  /PHKINP/ 
     $ FI0,FI1,FI2,FI3,FI4,FI5,TH0,TH1,TH3,TH4
     $,PARNEU,PARCH,BPAR,BSTA,BSTB
      REAL*8 PNEUTR(4),PAA(4),PHOT(4),PA(4),PB(4),PE(4),PP(4),PSUM(4)
      REAL*8 VEC(4)                                                
      REAL*8 RRR(8) 
      REAL*8 STRENG
      REAL*8  PI,AMEL,XK0,ALFINV,XLAM,A,B,C,X1,X2,FI0,TH0
      REAL*8  AMNE,AMCH,PRHARD,PHANGF,PHANXY,PHMAS,AMTO
      REAL*8  XMP,XP,C1,FIX1,C2,FIX2,WTA,XMK2,XPMAX,YOT3,YOT2,YOT1,WT
      REAL*8  AMCH2,AMNE2,AMTOST,QNEW,QOLD,GCHAR,GNEU,XK,XKM,XK0DEC,AXK
      REAL*8 PMOD,S2,PENE,PPED,COSTHG,SINTHG,PAIRB,GAMM
      INTEGER K
      LOGICAL JESLI,JESLIK                                                 
      DATA PI /3.141592653589793238462643D0/
      DATA AMEL,XK0 /.511D-3,1.0D-3/
      DATA ALFINV /137.01/                       
C                                                                  
      XLAM(A,B,C) = SQRT((A-B-C)**2-4.0*B*C)
C         
      PA(4)=MAX(PA(4),SQRT(PA(1)**2+PA(2)**2+PA(3)**2))
      PB(4)=MAX(PB(4),SQRT(PB(1)**2+PB(2)**2+PB(3)**2))
C 4-MOMENTUM OF THE NEUTRAL SYSTEM                                 
      DO 13 K=1,4
         PE(K)    =0D0
         PP(K)    =0D0
         PSUM(K)  =PA(K)+PB(K)
         PAA(K)   =PA(K)
 13      PNEUTR(K)=PB(K)
      IF ((PAA(4)+PNEUTR(4)). LT. .01D0 ) THEN
        JESLI=.FALSE.                                      
       RETURN
      ENDIF 
C MASSES OF THE NEUTRAL AND CHARGED SYSTEMS AND OVERALL MASS
C FIRST WE HAVE TO GO TO THE RESTFRAME TO GET RID OF INSTABILITIES 
C FROM BHLUMI OR ANYTHING ELSE            
C THIRD AXIS ALONG PNEUTR                                             
      X1 = PSUM(1)                                                  
      X2 = PSUM(2) 
      FI0  =PHANGF(X1,X2)           
      X1 = PSUM(3)                                                 
      X2 = SQRT(PSUM(1)**2+PSUM(2)**2)                            
      TH0  =PHANXY(X1,X2) 
      CALL PHSPAJ(-2,PNEUTR,PAA,PP,PE)    
      CALL PHLORT(3,-FI0,PNEUTR,VEC,PAA,PP,PE)
      CALL PHLORT(2,-TH0,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT3(-FI0,PSUM,PSUM)
      CALL PHROT2(-TH0,PSUM,PSUM)
      BSTA=(PSUM(4)-PSUM(3))/SQRT(PSUM(4)**2-PSUM(3)**2)
      BSTB=(PSUM(4)+PSUM(3))/SQRT(PSUM(4)**2-PSUM(3)**2)
      CALL PHLORT(1,BSTA,PNEUTR,VEC,PAA,PP,PE)
      CALL PHSPAJ(-1,PNEUTR,PAA,PP,PE)                                   
      AMNE=PHMAS(PNEUTR)                                              
      AMCH=PHMAS(PAA)                                                 
      IF (AMCH.LT.0.0D0) AMCH=AMEL                                   
      IF (AMNE.LT.0.D0) AMNE=0.0
      AMTO =PAA(4)+PNEUTR(4)
      CALL PHRAND(RRR,8)
      PRHARD= (1D0/PI/ALFINV)**2* (2.D0*LOG(AMTO/AMEL/2D0)) *
     $        LOG(1D0/XK0) * LOG(AMTO**2/2D0/AMEL**2)
C this just enforces hard pairs to be generated 'always'
C this is for the sake of tests only.
C       PRHARD=0.99  
C
         XMP=2D0*AMEL*EXP(RRR(1)*LOG(AMTO/2D0/AMEL)) 
         XP =AMTO*XK0*EXP(RRR(2)*LOG(1D0/XK0))    
         C1 =1.0-4D0*AMEL**2/AMTO**2*
     $       EXP(RRR(3)*LOG(AMTO**2/2D0/AMEL**2))
         FIX1=2D0*PI*RRR(4)
         C2  =1D0-2D0*RRR(5) 
         FIX2=2D0*PI*RRR(6) 
         JESLI=(RRR(7).LT.PRHARD).AND.
     $         (XMP.LT.(AMTO-AMNE-AMCH)).AND.
     $         (XP .GT.XMP)             .AND.
     $         (XP .LT.((AMTO**2+XMP**2-(AMCH+AMNE)**2)/2D0/AMTO))
C histograming .......................
         JESLIK=
     $         (XP .LT.((AMTO**2+XMP**2-(AMCH+AMNE)**2)/2D0/AMTO))
      WTA=0D0
      IF (JESLIK) WTA=1D0
!      CALL GMONIT( 0,101   ,WTA,1D0,0D0)
         JESLIK=
     $         (XMP.LT.(AMTO-AMNE-AMCH))
      WTA=0D0
      IF (JESLIK) WTA=1D0
!      CALL GMONIT( 0,102   ,WTA,1D0,0D0)
         JESLIK=
     $         (XMP.LT.(AMTO-AMNE-AMCH)).AND.
     $         (XP .GT.XMP)

      WTA=0D0
      IF (JESLIK) WTA=1D0
!      CALL GMONIT( 0,103   ,WTA,1D0,0D0)
         JESLIK=
     $         (XMP.LT.(AMTO-AMNE-AMCH)).AND.
     $         (XP .GT.XMP)             .AND.
     $         (XP .LT.((AMTO**2+XMP**2-(AMCH+AMNE)**2)/2D0/AMTO))
      WTA=0D0
      IF (JESLIK) WTA=1D0
!      CALL GMONIT( 0,104   ,WTA,1D0,0D0)
C end of histograming ................  
C ... rejection due to parameters out of phase space
      IF (.NOT.JESLI) RETURN
C ... jacobians weights etc. 
      XMK2=AMTO**2+XMP**2-2D0*XP*AMTO
      XPMAX=(AMTO**2+XMP**2-(AMCH+AMNE)**2)/2D0/AMTO
      YOT3=(1-C1+4D0*AMEL**2/AMTO**2)/(1-C1+XMP**2/AMTO**2)
      YOT2=(1-C1)*(1+C1)/2D0/(1-C1+XMP**2/AMTO**2)*
     $     XLAM(AMTO**2,XMK2,XMP**2)/(AMTO**2+XMP**2-XMK2)
      YOT1=XLAM(1D0,AMEL**2/XMP**2, AMEL**2/XMP**2)*
     $     ((1-C1+XMP**2/AMTO**2)/(1-C1+XMP**2/XP**2))**2*
     $     (1-XP/XPMAX+0.5D0*(XP/XPMAX)**2)*2D0/3D0
C note that the factor 2/3 in YOT1 above should be replaced by the 
C appropriate A-P kernel for gamma splitting to e+e- !!!!!!!
C the part of the weight below, should have average 1, but fluctuates 
C wildly. This cancelation is important ingredient of the leading logs.
C      YOT1=YOT1*
C     $     (1D0-XP/(0.5D0*AMTO))/(1D0-XP/XPMAX)*
C     $     XLAM(1D0,AMCH**2/XMK2,AMNEU**2/XMK2)/
C     $     XLAM(1D0,AMCH**2/AMTO**2,AMNEU**2/AMTO**2)
      WT=YOT1*YOT2*YOT3
C histograming .......................
C      CALL GMONIT( 0,105   ,WT  ,1D0,0D0) 
C      CALL GMONIT( 0,106   ,YOT1,1D0,0D0) 
C      CALL GMONIT( 0,107   ,YOT2,1D0,0D0) 
C      CALL GMONIT( 0,108   ,YOT3,1D0,0D0)
C      CALL GMONIT( 0,109   ,YOT4,1D0,0D0)
C end of histograming ................ 
C ... rejection due to weight
      IF (RRR(8).GT.WT) THEN
        JESLI=.FALSE.
        RETURN
      ENDIF
C                                                                     
C                                                                     
C FRAGMENTATION COMES !!                                              
C                                                                     
C THIRD AXIS ALONG PNEUTR                                             
      X1 = PNEUTR(1)                                                  
      X2 = PNEUTR(2)                                                 
      FI1  =PHANGF(X1,X2)                                             
      X1 = PNEUTR(3)                                                 
      X2 = SQRT(PNEUTR(1)**2+PNEUTR(2)**2)                            
      TH1  =PHANXY(X1,X2) 
      CALL PHSPAJ(0,PNEUTR,PAA,PP,PE)                                             
      CALL PHLORT(3,-FI1,PNEUTR,VEC,PAA,PP,PE)
      CALL PHLORT(2,-TH1,PNEUTR,VEC,PAA,PP,PE)
        VEC(4)=0D0                                                  
        VEC(3)=0D0                                                  
        VEC(2)=0D0                                                  
        VEC(1)=1D0                                                  
        FI2=PHANGF(VEC(1),VEC(2))                                              
        CALL PHLORT(3,-FI2,PNEUTR,VEC,PAA,PP,PE)
      CALL PHSPAJ(1,PNEUTR,PAA,PP,PE)                                             
C STEALING FROM PAA AND PNEUTR ENERGY FOR THE pair
C ====================================================  
C NEW MOMENTUM OF PAA AND PNEUTR (IN THEIR VERY REST FRAME)
C 1) PARAMETERS..... 
      AMCH2=AMCH**2
      AMNE2=AMNE**2
      AMTOST=XMK2
      QNEW=XLAM(AMTOST,AMNE2,AMCH2)/SQRT(AMTOST)/2D0
      QOLD=PNEUTR(3)
      GCHAR=(QNEW**2+QOLD**2+AMCH**2)/
     &      (QNEW*QOLD+SQRT((QNEW**2+AMCH**2)*(QOLD**2+AMCH**2)))
      GNEU= (QNEW**2+QOLD**2+AMNE**2)/
     &      (QNEW*QOLD+SQRT((QNEW**2+AMNE**2)*(QOLD**2+AMNE**2)))
C      GCHAR=(QOLD**2-QNEW**2)/(
C     &       QNEW*SQRT(QOLD**2+AMCH2)+QOLD*SQRT(QNEW**2+AMCH2)
C     &                        )
C      GCHAR=SQRT(1D0+GCHAR**2) 
C      GNEU=(QOLD**2-QNEW**2)/(
C     &       QNEW*SQRT(QOLD**2+AMNE2)+QOLD*SQRT(QNEW**2+AMNE2)
C     &                        )
C      GNEU=SQRT(1D0+GNEU**2) 
      IF(GNEU.LT.1..OR.GCHAR.LT.1.) THEN
        PRINT *,' PHPAR GBOOST LT 1., LIMIT OF PHASE SPACE '
     &         ,GNEU,GCHAR,QNEW,QOLD,AMTO,AMTOST,AMNE,AMCH
        PRINT *,XK,XKM,XK0DEC,AXK
        RETURN
      ENDIF
      PARCH =GCHAR+SQRT(GCHAR**2-1.0D0)
      PARNEU=GNEU -SQRT(GNEU**2 -1.0D0)
C 2) REDUCTIVE BOOSTS
      CALL PHBOST(PARNEU,VEC ,VEC )
      CALL PHBOST(PARNEU,PNEUTR,PNEUTR)
      CALL PHBOST(PARCH,PAA ,PAA )
      CALL PHSPAJ(2,PNEUTR,PAA,PP,PE)                                             
C TIME FOR THE PHOTON that is electron pair
      PMOD=XLAM(XMP**2,AMEL**2,AMEL**2)/XMP/2D0
      S2=SQRT(1D0-C2**2)
      PP(4)=XMP/2D0
      PP(3)=PMOD*C2
      PP(2)=PMOD*S2*COS(FIX2)
      PP(1)=PMOD*S2*SIN(FIX2)
      PE(4)= PP(4)
      PE(3)=-PP(3)
      PE(2)=-PP(2)
      PE(1)=-PP(1)
C PHOTON ENERGY and momentum IN THE REDUCED SYSTEM
      PENE=(AMTO**2-XMP**2-XMK2)/2D0/SQRT(XMK2)
      PPED=SQRT(PENE**2-XMP**2)
      FI3=FIX1
      COSTHG=C1
      SINTHG=SQRT(1D0-C1**2)
      X1 = -COSTHG
      X2 =  SINTHG
      TH3  =PHANXY(X1,X2)
      PHOT(1)=PMOD*SINTHG*COS(FI3) 
      PHOT(2)=PMOD*SINTHG*SIN(FI3)
C MINUS BECAUSE AXIS OPPOSITE TO PAA
      PHOT(3)=-PMOD*COSTHG
      PHOT(4)=PENE
C ROTATE TO PUT PHOTON ALONG THIRD AXIS
      X1 = PHOT(1)
      X2 = PHOT(2) 
      CALL PHLORT(3,-FI3,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT3(-FI3,PHOT,PHOT) 
      CALL PHLORT(2,-TH3,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT2(-TH3,PHOT,PHOT)
      CALL PHSPAJ(21,PNEUTR,PAA,PP,PE)                                             
C ... now get the pair !
      PAIRB=PENE/XMP+PPED/XMP
      CALL PHBOST(PAIRB,PE,PE)   
      CALL PHBOST(PAIRB,PP,PP)   
      CALL PHSPAJ(3,PNEUTR,PAA,PP,PE) 
      GAMM=(PNEUTR(4)+PAA(4)+PP(4)+PE(4))/AMTO
      BPAR=GAMM-SQRT(GAMM**2-1D0)
      CALL PHLORT(1, BPAR,PNEUTR,VEC,PAA,PP,PE)
      CALL PHBOST( BPAR,PHOT,PHOT)
      CALL PHSPAJ(4,PNEUTR,PAA,PP,PE)                                             
C BACK IN THE TAU REST FRAME BUT PNEUTR NOT YET ORIENTED.
      X1 = PNEUTR(1)
      X2 = PNEUTR(2)
      FI4  =PHANGF(X1,X2)
      X1 = PNEUTR(3)
      X2 = SQRT(PNEUTR(1)**2+PNEUTR(2)**2)
      TH4  =PHANXY(X1,X2)
      CALL PHLORT(3, FI4,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT3( FI4,PHOT,PHOT)
      CALL PHLORT(2,-TH4,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT2(-TH4,PHOT,PHOT)
      X1 = VEC(1)
      X2 = VEC(2)
      FI5=PHANGF(X1,X2)
      CALL PHLORT(3,-FI5,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT3(-FI5,PHOT,PHOT)
C PAA RESTORES ORIGINAL DIRECTION 
      CALL PHLORT(3, FI2,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT3( FI2,PHOT,PHOT)
      CALL PHLORT(2, TH1,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT2( TH1,PHOT,PHOT) 
      CALL PHLORT(3, FI1,PNEUTR,VEC,PAA,PP,PE)
      CALL PHROT3( FI1,PHOT,PHOT)
      CALL PHSPAJ(10,PNEUTR,PAA,PP,PE)
      CALL PHLORT(1,BSTB,PNEUTR,VEC,PAA,PP,PE)
      CALL PHLORT(2,TH0,PNEUTR,VEC,PAA,PP,PE)
      CALL PHLORT(3,FI0,PNEUTR,VEC,PAA,PP,PE)
      CALL PHSPAJ(11,PNEUTR,PAA,PP,PE)                                             
C      STOP                                             
      END
      SUBROUTINE PHPRTR(IBRAN,PHOT)       
C.----------------------------------------------------------------------
C.
C.    PHPRTR:   PHotos PaRticle TRansformation
C.
C.    Purpose:  THIS ROUTINE PERFORMS LORENTZ 
C.              TRANSFORMATIONS ON PHOT IN CONSTRUCTING
C.              Elementary branch after pair generation
C.              See routine PHOPAR for details of angles 
C.              and boost parameters
C.
C.    Input Parameter: 
C.                     IBRAN   = 1 for charged particle   
C.                             =-1 for neutral particle(s)
C.     
C.    Input/Output Parameters: PHOT
C.                                  Four-vector transformed
C.
C.    Author(s):  Z. Was,                         Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IBRAN
      REAL*8                                   
     $ FI0,FI1,FI2,FI3,FI4,FI5,TH0,TH1,TH3,TH4
     $,PARNEU,PARCH,BPAR,BSTA,BSTB
      COMMON /PHKINP/ 
     $ FI0,FI1,FI2,FI3,FI4,FI5,TH0,TH1,TH3,TH4
     $,PARNEU,PARCH,BPAR,BSTA,BSTB
      REAL*8 PHOT(4)
      CALL PHROT3(-FI0,PHOT,PHOT)
      CALL PHROT2(-TH0,PHOT,PHOT)
      CALL PHBOST(BSTA,PHOT,PHOT)
      CALL PHROT3(-FI1,PHOT,PHOT)
      CALL PHROT2(-TH1,PHOT,PHOT)
      CALL PHROT3(-FI2,PHOT,PHOT)
      IF(IBRAN.EQ.-1) THEN
        CALL PHBOST(PARNEU,PHOT,PHOT)
      ELSE
        CALL PHBOST(PARCH,PHOT,PHOT)
      ENDIF
      CALL PHROT3(-FI3,PHOT,PHOT) 
      CALL PHROT2(-TH3,PHOT,PHOT)
      CALL PHBOST(BPAR,PHOT,PHOT)
      CALL PHROT3( FI4,PHOT,PHOT)
      CALL PHROT2(-TH4,PHOT,PHOT)
      CALL PHROT3(-FI5,PHOT,PHOT)
      CALL PHROT3( FI2,PHOT,PHOT)
      CALL PHROT2( TH1,PHOT,PHOT) 
      CALL PHROT3( FI1,PHOT,PHOT)
      CALL PHBOST(BSTB,PHOT,PHOT)
      CALL PHROT2( TH0,PHOT,PHOT)
      CALL PHROT3( FI0,PHOT,PHOT)
      END
      FUNCTION PHMAS(VEC)
C.----------------------------------------------------------------------
C.
C.    PHMAS:   PHotos MASs
C.
C.    Purpose:  calculates mass of the four-vector
C.
C.    Input Parameter:    VEC
C.    Output Parameters:  PHMAS
C.                    
C.
C.    Author(s):  S.Jadach Z. Was,                Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 PHMAS
      REAL*8 VEC(4)
      PHMAS=SQRT(ABS(VEC(4)**2-VEC(1)**2-VEC(2)**2-VEC(3)**2)) 
      END
      SUBROUTINE PHLORT(KEY,PRM,PNEUTR,PNU,PAA,PP,PE)
C.----------------------------------------------------------------------
C.
C.    PHLORT:   PHotos LORentz Transformation
C.
C.    Purpose:  THIS ROUTINE PERFORMS LORENTZ 
C.              TRANSFORMATION ON MANY 4-VECTORS        
C.
C.    Input Parameter: 
C.                     KEY   =1    BOOST    ALONG   3RD AXIS 
C.                           =2    ROTATION AROUND 2ND AXIS   
C.                           =3    ROTATION AROUND 3RD AXIS 
C                      PRM         TRANSFORMATION PARAMETER 
C.                               - ANGLE OR EXP(HIPERANGLE).    
C.     
C.    Input/Output Parameters: PNEUTR,PNU,PAA,PP,PE
C.       Four-vectors transformed
C.
C.    Author(s):  Z. Was,                         Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE 
      REAL*8 PRM                 
      REAL*8 PNU(4),PAA(4),PE(4),PP(4)
      REAL*8 PNEUTR(4)
      INTEGER KEY                                                 
C                                                                     
      IF     (KEY.EQ.1) THEN                                           
        CALL PHBOST(PRM,PNEUTR,PNEUTR)                                 
        CALL PHBOST(PRM,PNU ,PNU )                                     
        CALL PHBOST(PRM,PAA ,PAA )                                     
        CALL PHBOST(PRM,PE ,PE )                                     
        CALL PHBOST(PRM,PP ,PP )                                     
      ELSEIF (KEY.EQ.2) THEN                                           
        CALL PHROT2(PRM,PNEUTR,PNEUTR)                                 
        CALL PHROT2(PRM,PNU ,PNU )                                     
        CALL PHROT2(PRM,PAA ,PAA )                                    
        CALL PHROT2(PRM,PE  ,PE  )                                     
        CALL PHROT2(PRM,PP  ,PP )                                    
      ELSEIF (KEY.EQ.3) THEN                                           
        CALL PHROT3(PRM,PNEUTR,PNEUTR)                                
        CALL PHROT3(PRM,PNU ,PNU )                                     
        CALL PHROT3(PRM,PAA ,PAA )                                    
        CALL PHROT3(PRM,PE  ,PE  )                                     
        CALL PHROT3(PRM,PP  ,PP )                                    
      ELSE                                                            
        PRINT *, 'STOP IN PHLORT. WRONG KEY'                        
        STOP                                                           
      ENDIF                                                            
      END                                                            

      SUBROUTINE PHROT2(PH1,PVEC,QVEC)
C.----------------------------------------------------------------------
C.
C.    PHROT2:   PHotos Rotation axis 2
C.
C.    Purpose:  Rotates vector PVEC into QVEC by angle PH1 around
C.              second axis
C.
C.    Input Parameter:    PH1,PVEC 
C.    Output Parameters:  QVEC
C.                    
C.
C.    Author(s):  S.Jadach Z. Was,                Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 PVEC(4),QVEC(4),RVEC(4)
      REAL*8 PHI,PH1,CS,SN
      INTEGER I
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)+SN*RVEC(3)
      QVEC(2)=RVEC(2)
      QVEC(3)=-SN*RVEC(1)+CS*RVEC(3)
      QVEC(4)=RVEC(4)
      RETURN
      END

      SUBROUTINE PHROT3(PH1,PVEC,QVEC)
C.----------------------------------------------------------------------
C.
C.    PHROT3:   PHotos Rotation axis 3
C.
C.    Purpose:  Rotates vector PVEC into QVEC by angle PH1 around
C.              third axis
C.
C.    Input Parameter:    PH1,PVEC 
C.    Output Parameters:  QVEC
C.                    
C.
C.    Author(s):  S.Jadach Z. Was,                Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 PVEC(4),QVEC(4),RVEC(4)
      REAL*8 PHI,PH1,CS,SN
      INTEGER I
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)-SN*RVEC(2)
      QVEC(2)= SN*RVEC(1)+CS*RVEC(2)
      QVEC(3)=RVEC(3)
      QVEC(4)=RVEC(4)
      END

      SUBROUTINE PHBOST(EXE,PVEC,QVEC)
C.----------------------------------------------------------------------
C.
C.    PHBOST:   PHotos Boost axis 3
C.
C.    Purpose:  Boost vector PVEC into QVEC by parameter EXE along
C.              third axis
C.
C.    Input Parameter:    EXE,PVEC 
C.    Output Parameters:  QVEC
C.                    
C.
C.    Author(s):  S.Jadach Z. Was,                Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 PVEC(4),QVEC(4),RVEC(4)
      REAL*8 EXE,RPL,RMI,QPL,QMI
      INTEGER I
C
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      RPL=RVEC(4)+RVEC(3)
      RMI=RVEC(4)-RVEC(3)
      QPL=RPL*EXE
      QMI=RMI/EXE
      QVEC(1)=RVEC(1)
      QVEC(2)=RVEC(2)
      QVEC(3)=(QPL-QMI)/2
      QVEC(4)=(QPL+QMI)/2
      RETURN
      END

      DOUBLE PRECISION FUNCTION PHANXY(X,Y) 
C.----------------------------------------------------------------------
C.
C.    PHANXY:   PHotos ANgle X Y
C.
C.    Purpose:  CALCULATES ANGLE IN (0,PI) RANGE OUT OF X-Y
C.
C.    Input Parameter:    X,Y 
C.    Output Parameters:  PHANXY Function name
C.                    
C.
C.    Author(s):  S.Jadach Z. Was,                Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE      
      REAL*8 PI,X,Y,PHANXY,THE                          
      DATA PI /3.141592653589793238462643D0/                            
C     
      IF(X.EQ.0D0.AND.Y.EQ.0.D0) THEN
       PHANXY=0D0
       RETURN
      ENDIF                                                      
      IF(ABS(Y).LT.ABS(X)) THEN                                         
        THE=ATAN(ABS(Y/X))                                              
        IF(X.LE.0D0) THE=PI-THE                                         
      ELSE                                                              
        THE=ACOS(X/SQRT(X**2+Y**2))                                     
      ENDIF                                                             
      PHANXY=THE                                                         
      RETURN                                                            
      END                                                               

      DOUBLE PRECISION FUNCTION PHANGF(X,Y)
C.----------------------------------------------------------------------
C.
C.    PHANGF:   PHotos ANGle Fi
C.
C.    Purpose:  CALCULATES ANGLE IN (0,2*PI) RANGE OUT OF X-Y
C.
C.    Input Parameter:    X,Y 
C.    Output Parameters:  PHANGF Function name
C.                    
C.
C.    Author(s):  S.Jadach Z. Was,                Created at:  01/06/93
C.                                                Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 PI,X,Y,THE
      DATA PI /3.1415926535897932D0/
      IF(X.EQ.0D0.AND.Y.EQ.0.D0) THEN
       PHANGF=0D0
       RETURN
      ENDIF
      IF(ABS(Y).LT.ABS(X)) THEN
        THE=ATAN(ABS(Y/X))
        IF(X.LE.0D0) THE=PI-THE
      ELSE
        THE=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      IF(Y.LT.0D0) THE=2D0*PI-THE
      PHANGF=THE
      END


      SUBROUTINE PHRAND(DRVEC,LEN)
C.----------------------------------------------------------------------
C.
C.    PHRAND:   Random Number interface 
C.
C.    Purpose:  One can make a choice of random number generator
C.
C.    Input Parameter:    LEN:    Number or random numbers to be generated
C.
C.    Output Parameters:  DRVEC vector of random numbers
C.
C.    Author(s):  S. Jadach, Z. Was,               Created at:  01/06/93
C.                                                 Last Update: 11/11/98
C.
C.----------------------------------------------------------------------
C     ***************************
C Switchable random number generator
C Translation to double precision
C     ***************************
      COMMON / PHRNDP / KEYRND
      save   / PHRNDP /
      DOUBLE PRECISION DRVEC(*)
      DIMENSION RVEC(1000)
      IF(LEN.LT.1.OR.LEN.GT.1000) GOTO 901
   10 CONTINUE
      IF(KEYRND.EQ.1) THEN
         CALL MARRAN(RVEC,LEN)
      ELSEIF(KEYRND.EQ.2) THEN
         CALL ECURAN(RVEC,LEN)
      ELSEIF(KEYRND.EQ.3) THEN
         CALL CARRAN(RVEC,LEN)
      ELSE
         GOTO 902
      ENDIF
C random numbers 0 and 1 not accepted
      DO 30 I=1,LEN
      IF(RVEC(I).LE.0E0.OR.RVEC(I).GE.1E0) THEN
        WRITE(6,*) ' +++++ PHRAND: RVEC=',RVEC(I)
        GOTO 10
      ENDIF
      DRVEC(I)=RVEC(I)
   30 CONTINUE
      RETURN
  901 WRITE(6,*) ' +++++ STOP IN PHRAND: LEN=',LEN
      STOP
  902 WRITE(6,*) ' +++++ STOP IN PHRAND: WRONG KEYRND',KEYRND
      STOP
      END


C=======================================================================
C=======================================================================
C=======================================================================
C==Received: by dxmint.cern.ch (cernvax) (5.57/3.14)
C== id AA13405; Wed, 23 Jan 91 17:19:06 +0100
C==Message-Id: <9101231619.AA13405@dxmint.cern.ch>
C==Received: by cernapo; Wed, 23 Jan 91 17:23:40 +0100
C==Received: by apojames.cern.ch; Wed, 23 Jan 91 17:05:23 CET
C==Date: Wed, 23 Jan 91 17:05:23 CET
C==From: james@cernapo.cern.ch (Frederick James)
C==To: jadach@cernvm
C==Subject: Random generators
C==
C==      PROGRAM PSEUDORAN
C==C  CPC # ABTK                                           CPC # ABTK
C==C         Pseudorandom generator demonstration (test case)
C==      DIMENSION RVEC(1000)
C==      DIMENSION VERI(5), ISD25(25)
C==C
C==C
C==C   ................................................
C==      WRITE(6,'(20X,A)') 'DEMONSTRATION OF PSEUDORANDOM GENERATORS'
C==      WRITE(6,'(20X,A)') 'MACHINE/SYSTEM: date:'
C==      WRITE(6,'(/20X,A/)') 'INITIALIZATION AND TEST OF PORTABILITY'
C==C   ................................................
C==C
C==C                   initialization and verification  RANMAR
C==        DO 40 I9= 1, 20
C==   40   CALL RANMAR(RVEC,1000)
C==      CALL RANMAR(RVEC,5)
C==      DO 41 I= 1 ,5
C==   41 VERI(I) = (4096.*RVEC(I))*(4096.)
C==      WRITE(6,'(A,5F12.1/)') '  RANMAR 20001  ',VERI
C==C
C==C                   initialization and verification  RANECU
C==      CALL RANECU(RVEC,1000)
C==      CALL RANECU(VERI,5)
C==      DO 52 I= 1 ,5
C==   52 VERI(I) = 4096.*(4096.*VERI(I))
C==      WRITE(6,'(A,5F12.1/)') '  RANECU 1001   ',VERI
C==C
C==C                   initialization and verification  RCARRY
C==      CALL RCARRY(RVEC,1000)
C==      CALL RCARRY(VERI,5)
C==      DO 62 I= 1 ,5
C==   62 VERI(I) = 4096.*(4096.*VERI(I))
C==      WRITE(6,'(A,5F12.1/)') '  RCARRY 1001   ',VERI
C==C
C==      WRITE(6,'(//20X,A/)') 'TEST OF REPEATABILITY'
C==C  .................................................
C==C                  verify restarting      RANMAR
C==      WRITE(6,'(/A)') '   THE NEXT LINE SHOULD BE REPEATED:'
C==      CALL RMARUT(IMAR1,IMAR2,IMAR3)
C==      CALL RANMAR(RVEC,777)
C==      CALL RANMAR(VERI,5)
C==      WRITE(6,'(A,5F12.9)') '       RANMAR 1 ',VERI
C==      CALL RMARIN(IMAR1,IMAR2,IMAR3)
C==      CALL RANMAR(RVEC,777)
C==      CALL RANMAR(VERI,5)
C==      WRITE(6,'(A,5F12.9)') '       RANMAR 2 ',VERI
C==C
C==C                  verify restarting      RANECU
C==      WRITE(6,'(/A)') '   THE NEXT LINE SHOULD BE REPEATED:'
C==      CALL RECUUT(IS1,IS2)
C==      CALL RANECU(RVEC,777)
C==      CALL RANECU(VERI,5)
C==      WRITE(6,'(A,5F12.9)') '       RANECU 1 ',VERI
C==      CALL RECUIN(IS1,IS2)
C==      CALL RANECU(RVEC,777)
C==      CALL RANECU(VERI,5)
C==      WRITE(6,'(A,5F12.9)') '       RANECU 2 ',VERI
C==C
C==C                  verify restarting      RCARRY
C==      WRITE(6,'(/A)') '   THE NEXT LINE SHOULD BE REPEATED:'
C==      CALL RCARUT(ISD25)
C==      CALL RCARRY(RVEC,777)
C==      CALL RCARRY(VERI,5)
C==      WRITE(6,'(A,5F12.9)') '       RCARRY 1 ',VERI
C==      CALL RCARIN(ISD25)
C==      CALL RCARRY(RVEC,777)
C==      CALL RCARRY(VERI,5)
C==      WRITE(6,'(A,5F12.9)') '       RCARRY 2 ',VERI
C==C
C==      STOP
C==      END
C=======================================================================
C=======================================================================
C=======================================================================

      SUBROUTINE MARRAN(RVEC,LENV)
C =======================S. JADACH===================================
C == This commes from F. James, The name of RANMAR is changed to   ==
C == MARRAN in order to avoid interference with the version        ==
C == already in use and the public library version (if present).   ==
C ==      THIS IS THE ONLY MODIFICATION !!!!                       ==
C ========================S. JADACH==================================
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        modified by F. James, 1988 and 1989, to generate a vector
C        of pseudorandom numbers RVEC of length LENV, and to put in
C        the COMMON block everything needed to specify currrent state,
C        and to add input and output entry points MARINI, MAROUT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANMAR:                                  ++
C!!!      CALL RANMAR (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL MARINI(I1,N1,N2)   initializes the generator from one ++
C!!!                   32-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++
C!!!                    output by MAROUT)                            ++
C!!!      CALL MAROUT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(*)
      COMMON/RASET1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RANMAR without MARINI.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      MARINI(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RANMAR, may be called before
C         generating pseudorandom numbers with RANMAR. The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      WRITE(6,'(A,5I10)')
     $' MARran INITIALIZED: IJ,KL,IJKL,NTOT,NTOT2=',IJ,KL,IJKL,NTOT,NTOT2
      DO 2 II= 1, 97
      S = 0.
      T = .5
      DO 3 JJ= 1, 24
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = 0.5*T
    2 U(II) = S
      TWOM24 = 1.0
      DO 4 I24= 1, 24
    4 TWOM24 = 0.5*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
        WRITE(6,'(A,I15)') ' MARINI SKIPPING OVER ',NOW
       DO 40 IDUM = 1, NTOT
       UNI = U(I97)-U(J97)
       IF (UNI .LT. 0.)  UNI=UNI+1.
       U(I97) = UNI
       I97 = I97-1
       IF (I97 .EQ. 0)  I97=97
       J97 = J97-1
       IF (J97 .EQ. 0)  J97=97
       C = C - CD
       IF (C .LT. 0.)  C=C+CM
   40  CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. 0.)  UNI=UNI+1.
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. 0.)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. 0.) UNI=UNI+1.
      RVEC(IVEC) = UNI
C             Replace exact zeros by uniform distr. *2**-24
         IF (UNI .EQ. 0.)  THEN
         ZUNI = TWOM24*U(2)
C             An exact zero here is very unlikely, but let's be safe.
         IF (ZUNI .EQ. 0.) ZUNI= TWOM24*TWOM24
         RVEC(IVEC) = ZUNI
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY MAROUT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END

      SUBROUTINE CARRAN(RVEC,LENV)
C         Add-and-carry random number generator proposed by
C         Marsaglia and Zaman in SIAM J. Scientific and Statistical
C             Computing, to appear probably 1990.
C         modified with enhanced initialization by F. James, 1990
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for CARRAN:                                  ++
C!!!      CALL CARRAN (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL CARINI(INT)     initializes the generator from one    ++
C!!!                   32-bit integer INT                            ++
C!!!      CALL CARRES(IVEC)    restarts the generator from vector    ++
C!!!                   IVEC of 25 32-bit integers (see CAROUT)       ++
C!!!      CALL CAROUT(IVEC)    outputs the current values of the 25  ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (TWOP12=4096.)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24
      LOGICAL NOTYET
      DATA NOTYET/.TRUE./
      DATA I24,J24,CARRY/24,10,0./
C
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = 314159265
         WRITE(6,'(A,I12)') ' CARRAN DEFAULT INITIALIZATION: ',JSEED
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
   50    CONTINUE
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .LT. SEEDS(14)) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(I24) - SEEDS(J24) - CARRY
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = I24 - 1
      IF (I24 .EQ. 0)  I24 = 24
      J24 = J24 - 1
      IF (J24 .EQ. 0)  J24 = 24
      RVEC(IVEC) = UNI
  100 CONTINUE
      RETURN
C           Entry to input and float integer seeds from previous run
      ENTRY CARRES(ISDEXT)
         TWOM24 = 1.
         DO 195 I= 1, 24
  195    TWOM24 = TWOM24 * 0.5
      WRITE(6,'(A)') ' FULL INITIALIZATION OF CARRAN WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = REAL(MOD(ISDEXT(25),10))*TWOM24
      ISD = ISDEXT(25)/10
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = ISD
      RETURN
C                    Entry to ouput seeds as integers
      ENTRY CAROUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ICARRY = 0
      IF (CARRY .GT. 0.)  ICARRY = 1
      ISDEXT(25) = 1000*J24 + 10*I24 + ICARRY
      RETURN
C                    Entry to initialize from one integer
      ENTRY CARINI(INSEED)
      JSEED = INSEED
      WRITE(6,'(A,I12)') ' CARRAN INITIALIZED FROM SEED ',INSEED
C      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
  350    CONTINUE
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .LT. SEEDS(14)) CARRY = TWOM24
      RETURN
      END

      SUBROUTINE ECURAN(RVEC,LEN)
C         Random number generator given by L'Ecuyer in
C            Comm. ACM Vol 31, p.742, 1988
C            modified by F. James to return a vector of numbers
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for ECURAN:                                  ++
C!!!      CALL ECURAN (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL ECUINI(I1,I2)    initializes the generator from two   ++
C!!!                   32-bit integers I1 and I2                     ++
C!!!      CALL ECUOUT(I1,I2)    outputs the current values of the    ++
C!!!                   two integer seeds, to be used for restarting  ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(*)
      SAVE ISEED1,ISEED2
      DATA ISEED1,ISEED2 /12345,67890/
C
      DO 100 I= 1, LEN
      K = ISEED1/53668
      ISEED1 = 40014*(ISEED1 - K*53668) - K*12211
      IF (ISEED1 .LT. 0) ISEED1=ISEED1+2147483563
C
      K = ISEED2/52774
      ISEED2 = 40692*(ISEED2 - K*52774) - K* 3791
      IF (ISEED2 .LT. 0) ISEED2=ISEED2+2147483399
C
      IZ = ISEED1 - ISEED2
      IF (IZ .LT. 1)  IZ = IZ + 2147483562
C
      RVEC(I) = REAL(IZ) * 4.656613E-10
  100 CONTINUE
      RETURN
C
      ENTRY ECUINI(IS1,IS2)
      ISEED1 = IS1
      ISEED2 = IS2
      RETURN
C
      ENTRY ECUOUT(IS1,IS2)
      IS1 = ISEED1
      IS2 = ISEED2
      RETURN
      END





















