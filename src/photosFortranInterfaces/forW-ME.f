
c ###### 
c  standard routine !!!
c ######
c$$$      FUNCTION PHOCOR(MPASQR,MCHREN,ME)
c$$$C.----------------------------------------------------------------------
c$$$C.
c$$$C.    PHOTOS:   PHOton radiation in decays CORrection weight from
c$$$C.              matrix elements
c$$$C.
c$$$C.    Purpose:  Calculate  photon  angle.  The reshaping functions  will
c$$$C.              have  to  depend  on the spin S of the charged particle.
c$$$C.              We define:  ME = 2 * S + 1 !
c$$$C.
c$$$C.    Input Parameters:  MPASQR:  Parent mass squared,
c$$$C.                       MCHREN:  Renormalised mass of charged system,
c$$$C.                       ME:      2 * spin + 1 determines matrix element
c$$$C.
c$$$C.    Output Parameter:  Function value.
c$$$C.
c$$$C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
c$$$C.                                                Last Update: 21/03/93
c$$$C.
c$$$C.----------------------------------------------------------------------
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION MPASQR,MCHREN,BETA,XX,YY,DATA
c$$$      INTEGER ME
c$$$      REAL*8 PHOCOR,PHOFAC,WT1,WT2,WT3
c$$$      DOUBLE PRECISION MCHSQR,MNESQR
c$$$      REAL*8 PNEUTR
c$$$      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
c$$$      DOUBLE PRECISION COSTHG,SINTHG
c$$$      REAL*8 XPHMAX,XPHOTO
c$$$      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
c$$$      INTEGER IREP
c$$$      REAL*8 PROBH,CORWT,XF
c$$$      COMMON/PHOPRO/PROBH,CORWT,XF,IREP
c$$$C--
c$$$C--   Shaping (modified by ZW)...
c$$$      XX=4.D0*MCHSQR/MPASQR*(1.D0-XPHOTO)/(1.D0-XPHOTO+(MCHSQR-MNESQR)/
c$$$     &MPASQR)**2
c$$$      IF (ME.EQ.1) THEN
c$$$        YY=1.D0
c$$$        WT3=(1.D0-XPHOTO/XPHMAX)/((1.D0+(1.D0-XPHOTO/XPHMAX)**2)/2.D0)
c$$$      ELSEIF (ME.EQ.2) THEN
c$$$        YY=0.5D0*(1.D0-XPHOTO/XPHMAX+1.D0/(1.D0-XPHOTO/XPHMAX))
c$$$        WT3=1.D0
c$$$      ELSEIF ((ME.EQ.3).OR.(ME.EQ.4).OR.(ME.EQ.5)) THEN
c$$$        YY=1.D0
c$$$        WT3=(1.D0+(1.D0-XPHOTO/XPHMAX)**2-(XPHOTO/XPHMAX)**3)/
c$$$     &  (1.D0+(1.D0-XPHOTO/XPHMAX)** 2)
c$$$      ELSE
c$$$        DATA=(ME-1.D0)/2.D0
c$$$        CALL PHOERR(6,'PHOCOR',DATA)
c$$$        YY=1.D0
c$$$        WT3=1.D0
c$$$      ENDIF
c$$$      BETA=SQRT(1.D0-XX)
c$$$      WT1=(1.D0-COSTHG*SQRT(1.D0-MCHREN))/(1.D0-COSTHG*BETA)
c$$$      WT2=(1.D0-XX/YY/(1.D0-BETA**2*COSTHG**2))*(1.D0+COSTHG*BETA)/2.D0
c$$$      WT2=WT2*PHOFAC(1)
c$$$      PHOCOR=WT1*WT2*WT3
c$$$      CORWT=PHOCOR
c$$$      IF (PHOCOR.GT.1.D0) THEN
c$$$        DATA=PHOCOR
c$$$        CALL PHOERR(3,'PHOCOR',DATA)
c$$$      ENDIF
c$$$      RETURN
c$$$      END
c$$$


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

      SUBROUTINE PHOBWnlo(WT)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOtos Boson W correction weight
C.
C.    Purpose:  calculates correction weight due to amplitudes of 
C.              emission from W boson. It is ecact, but not verified
C.              for exponentiation yet.
C.              
C.              
C.              
C.
C.    Input Parameters:  Common /PHOEVT/, with photon added.
C.                       wt  to be corrected
C.                       
C.                       
C.                       
C.    Output Parameters: wt
C.
C.    Author(s):  G. Nanava, Z. Was               Created at:  13/03/03
C.                                                Last Update: 01/11/12
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION WT
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      REAL*8 ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),JMOPHO(2,NMXPHO),
     &              JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER I,JJ,II,I3,I4,IJ
      DOUBLE PRECISION EMU,MCHREN,BETA,COSTHG,MPASQR,XPH,
     &                 PW(0:3),PMU(0:3),PPHOT(0:3),PNE(0:3),
     &                 B_PW(0:3),B_PNE(0:3),B_PMU(0:3),AMPSQR,SANC_MC_INIT

      DOUBLE PRECISION SANC_WT
      EXTERNAL SANC_WT 
      INTEGER  d_h_nmxhep         ! maximum number of particles
      PARAMETER (d_h_NMXHEP=10000)
      REAL*8  d_h_phep,  d_h_vhep ! to be real*4/ *8  depending on host
      INTEGER d_h_nevhep,d_h_nhep,d_h_isthep,d_h_idhep,d_h_jmohep,
     $        d_h_jdahep
      COMMON /ph_hepevt/
     $      d_h_nevhep,               ! serial number
     $      d_h_nhep,                 ! number of particles
     $      d_h_isthep(d_h_nmxhep),   ! status code
     $      d_h_idhep(d_h_nmxhep),    ! particle ident KF
     $      d_h_jmohep(2,d_h_nmxhep), ! parent particles
     $      d_h_jdahep(2,d_h_nmxhep), ! childreen particles
     $      d_h_phep(5,d_h_nmxhep),   ! four-momentum, mass [GeV]
     $      d_h_vhep(4,d_h_nmxhep)    ! vertex [mm]



      DOUBLE PRECISION  pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af

!      write(*,*) 'IDPHOs=',IDPHO(1),IDPHO(2),IDPHO(3),IDPHO(4),IDPHO(5)
!      write(*,*) 'IDPHOs=',JDAPHO(1,1),npho
!      write(*,*) 'd_h_IDPHOs=',d_h_IDhep(1),d_h_IDhep(2),d_h_IDhep(3),d_h_IDhep(4),d_h_IDhep(5)

C--
        IF (ABS(IDPHO(1)).EQ.24.AND.
     $     ABS(IDPHO(JDAPHO(1,1)  )).GE.11.AND.
     $     ABS(IDPHO(JDAPHO(1,1)  )).LE.16.AND.
     $     ABS(IDPHO(JDAPHO(1,1)+1)).GE.11.AND.
     $     ABS(IDPHO(JDAPHO(1,1)+1)).LE.16     ) THEN

           IF(
     $      ABS(IDPHO(JDAPHO(1,1)  )).EQ.11.OR.
     $      ABS(IDPHO(JDAPHO(1,1)  )).EQ.13.OR.
     $      ABS(IDPHO(JDAPHO(1,1)  )).EQ.15    ) THEN 
              I=JDAPHO(1,1)
           ELSE
              I=JDAPHO(1,1)+1
           ENDIF
c..        muon energy   
           EMU=PPHO(4,I)
c..        muon mass square
           MCHREN=ABS(PPHO(4,I)**2-PPHO(3,I)**2
     $               -PPHO(2,I)**2-PPHO(1,I)**2)
           BETA=SQRT(1- MCHREN/ PPHO(4,I)**2)
           COSTHG=(PPHO(3,I)*PPHO(3,NPHO)+PPHO(2,I)*PPHO(2,NPHO)
     $                                   +PPHO(1,I)*PPHO(1,NPHO))/
     $            SQRT(PPHO(3,I)**2+PPHO(2,I)**2+PPHO(1,I)**2)   /
     $            SQRT(PPHO(3,NPHO)**2+PPHO(2,NPHO)**2+PPHO(1,NPHO)**2)
           MPASQR=PPHO(4,1)**2
           XPH=PPHO(4,NPHO)

c...       Initialization of the W->l\nu\gamma 
c...       decay Matrix Element parameters 
           CALL SANC_INIT(ALPHA,PHLUN)


           MB=PPHO(4,1)                      ! W boson mass
           MF2=dsqrt(MCHREN)                 ! muon mass

           DO IJ=1,d_h_nhep
            IF(ABS(d_h_idhep(IJ)).EQ.24) I3=IJ ! position of W 
           ENDDO
           IF(
     $      ABS(d_h_idhep(d_h_jdahep(1,I3)  )).EQ.11.OR.
     $      ABS(d_h_idhep(d_h_jdahep(1,I3)  )).EQ.13.OR.
     $      ABS(d_h_idhep(d_h_jdahep(1,I3)  )).EQ.15    ) THEN 
              I4=d_h_jdahep(1,I3)              ! position of lepton
           ELSE
              I4=d_h_jdahep(1,I3)+1            ! position of lepton
           ENDIF


              IF (d_h_idhep(I3).EQ.-24) QB=-1D0  ! W boson charge
              IF (d_h_idhep(I3).EQ.+24) QB=+1D0    
              IF (d_h_idhep(I4).GT.0d0) QF2=-1D0 ! lepton charge
              IF (d_h_idhep(I4).LT.0d0) QF2=+1D0


c...          Particle momenta before foton radiation; effective Born level
              DO JJ=1,4
                B_PW(mod(JJ,4))=d_h_phep(JJ,I3)  ! W boson
                B_PNE(mod(JJ,4))=d_h_phep(JJ,I3)-d_h_phep(JJ,I4) ! neutrino
                B_PMU(mod(JJ,4))=d_h_phep(JJ,I4) ! muon
              ENDDO

c..        Particle monenta after photon radiation
           DO JJ=1,4
             PW(mod(JJ,4))=PPHO(JJ,1)
             PMU(mod(JJ,4))=PPHO(JJ,I)
             PPHOT(mod(JJ,4))=PPHO(JJ,NPHO)
             PNE(mod(JJ,4))=PPHO(JJ,1)-PPHO(JJ,I)-PPHO(JJ,NPHO)
           ENDDO
C two options of calculating neutrino (spectator) mass
           MF1=SQRT(ABS(B_PNE(0)**2-B_PNE(1)**2-B_PNE(2)**2-B_PNE(3)**2))
           MF1=SQRT(ABS(  PNE(0)**2-  PNE(1)**2-  PNE(2)**2-  PNE(3)**2))

          CALL SANC_INIT1(QB,QF2,MF2,MB)
          WT=WT*SANC_WT(PW,PNE,PMU,PPHOT,B_PW,B_PNE,B_PMU)
        ENDIF
!      write(*,*)   'AMPSQR/EIKONALFACTOR= ',   AMPSQR/EIKONALFACTOR
      END


C========================================================== ==
C========================================================== ==
C these will be public for PHOTOS functions of W_ME class   ==
C========================================================== ==
C========================================================== ==

      DOUBLE PRECISION FUNCTION SANC_WT(PW,PNE,PMU,PPHOT,B_PW,B_PNE,B_PMU)
      IMPLICIT NONE

      DOUBLE PRECISION WDecayAmplitudeSqrKS_1ph,WDecayBornAmpSqrKS_1ph
      DOUBLE PRECISION EIKONALFACTOR,WDecayEikonalSqrKS_1ph,AMPSQR,
     &                 PW(0:3),PMU(0:3),PPHOT(0:3),PNE(0:3),
     &                 B_PW(0:3),B_PNE(0:3),B_PMU(0:3)
      EXTERNAL WDecayAmplitudeSqrKS_1ph,WDecayEikonalSqrKS_1ph,WDecayBornAmpSqrKS_1ph
c..        Exact amplitude square      
           AMPSQR=WDecayAmplitudeSqrKS_1ph(PW,PNE,PMU,PPHOT)

           EIKONALFACTOR=WDecayBornAmpSqrKS_1ph(B_PW,B_PNE,B_PMU)
     &                  *WDecayEikonalSqrKS_1ph(PW,PNE,PMU,PPHOT)
      
c..        New weight
!           write(*,*) 'B_pne=',B_PNE
!           write(*,*) 'B_PMU=',B_PMU
!           write(*,*) 'bornie=',WDecayBornAmpSqrKS_1ph(B_PW,B_PNE,B_PMU)

!           write(*,*) ' '
!           write(*,*) '  pne=',pne
!           write(*,*) '  pmu=',pmu
!           write(*,*) 'pphot=',pphot
!           write(*,*) ' '
!           write(*,*) '  b_pw=',B_PW
!           write(*,*) '  b_pne=',B_PNE
!           write(*,*) 'b_pmu=',B_PMU
 
 !          write(*,*) 'cori=',AMPSQR/EIKONALFACTOR,AMPSQR,EIKONALFACTOR
 
      SANC_WT=AMPSQR/EIKONALFACTOR
c           
c          SANC_WT=(1-8*EMU*XPH*(1-COSTHG*BETA)*     
c     $           (MCHREN+2*XPH*SQRT(MPASQR))/
c     $            MPASQR**2/(1-MCHREN/MPASQR)/(4-MCHREN/MPASQR)) 

      RETURN
      END


      subroutine SANC_INIT1(QB0,QF20,MF20,MB0)
      IMPLICIT NONE
      DOUBLE PRECISION QB0,QF20,MF10,MF20,MB0 
      INTEGER           mcLUN
      DOUBLE PRECISION  pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af
      COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN
      QB =QB0
      QF2=QF20
      MF2=MF20
      MB =MB0
      RETURN
      END

      subroutine SANC_INIT(ALPHA,PHLUN)
      implicit none
      DOUBLE PRECISION  ALPHA
      DOUBLE PRECISION  spV(0:3),bet(0:3)
      INTEGER           mcLUN
      DOUBLE PRECISION  pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af
      COMMON /Kleiss_Stirling/spV,bet
      COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN   
      INTEGER PHLUN

      DOUBLE PRECISION SANC_MC_INIT
      SAVE SANC_MC_INIT
      DATA SANC_MC_INIT /-123456789D0/
c...       Initialization of the W->l\nu\gamma 
c...       decay Matrix Element parameters 
           IF (SANC_MC_INIT.EQ.-123456789D0) THEN
              SANC_MC_INIT=1D0

              PI=4*datan(1d0)
              QF1=0d0                           ! neutrino charge
              MF1=1d-10                         ! newutrino mass
              VF=1d0                            ! V&A couplings
              AF=1d0
              alphaI=1d0/ALPHA
              CW=0.881731727d0                  ! Weak Weinberg angle
              SW=0.471751166d0
           

c...          An auxilary K&S vectors
              bet(0)= 1d0
              bet(1)= 0.0722794881816159d0
              bet(2)=-0.994200045099866d0
              bet(3)= 0.0796363353729248d0 

              spV(0)= 0d0 
              spV(1)= 7.22794881816159d-002
              spV(2)=-0.994200045099866d0     
              spV(3)= 7.96363353729248d-002

              mcLUN = PHLUN
           ENDIF 

      return
      end
