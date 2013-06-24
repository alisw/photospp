
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
