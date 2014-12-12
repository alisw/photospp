

CC======================================================================
CC                                                                    CC
CC WRITTEN BY E. RICHTER- WAS                                         CC
CC                                                                    CC
CC    NOVEMBER 1992                                                   CC 
CC                                                                    CC  
CC======================================================================
C THIS IS UNPUBLISHED VERSION OF PROGRAM WHICH COULD BE APPLIED TO TEST
C LEPTONIC SPECTRUM FOR JETSET7.3+PHOTOS WITH SEMIANALYTICAL FORMULA 
C FOR HISTOGRAMING PRIVATE LIBBRARY GLIBK.FOR IS USED 
C (can be replaced by HBOOK)
C TWO DECAY CHANNELS COULD BE TESTED
C    -----TAU --> ELECTRON, NEUTRINO, NEUTRINO
C    ----- B  --> D0, ELECTRON, NEUTRINO 
C IN DECAYING PARTICLE REST FRAME AND LABORATORY FRAME
C VARIABLE IN COMMON /KEYS/ SWITCH OFF/ON CHANNELS
C    -- KEYTAU = 1   TAU DECAY IN REST FRAME
C    -- KEYTAU = 2   TAU DECAY IN LABORATORY FRAME    
C    -- KEYB   = 1   B   DECAY IN REST FRAME
C    -- KEYB   = 2   B   DECAY IN LABORATORY FRAME  


      PROGRAM PHOSEM
C     ***************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON  /    / BLAN(100000)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /KEYS/ KEYTAU,KEYB
      REAL*4  EBEAM
      COMMON /BEAM/ EBEAM
      COMMON /PHNUM/ IEV


c.... defining matrix for GBOOK
      CALL GLIMIT(100000)
c....define option
      KEYB=0 
      KEYTAU=3 
c.... definning file for printed histograms 
      IF(KEYB.EQ.1) NOUT=16
      IF(KEYB.EQ.2) NOUT=26
      IF(KEYTAU.EQ.1) NOUT=36
      IF(KEYTAU.EQ.2) NOUT=46
      IF(KEYTAU.EQ.3) NOUT=56
      CALL GOUTPU(NOUT)
C....initialization of main subroutine
      CALL RADIATE(-1)
c.... preaparing LUND for  desired event
      CALL PRELUND
C    
C....generating mode
c
      DO 10 IEV=1,15 00
      if(mod(iev,50).eq.1) write(6,*) 'event no=',iev
      if(mod(iev,10).eq.1) write(56,*) 'event no=',iev

      EBEAM = 45.0
cccc for B+- decays
      IF(KEYTAU.EQ.1) CALL LU1ENT(0,15,0.,0.,0.)
      IF(KEYTAU.EQ.2) CALL LU1ENT(0,15,EBEAM,0.,0.)
      IF(KEYTAU.EQ.3) CALL LU2ENT(0,5,-5,45.)

cccc for tau decays 
      IF(KEYB.EQ.1) CALL LU1ENT(0,521,0.,0.,0.)
      IF(KEYB.EQ.2) CALL LU1ENT(0,521,EBEAM,0.,0.)

cc      IF(IEV.LT.5) CALL LULIST(1)
      CALL RADIATE(0)
   10 CONTINUE
C
C....post generation mode
      CALL RADIATE(1)

C ------------WRITING HISTOS ON THE DISK ------------------------
      IF(KEYB.EQ.1) NOUTH=10
      IF(KEYB.EQ.2) NOUTH=20
      IF(KEYTAU.EQ.1) NOUTH=30
      IF(KEYTAU.EQ.2) NOUTH=40

      CALL GRFILE(NOUTH,' ','N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)
C ------------THE END OF HISTO WRITING -------------------------

      END
      
      SUBROUTINE PRELUND
C     *********************
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON /KEYS/ KEYTAU,KEYB

c....uotput will be writen on NOUT
      MSTU(11) = NOUT
c....
cccc      CALL LULIST(12)
c....
      IF(KEYB.NE.0) THEN
cccc      for B+- decays
        DO 100 J=309,316
 100    MDME(J,1)=0
        MDCY(105,1)=0
        PARJ(11) = 0.
      ENDIF
      IF(KEYTAU.NE.0.AND.KEYTAU.NE.3) THEN     
cccc      for tau decays
        DO 200 J=79,117
 200    MDME(J,1)=0

      ENDIF
      END
      SUBROUTINE RADIATE(MODE)
C     *************************
C     this is main subroutine, where histograming routines are called
C     and photons are generated throuhg CALL PHOTOS
C
      PARAMETER (LJNPAR=4000)                                       
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /LUJETS/ N7LU,K7LU(LJNPAR,5),P7LU(LJNPAR,5),V7LU(LJNPAR,5)
      PARAMETER(NMXHEP=2000)  
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),        
     &       JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),         
     &       VHEP(4,NMXHEP) 
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      COMMON /KEYS/ KEYTAU,KEYB
C
      DOUBLE PRECISION XRIAD,TMIN,TMAX

      IF(MODE.EQ.-1) THEN
C     =======================
c histograms initialization for photon number
      TMAX = 10.D0
      TMIN = 0.D0
      NBIN = 10
      CALL GBOOK1(101,'NUMBER OF RADIATED PHOTONS$ ',NBIN,TMIN,TMAX)       
c
C initialization of photos
      CALL PHOINI
C initialization of histograms for leptons spectrum
      IF(KEYB.EQ.1) CALL DISTRY(-1)
      IF(KEYB.EQ.2) CALL DISLAB(-1)
      IF(KEYTAU.EQ.1) CALL DTAU(-1)
      IF(KEYTAU.EQ.2) CALL DTLAB(-1)
C
      IEV =0 
      ELSEIF(MODE.EQ.0) THEN
C     =========================
      IEV = IEV+1
C histograming lepton spectrum before PHOTOS
      IF(KEYB.EQ.1) CALL DISTRY(0)
      IF(KEYB.EQ.2) CALL DISLAB(0)
      IF(KEYTAU.EQ.1) CALL DTAU(0)
      IF(KEYTAU.EQ.2) CALL DTLAB(0)
      NPOZ = N7LU
C conversion from common /LUJETS/ to common /HEPEVT/
      CALL LUHEPC(1)
c calling photos
      CALL PHOTOS(1)
ccc      IF(IEV.LT.2)CALL LULIST(1)
C conversion from common /HEPETV/ to common /LUJETS/
      CALL LUHEPC(2)
C store number of radiated photons
      XRIAD = (N7LU-NPOZ)*1D0
      if(xriad.gt.6.2d0.AND.IEV.LT.100101) call lulist(1)
      CALL GF1(101,XRIAD,1.D0)                              
C histograming lepton spectrum after PHOTOS
      IF(KEYB.EQ.1) CALL DISTRY(1)
      IF(KEYB.EQ.2) CALL DISLAB(1)
      IF(KEYTAU.EQ.1) CALL DTAU(1)
      IF(KEYTAU.EQ.2) CALL DTLAB(1)
c
      ELSEIF(MODE.EQ.1) THEN
C     =========================
c  printing few histograms
      CALL GPRINT(101)
      IF(KEYB.EQ.1) CALL DISTRY(2)
      IF(KEYB.EQ.2) CALL DISLAB(2)
      IF(KEYTAU.EQ.1) CALL DTAU(2)
      IF(KEYTAU.EQ.2) CALL DTLAB(2)
C 

      ENDIF

      END
C
    
      SUBROUTINE DTAU(MODE)
C     *************************
C      tau decays in tau rest frame
C
      PARAMETER (LJNPAR=4000)                                        
      COMMON /LUJETS/ N7LU,K7LU(LJNPAR,5),P7LU(LJNPAR,5),V7LU(LJNPAR,5) 
C
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /KEYS/ KEYTAU,KEYB
      DOUBLE PRECISION  TMAX, TMIN, FACT
      DOUBLE PRECISION XLEPT,YLEPT
      DIMENSION XLEPT(20),YLEPT(20)
      LOGICAL  NEUTRT, ELECTRO, NEUTRE
      DOUBLE PRECISION WT,SPEC1

      IF(MODE.EQ.-1) THEN
C     =======================
c histograms initialization for leptons momentum spectrum
      TMAX = 1.D0
      TMIN = 0.D0
      NBIN = 50
      CALL GBOOK1(10,'EL/EMAX BEFORE PHOTOS $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(10,'ERRO')
      CALL GBOOK1(20,'EL/EMAX AFTER PHOTOS  $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(20,'ERRO')
      CALL GBOOK1(30,'EL/EMAX CORRECTION    $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(30,'ERRO')
      CALL GBOOK1(50,'EL/EMAX alpha/analit    $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(50,'ERRO')
      NENTRY = 0
C
      ELSEIF(MODE.EQ.0) THEN
C     =======================
      ILEPT = 0
      DO 101 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 101
C PARENT= TAU                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.15) GO TO 101 
C CHECK DECAY TREE (TAU +-  --> neutrino  electron  neutrino)
      NEUTRE =.FALSE.
      NEUTRT =.FALSE.
      DO 111 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11) THEN
        ELECTRO=.TRUE.
        IELE = ND
      ENDIF
      IF(IDD.EQ.12) THEN
        NEUTRE =.TRUE.
        INEUE = ND
      ENDIF
      IF(IDD.EQ.16) THEN
        NEUTRT =.TRUE.
        INEUT    = ND
      ENDIF
 111  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRE.AND.NEUTRT) THEN
        NENTRY = NENTRY + 1
        ILEPT = ILEPT+1
        XMASST = SQRT(abs( P7LU(IPAR,4)**2-P7LU(IPAR,3)**2
     &                 -P7LU(IPAR,2)**2-P7LU(IPAR,1)**2 ))
        ELEPT = (P7LU(I,4)*P7LU(IPAR,4)-P7LU(I,1)*P7LU(IPAR,1) 
     &         - P7LU(I,2)*P7LU(IPAR,2)-P7LU(I,3)*P7LU(IPAR,3))
     &         /XMASST
        EMAX  = XMASST/2.
        X  = ELEPT/EMAX
C store lepton momenta and its order in jets tree
        XLEPT(ILEPT) = DBLE(X)
      ENDIF
  101 CONTINUE
C
      ELSEIF(MODE.EQ.1) THEN
C     =========================
      ILEPT = 0
      DO 201 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 201
C PARENT= TAU                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.15) GO TO 201 
C CHECK DECAY TREE (TAU +-  --> neutrino electron neutrino)                                
      NEUTRE=.FALSE.
      NEUTRT =.FALSE.
      IMES  = 0
      DO 211 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11)  ELECTRO=.TRUE.
      IF(IDD.EQ.12)  NEUTRE =.TRUE.
      IF(IDD.EQ.16) THEN
        NEUTRT =.TRUE.
        INEUT    = ND
      ENDIF
 211  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRT.AND.NEUTRE) THEN
        ILEPT = ILEPT+1
        XMASST = SQRT(abs( P7LU(IPAR,4)**2-P7LU(IPAR,3)**2
     &                 -P7LU(IPAR,2)**2-P7LU(IPAR,1)**2 ))
        ELEPT = (P7LU(I,4)*P7LU(IPAR,4)-P7LU(I,1)*P7LU(IPAR,1) 
     &         - P7LU(I,2)*P7LU(IPAR,2)-P7LU(I,3)*P7LU(IPAR,3))
     &         /XMASST
        EMAX  = XMASST/2.
        Y  = ELEPT/EMAX
C store lepton momenta and its order in jets tree
        YLEPT(ILEPT) = DBLE(Y)
C store the effect of x distribution orginated from radiation
C distribution without radiation
        CALL GF1(10,XLEPT(ILEPT),1D0)
C distribution with radiation 
        CALL GF1(20,YLEPT(ILEPT),1D0) 
C difference, parralel events used
        IF((XLEPT(ILEPT)-YLEPT(ILEPT)).LT.1D-3) THEN
          CALL GF1(30,XLEPT(ILEPT),0D0) 
          CALL GF1(30,YLEPT(ILEPT),0D0) 
        ELSE
          CALL GF1(30,XLEPT(ILEPT),-1D0) 
          CALL GF1(30,YLEPT(ILEPT), 1D0) 
        ENDIF
C event weighted by exact order alpha formula, histogram should be flat
        WT = 1D0/SPEC1(YLEPT(ILEPT))
        CALL GF1(50,YLEPT(ILEPT),WT) 
      ENDIF
  201 CONTINUE
C
      ELSEIF(MODE.EQ.2) THEN
C     =========================
      CALL GPRINT(10)
      CALL GPRINT(20)
      CALL GPRINT(30)
c.... normalize area under histogram to unity    
      FACT = DBLE(1.0/NENTRY)*NBIN
      CALL GOPERA(10,'+',10,10,FACT,0D0)
      CALL GOPERA(20,'+',20,20,FACT,0D0)
      CALL GOPERA(30,'+',30,30,FACT,0D0) 
      CALL GOPERA(50,'+',50,50,FACT,0D0)
c.... normalize ratio to 1 
      FACT = DBLE(1.D0/NENTRY)
      CALL GOPERA(50,'+',50,50,FACT,0D0)
      CALL GPRINT(50)
C
      ENDIF

      END
      FUNCTION SPEC1 (XX) 
C     **********************   
C     analytical formula for Order alpha lepton spectrum                             C                      
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DATA PI /3.141592653589793238462643D0/
      DATA ALFINV /137.03604D0/


      AMTAU = 1.7841D0
      AMEL = 0.511D-3
      X=XX                                                                      
      L=LOG(AMTAU/AMEL)                                                         
      BORN  = 2D0*X**2*(3-2*X)                                                      
      G=4*DILOG(X)-2D0/3D0*PI**2+3*L-4             +                            
     $  2*LOG(X)       *(3*LOG(1-X)-2*LOG(X)-2*L+1)+                            
     $  2*(2*L-1-1/X)  *LOG(1-X)                   +                            
     $  6*(1-X)/(3-2*X)*LOG(X)                     +                            
     $  (1-X)/3/X**2/(3-2*X)*                                                   
     $  ((5+17*X-34*X**2)*(L+LOG(X))-22*X+34*X**2 )                             
      SPEC1=BORN*(1D0+1D0/2D0/PI/ALFINV*G)
                                      
      END
      SUBROUTINE DTLAB(MODE)
C     *************************
C      tau decays in laboratory frame
C
      PARAMETER (LJNPAR=4000)                                        
      COMMON /LUJETS/ N7LU,K7LU(LJNPAR,5),P7LU(LJNPAR,5),V7LU(LJNPAR,5) 
C
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /KEYS/ KEYTAU,KEYB
      COMMON /BEAM/ EBEAM   
      DOUBLE PRECISION  TMAX, TMIN, FACT
      DOUBLE PRECISION XLEPT,YLEPT
      DIMENSION XLEPT(20),YLEPT(20)
      LOGICAL  NEUTRT, ELECTRO, NEUTRE
      DOUBLE PRECISION WT,SPEC1

      IF(MODE.EQ.-1) THEN
C     =======================
c histograms initialization for leptons momentum spectrum
      TMAX = 1.D0
      TMIN = 0.D0
      NBIN = 50
      CALL GBOOK1(510,'EL/EMAX BEFORE PHOTOS $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(510,'ERRO')
      CALL GBOOK1(520,'EL/EMAX AFTER PHOTOS  $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(520,'ERRO')
      CALL GBOOK1(530,'EL/EMAX CORRECTION    $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(530,'ERRO')
      NENTRY = 0
C
      ELSEIF(MODE.EQ.0) THEN
C     =======================
      ILEPT = 0
      DO 101 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 101
C PARENT= TAU                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.15) GO TO 101 
C CHECK DECAY TREE (TAU +-  --> neutrino  electron  neutrino)
      NEUTRE =.FALSE.
      NEUTRT =.FALSE.
      DO 111 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11) THEN
        ELECTRO=.TRUE.
        IELE = ND
      ENDIF
      IF(IDD.EQ.12) THEN
        NEUTRE =.TRUE.
        INEUE = ND
      ENDIF
      IF(IDD.EQ.16) THEN
        NEUTRT =.TRUE.
        INEUT    = ND
      ENDIF
 111  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRE.AND.NEUTRT) THEN
        NENTRY = NENTRY + 1
        ILEPT = ILEPT+1
        ELEPT =  P7LU(I,4) 
        X  = ELEPT/EBEAM
C store lepton momenta and its order in jets tree
        XLEPT(ILEPT) = DBLE(X)
      ENDIF
  101 CONTINUE
C
      ELSEIF(MODE.EQ.1) THEN
C     =========================
      ILEPT = 0
      DO 201 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 201
C PARENT= TAU                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.15) GO TO 201 
C CHECK DECAY TREE (TAU +-  --> neutrino electron neutrino)                                
      NEUTRE=.FALSE.
      NEUTRT =.FALSE.
      IMES  = 0
      DO 211 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11)  ELECTRO=.TRUE.
      IF(IDD.EQ.12)  NEUTRE =.TRUE.
      IF(IDD.EQ.16) THEN
        NEUTRT =.TRUE.
        INEUT    = ND
      ENDIF
 211  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRT.AND.NEUTRE) THEN
        ILEPT = ILEPT+1
        ELEPT = P7LU(I,4)
        Y  = ELEPT/EBEAM
C store lepton momenta and its order in jets tree
        YLEPT(ILEPT) = DBLE(Y)
C store the effect of x distribution orginated from radiation
C lepton distribution without radiation
        CALL GF1(510,XLEPT(ILEPT),1D0) 
C lepton distribution with radiation
        CALL GF1(520,YLEPT(ILEPT),1D0)
C difference parralel events used 
        IF((XLEPT(ILEPT)-YLEPT(ILEPT)).LT.1D-3) THEN
          CALL GF1(530,XLEPT(ILEPT),0D0) 
          CALL GF1(530,YLEPT(ILEPT),0D0) 
        ELSE
          CALL GF1(530,XLEPT(ILEPT),-1D0) 
          CALL GF1(530,YLEPT(ILEPT), 1D0) 
        ENDIF
      ENDIF
  201 CONTINUE
C
      ELSEIF(MODE.EQ.2) THEN
C     =========================
      CALL GPRINT(510)
      CALL GPRINT(520)
      CALL GPRINT(530)
C   normalize area under histogram to unity   
      FACT = DBLE(1.0/NENTRY)*NBIN
      CALL GOPERA(510,'+',510,510,FACT,0D0)
      CALL GOPERA(520,'+',520,520,FACT,0D0)
      CALL GOPERA(530,'+',530,530,FACT,0D0) 
C
      ENDIF

      END

      SUBROUTINE DISTRY(MODE)
C     *************************
C      B decay in B rest frame
C
      PARAMETER (LJNPAR=4000)                                        
      COMMON /LUJETS/ N7LU,K7LU(LJNPAR,5),P7LU(LJNPAR,5),V7LU(LJNPAR,5) 
C
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /KEYS/ KEYTAU,KEYB
      DOUBLE PRECISION  TMAX, TMIN, FACT,FACC,WMAT,XANA,XLUND,WMAT1
      DOUBLE PRECISION XLEPT,YLEPT
      DIMENSION XLEPT(20),YLEPT(20)
      LOGICAL  MESON, ELECTRO, NEUTRO
c function
      FOUR(I,J)=P7LU(I,4)*P7LU(J,4)-P7LU(I,3)*P7LU(J,3)
     $         -P7LU(I,2)*P7LU(J,2)-P7LU(I,1)*P7LU(J,1)

      IF(MODE.EQ.-1) THEN
C     =======================
c histograms initialization for leptons momentum spectrum
      TMAX = 1.D0
      TMIN = 0.D0
      NBIN = 50
      CALL GBOOK1(150,'EL/EMAX BEFORE PHOTOS W$',NBIN,TMIN,TMAX)       
      CALL GIDOPT(150,'ERRO')
      CALL GBOOK1(250,'EL/EMAX AFTER PHOTOS  W$',NBIN,TMIN,TMAX)       
      CALL GIDOPT(250,'ERRO')
      CALL GBOOK1(350,'EL/EMAX CORRECTION    W$',NBIN,TMIN,TMAX)       
      CALL GIDOPT(350,'ERRO')

      NENTRY = 0
C
      ELSEIF(MODE.EQ.0) THEN
C     =======================
      ILEPT = 0
      DO 101 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 101
C PARENT= B+- HADRON (MESON)                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.521) GO TO 101 
C CHECK DECAY TREE (B +-  --> D0  electron  neutrino)
      NEUTRO=.FALSE.
      MESON =.FALSE.
      IMES  = 0
      DO 111 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11) THEN
        ELECTRO=.TRUE.
        IELE = ND
      ENDIF
      IF(IDD.EQ.12) THEN
        NEUTRO =.TRUE.
        INEU = ND
      ENDIF
      IF(IDD.EQ.421) THEN
        MESON =.TRUE.
        IMES    = ND
      ENDIF
 111  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRO.AND.MESON) THEN
        NENTRY = NENTRY + 1
        ILEPT = ILEPT+1
        XMASSB = SQRT(abs( P7LU(IPAR,4)**2-P7LU(IPAR,3)**2
     &                 -P7LU(IPAR,2)**2-P7LU(IPAR,1)**2 )) 
        XMASSD = SQRT(abs( P7LU(IMES,4)**2-P7LU(IMES,3)**2
     &                 -P7LU(IMES,2)**2-P7LU(IMES,1)**2 ))
        ELEPT = (P7LU(I,4)*P7LU(IPAR,4)-P7LU(I,1)*P7LU(IPAR,1) 
     &         - P7LU(I,2)*P7LU(IPAR,2)-P7LU(I,3)*P7LU(IPAR,3))
     &         /XMASSB
        EMAX  = (XMASSB**2-XMASSD**2)/2./XMASSB
        X  = ELEPT/EMAX
C store lepton momenta and its order in jets tree
        XLEPT(ILEPT) = DBLE(X)
      ENDIF
c calculate weigh for matrix element to modified spectrum from JETSET
      XLUND = DBLE(FOUR(IPAR,INEU) * FOUR(IMES,IELE)) 
      XANA  = DBLE(
     $      2.* FOUR(IPAR,IELE )* FOUR(IPAR,INEU)  
     $        - FOUR(IPAR,IPAR )* FOUR(IELE,INEU) 
     $            )
      WMAT  = XANA/XLUND

  101 CONTINUE
C
      ELSEIF(MODE.EQ.1) THEN
C     =========================
      ILEPT = 0
      DO 201 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 201
C PARENT= B+- HADRON (MESON)                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.521) GO TO 201 
C CHECK DECAY TREE (B +-  --> D0 electron neutrino)                                
      NEUTRO=.FALSE.
      MESON =.FALSE.
      IMES  = 0
      DO 211 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11)  ELECTRO=.TRUE.
      IF(IDD.EQ.12)  NEUTRO =.TRUE.
      IF(IDD.EQ.421) THEN
        MESON =.TRUE.
        IMES    = ND
      ENDIF
 211  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRO.AND.MESON) THEN
        ILEPT = ILEPT+1
        XMASSB = SQRT(abs( P7LU(IPAR,4)**2-P7LU(IPAR,3)**2
     &                 -P7LU(IPAR,2)**2-P7LU(IPAR,1)**2 )) 
        XMASSD = SQRT(abs( P7LU(IMES,4)**2-P7LU(IMES,3)**2
     &                 -P7LU(IMES,2)**2-P7LU(IMES,1)**2 )) 
        ELEPT = (P7LU(I,4)*P7LU(IPAR,4)-P7LU(I,1)*P7LU(IPAR,1) 
     &         - P7LU(I,2)*P7LU(IPAR,2)-P7LU(I,3)*P7LU(IPAR,3))
     &         /XMASSB
        EMAX  = (XMASSB**2-XMASSD**2)/2./XMASSB
        Y  = ELEPT/EMAX
C store lepton momenta and its order in jets tree
        YLEPT(ILEPT) = DBLE(Y)
C store the effect of x distribution orginated from radiation
C distribution without radiation, JETSET matrix element corrected
        CALL GF1(150,XLEPT(ILEPT),WMAT) 
C distribution with radiative correction, JETSET matrix element corrected
        CALL GF1(250,YLEPT(ILEPT),WMAT)
C difference, parralel events used
        IF((XLEPT(ILEPT)-YLEPT(ILEPT)).LT.1D-3) THEN
          CALL GF1(350,XLEPT(ILEPT),0D0) 
          CALL GF1(350,YLEPT(ILEPT),0D0) 
        ELSE
          CALL GF1(350,XLEPT(ILEPT),-WMAT) 
          CALL GF1(350,YLEPT(ILEPT), WMAT) 
        ENDIF
      ENDIF

  201 CONTINUE
C
      ELSEIF(MODE.EQ.2) THEN
C     =========================
      CALL GPRINT(150)
      CALL GPRINT(250)
      CALL GPRINT(350)
C weighet event used, only histogram 150 normalized exactly to unity
      CALL LICFAC(150,FACC,NBIN)
      FACT=FACC 
      CALL RENHIE(150,FACT,NBIN)
      CALL RENHIE(250,FACT,NBIN)
      CALL RENHIE(350,FACT,NBIN)   
      FACT = 1D0*NBIN
      CALL GOPERA(150,'+',150,150,FACT,0D0)
      CALL GOPERA(250,'+',250,250,FACT,0D0)
      CALL GOPERA(350,'+',350,350,FACT,0D0)
C
      ENDIF

      END
      SUBROUTINE DISLAB(MODE)
C     *************************
C       B decay in laboratory frame
C
      PARAMETER (LJNPAR=4000)                                        
      COMMON /LUJETS/ N7LU,K7LU(LJNPAR,5),P7LU(LJNPAR,5),V7LU(LJNPAR,5) 
C
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /KEYS/ KEYTAU,KEYB
      COMMON /BEAM/ EBEAM
      DOUBLE PRECISION  TMAX, TMIN, FACT,FACC,WMAT,XANA,XLUND
      DOUBLE PRECISION XLEPT,YLEPT
      DIMENSION XLEPT(20),YLEPT(20)
      LOGICAL  MESON, ELECTRO, NEUTRO
c function
      FOUR(I,J)=P7LU(I,4)*P7LU(J,4)-P7LU(I,3)*P7LU(J,3)
     $         -P7LU(I,2)*P7LU(J,2)-P7LU(I,1)*P7LU(J,1)

      IF(MODE.EQ.-1) THEN
C     =======================
c histograms initialization for leptons momentum spectrum
      TMAX = 1.D0
      TMIN = 0.D0
      NBIN = 50
      CALL GBOOK1(10100,'EL/EMAX BEFORE PHOTOS $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(10100,'ERRO')
      CALL GBOOK1(10200,'EL/EMAX AFTER PHOTOS  $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(10200,'ERRO')
      CALL GBOOK1(10300,'EL/EMAX CORRECTION    $',NBIN,TMIN,TMAX)       
      CALL GIDOPT(10300,'ERRO')
      NENTRY = 0
C
      ELSEIF(MODE.EQ.0) THEN
C     =======================
      ILEPT = 0
      DO 101 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 101
C PARENT= B+- HADRON (MESON)                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.521) GO TO 101 
C CHECK DECAY TREE (B +-  --> D0  electron  neutrino)
      NEUTRO=.FALSE.
      MESON =.FALSE.
      IMES  = 0
      DO 111 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11) THEN
        ELECTRO=.TRUE.
        IELE = ND
      ENDIF
      IF(IDD.EQ.12) THEN
        NEUTRO =.TRUE.
        INEU = ND
      ENDIF
      IF(IDD.EQ.421) THEN
        MESON =.TRUE.
        IMES    = ND
      ENDIF
 111  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRO.AND.MESON) THEN
        NENTRY = NENTRY + 1
        ILEPT = ILEPT+1
        ELEPT = P7LU(I,4)
        X  = ELEPT/EBEAM
C store lepton momenta and its order in jets tree
        XLEPT(ILEPT) = DBLE(X)
      ENDIF
c calculate weigh for matrix element--> JETSET matrix element corrected
      XLUND = DBLE(FOUR(IPAR,INEU) * FOUR(IMES,IELE))
      XANA  = DBLE(
     $      2.* FOUR(IPAR,IELE )* FOUR(IPAR,INEU)  
     $        - FOUR(IPAR,IPAR )* FOUR(IELE,INEU) 
     $            )
      WMAT = XANA/XLUND

  101 CONTINUE
C
      ELSEIF(MODE.EQ.1) THEN
C     =========================
      ILEPT = 0
      DO 201 I=3,N7LU                                               
      IDENT=IABS(K7LU(I,2))                                         
C LOOK FOR ELECTRON                                           
      IF(IDENT.NE.11) GO TO 201
C PARENT= B+- HADRON (MESON)                    
      IPAR=K7LU(I,3)                                                
      IDPAR=IABS(K7LU(IPAR,2))                                      
      IF(IDPAR.NE.521) GO TO 201 
C CHECK DECAY TREE (B +-  --> D0 electron neutrino)                                
      NEUTRO=.FALSE.
      MESON =.FALSE.
      IMES  = 0
      DO 211 ND = K7LU(IPAR,4),K7LU(IPAR,5)
      IDD = IABS(K7LU(ND,2))
      IF(IDD.EQ.11)  ELECTRO=.TRUE.
      IF(IDD.EQ.12)  NEUTRO =.TRUE.
      IF(IDD.EQ.421) THEN
        MESON =.TRUE.
        IMES    = ND
      ENDIF
 211  CONTINUE
C if decay tree accepted 
      IF(ELECTRO.AND.NEUTRO.AND.MESON) THEN
        ILEPT = ILEPT+1
        ELEPT = P7LU(I,4)
        Y  = ELEPT/EBEAM
C store lepton momenta and its order in jets tree
        YLEPT(ILEPT) = DBLE(Y)
C store the effect of x distribution orginated from radiation
C no radiative correction to distribution
        CALL GF1(10100,XLEPT(ILEPT),WMAT)
C radiative correction to distribution included,
C  JETSET matrix element corrected 
        CALL GF1(10200,YLEPT(ILEPT),WMAT) 
C difference, parralel events used
        IF((XLEPT(ILEPT)-YLEPT(ILEPT)).LT.1D-3) THEN
          CALL GF1(10300,XLEPT(ILEPT),0D0) 
          CALL GF1(10300,YLEPT(ILEPT),0D0) 
        ELSE
          CALL GF1(10300,XLEPT(ILEPT),-WMAT) 
          CALL GF1(10300,YLEPT(ILEPT), WMAT) 
        ENDIF
      ENDIF
  201 CONTINUE
C
      ELSEIF(MODE.EQ.2) THEN
C     =========================
      CALL GPRINT(10100)
      CALL GPRINT(10200)
      CALL GPRINT(10300)
      CALL LICFAC(10100,FACC,NBIN)
      FACT=FACC 
      CALL RENHIE(10100,FACT,NBIN)
      CALL RENHIE(10200,FACT,NBIN)
      CALL RENHIE(10300,FACT,NBIN)
      FACT = 1D0*NBIN
      CALL GOPERA(10100,'+',10100,10100,FACT,0D0)
      CALL GOPERA(10200,'+',10200,10200,FACT,0D0)
      CALL GOPERA(10300,'+',10300,10300,FACT,0D0)
C
      ENDIF

      END

      SUBROUTINE RENHIE(ID,FACT,NB)
C     ****************************
C errors taken into account
C     INTRODUCES HISTOGRAM NORMALISATION
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(400),XE(400)
      FAC=FACT
      CALL GUNPAK(ID, X,'    ',1)
      CALL GUNPAK(ID,XE,'ERRO',1)
      CALL GRESET(ID,' ')
      SUM=0D0
      DO 10 I=1,NB
   10 SUM=SUM+X(I)
      IF(SUM.EQ.0D0) RETURN
      IF(FAC.EQ.0D0) FAC=1D0/SUM
      DO 20 I=1,NB
      X(I) = X(I)*FAC
      XE(I)=XE(I)*FAC
   20 CONTINUE
      CALL GPAK( ID,X )
      CALL GPAKE(ID,XE)
      
      END
      SUBROUTINE LICFAC(ID,FACT,NB)
C     ****************************
C errors taken into account
C     INTRODUCES HISTOGRAM NORMALISATION
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(400),XE(400)

      CALL GUNPAK(ID, X,'    ',1)
      CALL GUNPAK(ID,XE,'ERRO',1)
      CALL GRESET(ID,' ')
      SUM=0D0
      DO 10 I=1,NB
   10 SUM=SUM+X(I)
      IF(SUM.EQ.0D0) RETURN
      FACT=1D0/SUM

      CALL GPAK( ID,X )
      CALL GPAKE(ID,XE)
      
      END

      DOUBLE PRECISION FUNCTION DILOG(X)
C-------------------------------------------- REMARKS ---------------
C DILOGARITHM FUNCTION: DILOG(X)=INT( -LN(1-Z)/Z ) , 0 < Z < X .
C THIS IS THE CERNLIB VERSION.
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Z=-1.644934066848226D0
      IF(X .LT.-1.D0) GO TO 1
      IF(X .LE. 0.5D0) GO TO 2
      IF(X .EQ. 1.D0) GO TO 3
      IF(X .LE. 2.D0) GO TO 4
      Z=3.289868133696453D0
    1 T=1.D0/X
      S=-0.5D0
      Z=Z-0.5D0*DLOG(DABS(X))**2
      GO TO 5
    2 T=X
      S=0.5D0
      Z=0.D0
      GO TO 5
    3 DILOG=1.644934066848226D0
      RETURN
    4 T=1.D0-X
      S=-0.5D0
      Z=1.644934066848226D0-DLOG(X)*DLOG(DABS(T))
    5 Y=2.666666666666667D0*T+0.666666666666667D0
      B=      0.000000000000001D0
      A=Y*B  +0.000000000000004D0
      B=Y*A-B+0.000000000000011D0
      A=Y*B-A+0.000000000000037D0
      B=Y*A-B+0.000000000000121D0
      A=Y*B-A+0.000000000000398D0
      B=Y*A-B+0.000000000001312D0
      A=Y*B-A+0.000000000004342D0
      B=Y*A-B+0.000000000014437D0
      A=Y*B-A+0.000000000048274D0
      B=Y*A-B+0.000000000162421D0
      A=Y*B-A+0.000000000550291D0
      B=Y*A-B+0.000000001879117D0
      A=Y*B-A+0.000000006474338D0
      B=Y*A-B+0.000000022536705D0
      A=Y*B-A+0.000000079387055D0
      B=Y*A-B+0.000000283575385D0
      A=Y*B-A+0.000001029904264D0
      B=Y*A-B+0.000003816329463D0
      A=Y*B-A+0.000014496300557D0
      B=Y*A-B+0.000056817822718D0
      A=Y*B-A+0.000232002196094D0
      B=Y*A-B+0.001001627496164D0
      A=Y*B-A+0.004686361959447D0
      B=Y*A-B+0.024879322924228D0
      A=Y*B-A+0.166073032927855D0
      A=Y*A-B+1.935064300869969D0
      DILOG=S*T*(A-B)+Z
      END
     

C===========================================================
C implementation of electron positron pairs
C===========================================================
      SUBROUTINE PHOPAR(IPARR,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   Photon radiation in decays
C.
C.    Purpose:  e+e- pairs  are  generated  in
C.              the decay of the IPPAR-th particle in the HEP-like
C.              common /PHOEVT/.  Radiation takes place from one
C.              of the charged daughters of the decaying particle IPPAR
C.           
C.           
C.
C.    Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
C.                                /PHOEVT/ and the common itself,
C.                                NHEP0 length of the /HEPEVT/ entry
C.                                before starting any activity on this
C.                                IPPAR decay.
C.    Output Parameters:  Common  /HEPEVT/, either  with  or  without a
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
      REAL PHOCHA,MASSUM
      INTEGER IP,IPARR,IPPAR,I,J,NCHARG,NLAST,K
      INTEGER IDABS,IDUM
      INTEGER NMXPHO
      PARAMETER (NMXPHO=100)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER PDMMAX
      PARAMETER (PDMMAX=100)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(PDMMAX)
      INTEGER CHAPOI(PDMMAX)
      INTEGER IREP
      INTEGER NHEP0
      LOGICAL BOOST,JESLI
      INTEGER KEYRND
c-- type of RN generator
      COMMON / RANPAR / KEYRND
      save   / RANPAR /
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
            CALL TRYPAR(JESLI,STRENG,PCHAR,PNEU,PELE,PPOZ)
          IF (JESLI) THEN
           DO I=JDAPHO(1,IP),JDAPHO(2,IP)
            DO K=1,4
              BUF(K)=PPHO(K,I)
            ENDDO
            IF (I.EQ.CHAPOI(NCHARG)) THEN
               CALL PARTRA( 1,BUF)
            ELSE
               CALL PARTRA(-1,BUF)
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
      SUBROUTINE DUMPZ(NOUT)     
*     **********************     
* THIS PRINTS OUT FOUR MOMENTA OF PHOTONS 
* ON OUTPUT UNIT NOUT
      IMPLICIT REAL*8(A-H,O-Z)   
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      COMMON /PAIRST/   PELE(4),PPOZ(4),NPAIRS  
      REAL*8 SUM(4)     
      WRITE(NOUT,*) '=====================DUMPS====================' 
      WRITE(NOUT,3100) ' P2',(P2(K),K=1,4)   
      WRITE(NOUT,3100) ' Q2',(Q2(K),K=1,4)   
      WRITE(NOUT,3100) ' PE',(PELE(K),K=1,4)   
      WRITE(NOUT,3100) ' PP',(PPOZ(K),K=1,4)   
      DO 100 I=1,NPHOT  
  100 WRITE(NOUT,3100) 'PHO',(PHOT(I,K),K=1,4)        
      DO 200 K=1,4      
  200 SUM(K)=P2(K)+Q2(K)+PELE(K)+PPOZ(K)         
      DO 210 I=1,NPHOT  
      DO 210 K=1,4      
  210 SUM(K)=SUM(K)+PHOT(I,K)    
      WRITE(NOUT,3100) 'SUM',(SUM(K),K=1,4)           
 3100 FORMAT(1X,A3,1X,5F18.13)   
      END   
      SUBROUTINE SPAJ(KUDA,PA,PB,PP,PE)     
*     **********************     
* THIS PRINTS OUT FOUR MOMENTA OF PHOTONS 
* ON OUTPUT UNIT NOUT
      IMPLICIT REAL*8(A-H,O-Z)   
      REAL*8 SUM(4),PA(4),PB(4),PP(4),PE(4)
      DATA KLUCZ /0/
      IF (KLUCZ.EQ.0) RETURN
      NOUT=56     
      WRITE(NOUT,*) KUDA,'=====================SPAJ====================' 
      WRITE(NOUT,3100) ' P2',(PA(K),K=1,4)   
      WRITE(NOUT,3100) ' Q2',(PB(K),K=1,4)   
      WRITE(NOUT,3100) ' PE',(PE(K),K=1,4)   
      WRITE(NOUT,3100) ' PP',(PP(K),K=1,4)   
      DO 200 K=1,4      
  200 SUM(K)=PA(K)+PB(K)+PE(K)+PP(K)         
      WRITE(NOUT,3100) 'SUM',(SUM(K),K=1,4)           
      NOUT=16     
      WRITE(NOUT,*) KUDA,'=====================SPAJ====================' 
      WRITE(NOUT,3100) ' P2',(PA(K),K=1,4)   
      WRITE(NOUT,3100) ' Q2',(PB(K),K=1,4)   
      WRITE(NOUT,3100) ' PE',(PE(K),K=1,4)   
      WRITE(NOUT,3100) ' PP',(PP(K),K=1,4)   
      WRITE(NOUT,3100) 'SUM',(SUM(K),K=1,4)           
 3100 FORMAT(1X,A3,1X,5F18.13)   
      END   
   
      SUBROUTINE TRYPAR(JESLI,STRENG,PA,PB,PE,PP)       
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      COMMON  /PARKIN/ 
     $ FI0,FI1,FI2,FI3,FI4,FI5,TH0,TH1,TH3,TH4
     $,PARNEU,PARCH,BPAR,BSTA,BSTB
      REAL*8 PNEUTR(4),PAA(4),PHOT(4),PA(4),PB(4),PE(4),PP(4),PSUM(4)
      REAL*8 VEC(4)                                                
      REAL*8 RRR(8) 
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
      FI0  =ANGFI(X1,X2)           
      X1 = PSUM(3)                                                 
      X2 = SQRT(PSUM(1)**2+PSUM(2)**2)                            
      TH0  =ANGXY(X1,X2) 
      CALL SPAJ(-2,PNEUTR,PAA,PP,PE)    
      CALL LORTRA(3,-FI0,PNEUTR,VEC,PAA,PP,PE)
      CALL LORTRA(2,-TH0,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD3(-FI0,PSUM,PSUM)
      CALL ROTOD2(-TH0,PSUM,PSUM)
      BSTA=(PSUM(4)-PSUM(3))/SQRT(PSUM(4)**2-PSUM(3)**2)
      BSTB=(PSUM(4)+PSUM(3))/SQRT(PSUM(4)**2-PSUM(3)**2)
      CALL LORTRA(1,BSTA,PNEUTR,VEC,PAA,PP,PE)
      CALL SPAJ(-1,PNEUTR,PAA,PP,PE)                                   
      AMNE=AMAST(PNEUTR)                                              
      AMCH=AMAST(PAA)                                                 
      IF (AMCH.LT.0.0D0) AMCH=AMEL                                   
      IF (AMNE.LT.0.D0) AMNE=0.0
      AMTO =PAA(4)+PNEUTR(4)
      CALL VARRAN(RRR,8)           
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
      CALL GMONIT( 0,101   ,WTA,1D0,0D0)
         JESLIK=
     $         (XMP.LT.(AMTO-AMNE-AMCH))
      WTA=0D0
      IF (JESLIK) WTA=1D0
      CALL GMONIT( 0,102   ,WTA,1D0,0D0)
         JESLIK=
     $         (XMP.LT.(AMTO-AMNE-AMCH)).AND.
     $         (XP .GT.XMP)

      WTA=0D0
      IF (JESLIK) WTA=1D0
      CALL GMONIT( 0,103   ,WTA,1D0,0D0)
         JESLIK=
     $         (XMP.LT.(AMTO-AMNE-AMCH)).AND.
     $         (XP .GT.XMP)             .AND.
     $         (XP .LT.((AMTO**2+XMP**2-(AMCH+AMNE)**2)/2D0/AMTO))
      WTA=0D0
      IF (JESLIK) WTA=1D0
      CALL GMONIT( 0,104   ,WTA,1D0,0D0)
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
      FI1  =ANGFI(X1,X2)                                             
      X1 = PNEUTR(3)                                                 
      X2 = SQRT(PNEUTR(1)**2+PNEUTR(2)**2)                            
      TH1  =ANGXY(X1,X2) 
      CALL SPAJ(0,PNEUTR,PAA,PP,PE)                                             
      CALL LORTRA(3,-FI1,PNEUTR,VEC,PAA,PP,PE)
      CALL LORTRA(2,-TH1,PNEUTR,VEC,PAA,PP,PE)
        VEC(4)=0D0                                                  
        VEC(3)=0D0                                                  
        VEC(2)=0D0                                                  
        VEC(1)=1D0                                                  
        FI2=ANGFI(VEC(1),VEC(2))                                              
        CALL LORTRA(3,-FI2,PNEUTR,VEC,PAA,PP,PE)
      CALL SPAJ(1,PNEUTR,PAA,PP,PE)                                             
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
        PRINT *,' TRYPAR GBOOST LT 1., LIMIT OF PHASE SPACE '
     &         ,GNEU,GCHAR,QNEW,QOLD,AMTO,AMTOST,AMNE,AMCH
        PRINT *,XK,XKM,XK0DEC,AXK
        RETURN
      ENDIF
      PARCH =GCHAR+SQRT(GCHAR**2-1.0D0)
      PARNEU=GNEU -SQRT(GNEU**2 -1.0D0)
C 2) REDUCTIVE BOOSTS
      CALL BOSTD3(PARNEU,VEC ,VEC )
      CALL BOSTD3(PARNEU,PNEUTR,PNEUTR)
      CALL BOSTD3(PARCH,PAA ,PAA )
      CALL SPAJ(2,PNEUTR,PAA,PP,PE)                                             
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
      TH3  =ANGXY(X1,X2)
      PHOT(1)=PMOD*SINTHG*COS(FI3) 
      PHOT(2)=PMOD*SINTHG*SIN(FI3)
C MINUS BECAUSE AXIS OPPOSITE TO PAA
      PHOT(3)=-PMOD*COSTHG
      PHOT(4)=PENE
C ROTATE TO PUT PHOTON ALONG THIRD AXIS
      X1 = PHOT(1)
      X2 = PHOT(2) 
      CALL LORTRA(3,-FI3,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD3(-FI3,PHOT,PHOT) 
      CALL LORTRA(2,-TH3,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD2(-TH3,PHOT,PHOT)
      CALL SPAJ(21,PNEUTR,PAA,PP,PE)                                             
C ... now get the pair !
      PAIRB=PENE/XMP+PPED/XMP
      CALL BOSTD3(PAIRB,PE,PE)   
      CALL BOSTD3(PAIRB,PP,PP)   
      CALL SPAJ(3,PNEUTR,PAA,PP,PE) 
      GAMM=(PNEUTR(4)+PAA(4)+PP(4)+PE(4))/AMTO
      BPAR=GAMM-SQRT(GAMM**2-1D0)
      CALL LORTRA(1, BPAR,PNEUTR,VEC,PAA,PP,PE)
      CALL BOSTD3( BPAR,PHOT,PHOT)
      CALL SPAJ(4,PNEUTR,PAA,PP,PE)                                             
C BACK IN THE TAU REST FRAME BUT PNEUTR NOT YET ORIENTED.
      X1 = PNEUTR(1)
      X2 = PNEUTR(2)
      FI4  =ANGFI(X1,X2)
      X1 = PNEUTR(3)
      X2 = SQRT(PNEUTR(1)**2+PNEUTR(2)**2)
      TH4  =ANGXY(X1,X2)
      CALL LORTRA(3, FI4,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD3( FI4,PHOT,PHOT)
      CALL LORTRA(2,-TH4,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD2(-TH4,PHOT,PHOT)
      X1 = VEC(1)
      X2 = VEC(2)
      FI5=ANGFI(X1,X2)
      CALL LORTRA(3,-FI5,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD3(-FI5,PHOT,PHOT)
C PAA RESTORES ORIGINAL DIRECTION 
      CALL LORTRA(3, FI2,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD3( FI2,PHOT,PHOT)
      CALL LORTRA(2, TH1,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD2( TH1,PHOT,PHOT) 
      CALL LORTRA(3, FI1,PNEUTR,VEC,PAA,PP,PE)
      CALL ROTOD3( FI1,PHOT,PHOT)
      CALL SPAJ(10,PNEUTR,PAA,PP,PE)
      CALL LORTRA(1,BSTB,PNEUTR,VEC,PAA,PP,PE)
      CALL LORTRA(2,TH0,PNEUTR,VEC,PAA,PP,PE)
      CALL LORTRA(3,FI0,PNEUTR,VEC,PAA,PP,PE)
      CALL SPAJ(11,PNEUTR,PAA,PP,PE)                                             
C      STOP                                             
      END
      SUBROUTINE PARTRA(IBRAN,PHOT)       
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      COMMON /PARKIN/ 
     $ FI0,FI1,FI2,FI3,FI4,FI5,TH0,TH1,TH3,TH4
     $,PARNEU,PARCH,BPAR,BSTA,BSTB
      REAL*8 PHOT(4)
      CALL ROTOD3(-FI0,PHOT,PHOT)
      CALL ROTOD2(-TH0,PHOT,PHOT)
      CALL BOSTD3(BSTA,PHOT,PHOT)
      CALL ROTOD3(-FI1,PHOT,PHOT)
      CALL ROTOD2(-TH1,PHOT,PHOT)
      CALL ROTOD3(-FI2,PHOT,PHOT)
      IF(IBRAN.EQ.-1) THEN
        CALL BOSTD3(PARNEU,PHOT,PHOT)
      ELSE
        CALL BOSTD3(PARCH,PHOT,PHOT)
      ENDIF
      CALL ROTOD3(-FI3,PHOT,PHOT) 
      CALL ROTOD2(-TH3,PHOT,PHOT)
      CALL BOSTD3(BPAR,PHOT,PHOT)
      CALL ROTOD3( FI4,PHOT,PHOT)
      CALL ROTOD2(-TH4,PHOT,PHOT)
      CALL ROTOD3(-FI5,PHOT,PHOT)
      CALL ROTOD3( FI2,PHOT,PHOT)
      CALL ROTOD2( TH1,PHOT,PHOT) 
      CALL ROTOD3( FI1,PHOT,PHOT)
      CALL BOSTD3(BSTB,PHOT,PHOT)
      CALL ROTOD2( TH0,PHOT,PHOT)
      CALL ROTOD3( FI0,PHOT,PHOT)
      END
      FUNCTION AMAST(VEC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VEC(4)
      AMAST=SQRT(ABS(VEC(4)**2-VEC(1)**2-VEC(2)**2-VEC(3)**2)) 
      END
      SUBROUTINE LORTRA(KEY,PRM,PNEUTR,PNU,PAA,PP,PE)
C --------------------------------------------------------------------- 
C THIS ROUTINE PERFORMS LORENTZ TRANSFORMATION ON MANY 4-VECTORS        
C KEY   =1    BOOST    ALONG   3RD AXIS                                 
C       =2    ROTATION AROUND 2ND AXIS                                  
C       =3    ROTATION AROUND 3RD AXIS                                 
C PRM         TRANSFORMATION PARAMETER - ANGLE OR EXP(HIPERANGLE).     
C                                                                      
C     called by : RADCOR                                                
C --------------------------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      REAL*8 PNU(4),PAA(4)
      REAL*8 PNEUTR(4)                                                 
C                                                                     
      IF     (KEY.EQ.1) THEN                                           
        CALL BOSTD3(PRM,PNEUTR,PNEUTR)                                 
        CALL BOSTD3(PRM,PNU ,PNU )                                     
        CALL BOSTD3(PRM,PAA ,PAA )                                     
        CALL BOSTD3(PRM,PE ,PE )                                     
        CALL BOSTD3(PRM,PP ,PP )                                     
      ELSEIF (KEY.EQ.2) THEN                                           
        CALL ROTOD2(PRM,PNEUTR,PNEUTR)                                 
        CALL ROTOD2(PRM,PNU ,PNU )                                     
        CALL ROTOD2(PRM,PAA ,PAA )                                    
        CALL ROTOD2(PRM,PE  ,PE  )                                     
        CALL ROTOD2(PRM,PP  ,PP )                                    
      ELSEIF (KEY.EQ.3) THEN                                           
        CALL ROTOD3(PRM,PNEUTR,PNEUTR)                                
        CALL ROTOD3(PRM,PNU ,PNU )                                     
        CALL ROTOD3(PRM,PAA ,PAA )                                    
        CALL ROTOD3(PRM,PE  ,PE  )                                     
        CALL ROTOD3(PRM,PP  ,PP )                                    
      ELSE                                                            
        PRINT *, 'STOP IN LOTRA. WRONG KEYTRA'                        
        STOP                                                           
      ENDIF                                                            
      END                                                            
      DOUBLE PRECISION FUNCTION ANGXY(X,Y)                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DATA PI /3.141592653589793238462643D0/                            
C     
      IF(X.EQ.0D0.AND.Y.EQ.0.D0) THEN
       ANGXY=0D0
       RETURN
      ENDIF                                                      
      IF(ABS(Y).LT.ABS(X)) THEN                                         
        THE=ATAN(ABS(Y/X))                                              
        IF(X.LE.0D0) THE=PI-THE                                         
      ELSE                                                              
        THE=ACOS(X/SQRT(X**2+Y**2))                                     
      ENDIF                                                             
      ANGXY=THE                                                         
      RETURN                                                            
      END                                                               


      SUBROUTINE ROTOD2(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
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
      SUBROUTINE ROTOD3(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
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
      SUBROUTINE BOSTD3(EXE,PVEC,QVEC)
C ----------------------------------------------------------------------
C BOOST ALONG Z AXIS, EXE=EXP(ETA), ETA= HIPERBOLIC VELOCITY.
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
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

      DOUBLE PRECISION FUNCTION ANGFI(X,Y)
C     *******************
* CALCULATES ANGLE IN (0,2*PI) RANGE OUT OF X-Y
*     ***********************
      IMPLICIT REAL*8(A-H,O-Z)
      DATA PI /3.1415926535897932D0/
      IF(X.EQ.0D0.AND.Y.EQ.0.D0) THEN
       ANGFI=0D0
       RETURN
      ENDIF
      IF(ABS(Y).LT.ABS(X)) THEN
        THE=ATAN(ABS(Y/X))
        IF(X.LE.0D0) THE=PI-THE
      ELSE
        THE=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      IF(Y.LT.0D0) THE=2D0*PI-THE
      ANGFI=THE
      END


      SUBROUTINE VARRAN(DRVEC,LEN)
C     ***************************
C Switchable random number generator
C Translation to double precision
C     ***************************
      COMMON / RANPAR / KEYRND
      save   / RANPAR /
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
        WRITE(6,*) ' +++++ VARRAN: RVEC=',RVEC(I)
        GOTO 10
      ENDIF
      DRVEC(I)=RVEC(I)
   30 CONTINUE
      RETURN
  901 WRITE(6,*) ' +++++ STOP IN VARRAN: LEN=',LEN
      STOP
  902 WRITE(6,*) ' +++++ STOP IN VARRAN: WRONG KEYRND',KEYRND
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





















