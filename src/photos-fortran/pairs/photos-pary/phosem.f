

CC======================================================================
CC                                                                    CC
CC WRITTEN BY E. RICHTER- WAS                                         CC
CC                                                                    CC
CC    NOVEMBER 1992                                                   CC 
CC                                                                    CC  
CC======================================================================
C THIS IS UNPUBLISHED VERSION OF PROGRAM WHICH COULD BE APPLIED TO TEST
C LEPTONIC SPECTRUM FOR JETSET7.3+PHOTOS WITH SEMIANALYTICAL FORMULA 
C FOR HISTOGRAMING PRIVATE LIBBRARY GLIBK.FOR IS USED (can be replaced by HBOOK)
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
      DO 10 IEV=1,1 00
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
      if(xriad.gt.0.2d0.AND.IEV.LT.100101) call lulist(1)
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
     





















