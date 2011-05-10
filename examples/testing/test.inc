EVENTS = 1000

VIEWER=gv
COMPARE=../compare.sh

# Programs
PHOTOS_ONLY   =../photos_test.exe
PHOTOS_TAUOLA =../photos_tauola_test.exe

# PHOTOS_ONLY -------------------------------------------

# Pythia modes
DEFAULT_MODE =0 # ee@200GeV    (for Wenu, Wmunu)
Z_MU_MODE    =1 # ee@91.187GeV (for Zmumu)
PP_MODE      =2 # pp@1.4TeV    (for Zee, ttbar)

# Special modes
TOP_MODE     =1 # Turn on MC-Tester fix for Top decays
NLO_MODE     =2 # Turn on Photos NLO corrections

# PHOTOS_TAUOLA -----------------------------------------

# Tauola generation modes
PIPI_MODE    =3 # ee@91.187GeV (for Ztautau)
RHORHO_MODE  =4 # ee@500GeV    (for Htautau)
K0S_MODE     =7 # ee@91.187GeV (for ScalarNLO)
