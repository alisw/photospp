############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.Zmumu.root-O1-NLO-100M
#NAME_BENCHMARK=mc-tester.Zmumu.root-NLO-100M
#NAME_BENCHMARK=mc-tester.Zmumu.root-attempt
#NAME_BENCHMARK=mc-tester.root-Zmumu-photos-F77-100M-sym-10MeV
# Or use other benchmark reference files:
# NAME_BENCHMARK=mc-tester.root-Zmumu-kkmc-F77-10M-sym-1GeV
# NAME_BENCHMARK=mc-tester.root-Zmumu-CEEX1-100M-nosym-1GeV
# NAME_BENCHAMRK=mc-tester.root-Zmumu-CEEX2-100M-nosym-1GeV

CONF=pythia_Z.conf
NAME=ZmumuNLO

PROGRAM=$(PHOTOS_ONLY)
PYTHIA_MODE=$(Z_MU_MODE)
SPECIAL_MODE=$(NLO_MODE)
ONE_PHOTON_MODE=1

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(PYTHIA_MODE) $(EVENTS) $(SPECIAL_MODE) $(ONE_PHOTON_MODE)
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
