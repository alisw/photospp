############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-Wmunu-photos-100M-1-phot-NLO
# Or use other benchmark reference files:
#NAME_BENCHMARK=mc-tester.root-Wmunu-photos-100M-1-phot
#NAME_BENCHMARK=mc-tester.root-Wmunu-photos-100M-nosym-10MeV-1-phot-NLO
#NAME_BENCHMARK=mc-tester.root-Wmunu-photos-100M-sym-10MeV

CONF=pythia_W.conf
NAME=WmunuNLO

PROGRAM=$(PHOTOS_ONLY)
PYTHIA_MODE=$(DEFAULT_MODE)
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
