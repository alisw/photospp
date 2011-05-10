############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-Zmumu-photos-F77-100M-sym-10MeV
# Or use other benchmark reference files:
# NAME_BENCHMARK=mc-tester.root-Zmumu-kkmc-F77-10M-sym-1GeV
# NAME_BENCHMARK=mc-tester.root-Zmumu-CEEX1-100M-nosym-1GeV
# NAME_BENCHAMRK=mc-tester.root-Zmumu-CEEX2-100M-nosym-1GeV

CONF=pythia_Z.conf
NAME=Zmumu

PROGRAM=$(PHOTOS_ONLY)
MODE=   $(Z_MU_MODE)

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(MODE) $(EVENTS)
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
