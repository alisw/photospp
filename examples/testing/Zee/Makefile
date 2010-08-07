############## Make tests ###########################
include ../test.inc

#NAME_BENCHMARK=NEW-Zee-100M-10MeV.root
NAME_BENCHMARK=mc-tester.root-Zee-1M-good #photos-100M-sym-10MeV #or use other benchmark

CONF=pythia_Z.conf
NAME=Zee

PROGRAM=$(PHOTOS_ONLY)
MODE=   $(PP_MODE)

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(MODE) 2>/dev/null
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
