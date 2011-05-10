############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-Wmunu-photos-100M-sym-10MeV


CONF=pythia_W.conf
NAME=Wmunu

PROGRAM=$(PHOTOS_ONLY)
MODE=   $(DEFAULT_MODE)

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
