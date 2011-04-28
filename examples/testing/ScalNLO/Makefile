############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-pispectr-RCall
# Or use other benchmark reference files:
# NAME_BENCHMARK=mc-tester.root-pispectr-RCtau
# NAME_BENCHMARK=mc-tester.root-pispectr-RCno

CONF=pythia_Z.conf
NAME=Ztautau

PROGRAM=$(PHOTOS_TAUOLA)
#MODE=   $(PIPI_MODE)
MODE=   7

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(MODE) #2>/dev/null
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
