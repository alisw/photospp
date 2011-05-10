############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-1000M-NLO-EXP
# Or use other benchmark reference files:
#NAME_BENCHMARK=mc-tester.root-1000M-NLO-EXP
#NAME_BENCHMARK=mc-tester.root-1000M-LO-ALPHA
#NAME_BENCHMARK=mc-tester.root-1000M-NLO-ALPHA

CONF=pythia_Z.conf
NAME=ScalarNLO

PROGRAM=$(PHOTOS_TAUOLA)
MODE=   $(K0S_MODE)

USE_ALPHA = 0
USE_NLO   = 1

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(MODE) $(USE_ALPHA) $(USE_NLO) #2>/dev/null
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
