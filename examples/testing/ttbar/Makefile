############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-ttbar-photos

CONF=pythia_ttbar.conf
NAME=ttbar

PROGRAM=$(PHOTOS_ONLY)
MODE=   $(TOP_MODE)

tester.$(NAME).pdf: mc-tester.$(NAME).root
	@echo
	@echo "IMPORTANT:"
	@echo "Following example is missing hadronization of b b_bar pair"
	@echo
#	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
#	mv tester.pdf tester.$(NAME).pdf
#	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(MODE) 2>/dev/null
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root
