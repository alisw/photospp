SUBDIRS = photos-fortran src

# Post-installation: move all files from 'include' to 'include/Photos'
install-data-hook:
	mkdir -p $(prefix)/include/Photos && \
	mv $(prefix)/include/*.h $(prefix)/include/Photos/.
