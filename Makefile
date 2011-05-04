TRIANGLE_VERSION=$(shell grep Version pkg/DESCRIPTION | perl -p -e "s/Version: //;")
TRIANGLE_SVN_REVISION=$(shell svn info -R | grep "Revision:" | perl -p -e 's/Revision: //;' | sort -n -r | head -1)
TRIANGLE_SVN_REVISION1=$(shell echo $(TRIANGLE_SVN_REVISION) + 1 | bc) 
PACKAGE=Triangle_$(TRIANGLE_VERSION).tar.gz

roxygen:
	rm -f pkg/man/*
	R CMD roxygen -d pkg

package: roxygen 
	rm -f pkg/R/*~
	R CMD build pkg

install: package
	R CMD INSTALL --latex $(PACKAGE) 

revision:
	@echo $(TRIANGLE_SVN_REVISION)
	@echo $(TRIANGLE_SVN_REVISION1)
