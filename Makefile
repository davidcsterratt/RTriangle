TRIANGLE_VERSION=$(shell grep Version pkg/DESCRIPTION | perl -p -e "s/Version: //;")
TRIANGLE_SVN_REVISION=$(shell svn info -R | grep "Revision:" | perl -p -e 's/Revision: //;' | sort -n -r | head -1)
TRIANGLE_SVN_REVISION1=$(shell echo $(TRIANGLE_SVN_REVISION) + 1 | bc) 
PACKAGE=RTriangle_$(TRIANGLE_VERSION).tar.gz

roxygen:
	rm -f pkg/man/*
	echo "library(roxygen2) ; roxygenize(\"pkg\")" |	R --no-restore --slave

package: roxygen 
	rm -f pkg/R/*~
	R CMD build pkg

install: package
	R CMD INSTALL --latex $(PACKAGE) 

doc: package
	rm -f RTriangle.pdf
	R CMD Rd2dvi --pdf --output=RTriangle.pdf pkg 

check:
	R CMD check $(PACKAGE)

revision:
	@echo $(TRIANGLE_SVN_REVISION)
	@echo $(TRIANGLE_SVN_REVISION1)

changelog:
	cd pkg &&	svn2cl -i -r HEAD:{2010-01-01}
