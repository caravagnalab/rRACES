# prepare the package for release
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

.PHONY: all docs build build-cran install check travis clean

docs:
	R --vanilla --silent -e 'Rcpp::compileAttributes()'
	R --vanilla --silent -e 'roxygen2::roxygenise(roclets=c("rd"))'

build: 
	(cd .. && R CMD build --no-manual $(PKGSRC))

build-cran: 
	(cd .. && R CMD build $(PKGSRC))

install: build
	(cd .. && R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz)

check: build-cran
	(cd .. && R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran)

travis: build
	(cd ..; && R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --no-manual)

clean:
	rm -rf _install _builds RACES 
	rm src/RcppExports.cpp src/*.o src/*.so R/RcppExports.R
	rm -rf man/*.Rd
	rm .roxygen.lock
	(cd .. && rm -rf $(PKGNAME).Rcheck)
