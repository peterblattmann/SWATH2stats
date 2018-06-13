export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean roxygen reference check build test

build:
	@echo "Performing build with R CMD build hpgltools"
	R CMD build .

check:
	@echo "Performing check with R CMD check hpgltools"
	export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check . --no-build-vignettes

clean:
	@echo "Cleaning up"
	@echo "Not currently doing anything here."

clean_vignette:
	rm -f vignettes/*.rda vignettes/*.map vignettes/*.Rdata

dep_push: deps snap
	echo "Setting default commit message and pushing."
	git commit -a -m 'packrat modification.'
	git push

deps:
	@echo "Invoking devtools::install_dev_deps()"
	R -e "suppressPackageStartupMessages(suppressMessages(source('http://bioconductor.org/biocLite.R')));\
all = as.data.frame(devtools::dev_package_deps('.', dependencies=TRUE)); needed = all[['diff']] < 0; needed = all[needed, 'package']; for (t in needed) { biocLite(t) }"

document: roxygen vignette reference

inst: roxygen install test
	@echo "Regenerated documentation, installed, and tested."

install:
	@echo "Performing R CMD INSTALL hpgltools globally."
	R CMD INSTALL .

install_bioconductor:
	R -e "library(hpgltools); bioc_all()"

push:
	echo "Pushing to github."
	git commit -a && git push

reference:
	@echo "Generating reference manual with R CMD Rd2pdf"
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf --no-preview

roxygen:
	@echo "Generating documentation with devtools::document()"
	R -e "suppressPackageStartupMessages(devtools::document()); warnings()"

suggests:
	@echo "Installing suggested packages."
	R -e "source('http://bioconductor.org/biocLite.R');\
library(desc);\
d = description\$$new(); suggests = d\$$get('Suggests');\
 suggests = gsub(pattern='\\n', replacement='', x=suggests);\
 suggests = gsub(pattern=' ', replacement='', x=suggests);\
 suggests = strsplit(x=suggests, split=',');\
 for (pkg in suggests[[1]]) { if (! pkg %in% installed.packages()) { biocLite(pkg); } else { message(paste0(pkg, ' is already installed.')) } };"

test: roxygen
	@echo "Running run_tests.R"
	tests/test-all.R

update:
	R -e "source('http://bioconductor.org/biocLite.R'); biocLite(); library(BiocInstaller); biocValid()"

update_bioc:
	R -e "source('http://bioconductor.org/biocLite.R'); biocLite(); biocLite('BiocUpgrade');"

vignette:
	@echo "Building vignettes with devtools::build_vignettes()"
	R -e "devtools::build_vignettes()"

vt:	clean_vignette vignette reference install
