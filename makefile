# Makefile for preparing a release version

git_version_string:=$(shell git version 2>/dev/null)
ifneq ($(findstring git version, $(git_version_string)), git version)
 $(error Makefile needs git)
endif

pydsm_version:=$(shell git describe --tags)
ifeq ($(pydsm_version),)
 $(error Makefile should be run in git workdir)
endif

dist:	source_dist html_doc_dist

source_dist:
	git archive --format=zip -9 --prefix=pydsm-$(pydsm_version)/ \
	-o pydsm-$(pydsm_version).zip HEAD

html_doc_dist:
	python setup.py build
	rm -fr doc/build/pydsm-doc-html-$(pydsm_version)
	make -C doc html
	cp -a doc/build/html doc/build/pydsm-doc-html-$(pydsm_version)
	rm -fr doc/build/pydsm-doc-html-$(pydsm_version)/_sources
	rm -f pydsm-doc-html-$(pydsm_version).zip
	cd doc/build ; \
	zip -r ../../pydsm-doc-html-$(pydsm_version).zip \
	  pydsm-doc-html-$(pydsm_version)
	rm -fr doc/build/pydsm-doc-html*

.phony:	dist source_dist html_doc_dist