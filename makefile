# Copyright (c) 2013-2014, Sergio Callegari
# All rights reserved.

# This file is part of PyDSM.

# PyDSM is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# PyDSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with PyDSM.  If not, see <http://www.gnu.org/licenses/>.

# Deprecated makefile for preparing a built version
# This is still useful for marking the package with a git description

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
