# -*- coding: utf-8 -*-

# Copyright (c) 2014, Sergio Callegari
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

"""docdist

Implements a Distutils 'docdist' subcommand (create zipfile with html
documentation).
"""

import os
from distutils import dir_util
from distutils.cmd import Command

__all__ = ["docdist"]


class docdist(Command):

    description = "Create zip file containing standalone html docs"

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        cmd = self.get_finalized_command('build')
        self.build_dir = os.path.join(cmd.build_base, 'docdist')
        self.dist_dir = "dist"
        self.format = "zip"

    def run(self):
        # call build sphinx to build docs
        self.run_command("build_sphinx")
        cmd = self.get_finalized_command("build_sphinx")
        source_dir = cmd.builder_target_dir

        # copy to directory with appropriate name
        dist = self.distribution
        arc_name = "%s-doc-html-%s" % (dist.get_name(), dist.get_version())
        tmp_dir = os.path.join(self.build_dir, arc_name)
        if os.path.exists(tmp_dir):
            dir_util.remove_tree(tmp_dir, dry_run=self.dry_run)
        self.copy_tree(source_dir, tmp_dir, preserve_symlinks=True)
        dir_util.remove_tree(os.path.join(tmp_dir, '_sources'))

        # make archive from dir
        arc_base = os.path.join(self.dist_dir, arc_name)
        self.arc_filename = self.make_archive(arc_base, self.format,
                                              self.build_dir)

        dir_util.remove_tree(tmp_dir, dry_run=self.dry_run)
