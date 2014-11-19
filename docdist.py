# -*- coding: utf-8 -*-
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
