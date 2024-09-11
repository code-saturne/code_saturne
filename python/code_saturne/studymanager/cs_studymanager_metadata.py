# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, logging
from code_saturne.base.cs_package import package as cs_pkg

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# class handling the metadata
#-------------------------------------------------------------------------------

class study_metadata:
    """Class handling metadata of a studymanager study"""

    #---------------------------------------------------------------------------

    def __init__(self, pkg: cs_pkg, filename: str):
        """
        Constructor for the class.
        @param pkg      cs_package
        @param filename name of the study xml file
        """
        from code_saturne.model import XMLengine
        from code_saturne.studymanager.cs_studymanager_parser import Parser

        _smgr = XMLengine.Case(package=pkg,
                               file_name=filename,
                               studymanager=True)
        _smgr['xmlfile'] = filename
        parser = Parser(filename, doc=_smgr.doc)

        sname = parser.getStudiesLabel()[0]

        md = parser.getStudyMetadata(sname)

        self.study = md['study']
        self.cases = md['cases']
        self.n_cases = len(self.cases)

    #---------------------------------------------------------------------------

    def __repr__(self):
        """dunder method"""

        return f'study_metada("study":{self.study}, "cases":{self.cases})'

    #---------------------------------------------------------------------------

    def __str__(self):
        """dunder method"""

        _str = "STUDY '" + self.study['name'] + "'\n"
        _str += "Study tags:"
        if self.study['tags::study']:
            for t in self.study['tags::study']:
                _str += " " + t
        else:
            _str += " N/A"
        _str += "\n"
        _str += "-----\n\n"

        _str += "Cases tags:"
        if self.study['tags::cases']:
            for t in self.study['tags::cases']:
                _str += " " + t
        else:
            _str += " N/A"
        _str += "\n"
        _str += "-----\n\n"

        _str += "keywords:\n"
        _str += "---------\n"
        if self.study['keywords']:
            for kw in self.study['keywords']:
                _str += " - " + kw + "\n"
        else:
            _str += " -> N/A\n"
        _str += "\n"

        _str += "Cases:\n"
        _str += "------\n"

        if not self.cases:
            _str += " No cases listed in this study.\n\n"

        for c in self.cases:
            _str += " * CASE '"+ c['name'] + "'\n"

            if c['tags']:
                _str += "   # tags:"
                for t in c['tags']:
                    _str += " " + t
                _str += "\n"
                _str += "   -------\n"
                _str += "\n"

            if c['vnvitem']:
                _str += "   # description:\n"
                _str += "   --------------\n"

                for itm in c['vnvitem']:
                    _str += "   - " + itm + "\n"

                _str += "\n"

        return _str

    #---------------------------------------------------------------------------

    def log(self, fname: str = None):
        """
        Log metadata to stdout or file.
        @param fname Name of file to log metadata. If None logged to stdout.
        """

        if fname:
            with open(fname, "w") as logfile:
                logfile.write(self.__str__())

        else:
            print(self)

    #---------------------------------------------------------------------------

    def dump_keywords(self, path: str = None):
        """
        Dump keywords to file <path>/study_keywords.tex
        @param path path to write file. If None file is written in current location.
        """

        f2open = 'study_keywords.tex'
        if path and os.path.isdir(path):
            f2open = os.path.join(path, f2open)

        with open(f2open, "w") as kwfile:
            kw_str = '\\keywords{' + ', '.join(self.study['keywords']) + '}'
            kwfile.write(kw_str)

    #---------------------------------------------------------------------------

    def dump_readme(self, path: str = None):
        """
        Dump readme data to file <path>/study_readme.tex
        @param path path to write file. If None file is written in current location.
        """

        f2open = 'study_readme.tex'
        if path and os.path.isdir(path):
            f2open = os.path.join(path, f2open)

        with open(f2open, 'w') as readmefile:
            for c in self.cases:
                readmefile.write(f'\\vnvcase[{c["name"]}:]\n')
                for itm in c['vnvitem']:
                    readmefile.write(f'\\vnvitem {itm}\n')
                readmefile.write('\n')

#-------------------------------------------------------------------------------
