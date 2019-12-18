# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

import os, sys
import shutil, re
import subprocess
import threading
import string
import time
import logging
import fnmatch

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.cs_exec_environment import get_shell_type, enquote_arg
from code_saturne.cs_compile import files_to_compile, compile_and_link
from code_saturne import cs_create
from code_saturne.cs_create import set_executable
from code_saturne import cs_runcase

from code_saturne.model import XMLengine
from code_saturne.studymanager.cs_studymanager_pathes_model import PathesModel

from code_saturne.studymanager.cs_studymanager_parser import Parser
from code_saturne.studymanager.cs_studymanager_texmaker import Report1, Report2

try:
    from code_saturne.studymanager.cs_studymanager_drawing import Plotter
except Exception:
    print("Warning: import studymanager Plotter failed. Plotting disabled.\n")
    pass

from code_saturne.studymanager.cs_studymanager_run import run_studymanager_command
from code_saturne.studymanager.cs_studymanager_xml_init import smgr_xml_init

#-------------------------------------------------------------------------------
# log config.
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
#log.setLevel(logging.DEBUG)
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------

def nodot(item):
    return item[0] != '.'

#-------------------------------------------------------------------------------

def create_base_xml_file(filepath, pkg):
    """Create studymanager XML file.
    """

    filename = os.path.basename(filepath)
    if os.path.isfile(filepath):
        print("Can not create XML file of parameter:\n" \
              + filepath + " already exists.")
        sys.exit(1)

    # using xml engine from Code_Saturne GUI
    smgr = XMLengine.Case(package=pkg, studymanager=True)
    smgr['xmlfile'] = filename
    pm = PathesModel(smgr)

    # empty repo and dest
    pm.setRepositoryPath('')
    pm.setDestinationPath('')

    return smgr

#-------------------------------------------------------------------------------

def init_xml_file_with_study(smgr, studyp):
    """Initialize XML file with study content and save it.
    """

    smgr_node = smgr.xmlGetNode('studymanager')

    studyd = os.path.basename(studyp)
    st_node = smgr_node.xmlInitChildNode('study', label = studyd)
    st_node['status'] = "on"

    cases = []
    for elt in os.listdir(studyp):
        eltd = os.path.join(studyp, elt)
        if os.path.isdir(eltd):
            if isCase(eltd):
                cases.append(elt)

    cases.sort()
    for case in cases:
        c_node = st_node.xmlInitChildNode("case", label = case)
        c_node['status']  = "on"
        c_node['compute'] = "on"
        c_node['post']    = "on"

    smgr.xmlSaveDocument()

#-------------------------------------------------------------------------------

def isStudy(dirpath):
    """Try to determine if dirpath is a Code_Saturne study directory.
    """

    meshd = os.path.join(dirpath, 'MESH')
    is_study = os.path.isdir(meshd)

    return is_study

#-------------------------------------------------------------------------------

def isCase(dirpath):
    """Try to determine if dirpath is a Code_Saturne case directory.
    """

    datad = os.path.join(dirpath, 'DATA')
    scriptd = os.path.join(dirpath, 'SCRIPTS')
    is_case = os.path.isdir(datad) and os.path.isdir(scriptd)

    return is_case

#===============================================================================
# Case class
#===============================================================================

class Case(object):
    def __init__(self, pkg, rlog, diff, parser, study, data, repo, dest):
        """
        @type data: C{Dictionary}
        @param data: contains all keyword and value read in the parameters file
        """
        self.__log      = rlog
        self.__diff     = diff
        self.__parser   = parser
        self.__study    = study
        self.__data     = data
        self.__repo     = repo
        self.__dest     = dest

        self.node       = data['node']
        self.label      = data['label']
        self.compute    = data['compute']
        self.plot       = data['post']
        self.run_id     = data['run_id']
        self.tags       = data['tags']
        self.compare    = data['compare']

        self.is_compiled= "not done"
        self.is_run     = "not done"
        self.is_time    = None
        self.is_plot    = "not done"
        self.is_compare = "not done"
        self.threshold  = "default"
        self.diff_value = [] # list of differences (in case of comparison)
        self.m_size_eq  = True # mesh sizes equal (in case of comparison)
        self.subdomains = None
        self.run_dir    = ""

        self.resu = 'RESU'

        # Specific case for coupling

        coupling = os.path.join(self.__repo, self.label, "coupling_parameters.py")
        if os.path.isfile(coupling):

            from code_saturne import cs_case_coupling
            try:
                exec(compile(open(coupling).read(), '<string>', 'exec'))
            except Exception:
                execfile(coupling)

            run_ref = os.path.join(self.__repo, self.label, "runcase")
            self.exe, self.pkg = self.__get_exe(pkg, run_ref)
            self.subdomains = []
            for d in locals()['domains']:
                if d['solver'].lower() in ('code_saturne', 'neptune_cfd'):
                    self.subdomains.append(d['domain'])
                elif d['solver'].lower() == "syrthes":
                    syrthes = True

            self.resu = 'RESU_COUPLING'

            # insert syrthes path
            if syrthes:
                syrthes_insert = cs_create.syrthes_path_line(pkg)
                if syrthes_insert:
                    fd = open(coupling)
                    fd_lines = fd.readlines()
                    fd.close()
                    fd = open(coupling, 'w')
                    for line in fd_lines:
                        if "sys.path.insert" in line:
                            fd.write(syrthes_insert)
                        else:
                            fd.write(line)
                    fd.close()

        else:
            run_ref = os.path.join(self.__repo, self.label, "SCRIPTS", "runcase")
            self.exe, self.pkg = self.__get_exe(pkg, run_ref)

    #---------------------------------------------------------------------------

    def __get_exe(self, old_pkg, run_ref):
        """
        Return the name of the exe of the case, in order to mix
        Code_Saturne and NEPTUNE_CFD test cases in the same study.
        """

        # Read the runcase script from the Repository

        runcase = cs_runcase.runcase(run_ref)

        if runcase.cmd_name == "code_saturne":
            from code_saturne.cs_package import package
        elif runcase.cmd_name == "neptune_cfd":
            from neptune_cfd.nc_package import package
        pkg = package()

        return runcase.cmd_name, pkg

    #---------------------------------------------------------------------------

    def __update_domain(self, subdir, xmlonly=False):
        """
        Update path for the script in the Repository.
        """
        # 1) Load the xml file of parameters in order to update it
        #    with the __backwardCompatibility method.

        from code_saturne.model.XMLengine import Case

        if self.exe == "code_saturne":
            from code_saturne.model.XMLinitialize import XMLinit as solver_xml_init
        elif self.exe == "neptune_cfd":
            from code_saturne.model.XMLinitializeNeptune import XMLinit as solver_xml_init

        for fn in os.listdir(os.path.join(self.__repo, subdir, "DATA")):
            fp = os.path.join(self.__repo, subdir, "DATA", fn)
            if os.path.isfile(fp):
                fd = os.open(fp , os.O_RDONLY)
                f = os.fdopen(fd)
                l = f.readline()
                f.close()
                if l.startswith('''<?xml version="1.0" encoding="utf-8"?><Code_Saturne_GUI''') or \
                   l.startswith('''<?xml version="1.0" encoding="utf-8"?><NEPTUNE_CFD_GUI'''):
                    try:
                        case = Case(package = self.pkg, file_name = fp)
                    except:
                        print("Parameters file reading error.\n")
                        print("This file is not in accordance with XML specifications.\n")
                        sys.exit(1)

                    case['xmlfile'] = fp
                    case.xmlCleanAllBlank(case.xmlRootNode())
                    solver_xml_init(case).initialize()
                    case.xmlSaveDocument()


        # 2) Recreate missing directories which are compulsory
        # git removes these usually empty directories
        if not xmlonly:
            git_rm_dirs = ["RESU", "SRC"]
            for gd in git_rm_dirs:
                r = os.path.join(self.__repo, subdir, gd)
                if not os.path.isdir(r):
                    os.makedirs(r)

        # 3) Update the GUI script from the Repository
        data_subdir = os.path.join(self.__repo, subdir, "DATA")
        self.update_gui_script_path(data_subdir, None, xmlonly)

        # 4) Update the runcase script from the Repository
        scripts_subdir = os.path.join(self.__repo, subdir, "SCRIPTS")
        self.update_runcase_path(scripts_subdir, None, xmlonly)

    #---------------------------------------------------------------------------

    def update_gui_script_path(self, subdir, destdir=None, xmlonly=False):
        """
        Update path for the script in the Repository.
        """
        if not destdir:
            gui_script = os.path.join(self.__repo, subdir, self.pkg.guiname)
        else:
            gui_script = os.path.join(destdir, subdir, self.pkg.guiname)

        have_gui = 1
        try:
            f = open(gui_script, mode = 'r')
        except IOError:
           print("Warning SaturneGUI does not exist: %s\n" % gui_script)
           have_gui = 0

        if have_gui:
           lines = f.readlines()
           f.close()

           for i in range(len(lines)):
               if re.search(r'^export PATH=', lines[i]):
                   if xmlonly:
                       lines[i] = 'export PATH="":$PATH\n'
                   else:
                       lines[i] = 'export PATH="' + self.pkg.get_dir('bindir')\
                                                  +'":$PATH\n'
           f = open(gui_script, mode = 'w')
           f.writelines(lines)
           f.close()

           set_executable(gui_script)

    #---------------------------------------------------------------------------

    def update_runcase_path(self, subdir, destdir=None, xmlonly=False):
        """
        Update path for the script in the Repository.
        """
        if not destdir:
            batch_file = os.path.join(self.__repo, subdir, "runcase")
        else:
            batch_file = os.path.join(destdir, subdir, "runcase")

        try:
            f = open(batch_file, mode = 'r')
        except IOError:
            print("Error: can not open %s\n" % batch_file)
            sys.exit(1)

        lines = f.readlines()
        f.close()

        for i in range(len(lines)):
            if lines[i].strip()[0:1] == '#':
                continue
            if xmlonly:
                if re.search(r'^export PATH=', lines[i]):
                    lines[i] = 'export PATH="":$PATH\n'
            else:
                if re.search(r'^export PATH=', lines[i]):
                    lines[i] = 'export PATH="' + self.pkg.get_dir('bindir') +'":$PATH\n'

        f = open(batch_file, mode = 'w')
        f.writelines(lines)
        f.close()

        set_executable(batch_file)

    #---------------------------------------------------------------------------

    def update(self, xmlonly=False):
        """
        Update path for the script in the Repository.
        """
        # 1) Load the xml file of parameters in order to update it
        #    with the __backwardCompatibility method.

        if self.subdomains:
            cdirs = []
            for d in self.subdomains:
                cdirs.append(os.path.join(self.label, d))
        else:
            cdirs = (self.label,)

        for d in cdirs:
            self.__update_domain(d, xmlonly)

        # Update the runcase script from the repository in case of coupling

        if self.subdomains:
            case_dir = os.path.join(self.__repo, self.label)
            self.update_runcase_path(case_dir, None, xmlonly)
            # recreate possibly missing but compulsory directory
            r = os.path.join(case_dir, "RESU_COUPLING")
            if not os.path.isdir(r):
                os.makedirs(r)


    #---------------------------------------------------------------------------

    def test_compilation(self, study_path, log):
        """
        Test compilation of sources for current case (if some exist).
        @rtype: C{String}
        @return: compilation test status (None if no files to compile).
        """

        if self.subdomains:
            sdirs = []
            for sd in self.subdomains:
                sdirs.append(os.path.join(study_path, self.label, sd, 'SRC'))
        else:
            sdirs = (os.path.join(study_path, self.label, 'SRC'),)

        # compilation test mode
        dest_dir = None

        self.is_compiled = None
        retcode = 0

        # loop over subdomains
        for s in sdirs:
            src_files = files_to_compile(s)

            if len(src_files) > 0:
                self.is_compiled = "OK"
                retcode += compile_and_link(self.pkg, s, dest_dir,
                                            stdout=log, stderr=log)

        if retcode > 0:
            self.is_compiled = "KO"

        return self.is_compiled

    #---------------------------------------------------------------------------

    def __suggest_run_id(self):

        cmd = enquote_arg(os.path.join(self.pkg.get_dir('bindir'), self.exe)) + " run --suggest-id"
        if self.subdomains:
            cmd += " --coupling=coupling_parameters.py"
        p = subprocess.Popen(cmd,
                             shell=True,
                             executable=get_shell_type(),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        i = p.communicate()[0]
        run_id = " ".join(i.split())

        return run_id, os.path.join(self.__dest, self.label, self.resu, run_id)

    #---------------------------------------------------------------------------

    def __updateRuncase(self, run_id):
        """
        Update the command line in the launcher C{runcase}.
        """
        from code_saturne.cs_exec_environment import separate_args, \
            get_command_single_value

        # Prepare runcase path
        scripts_repo = os.path.join(self.__repo, self.label)
        scripts_dest = os.path.join(self.__dest, self.label)
        if not self.subdomains:
            scripts_repo = os.path.join(scripts_repo, "SCRIPTS")
            scripts_dest = os.path.join(scripts_dest, "SCRIPTS")

        run_ref = os.path.join(scripts_repo, "runcase")
        run_ref_win = os.path.join(scripts_repo, "runcase.bat")
        run_new = os.path.join(scripts_dest, "runcase")

        if sys.platform.startswith('win'):
            run_new = os.path.join(scripts_dest, "runcase.bat")

        # Read runcase from repo
        path = run_ref
        if not os.path.isfile(path):
            path = run_ref_win
        if not os.path.isfile(path):
            print("Error: could not find %s (or %s)\n" % run_ref, run_ref_win)
            sys.exit(1)

        runcase_repo = cs_runcase.runcase(path, create_if_missing=False)

        # Read runcase from dest
        path = run_new
        runcase_dest = cs_runcase.runcase(path, create_if_missing=False,
                                          ignore_batch=True)

        # Assign run command from repo in dest
        runcase_dest.set_run_args(runcase_repo.get_run_args())

        # set run_id in dest
        runcase_dest.set_run_id(run_id=run_id)

        # Set number of processors if provided
        n_procs = self.__data['n_procs']
        if n_procs:
            runcase_dest.set_nprocs(n_procs)

        # Write runcase
        runcase_dest.save()

    #---------------------------------------------------------------------------

    def run(self):
        """
        Check if a run with same result subdirectory name exists
        and launch run if not.
        """
        home = os.getcwd()
        if self.subdomains:
            os.chdir(os.path.join(self.__dest, self.label))
        else:
            os.chdir(os.path.join(self.__dest, self.label, 'SCRIPTS'))

        if self.run_id:
            run_id = self.run_id
            run_dir = os.path.join(self.__dest, self.label, self.resu, run_id)

            if os.path.isdir(run_dir):
                if os.path.isfile(os.path.join(run_dir, "error")):
                    self.is_run = "KO"
                    error = 1
                else:
                    self.is_run = "OK"
                    error = 0
                os.chdir(home)

                return error

        else:
            run_id, run_dir = self.__suggest_run_id()

            while os.path.isdir(run_dir):
                run_id, run_dir = self.__suggest_run_id()

        self.run_id  = run_id
        self.run_dir = run_dir

        self.__updateRuncase(run_id)

        if sys.platform.startswith('win'):
            error, self.is_time = run_studymanager_command("runcase.bat", self.__log)
        else:
            error, self.is_time = run_studymanager_command("./runcase", self.__log)

        if not error:
            self.is_run = "OK"
        else:
            self.is_run = "KO"

        os.chdir(home)

        return error

    #---------------------------------------------------------------------------

    def runCompare(self, studies, r, d, threshold, args, reference=None):
        home = os.getcwd()

        node = None

        if reference:
            result = os.path.join(reference, self.label, self.resu)
        else:
            result = os.path.join(self.__repo, self.label, self.resu)
        # check_dir called again here to get run_id (possibly date-hour)
        repo, msg = self.check_dir(node, result, r, "repo")
        if msg:
            studies.reporting(msg)
        repo = os.path.join(result, repo, 'checkpoint', 'main')

        result = os.path.join(self.__dest, self.label, self.resu)
        # check_dir called again here to get run_id (possibly date-hour)
        dest, msg = self.check_dir(node, result, d, "dest")
        if msg:
            studies.reporting(msg)
        dest = os.path.join(result, dest, 'checkpoint', 'main')

        cmd = self.__diff + ' ' + repo + ' ' + dest

        self.threshold = "default"
        if threshold != None:
            cmd += ' --threshold ' + threshold
            self.threshold = threshold

        if args != None:
            cmd += (" " + args)
            l = args.split()
            try:
                i = l.index('--threshold')
                self.threshold = l[i+1]
            except:
                pass

        l = subprocess.Popen(cmd,
                             shell=True,
                             executable=get_shell_type(),
                             stdout=subprocess.PIPE,
                             universal_newlines=True).stdout
        lines = l.readlines()

        # list of field differences
        tab = []
        # meshes have same sizes
        m_size_eq = True

        # studymanager compare log only for field of real values
        for i in range(len(lines)):
            # only select line with "Type" (english and french) and ';'
            # since this should be only true for heads of section
            if lines[i].find("Type") != -1 and lines[i].find(";") != -1:
                line = [x.replace("\""," ").strip() for x in lines[i].split(";")]
                name = line[0]
                info = [x.split(":") for x in line[1:]]
                info = [[x[0].strip(),x[1].strip()] for x in info]

                # section with at least 2 informations (location, type) after
                # their name, and of type r (real)
                if len(info) >= 2 and info[1][1] in ['r4', 'r8']:
                    # if next line contains size, this means sizes are different
                    if lines[i+1].find("Taille") != -1 or lines[i+1].find("Size") != -1:
                        m_size_eq = False
                        break
                    else:
                        line = [x.strip() for x in lines[i+1].split(";")]
                        vals = [x.split(":") for x in line]
                        vals = [[x[0].strip(),x[1].strip()] for x in vals]
                        tab.append([name.replace("_", "\_"),
                                    vals[1][1],
                                    vals[2][1],
                                    self.threshold])

        os.chdir(home)

        return tab, m_size_eq

    #---------------------------------------------------------------------------

    def run_ok(self, run_dir):
        """
        Check if a result directory contains an error file
        or if it doesn't contain a summary file
        """
        if not os.path.isdir(run_dir):
            print("Error: the result directory %s does not exist." % run_dir)
            sys.exit(1)

        msg = ""
        ok = True

        f_error = os.path.join(run_dir, 'error')
        if os.path.isfile(f_error):
            ok = False
            msg += "the result directory %s in case %s " \
                   "contains an error file." % (os.path.basename(run_dir),
                                                self.label)

        f_summary = os.path.join(run_dir, 'summary')
        if not os.path.isfile(f_summary):
            ok = False
            msg += "the result directory %s in case %s " \
                   "does not contain any summary file." \
                   % (os.path.basename(run_dir), self.label)

        return ok, msg

    #---------------------------------------------------------------------------

    def check_dir(self, node, result, rep, attr):
        """
        Check coherency between xml file of parameters and repository or
        destination.
        """
        msg = "Warning: "

        if not os.path.isdir(result):
            msg += "the directory %s " \
                   "does not exist." % (result)
            return None, msg

        # 1. The result directory is given
        if rep != "":
            # check if it exists
            rep_f = os.path.join(result, rep)
            if not os.path.isdir(rep_f):
                msg += "the result directory %s " \
                       "does not exist." % (rep_f)
                return None, msg

            run_ok = self.run_ok(rep_f)
            if not run_ok[0]:
                return None, msg+run_ok[1]

        # 2. The result directory must be found/read automatically;
        elif rep == "":
            # check if there is at least one result directory.
            if len(list(filter(nodot, os.listdir(result)))) == 0:
                msg += "there is no result directory in %s." % (result)
                return None, msg

            # if no run_id is specified in the xml file
            # only one result directory allowed in RESU
            if len(list(filter(nodot, os.listdir(result)))) > 1 \
               and self.run_id == "":
                msg += "there are several result directories in %s " \
                       "and no run id specified." % (result)
                return None, msg

            rep = self.run_id
            # if no run_id is specified in the xml file
            # the only result directory present in RESU is taken
            if rep == "":
                rep = list(filter(nodot, os.listdir(result)))[0]

            rep_f = os.path.join(result, rep)
            if not os.path.isdir(rep_f):
                msg += "the result directory %s " \
                       "does not exist." % (rep_f)
                return None, msg

            run_ok = self.run_ok(rep_f)
            if not run_ok[0]:
                return None, msg+run_ok[1]

            # 3. Update the file of parameters with the name of the result directory
            if node:
                self.__parser.setAttribute(node, attr, rep)

        return rep, None

    def check_dirs(self, node, repo, dest, reference=None):
        """
        Check coherency between xml file of parameters and repository and destination.
        """
        msg = None

        if repo != None:
            # build path to RESU directory with path to study and case label in the repo
            if reference:
                result = os.path.join(reference, self.label, self.resu)
            else:
                result = os.path.join(self.__repo, self.label, self.resu)
            rep, msg = self.check_dir(node, result, repo, "repo")

        if dest != None:
            # build path to RESU directory with path to study and case label in the dest
            result = os.path.join(self.__dest, self.label, self.resu)
            rep, msg = self.check_dir(node, result, dest, "dest")

        return msg

#===============================================================================
# Study class
#===============================================================================

class Study(object):
    """
    Create, run and compare all cases for a given study.
    """
    def __init__(self, pkg, parser, study, exe, dif, rlog, n_procs=None,
                 force_rm=False, force_overwrite=False, with_tags=None,
                 without_tags=None, debug=False):
        """
        Constructor.
          1. initialize attributes,
          2. build the list of the cases,
          3. build the list of the keywords from the runcase.
        @type parser: C{Parser}
        @param parser: instance of the parser
        @type study: C{String}
        @param study: label of the current study
        @type exe: C{String}
        @param exe: name of the solver executable: C{code_saturne} or C{neptune_cfd}.
        @type dif: C{String}
        @param dif: name of the diff executable: C{cs_io_dump -d}.
        @n_procs: C{int}
        @param n_procs: number of requested processors
        @type force_rm: C{True} or C{False}
        @param force_rm: remove always existing cases
        @type force_overwrite: C{True} or C{False}
        @param force_overwrite: overwrite files in dest by files in repo
        @type with_tags: C{List}
        @param with_tags: list of tags given at the command line
        @type without_tags: C{List}
        @param without_tags: list of tags given at the command line
        @type debug: C{True} or C{False}
        @param debug: if true, increase verbosity in stdout
        """
        # Initialize attributes
        self.__package  = pkg
        self.__parser   = parser
        self.__main_exe = exe
        self.__diff     = dif
        self.__log      = rlog
        self.__force_rm = force_rm
        self.__force_ow = force_overwrite
        self.__debug    = debug

        self.__repo = os.path.join(self.__parser.getRepository(),  study)
        self.__dest = os.path.join(self.__parser.getDestination(), study)

        if not os.path.isdir(self.__repo):
            print("Error: the directory %s does not exist" % self.__repo)
            sys.exit(1)

        self.label = study

        self.cases = []
        self.matplotlib_figures = []
        self.input_figures = []
        self.case_labels = []

        # get list of cases in study
        on_cases = parser.getStatusOnCasesLabels(study)
        if not on_cases:

            print("\n\n\nWarning: no case defined in %s study\n\n\n" % study)

        else:

            for data in self.__parser.getStatusOnCasesKeywords(self.label):

                # n_procs given in smgr command line overwrites n_procs by case
                if n_procs:
                    data['n_procs'] = str(n_procs)

                # check if every tag passed by option --with-tags belongs to
                # list of tags of the current case
                tagged = False
                if with_tags and data['tags']:
                    tagged = all(tag in data['tags'] for tag in with_tags)
                elif not with_tags:
                    tagged = True

                # check if none of tags passed by option --without-tags
                # belong to list of tags of the current case
                exclude = False
                if without_tags and data['tags']:
                    exclude = any(tag in data['tags'] for tag in without_tags)

                # do not append case if tags do not match
                if tagged and not exclude:
                    c = Case(pkg,
                             self.__log,
                             self.__diff,
                             self.__parser,
                             self.label,
                             data,
                             self.__repo,
                             self.__dest)
                    self.cases.append(c)
                    self.case_labels.append(c.label)

    #---------------------------------------------------------------------------

    def create_case(self, c, log_lines):
        """
        Create a case in a study
        """
        e = os.path.join(c.pkg.get_dir('bindir'), c.exe)
        if c.subdomains:
            os.mkdir(c.label)
            os.chdir(c.label)
            refdir = os.path.join(self.__repo, c.label)
            retval = 1
            for node in os.listdir(refdir):
                ref = os.path.join(self.__repo, c.label, node)
                if node in c.subdomains:
                    cmd = e + " create --case " + node \
                          + " --quiet --noref --copy-from " \
                          + ref
                    node_retval, t = run_studymanager_command(cmd, self.__log)
                    # negative retcode is kept
                    retval = min(node_retval,retval)
                elif os.path.isdir(ref):
                    shutil.copytree(ref, node, symlinks=True)
                else:
                    shutil.copy2(ref, node)
            c.update_runcase_path(c.label, destdir=self.__dest)
            os.chdir(self.__dest)
        else:
            cmd = e + " create --case " + c.label  \
                  + " --quiet --noref --copy-from "    \
                  + os.path.join(self.__repo, c.label)
            retval, t = run_studymanager_command(cmd, self.__log)
        if retval == 0:
            log_lines += ['    - create case: ' + c.label]
        else:
            log_lines += ['    - create case: %s --> FAILED' % c.label]

    #---------------------------------------------------------------------------

    def create_cases(self):
        """
        Create a single study with all its cases.
        """
        repbase = os.getcwd()

        # Create study if necessary
        if not os.path.isdir(self.__dest):
            # build instance of study class
            cr_study = cs_create.Study(self.__package,
                                       self.label,
                                       [],   # cases
                                       [],   # syrthes cases
                                       None, # cathare case
                                       None, # python case
                                       None, # copy
                                       False,# import_only
                                       False,# use ref
                                       0)    # quiet

            # TODO: copy-from for study. For now, an empty study
            # is created and cases are created one by one with
            # copy-from

            # create empty study
            cr_study.create()

            # Link meshes and copy other files
            ref = os.path.join(self.__repo, "MESH")
            if os.path.isdir(ref):
                l = os.listdir(ref)
                meshes = []
                for cpr in ["", ".gz"]:
                    for fmt in ["unv",
                                "med",
                                "ccm",
                                "cgns",
                                "neu",
                                "msh",
                                "des"]:
                        meshes += fnmatch.filter(l, "*." + fmt + cpr)
                des = os.path.join(self.__dest, "MESH")
                for m in l:
                    if m in meshes:
                        if sys.platform.startswith('win'):
                            shutil.copy2(os.path.join(ref, m), os.path.join(des, m))
                        else:
                            os.symlink(os.path.join(ref, m), os.path.join(des, m))
                    elif m != ".svn":
                        t = os.path.join(ref, m)
                        if os.path.isdir(t):
                            shutil.copytree(t, os.path.join(des, m))
                        elif os.path.isfile(t):
                            shutil.copy2(t, des)

            # Copy external scripts for post-processing
            ref = os.path.join(self.__repo, "POST")
            if os.path.isdir(ref):
                des = os.path.join(self.__dest, "POST")
                shutil.rmtree(des)
                shutil.copytree(ref, des, symlinks=True)

        # Change directory to destination directory
        os.chdir(self.__dest)

        log_lines = []
        for c in self.cases:
            if not os.path.isdir(c.label):
                self.create_case(c, log_lines);
            else:
                if self.__force_rm == True:
                    if self.__debug:
                        print("Warning: case %s exists in the destination "
                              "and will be overwritten." % c.label)
                    # Build short path to RESU dir. such as 'CASE1/RESU'
                    _dest_resu_dir = os.path.join(c.label, 'RESU')
                    if os.path.isdir(_dest_resu_dir):
                        shutil.rmtree(_dest_resu_dir)
                        os.makedirs(_dest_resu_dir)
                else:
                    if self.__debug:
                        print("Warning: case %s exists in the destination. "
                              "It won't be overwritten." % c.label)

                # if overwrite option enabled, overwrite content of DATA, SRC, SCRIPTS
                if self.__force_ow:
                    dirs_to_overwrite = ["DATA", "SRC", "SCRIPTS"]
                    self.overwriteDirectories(dirs_to_overwrite,
                                              case_label=c.label)
                    # update path in gui script
                    data_subdir = os.path.join(c.label, "DATA")
                    c.update_gui_script_path(data_subdir, self.__dest,
                                             xmlonly=False)

                    # update path in runcase script
                    scripts_subdir = os.path.join(c.label, "SCRIPTS")
                    c.update_runcase_path(scripts_subdir, self.__dest,
                                          xmlonly=False)

        os.chdir(repbase)

        if self.__force_ow:
            dirs_to_overwrite = ["POST", "MESH"]
            self.overwriteDirectories(dirs_to_overwrite)

        return log_lines

    #---------------------------------------------------------------------------

    def overwriteDirectories(self, dirs_to_overwrite, case_label=""):
        """
        Overwrite given directories in the Study tree.
        Label of case is an empty string by default.
        """
        for _dir in dirs_to_overwrite:
            ref = os.path.join(self.__repo, case_label, _dir)
            if os.path.isdir(ref):
                dest = os.path.join(self.__dest, case_label, _dir)
                for _ref_dir, _dirs, _files in os.walk(ref):
                    _dest_dir = _ref_dir.replace(ref, dest, 1)
                    for _file in _files:
                        _ref_file = os.path.join(_ref_dir, _file)
                        _dest_file = os.path.join(_dest_dir, _file)
                        if os.path.isfile(_dest_file):
                            os.remove(_dest_file)
                        shutil.copy2(_ref_file, _dest_dir)

    #---------------------------------------------------------------------------

    def getRunDirectories(self):
        list_cases = []
        list_dir   = []

        for case in self.cases:
            if case.is_run != "KO":
                list_cases.append(case.label)
                list_dir.append(case.run_dir)

        return " ".join(list_cases), " ".join(list_dir)

    #---------------------------------------------------------------------------

    def disable_case(self, case):
        try:
            self.cases.remove(case)
            self.case_labels.remove(case.label)
        except Exception:
            pass

        msg = "    - Case %s" % (case.label)
        if case.run_id != "":
            msg += ", run id %s" % (case.run_id)
        msg += " --> DISABLED"

        return msg

    #---------------------------------------------------------------------------

    def needs_report_detailed(self, postpro):
        """
        check if study needs a section in the detailed report
        (for figures, comparison or input)
        """
        # study has figures or input figures
        needs = self.matplotlib_figures or self.input_figures

        for case in self.cases:
            if case.is_compare == "done":
                needs = True
                break

            # handle the input nodes that are inside case nodes
            if case.plot == "on" and case.is_run != "KO":
                nodes = self.__parser.getChildren(case.node, "input")
                if nodes:
                    needs = True
                    break

        # handle the input nodes that are inside postpro nodes
        if postpro:
            script, label, nodes, args = self.__parser.getPostPro(self.label)
            for i in range(len(label)):
                if script[i]:
                    input_nodes = self.__parser.getChildren(nodes[i], "input")
                    if input_nodes:
                        needs = True
                        break

        return needs

#===============================================================================
# Studies class
#===============================================================================

class Studies(object):
    """
    Manage all Studies and all Cases described in the files of parameters.
    """
    def __init__(self, pkg, options, exe, dif):
        """
        Constructor.
          1. create if necessary the destination directory,
          2. initialize the parser and the plotter,
          3. build the list of the studies,
          4. start the report.
        @type options: C{Structure}
        @param options: structure the parameters options.
        @type exe: C{String}
        @param exe: name of the solver executable: C{code_saturne} or C{neptune_cfd}.
        @type dif: C{String}
        @param dif: name of the diff executable: C{cs_io_dump -d}.
        """

        # try to determine if current directory is a study one
        cwd = os.getcwd()
        is_study = isStudy(cwd)
        studyd = None
        studyp = None
        if is_study:
            # default study directory is current one
            studyp = cwd

        smgr = None

        # Create file of parameters

        filename = options.filename
        if options.create_xml and is_study:
            if filename == None:
                studyd = os.path.basename(studyp)
                filename = "smgr_" + studyd + ".xml"

            filepath = os.path.join(studyp, filename)
            smgr = create_base_xml_file(filepath, pkg)

            init_xml_file_with_study(smgr, studyp)

        elif options.create_xml and not is_study:
            msg =   "Can not create XML file of parameter:\n" \
                  + "current directory is apparently not a study (no MESH directory).\n"
            sys.exit(msg)

        if filename == None:
            msg =    "A file of parameters must be specified or created " \
                   + "for studymanager to run.\n" \
                   + "See help message and use '--file' or '--create-xml' option.\n"
            sys.exit(msg)

        # create a first xml parser only for
        #   the repository verification and
        #   the destination creation

        if os.path.isfile(filename):
            self.__parser = Parser(filename)
        else:
            msg = "Specified XML parameter file for studymanager does not exist.\n"
            sys.exit(msg)

        # call smgr xml backward compatibility

        if not smgr:
            smgr = XMLengine.Case(package=pkg, file_name=filename, studymanager=True)
            smgr['xmlfile'] = filename

            # minimal modification of xml for now
            smgr_xml_init(smgr).initialize(reinit_indices = False)
            smgr.xmlSaveDocument(prettyString=False)

        self.__xmlupdate = options.update_xml

        # set repository
        if len(options.repo_path) > 0:
            self.__parser.setRepository(options.repo_path)
        self.__repo = self.__parser.getRepository()
        if self.__repo:
            if not os.path.isdir(self.__repo):
                msg="Studies.__init__() >> self.__repo = {0}\n".format(self.__repo)
                sys.exit(msg+"Error: repository path is not valid.\n")
        else: # default value
            # if current directory is a study
            # set repository as directory containing the study
            if is_study:
                studyd = os.path.basename(studyp)
                self.__parser.setRepository(os.path.join(studyp,".."))
                self.__repo = self.__parser.getRepository()
            else:
                msg =   "Can not set a default repository directory:\n" \
                      + "current directory is apparently not a study (no MESH directory).\n" \
                      + "Add a repository path to the parameter file or use the command " \
                      + "line option (--repo=..).\n"
                sys.exit(msg)

        # set destination
        if self.__xmlupdate:
            self.__dest = self.__repo
        else:
            if len(options.dest_path) > 0:
                self.__parser.setDestination(options.dest_path)
            self.__dest = self.__parser.getDestination()
            if not self.__dest: # default value
                # if current directory is a study
                # set destination as a directory "../RUN_(study_name)
                if is_study and studyd != None:
                    self.__parser.setDestination(os.path.join(studyp,
                                                              "../RUN_"+studyd))
                    self.__dest = self.__parser.getDestination()
                else:
                  msg =   "Can not set a default destination directory:\n" \
                        + "current directory is apparently not a study (no MESH directory).\n" \
                        + "Add a destination path to the parameter file or use the command " \
                        + "line option (--dest=..).\n"
                  sys.exit(msg)

        # create if necessary the destination directory

        if not os.path.isdir(self.__dest):
            os.makedirs(self.__dest)

        # copy the xml file of parameters for update and restart

        file = os.path.join(self.__dest, os.path.basename(filename))
        try:
            shutil.copyfile(filename, file)
        except:
            pass

        # create a new parser, which is definitive and the plotter

        self.__parser  = Parser(file)
        self.__parser.setDestination(self.__dest)
        self.__parser.setRepository(self.__repo)
        if options.debug:
            print(" Studies >> Repository  >> ", self.__repo)
            print(" Studies >> Destination >> ", self.__dest)
        try:
            self.__plotter = Plotter(self.__parser)
        except Exception:
            self.__plotter = None

        # create list of restricting and excluding tags

        self.__with_tags = None
        if options.with_tags:
            with_tags = re.split(',', options.with_tags)
            self.__with_tags = [tag.strip() for tag in with_tags]
        self.__without_tags = None
        if options.without_tags:
            without_tags = re.split(',', options.without_tags)
            self.__without_tags = [tag.strip() for tag in without_tags]

        # build the list of the studies

        doc = os.path.join(self.__dest, options.log_file)
        self.__log = open(doc, "w")
        self.labels  = self.__parser.getStudiesLabel()
        self.studies = []
        for l in self.labels:
            self.studies.append( [l, Study(pkg, self.__parser, l, \
                                           exe, dif, self.__log, \
                                           options.n_procs, \
                                           options.remove_existing, \
                                           options.force_overwrite, \
                                           self.__with_tags, \
                                           self.__without_tags, \
                                           options.debug)] )
            if options.debug:
                print(" >> Append study ", l)

        # start the report

        self.report = os.path.join(self.__dest, "report.txt")
        self.reportFile = open(self.report, mode='w')
        self.reportFile.write('\n')

        # attributes

        self.__debug       = options.debug
        self.__quiet       = options.quiet
        self.__running     = options.runcase
        self.__n_iter      = options.n_iterations
        self.__compare     = options.compare
        self.__ref         = options.reference
        self.__postpro     = options.post
        self.__default_fmt = options.default_fmt
        # do not use tex in matplotlib (built-in mathtext is used instead)
        self.__dis_tex     = options.disable_tex
        # tex reports compilation with pdflatex
        self.__pdflatex    = not options.disable_pdflatex

        # in case of restart

        iok = 0
        for l, s in self.studies:
            for case in s.cases:
                if case.compute == 'on':
                   iok+=1
        if not iok:
            self.__running = False

        if self.__xmlupdate:
            os.remove(file)
            os.remove(doc)

    #---------------------------------------------------------------------------

    def getDestination(self):
        """
        @rtype: C{String}
        @return: destination directory of all studies.
        """
        if self.__dest == None:
            msg=" cs_studymanager_study.py >> Studies.getDestination()"
            msg+=" >> self.__dest is not set"
            sys.exit(msg)
        return self.__dest

    #---------------------------------------------------------------------------

    def getRepository(self):
        """
        @rtype: C{String}
        @return: repository directory of all studies.
        """
        if self.__repo == None:
            msg=" cs_studymanager_study.py >> Studies.getRepository()"
            msg+=" >> self.__repo is not set"
            sys.exit(msg)
        return self.__repo

    #---------------------------------------------------------------------------

    def reporting(self, msg, stdout=True, report=True, status=False):
        """
        Write message on standard output and/or in report.
        @type l: C{String}
        @param l: the sentence to be written.
        """
        s  = ""
        if not status:
            s = chr(10)

        if stdout and not self.__quiet:
            sys.stdout.write (msg + chr(13) + s)
            sys.stdout.flush()

        if report:
            self.reportFile.write(msg + '\n')
            self.reportFile.flush()

    #---------------------------------------------------------------------------

    def report_action_location(self, header_msg, destination=True):
        """
        Add a message to report/stdout at head of sections
        specifying if action is performed in dest or repo.
        """
        if destination:
            header_msg = header_msg + " (in destination)"
        else:
            header_msg = header_msg + " (in repository)"

        self.reporting(header_msg)

    #---------------------------------------------------------------------------

    def updateRepository(self, xml_only=False):
        """
        Update all studies and all cases.
        """
        for l, s in self.studies:
            self.reporting('  o Update repository: ' + l)
            for case in s.cases:
                self.reporting('    - update  %s' % case.label)
                case.update(xml_only)

        self.reporting('')

    #---------------------------------------------------------------------------

    def create_studies(self):
        """
        Create all studies and all cases.
        """
        for l, s in self.studies:
            create_msg = "  o Create study " + l
            dest = True
            self.report_action_location(create_msg, dest)
            log_lines = s.create_cases()
            for line in log_lines:
                self.reporting(line)

        self.reporting('')

    #---------------------------------------------------------------------------

    def test_compilation(self):
        """
        Compile sources of all runs with compute attribute at on.
        """
        iko = 0
        for l, s in self.studies:
            # build study dir. (in repo.)
            study_path = os.path.join(self.__repo, l)

            self.reporting('  o Compile study: ' + l)

            for case in s.cases:
                if case.compute == 'on':

                    # test compilation (logs are redirected to smgr log file)
                    is_compiled = case.test_compilation(study_path, self.__log)

                    # report
                    if is_compiled == "OK":
                        self.reporting('    - compile %s --> OK' % case.label)
                    elif is_compiled == "KO":
                        self.reporting('    - compile %s --> FAILED' % case.label)
                        iko+=1

        self.reporting('')

        if iko:
            self.reporting('Error: compilation failed for %s case(s).\n' % iko)
            sys.exit(1)

    #---------------------------------------------------------------------------

    def prepro(self, l, s, case):
        """
        Launch external additional scripts with arguments.
        """
        pre, label, nodes, args = self.__parser.getPrepro(case.node)
        if self.__debug:
            print(" >> prepro ", pre)
            print(" >> label ", label)
            print(" >> nodes ", nodes)
            print(" >> args  ", args)
        for i in range(len(label)):
            if pre[i]:
                # search if the script is in the MESH directory
                # if not, the script is searched in the directories
                # of the current case
                cmd = os.path.join(self.__dest, l, "MESH", label[i])
                if self.__debug:
                    print("Path to prepro script ", cmd)
                if not os.path.isfile(cmd):
                    filePath = ""
                    for root, dirs, fs in os.walk(os.path.join(self.__dest, l, case.label)):
                        if label[i] in fs:
                            filePath = root
                            break

                    cmd = os.path.join(filePath, label[i])

                if os.path.isfile(cmd):
                    sc_name = os.path.basename(cmd)
                    # ensure script is executable
                    set_executable(cmd)

                    cmd += " " + args[i]
                    cmd += " -c " + os.path.join(self.__dest, l, case.label)
                    repbase = os.getcwd()
                    os.chdir(os.path.join(self.__dest, l, "MESH"))

                    # Prepro external script often need install python directory
                    # and package python directory: code_saturne or neptune_cfd
                    p_dir = case.pkg.get_dir('pythondir')
                    pkg_dir = case.pkg.get_dir('pkgpythondir')
                    p_dirs = p_dir + ":" + pkg_dir

                    # if package is neptune_cfd, prepro script often needs
                    # code_saturne package python directory
                    cs_pkg_dir = None
                    if case.pkg.name == 'neptune_cfd':
                        cs_pkg_dir = os.path.join(pkg_dir, '../code_saturne')
                        cs_pkg_dir = os.path.normpath(cs_pkg_dir)
                        p_dirs = p_dirs + ":" + cs_pkg_dir

                    retcode, t = run_studymanager_command(cmd,
                                                          self.__log,
                                                          pythondir = p_dirs)
                    stat = "FAILED" if retcode != 0 else "OK"

                    os.chdir(repbase)

                    self.reporting('    - script %s --> %s (%s s)' % (stat, sc_name, t),
                                   stdout=True, report=False)

                    self.reporting('    - script %s --> %s (%s s)' % (stat, cmd, t),
                                   stdout=False, report=True)

                else:
                    self.reporting('    - script %s not found' % cmd)


    #---------------------------------------------------------------------------

    def run(self):
        """
        Update and run all cases.
        Warning, if the markup of the case is repeated in the xml file of parameters,
        the run of the case is also repeated.
        """
        for l, s in self.studies:
            self.reporting("  o Prepro scripts and runs for study: " + l)
            for case in s.cases:
                self.prepro(l, s, case)
                if self.__running:
                    if case.compute == 'on' and case.is_compiled != "KO":

                        if self.__n_iter is not None:
                            if case.subdomains:
                                case_dir = os.path.join(self.__dest, s.label, case.label,
                                                        case.subdomains[0], "DATA")
                            else:
                                case_dir = os.path.join(self.__dest, s.label, case.label, "DATA")
                            os.chdir(case_dir)
                            # Create a control_file in each case DATA
                            if not os.path.exists('control_file'):
                                control_file = open('control_file','w')
                                control_file.write("time_step_limit " + str(self.__n_iter) + "\n")
                                # Flush to ensure that control_file content is seen
                                # when control_file is copied to the run directory on all systems
                                control_file.flush()
                                control_file.close

                        self.reporting('    - running %s ...' % case.label,
                                       stdout=True, report=False, status=True)
                        error = case.run()
                        if case.is_time:
                            is_time = "%s s" % case.is_time
                        else:
                            is_time = "existed already"

                        if not error:
                            if not case.run_id:
                                self.reporting("    - run %s --> Warning suffix"
                                               " is not read" % case.label)

                            self.reporting('    - run %s --> OK (%s) in %s' \
                                           % (case.label, \
                                              is_time, \
                                              case.run_id))
                            self.__parser.setAttribute(case.node,
                                                       "compute",
                                                       "off")

                            # update dest="" attribute
                            n1 = self.__parser.getChildren(case.node, "compare")
                            n2 = self.__parser.getChildren(case.node, "script")
                            n3 = self.__parser.getChildren(case.node, "data")
                            n4 = self.__parser.getChildren(case.node, "probe")
                            n5 = self.__parser.getChildren(case.node, "resu")
                            n6 = self.__parser.getChildren(case.node, "input")
                            for n in n1 + n2 + n3 + n4 + n5 + n6:
                                if self.__parser.getAttribute(n, "dest") == "":
                                    self.__parser.setAttribute(n, "dest", case.run_id)
                        else:
                            if not case.run_id:
                                self.reporting('    - run %s --> FAILED (%s)' \
                                               % (case.label, is_time))
                            else:
                                self.reporting('    - run %s --> FAILED (%s) in %s' \
                                               % (case.label, \
                                                  is_time, \
                                                  case.run_id))

                        self.__log.flush()

        self.reporting('')

    #---------------------------------------------------------------------------

    def check_compare(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to make a comparison with a file which does not exist.
        """
        for l, s in self.studies:
            check_msg = "  o Check compare of study: " + l
            self.report_action_location(check_msg, destination)

            # reference directory passed in studymanager command line overwrites
            # destination in all cases (even if compare is defined by a compare
            # markup with a non empty destination)

            ref = None
            if self.__ref:
                ref = os.path.join(self.__ref, s.label)
            cases_to_disable = []
            for case in s.cases:
                if case.compare == 'on' and case.is_run != "KO":
                    compare, nodes, repo, dest, threshold, args = self.__parser.getCompare(case.node)
                    if compare:
                        is_checked = False
                        for i in range(len(nodes)):
                            if compare[i]:
                                is_checked = True
                                if destination == False:
                                    dest[i]= None
                                msg = case.check_dirs(nodes[i], repo[i], dest[i], reference=ref)
                                if msg:
                                    self.reporting(msg)
                                    cases_to_disable.append(case)

                    if not compare or not is_checked:
                        node = None
                        repo = ""
                        dest = ""
                        if destination == False:
                            dest = None
                        msg = case.check_dirs(node, repo, dest, reference=ref)
                        if msg:
                            self.reporting(msg)
                            cases_to_disable.append(case)

            for case in cases_to_disable:
                msg = s.disable_case(case)
                self.reporting(msg)

        self.reporting('')

    #---------------------------------------------------------------------------

    def compare_case_and_report(self, case, repo, dest, threshold, args, reference=None):
        """
        Compare the results for one computation and report
        """
        case.is_compare = "done"
        diff_value, m_size_eq = case.runCompare(self,
                                                repo, dest,
                                                threshold, args,
                                                reference=reference)

        case.diff_value += diff_value
        case.m_size_eq = case.m_size_eq and m_size_eq

        if args:
            s_args = 'with args: %s' % args
        else:
            s_args = 'default mode'

        if not m_size_eq:
            self.reporting('    - compare %s (%s) --> DIFFERENT MESH SIZES FOUND' % (case.label, s_args))
        elif diff_value:
            self.reporting('    - compare %s (%s) --> DIFFERENCES FOUND' % (case.label, s_args))
        else:
            self.reporting('    - compare %s (%s) --> NO DIFFERENCES FOUND' % (case.label, s_args))


    #---------------------------------------------------------------------------

    def compare(self):
        """
        Compare the results of the new computations with those from the Repository.
        """
        if self.__compare:
            for l, s in self.studies:
                self.reporting('  o Compare study: ' + l)
                # reference directory passed in studymanager command line overwrites
                # destination in all cases (even if compare is defined by a
                # compare markup with a non empty destination)
                ref = None
                if self.__ref:
                    ref = os.path.join(self.__ref, s.label)
                for case in s.cases:
                    if case.compare == 'on' and case.is_run != "KO":
                        is_compare, nodes, repo, dest, t, args = self.__parser.getCompare(case.node)
                        if is_compare:
                            for i in range(len(nodes)):
                                if is_compare[i]:
                                    self.compare_case_and_report(case,
                                                                 repo[i],
                                                                 dest[i],
                                                                 t[i],
                                                                 args[i],
                                                                 reference=ref)
                        if not is_compare or case.is_compare != "done":
                            repo = ""
                            dest = ""
                            t    = None
                            args = None
                            self.compare_case_and_report(case,
                                                         repo,
                                                         dest,
                                                         t,
                                                         args,
                                                         reference=ref)

        self.reporting('')

    #---------------------------------------------------------------------------

    def check_script(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to run a script with a file which does not exist.
        """
        scripts_checked = False
        for l, s in self.studies:
            # search for scripts to check before
            check_scripts = False
            for case in s.cases:
                script, label, nodes, args, repo, dest = \
                    self.__parser.getScript(case.node)
                if nodes:
                    check_scripts = True
                    break

            if not check_scripts:
                continue

            # if scripts have to be checked
            check_msg = "  o Check scripts of study: " + l
            self.report_action_location(check_msg, destination)

            scripts_checked = True

            cases_to_disable = []
            for case in s.cases:
                script, label, nodes, args, repo, dest = \
                    self.__parser.getScript(case.node)
                for i in range(len(nodes)):
                    if script[i] and case.is_run != "KO":
                        if destination == False:
                            dest[i] = None
                        msg = case.check_dirs(nodes[i], repo[i], dest[i])
                        if msg:
                            self.reporting(msg)
                            cases_to_disable.append(case)

            for case in cases_to_disable:
                msg = s.disable_case(case)
                self.reporting(msg)

        if scripts_checked:
            self.reporting('')

        return scripts_checked

    #---------------------------------------------------------------------------

    def scripts(self):
        """
        Launch external additional scripts with arguments.
        """
        for l, s in self.studies:
            self.reporting("  o Run scripts of study: " + l)
            for case in s.cases:
                script, label, nodes, args, repo, dest = self.__parser.getScript(case.node)
                for i in range(len(label)):
                    if script[i] and case.is_run != "KO":
                        cmd = os.path.join(self.__dest, l, "POST", label[i])
                        if os.path.isfile(cmd):
                            sc_name = os.path.basename(cmd)
                            # ensure script is executable
                            set_executable(cmd)

                            cmd += " " + args[i]
                            if repo[i]:
                                r = os.path.join(self.__repo,  l, case.label, "RESU", repo[i])
                                cmd += " -r " + r
                            if dest[i]:
                                d = os.path.join(self.__dest, l, case.label, "RESU", dest[i])
                                cmd += " -d " + d
                            retcode, t = run_studymanager_command(cmd, self.__log)
                            stat = "FAILED" if retcode != 0 else "OK"

                            self.reporting('    - script %s --> %s (%s s)' % (stat, sc_name, t),
                                           stdout=True, report=False)

                            self.reporting('    - script %s --> %s (%s s)' % (stat, cmd, t),
                                           stdout=True, report=False)
                        else:
                            self.reporting('    - script %s not found' % cmd)

        self.reporting('')

    #---------------------------------------------------------------------------

    def postpro(self):
        """
        Launch external additional scripts with arguments.
        """
        for l, s in self.studies:
            # fill results directories and ids for the cases of the current study
            # that were not run by the current studymanager command
            for case in s.cases:
                if case.is_run != "KO":
                    if case.run_dir == "":
                        resu = os.path.join(self.__dest, l, case.label, case.resu)
                        rep, msg = case.check_dir(None, resu, "", "dest")
                        if msg:
                            self.reporting(msg)
                            s.disable_case(case)
                        else:
                            case.run_id = rep
                            case.run_dir = os.path.join(resu, rep)

            script, label, nodes, args = self.__parser.getPostPro(l)
            if not label:
                continue

            self.reporting('  o Postprocessing cases of study: ' + l)
            for i in range(len(label)):
                if script[i]:
                    cmd = os.path.join(self.__dest, l, "POST", label[i])
                    if os.path.isfile(cmd):
                        sc_name = os.path.basename(cmd)
                        # ensure script is executable
                        set_executable(cmd)

                        list_cases, list_dir = s.getRunDirectories()
                        cmd += ' ' + args[i] + ' -c "' + list_cases + '" -d "' \
                               + list_dir + '" -s ' + l

                        self.reporting('    - running postpro %s' % sc_name,
                                       stdout=True, report=False, status=True)

                        retcode, t = run_studymanager_command(cmd, self.__log)
                        stat = "FAILED" if retcode != 0 else "OK"

                        self.reporting('    - postpro %s --> %s (%s s)' \
                                       % (stat, sc_name, t),
                                       stdout=True, report=False)

                        self.reporting('    - postpro %s --> %s (%s s)' \
                                       % (stat, cmd, t),
                                       stdout=False, report=True)
                    else:
                        self.reporting('    - postpro %s not found' % cmd)

        self.reporting('')

    #---------------------------------------------------------------------------

    def check_data(self, case, destination=True):
        """
        Check coherency between xml file of parameters and repository
        for data markups of a run.
        """
        for node in self.__parser.getChildren(case.node, "data"):
            plots, file, dest, repo = self.__parser.getResult(node)
            if destination == False:
                dest = None
            msg = case.check_dirs(node, repo, dest)
            if msg:
                self.reporting(msg)
                return False

        return True

    #---------------------------------------------------------------------------

    def check_probes(self, case, destination=True):
        """
        Check coherency between xml file of parameters and repository
        for probes markups of a run.
        """
        for node in self.__parser.getChildren(case.node, "probes"):
            file, dest, fig = self.__parser.getProbes(node)
            if destination == False:
                dest = None
            repo = None
            msg = case.check_dirs(node, repo, dest)
            if msg:
                self.reporting(msg)
                return False

        return True

    #---------------------------------------------------------------------------

    def check_input(self, case, destination=True):
        """
        Check coherency between xml file of parameters and repository
        for probes markups of a run.
        """
        for node in self.__parser.getChildren(case.node, "input"):
            file, dest, repo, tex = self.__parser.getInput(node)
            if destination == False:
                dest = None
            msg = case.check_dirs(node, repo, dest)
            if msg:
                self.reporting(msg)
                return False

        return True

    #---------------------------------------------------------------------------

    def check_plots_and_input(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to make a plot of a file which does not exist.
        """
        for l, s in self.studies:
            check_msg = "  o Check plots and input of study: " + l
            self.report_action_location(check_msg, destination)

            cases_to_disable = []
            for case in s.cases:
                if case.plot == "on" and case.is_run != "KO":

                    if not self.check_data(case, destination):
                        cases_to_disable.append(case)
                        continue

                    if not self.check_probes(case, destination):
                        cases_to_disable.append(case)
                        continue

                    if not self.check_input(case, destination):
                        cases_to_disable.append(case)
                        continue

            for case in cases_to_disable:
                msg = s.disable_case(case)
                self.reporting(msg)

        self.reporting('')

    #---------------------------------------------------------------------------

    def plot(self):
        """
        Plot data.
        """
        if self.__plotter:
            for l, s in self.studies:
                self.reporting('  o Plot study: ' + l)
                self.__plotter.plot_study(l, s,
                                          self.__dis_tex,
                                          self.__default_fmt)

        self.reporting('')

    #---------------------------------------------------------------------------

    def report_input(self, doc2, i_nodes, s_label, c_label=None):
        """
        Add input to report detailed.
        """
        for i_node in i_nodes:
            f, dest, repo, tex = self.__parser.getInput(i_node)
            doc2.appendLine("\\subsubsection{%s}" % f)

            if dest:
                d = dest
                dd = self.__dest
            elif repo:
                d = repo
                dd = self.__repo
            else:
                d = ""
                dd = ""

            if c_label:
                ff = os.path.join(dd, s_label, c_label, 'RESU', d, f)
            else:
                ff = os.path.join(dd, s_label, 'POST', d, f)

            if not os.path.isfile(ff):
                print("\n\nWarning: this file does not exist: %s\n\n" % ff)
            elif ff[-4:] in ('.png', '.jpg', '.pdf') or ff[-5:] == '.jpeg':
                doc2.addFigure(ff)
            elif tex == 'on':
                doc2.addTexInput(ff)
            else:
                doc2.addInput(ff)

    #---------------------------------------------------------------------------

    def build_reports(self, report1, report2):
        """
        @type report1: C{String}
        @param report1: name of the global report.
        @type report2: C{String}
        @param report2: name of the detailed report.
        @rtype: C{List} of C{String}
        @return: list of file to be attached to the report.
        """
        attached_files = []

        # First global report
        doc1 = Report1(self.__dest,
                       report1,
                       self.__log,
                       self.report,
                       self.__parser.write(),
                       self.__pdflatex)

        for l, s in self.studies:
            for case in s.cases:
                if case.diff_value or not case.m_size_eq:
                    is_nodiff = "KO"
                else:
                    is_nodiff = "OK"

                doc1.add_row(l,
                             case.label,
                             case.is_compiled,
                             case.is_run,
                             case.is_time,
                             case.is_compare,
                             is_nodiff)

        attached_files.append(doc1.close())

        # Second detailed report
        if self.__compare or self.__postpro:
            doc2 = Report2(self.__dest,
                           report2,
                           self.__log,
                           self.__pdflatex)

            for l, s in self.studies:
                if not s.needs_report_detailed(self.__postpro):
                    continue

                doc2.appendLine("\\section{%s}" % l)

                if s.matplotlib_figures or s.input_figures:
                    doc2.appendLine("\\subsection{Graphical results}")
                    for g in s.matplotlib_figures:
                        doc2.addFigure(g)
                    for g in s.input_figures:
                        doc2.addFigure(g)

                for case in s.cases:
                    if case.is_compare == "done":
                        run_id = None
                        if case.run_id != "":
                            run_id = case.run_id
                        doc2.appendLine("\\subsection{Comparison for case "
                                        "%s (run_id: %s)}"
                                        % (case.label, run_id))
                        if not case.m_size_eq:
                            doc2.appendLine("Repository and destination "
                                            "have apparently not been run "
                                            "with the same mesh (sizes do "
                                            "not match).")
                        elif case.diff_value:
                            doc2.add_row(case.diff_value, l, case.label)
                        elif self.__compare:
                            doc2.appendLine("No difference between the "
                                            "repository and the "
                                            "destination.")

                    # handle the input nodes that are inside case nodes
                    if case.plot == "on" and case.is_run != "KO":
                        nodes = self.__parser.getChildren(case.node, "input")
                        if nodes:
                            doc2.appendLine("\\subsection{Results for "
                                            "case %s}" % case.label)
                            self.report_input(doc2, nodes, l, case.label)

                # handle the input nodes that are inside postpro nodes
                if self.__postpro:
                    script, label, nodes, args = self.__parser.getPostPro(l)

                    needs_pp_input = False
                    for i in range(len(label)):
                        if script[i]:
                            input_nodes = \
                                self.__parser.getChildren(nodes[i], "input")
                            if input_nodes:
                                needs_pp_input = True
                                break

                    if needs_pp_input:
                        doc2.appendLine("\\subsection{Results for "
                                        "post-processing cases}")
                        for i in range(len(label)):
                            if script[i]:
                                input_nodes = \
                                    self.__parser.getChildren(nodes[i], "input")
                                if input_nodes:
                                    self.report_input(doc2, input_nodes, l)

            attached_files.append(doc2.close())

        return attached_files

    #---------------------------------------------------------------------------

    def getlabel(self):
        return self.labels

    #---------------------------------------------------------------------------

    def logs(self):
        try:
            self.reportFile.close()
        except:
            pass
        self.reportFile = open(self.report, mode='r')
        return self.reportFile.read()

    #---------------------------------------------------------------------------

    def __del__(self):
        try:
            self.__log.close()
        except:
            pass
        try:
            self.reportFile.close()
        except:
            pass

#-------------------------------------------------------------------------------
