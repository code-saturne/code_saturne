# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

from code_saturne.base.cs_exec_environment import get_shell_type, enquote_arg
from code_saturne.base.cs_compile import files_to_compile, compile_and_link
from code_saturne.base import cs_create
from code_saturne.base.cs_create import set_executable, create_local_launcher
from code_saturne.base import cs_run_conf

from code_saturne.model import XMLengine
from code_saturne.studymanager.cs_studymanager_pathes_model import PathesModel

from code_saturne.studymanager.cs_studymanager_parser import Parser
from code_saturne.studymanager.cs_studymanager_texmaker import Report

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
    msg = None
    if os.path.isfile(filepath):
        msg = "Can not create XML file of parameter:\n" \
            + filepath + " already exists."
        return None, msg

    # using xml engine from code_saturne GUI
    smgr = XMLengine.Case(package=pkg, studymanager=True)
    smgr['xmlfile'] = filename
    pm = PathesModel(smgr)

    # empty repo and dest
    pm.setRepositoryPath('')
    pm.setDestinationPath('')

    return smgr, msg

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
    """Try to determine if dirpath is a code_saturne study directory.
    """

    meshd = os.path.join(dirpath, 'MESH')
    is_study = os.path.isdir(meshd)

    return is_study

#-------------------------------------------------------------------------------

def isCase(dirpath):
    """Try to determine if dirpath is a code_saturne case directory.
    """

    # Verify that DATA folder exists with a xml file inside
    datad = os.path.join(dirpath, 'DATA')

    found_xml = False
    if os.path.isdir(datad):
        for elt in os.listdir(datad):
            if ".xml" in str(elt):
                found_xml = os.path.isfile(os.path.join(datad, elt))

    return found_xml

#===============================================================================
# Case class
#===============================================================================

class Case(object):
    def __init__(self, pkg, rlog, diff, parser, study, data, repo, dest):
        """
        @type data: C{Dictionary}
        @param data: contains all keyword and value read in the parameters file
        """
        self.__log_file  = rlog
        self.__diff      = diff
        self.__parser    = parser
        self.__data      = data
        self.__repo      = repo
        self.__dest      = dest

        self.pkg         = pkg
        self.study       = study

        self.node        = data['node']
        self.label       = data['label']
        self.compute     = data['compute']
        self.plot        = data['post']
        self.run_id      = data['run_id']
        self.compare     = data['compare']
        self.n_procs     = data['n_procs']
        self.depends     = data['depends']

        self.parametric  = data['parametric']
        self.notebook    = data['notebook']
        self.kw_args     = data['kw_args']

        self.is_compiled = "not done"
        self.is_run      = "not done"
        self.is_time     = None
        self.is_plot     = "not done"
        self.is_compare  = "not done"
        self.disabled    = False
        self.threshold   = "default"
        self.diff_value  = [] # list of differences (in case of comparison)
        self.m_size_eq   = True # mesh sizes equal (in case of comparison)
        self.subdomains  = None
        self.level       = None # level of the node in the dependency graph

        # Run_dir and Title are based on study, label and run_id
        self.resu = "RESU"
        self.run_dir = os.path.join(self.__dest, self.label, self.resu,
                                    self.run_id)
        self.title = study + "/" + self.label + "/" +  self.resu + "/" \
                   + self.run_id

        # Check for coupling
        # TODO: use run.cfg info, so as to allow another coupling parameters
        #       path ('run.cfg' is fixed, 'coupling_parameters.py' is not).

        coupling = False
        run_conf = None
        run_config_path = os.path.join(self.__repo, self.label, "run.cfg")
        if os.path.isfile(run_config_path):
            run_conf = cs_run_conf.run_conf(run_config_path, package=self.pkg)
            if run_conf.get("setup", "coupled_domains") != None:
                coupling = True

        if coupling:
            self.resu = "RESU_COUPLING"

            # Apply coupling parameters information

            from code_saturne.base import cs_case_coupling

            coupled_domains = run_conf.get_coupling_parameters()

            self.subdomains = []
            for d in coupled_domains:
                if d['solver'].lower() in ('code_saturne', 'neptune_cfd'):
                    self.subdomains.append(d['domain'])

            self.run_dir = os.path.join(self.__dest, self.label, self.resu,
                                        self.run_id, self.subdomains[0])

        self.exe = os.path.join(pkg.get_dir('bindir'),
                                pkg.name + pkg.config.shext)

    #---------------------------------------------------------------------------

    def prepare_run_folder(self):
        """
        Prepare a run folder in destination directory run_dir
        """
        log_lines = []
        e = os.path.join(self.pkg.get_dir('bindir'), self.exe)
        home = os.getcwd()
        os.chdir(self.__dest)
        if self.subdomains:
            if not os.path.isdir(self.label):
                os.mkdir(self.label)
            os.chdir(self.label)
            refdir = os.path.join(self.__repo, self.label)
            retval = 1
            resu_coupling = None
            for node in os.listdir(refdir):
                ref = os.path.join(self.__repo, self.label, node)

                # only loop on code_saturne subdomains
                if node in self.subdomains:

                    # generate folder in dest/STUDY/CASE/RESU_COUPLING/
                    cmd = e + " run --stage --case " + refdir \
                            + " --dest " + self.__dest \
                            + " --id " + self.run_id

                    if self.notebook:
                        cmd += " --notebook-args " + self.notebook

                    if self.parametric:
                        cmd += " --parametric-args " + '"' + self.parametric + '"'

                    if self.kw_args:
                        if self.kw_args.find(" ") < 0:
                            self.kw_args += " "  # workaround for arg-parser issue
                        cmd += " --kw-args " + '"' + self.kw_args + '"'

                    # create run_case.log in dest/STUDY (will be moved later)
                    file_name = os.path.join(self.__dest, "run_case.log")
                    log_run = open(file_name, mode='w')

                    node_retval, t = run_studymanager_command(cmd, log_run)

                    # move run_case.log in run_dir
                    if os.path.isdir(self.run_dir):
                        os.replace(file_name, os.path.join(self.run_dir,
                                                           "run_case.log"))
                    else:
                        err_file = os.path.join(self.__dest, "run_" + self.label
                                 + "_" + self.run_id + ".log")
                        os.replace(file_name, err_file)

                    # negative retcode is kept
                    retval = min(node_retval,retval)

                elif not os.path.isdir(ref):
                    shutil.copy2(ref, node)

            create_local_launcher(self.pkg, self.__dest)
            os.chdir(self.__dest)
        else:
            refdir = os.path.join(self.__repo, self.label)

            cmd = e + " run --stage --case " + refdir \
                    + " --dest " + self.__dest \
                    + " --id " + self.run_id

            if self.notebook:
                cmd += " --notebook-args " + self.notebook

            if self.parametric:
                cmd += " --parametric-args " + '"' + self.parametric + '"'

            if self.kw_args:
                if self.kw_args.find(" ") < 0:
                    self.kw_args += " "  # workaround for arg-parser issue
                cmd += " --kw-args " + '"' + self.kw_args + '"'

            # create run_case.log in dest/STUDY/CASE
            file_name = os.path.join(self.__dest, "run_case.log")
            log_run = open(file_name, mode='w')

            retval, t = run_studymanager_command(cmd, log_run)

            # move run_case.log in run_dir
            if os.path.isdir(self.run_dir):
                os.replace(file_name, os.path.join(self.run_dir,
                                                   "run_case.log"))
            else:
                err_file = os.path.join(self.__dest, "run_" + self.label
                         + "_" + self.run_id + ".log")
                os.replace(file_name, err_file)

        if retval == 0:
            log_lines += ['      * prepare run folder: ' + self.title]

        else:
            log_lines += ['      * prepare run folder: %s --> FAILED'
                          % self.title]
            log_lines += ['        - see run_CASE_run_id.log in dest/STUDY']
            self.compute = "off"
            self.post = "off"
            self.compare = "off"

        os.chdir(home)
        return log_lines

    #---------------------------------------------------------------------------

    def __update_setup(self, subdir):
        """
        Update setup file in the Repository.
        """
        # Load setup.xml file in order to update it
        # with the __backwardCompatibility method.

        from code_saturne.model.XMLengine import Case

        msg = None
        for fn in os.listdir(os.path.join(self.__repo, subdir, "DATA")):
            fp = os.path.join(self.__repo, subdir, "DATA", fn)
            if os.path.isfile(fp):
                try:
                    fd = os.open(fp , os.O_RDONLY)
                    f = os.fdopen(fd)
                    l = f.readline()
                    f.close()
                except Exception:
                    continue
                xml_type = None
                if l.startswith('''<?xml version="1.0" encoding="utf-8"?><Code_Saturne_GUI'''):
                    xml_type = 'code_saturne'
                elif l.startswith('''<?xml version="1.0" encoding="utf-8"?><NEPTUNE_CFD_GUI'''):
                    xml_type = 'neptune_cfd'
                else:
                    continue
                try:
                    case = Case(package = self.pkg, file_name = fp)
                except:
                    msg = "Parameters file reading error.\n" \
                        + "This file is not in accordance with XML specifications."
                    return msg

                case['xmlfile'] = fp
                case.xmlCleanAllBlank(case.xmlRootNode())

                if xml_type == 'code_saturne':
                    from code_saturne.model.XMLinitialize import XMLinit as cs_solver_xml_init
                    cs_solver_xml_init(case).initialize()
                elif xml_type == 'neptune_cfd':
                    try:
                        from code_saturne.model.XMLinitializeNeptune import XMLinitNeptune as nc_solver_xml_init
                        nc_solver_xml_init(case).initialize()
                    except ImportError:
                        # Avoid completely failing an update of cases with
                        # mixed solver types when neptune_cfd is not available
                        # (will fail if really trying to run those cases)
                        msg = "Failed updating a neptune_cfd XML file as " \
                            + "neptune_cfd is not available."
                        pass

                case.xmlSaveDocument()

        return msg

    #---------------------------------------------------------------------------

    def update(self):
        """
        Update path for the script in the Repository.
        """
        # Load the setup file in order to update it
        # with the __backwardCompatibility method.

        if self.subdomains:
            cdirs = []
            for d in self.subdomains:
                cdirs.append(os.path.join(self.label, d))
        else:
            cdirs = (self.label,)

        error = None
        for d in cdirs:
            error = self.__update_setup(d)

        return error

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
                retcode += compile_and_link(self.pkg, self.pkg.solver, s, dest_dir,
                                            stdout=log, stderr=log)

        if retcode > 0:
            self.is_compiled = "KO"

        return self.is_compiled

    #---------------------------------------------------------------------------

    def run(self, resource_name=None):
        """
        Check if a run with same result subdirectory name exists
        and launch run if not.
        """
        home = os.getcwd()
        os.chdir(os.path.join(self.__dest, self.label))

        e = os.path.join(self.pkg.get_dir('bindir'), self.exe)

        refdir = os.path.join(self.__repo, self.label)

        # After the stage within run_id folder in destination
        # do initialize, execute and finalize steps
        run_cmd = e + " run --no-stage" \
                + " --case " + refdir \
                + " --dest " + self.__dest \
                + " --id " + self.run_id

        if self.kw_args:
            if self.kw_args.find(" ") < 0:
                self.kw_args += " "  # workaround for arg-parser issue
            run_cmd += " --kw-args " + '"' + self.kw_args + '"'

        n_procs = self.__data['n_procs']
        if n_procs:
            run_cmd += " -n " + n_procs

        if resource_name:
            run_cmd += " --with-resource " + resource_name

        # append run_case.log in run_dir
        file_name = os.path.join(self.run_dir, "run_case.log")
        log_run = open(file_name, mode='a')

        error, self.is_time = run_studymanager_command(run_cmd, log_run)

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
        if not os.path.isfile(repo):
            repo += '.csc'

        result = os.path.join(self.__dest, self.label, self.resu)
        # check_dir called again here to get run_id (possibly date-hour)
        dest, msg = self.check_dir(node, result, d, "dest")
        if msg:
            studies.reporting(msg)
        dest = os.path.join(result, dest, 'checkpoint', 'main.csc')

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
            msg = "Error: the result directory %s does not exist." % run_dir
            return False, msg

        msg = ""
        ok = True

        f_error = os.path.join(run_dir, 'error')
        if os.path.isfile(f_error):
            ok = False
            msg += "the result directory %s in case %s " \
                   "contains an error file." % (os.path.basename(run_dir),
                                                self.title)

        f_summary = os.path.join(run_dir, 'summary')
        if not os.path.isfile(f_summary):
            ok = False
            msg += "the result directory %s in case %s " \
                   "does not contain any summary file." \
                   % (os.path.basename(run_dir), self.title)

        return ok, msg

    #---------------------------------------------------------------------------

    def check_dir(self, node, result, rep, attr):
        """
        Check coherency between xml file of parameters and repository or
        destination.
        """
        msg = "Warning: "

        if not os.path.isdir(result):
            msg += "the directory %s does not exist." %(result)
            return None, msg

        # 1. The result directory is given
        if rep != "":
            # check if it exists
            rep_f = os.path.join(result, rep)
            if not os.path.isdir(rep_f):
                msg += "the result directory %s does not exist." %(rep_f)
                return None, msg

            run_ok = self.run_ok(rep_f)
            if not run_ok[0]:
                return None, msg+run_ok[1]

        # 2. The result directory must be found/read automatically;
        elif rep == "":
            # check if there is at least one result directory.
            if len(list(filter(nodot, os.listdir(result)))) == 0:
                msg += "there is no result directory in %s." %(result)
                return None, msg

            # if no run_id is specified in the xml file
            # only one result directory allowed in RESU
            if len(list(filter(nodot, os.listdir(result)))) > 1 \
               and self.run_id == "run1":
                msg += "there are several result directories in %s " \
                       "and no run id specified." % (result)
                return None, msg

            rep = self.run_id
            # if no run_id is specified in the xml file
            # the only result directory present in RESU is taken
            if rep == "run1":
                rep = list(filter(nodot, os.listdir(result)))[0]

            rep_f = os.path.join(result, rep)
            if not os.path.isdir(rep_f):
                msg += "the result directory %s does not exist." %(rep_f)
                return None, msg

            run_ok = self.run_ok(rep_f)
            if not run_ok[0]:
                return None, msg+run_ok[1]

            # 3. Update the file of parameters with the name of the result directory
            if node:
                self.__parser.setAttribute(node, attr, rep)

        return rep, None

    #---------------------------------------------------------------------------

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

    #---------------------------------------------------------------------------

    def disable(self):

        self.disabled = True
        msg = "    - Case %s --> DISABLED" %(self.title)

        return msg

#===============================================================================
# Study class
#===============================================================================

class Study(object):
    """
    Create, run and compare all cases for a given study.
    """
    def __init__(self, pkg, parser, study, exe, dif, rlog, n_procs=None,
                 with_tags=None, without_tags=None):
        """
        Constructor.
          1. initialize attributes,
          2. build the list of the cases,
          3. build the list of the keywords from the case
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
        @type with_tags: C{List}
        @param with_tags: list of tags given at the command line
        @type without_tags: C{List}
        @param without_tags: list of tags given at the command line
        """
        # Initialize attributes
        self.__package  = pkg
        self.__parser   = parser
        self.__main_exe = exe
        self.__diff     = dif
        self.__log_file = rlog

        # read repository and destination from smgr file
        # based on working directory if information is not available
        try:
            self.__repo = os.path.join(self.__parser.getRepository(),  study)
        except:
            self.__repo = os.path.join(os.getcwd(),  study)
        try:
            self.__dest = os.path.join(self.__parser.getDestination(), study)
        except:
            self.__dest = self.__repo

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

                # TODO: move tag's filtering to the graph level
                # not done for now as the POST step needs tags also
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
                             self.__log_file,
                             self.__diff,
                             self.__parser,
                             self.label,
                             data,
                             self.__repo,
                             self.__dest)
                    self.cases.append(c)
                    self.case_labels.append(c.label)

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
        self.__log = None
        self.__log_compile = None
        self.__log_file = None

        # Store options

        self.__pkg            = pkg
        self.__create_xml     = options.create_xml
        self.__update_smgr    = options.update_smgr
        self.__update_setup   = options.update_setup
        self.__force_rm       = options.remove_existing
        self.__disable_ow     = options.disable_overwrite
        self.__debug          = options.debug
        self.__n_procs        = options.n_procs
        self.__filter_level   = options.filter_level
        self.__filter_n_procs = options.filter_n_procs
        # Use the provided resource name if forced
        self.__resource_name  = options.resource_name
        self.__quiet          = options.quiet
        self.__running        = options.runcase
        self.__n_iter         = options.n_iterations
        self.__compare        = options.compare
        self.__ref            = options.reference
        self.__postpro        = options.post
        self.__default_fmt    = options.default_fmt
        # do not use tex in matplotlib (built-in mathtext is used instead)
        self.__dis_tex        = options.disable_tex
        # tex reports compilation with pdflatex
        self.__pdflatex       = not options.disable_pdflatex

        # Create file of parameters

        filename = options.filename
        if self.__create_xml and is_study:
            if filename is None:
                studyd = os.path.basename(studyp)
                filename = "smgr.xml"

            filepath = os.path.join(studyp, filename)
            smgr, error = create_base_xml_file(filepath, self.__pkg)

            if error:
                self.reporting(error, report=False, exit=True)
            else:
                init_xml_file_with_study(smgr, studyp)
                self.reporting(" ", report=False)
                self.reporting(" New smgr.xml file was created succesfully",
                               report=False)
                self.reporting(" ", report=False)
                return

        elif self.__create_xml and not is_study:
            msg = "Error: can not create XML file of parameter:\n" \
                + "current directory is apparently not a study (no MESH " \
                + "directory)."
            self.reporting(msg, report=False, exit=True)

        if filename is None:
            msg = "Error: a file of parameters must be specified or created" \
                + " for studymanager to run.\n" \
                + "See help message and use '--file' or '--create-xml'" \
                + " option."
            self.reporting(msg, report=False, exit=True)

        # create a first smgr parser only for
        #   the repository verification and
        #   the destination creation

        if os.path.isfile(filename):
            self.__parser = Parser(filename)
        else:
            msg = "Error: specified XML parameter file for studymanager does" \
                + " not exist."
            self.reporting(msg, report=False, exit=True)

        # call smgr xml backward compatibility

        if not smgr:
            smgr = XMLengine.Case(package=self.__pkg, file_name=filename,
                                  studymanager=True)
            smgr['xmlfile'] = filename

            # minimal modification of xml for now
            smgr_xml_init(smgr).initialize()
            if self.__update_smgr:
                smgr.xmlSaveDocument(prettyString=False)

                self.reporting(" ", report=False)
                self.reporting(" SMGR parameter file was updated succesfully",
                               report=False)
                self.reporting(" Note that update is not necessary before the"\
                               + " run step", report=False)
                self.reporting(" ", report=False)

        # set repository
        if len(options.repo_path) > 0:
            self.__parser.setRepository(options.repo_path)
        self.__repo = self.__parser.getRepository()
        if self.__repo:
            if not os.path.isdir(self.__repo):
                msg = "Error: repository path is not valid: " + self.__repo
                self.reporting(msg, report=False, exit=True)
        else: # default value
            # if current directory is a study
            # set repository as directory containing the study
            if is_study:
                studyd = os.path.basename(studyp)
                self.__parser.setRepository(os.path.join(studyp,".."))
                self.__repo = self.__parser.getRepository()
            else:
                msg = "Can not set a default repository directory:\n" \
                    + "current directory is apparently not a study (no MESH" \
                    + " directory).\n" \
                    + "Add a repository path to the parameter file or use" \
                    + " the command line option (--repo=..)."
                self.reporting(msg, report=False, exit=True)

        # set destination
        self.__dest = None
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
                msg = "Can not set a default destination directory:\n" \
                    + "current directory is apparently not a study (no MESH" \
                    + " directory).\n" \
                    + "Add a destination path to the parameter file or use" \
                    + " the command line option (--dest=..).\n"
                self.reporting(msg, report=False, exit=True)

        if options.runcase or options.compare or options.post:

            # create if necessary the destination directory

            if not os.path.isdir(self.__dest):
                os.makedirs(self.__dest)

            # copy the smgr file in destination for update and restart

            file = os.path.join(self.__dest, os.path.basename(filename))
            try:
                shutil.copyfile(filename, file)
            except:
                pass

            # create definitive parser for smgr file in destination

            self.__parser  = Parser(file)
            self.__parser.setDestination(self.__dest)
            self.__parser.setRepository(self.__repo)
            if options.debug:
                msg = " Studies repository: " + self.__repo + "\n" \
                  + " Studies destination: " + self.__dest
                self.reporting(msg, report=False)

            # create studymanager log file

            self.__log_name = os.path.join(self.__dest, options.log_file)
            self.__log_file = open(self.__log_name, "w")

            # create plotter and post log file

            if options.post:
                try:
                    self.__plotter = Plotter(self.__parser)
                except Exception:
                    self.__plotter = None
                self.__log_post_name = os.path.join(self.__dest,
                                                    "smgr_post_pro.log")
                self.__log_post_file = open(self.__log_post_name, "w")

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

        self.labels  = self.__parser.getStudiesLabel()
        self.studies = []
        for l in self.labels:
            self.studies.append( [l, Study(self.__pkg, self.__parser, l, \
                                           exe, dif, self.__log_file, \
                                           options.n_procs, \
                                           self.__with_tags, \
                                           self.__without_tags,)] )
            if options.debug:
                self.reporting(" Append study:" + l, report=False)

        # in case of restart

        iok = 0
        for l, s in self.studies:
            for case in s.cases:
                if case.compute == 'on':
                   iok+=1
        if not iok:
            self.__running = False

        # Handle relative paths:
        if self.__ref:
            if not os.path.isabs(self.__ref):
                self.__ref = os.path.join(os.getcwd(), self.__ref)

    #---------------------------------------------------------------------------

    def getDestination(self):
        """
        @rtype: C{String}
        @return: destination directory of all studies.
        """
        if self.__dest is None:
            msg = " cs_studymanager_study.py >> Studies.getDestination()" \
                + " >> self.__dest is not set"
            self.reporting(msg, exit=True)
        return self.__dest

    #---------------------------------------------------------------------------

    def getRepository(self):
        """
        @rtype: C{String}
        @return: repository directory of all studies.
        """
        if self.__repo is None:
            msg = " cs_studymanager_study.py >> Studies.getRepository()" \
                + " >> self.__repo is not set"
            self.reporting(msg, exit=True)
        return self.__repo

    #---------------------------------------------------------------------------

    def reporting(self, msg, stdout=True, report=True, status=False,
                  exit=False):
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
            self.__log_file.write(msg + '\n')
            self.__log_file.flush()

        if exit:
            sys.exit(msg)

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

    def updateSetup(self):
        """
        Update setup files in all cases.
        """

        for l, s in self.studies:
            self.reporting('  o In repository: ' + l, report=False)
            for case in s.cases:
                self.reporting('    - update setup file in %s' % case.label,
                               report=False)
                error = case.update()
                if error:
                   self.reporting(error, report=False, exit=True)

        self.reporting('',report=False)

    #---------------------------------------------------------------------------

    def create_studies(self, run_step):
        """
        Create all studies and all cases.
        """

        self.reporting("  o Create all studies and run folders")
        if self.__force_rm:
            self.reporting("    /!\ All earlier run folders are erased (option"
                           " --rm activated)")
        else:
            self.reporting("    /!\ All earlier run folders will not be erased."
                           " Use --rm option to do so.")

        study_list = []
        case_list = []

        for case in self.graph.graph_dict:

            # first step: create study of the case if necessary
            study = case.study
            if study not in study_list:
                self.create_study(study)
                study_list.append(study)

            # only run step requires to create cases
            if run_step:
                # second step: clean RESU folder if necessary
                case_name = case.study + "/" + case.label
                if case_name not in case_list and self.__force_rm:
                    case_list.append(case_name)
                    # Build short path to RESU dir. such as 'CASE1/RESU'
                    _dest_resu_dir = os.path.join(self.__dest, case.study,
                                                  case.label, 'RESU')
                    if os.path.isdir(_dest_resu_dir):
                        if os.listdir(_dest_resu_dir):
                            shutil.rmtree(_dest_resu_dir)
                            os.makedirs(_dest_resu_dir)
                            self.reporting("  All earlier results in case %s/RESU "
                                           "are removed (option --rm activated)"
                                           %case.label)

                # thrid step: prepare run folder
                log_lines = case.prepare_run_folder()

                for line in log_lines:
                    self.reporting(line)

        self.reporting('')

    #---------------------------------------------------------------------------

    def create_study(self, study):

        dest_study = os.path.join(self.__dest, study)
        repo_study = os.path.join(self.__repo, study)

        new_study = False
        home = os.getcwd()
        os.chdir(self.__dest)
        # Create study if necessary
        if not os.path.isdir(dest_study):
            new_study = True
            # build instance of study class
            cr_study = cs_create.study(self.__pkg, study)    # quiet

            create_msg = "    - Create study " + study
            self.report_action_location(create_msg)

            # create empty study
            cr_study.create()

        # write content of MESH and POST is the study is new
        # overwrite them if not disabled
        if new_study or not self.__disable_ow:

            if not new_study:
                self.reporting("  /!\ POST folder is overwritten in %s"
                               " use option --dow to disable overwrite" %study)

            # Copy external scripts for post-processing
            ref = os.path.join(repo_study, "POST")
            if os.path.isdir(ref):
                des = os.path.join(dest_study, "POST")
                shutil.rmtree(des)
                shutil.copytree(ref, des, symlinks=True)
        os.chdir(home)

    #---------------------------------------------------------------------------

    def dump_graph(self):
        """
        Dump dependency graph based on all studies and all cases.
        Can be limited to a sub graph is filters an tags are given
        """

        filter_level   = self.__filter_level
        filter_n_procs = self.__filter_n_procs

        self.reporting("  o Dump dependency graph with option :")
        self.reporting("     - level=" + str(filter_level))
        self.reporting("     - n_procs=" + str(filter_n_procs))

        # create the global graph with all cases of all studies without filtering
        global_graph = dependency_graph()
        for l, s in self.studies:
            for case in s.cases:
                if case.compute == "on":
                    global_graph.add_node(case)

        # extract the sub graph based on filters and tags
        if filter_level is not None or filter_n_procs is not None:

            sub_graph = dependency_graph()
            for node in global_graph.graph_dict:

                # check if the level of the case is the targeted one
                # only effective if filter_level is prescribed
                target_level = True
                if filter_level is not None:
                    target_level = node.level is int(filter_level)

                # check if the number of procs of the case is the targeted one
                # only effective if filter_n_procs is prescribed
                target_n_procs = True
                if filter_n_procs is not None:
                    target_n_procs = node.n_procs is int(filter_n_procs)

                if target_level and target_n_procs:
                    sub_graph.add_node(node)

            self.graph = sub_graph

        else:
            self.graph = global_graph

        self.reporting('')

    #---------------------------------------------------------------------------

    def test_compilation(self):
        """
        Compile sources of all runs with compute attribute at on.
        """

        # create compilation log file in repository
        log_comp_name = os.path.join(self.__repo, "smgr_compilation.log")
        log_comp_file = open(log_comp_name, "w")

        iko = 0
        for l, s in self.studies:
            self.reporting('  o Compile study: ' + l + ' (in repository)',
                           report=False)
            for case in s.cases:

                # build case dir. (in repo.)
                study_path = os.path.join(self.__repo, case.study)

                if case.compute == 'on':

                    # test compilation
                    is_compiled = case.test_compilation(study_path,
                                                        log_comp_file)

                    # report
                    if is_compiled == "OK":
                        self.reporting('    - compile %s --> OK' %case.title,
                                       report=False)
                    elif is_compiled == "KO":
                        self.reporting('    - compile %s --> FAILED' %case.title,
                                       report=False)
                        iko+=1

        self.reporting('',report=False)

        log_comp_file.close()

        if iko:
            self.reporting('  Error: compilation failed for %s case(s).\n'
                           %iko, report=False, exit=False)
        self.reporting('  Compilation logs available in  %s \n'
                       %log_comp_name, report=False, exit=False)

    #---------------------------------------------------------------------------

    def check_prepro(self, case):
        """
        Launch external additional scripts with arguments.
        """
        pre, label, nodes, args = self.__parser.getPrepro(case.node)
        iko = 0
        for i in range(len(label)):
            if pre[i]:
                cmd = os.path.basename(label[i])
                self.reporting('    - script %s --> FAILED (%s)' % (cmd),
                               stdout=True, report=False)
                self.reporting('    - script %s --> FAILED (%s)' % (cmd),
                               stdout=False, report=True)
                iko += 1

        if iko:
            self.reporting('  Error: "prepro" tag present for %s case(s).\n'
                           %iko, report=False, exit=False)
            self.reporting('        "prepro" must be updated to one of\n' +
                           '           "notebook_args", "parametric_args", or "kw_args"\n',
                           report = False, exit=False)

    #---------------------------------------------------------------------------

    def run(self):
        """
        Update and run all cases.
        Warning, if the markup of the case is repeated in the xml file of parameters,
        the run of the case is also repeated.
        """

        home = os.getcwd()

        self.reporting("  o Run all cases")

        for case in self.graph.graph_dict:
            self.check_prepro(case)
            if self.__running:
                if case.compute == 'on' and case.is_compiled != "KO":

                    if self.__n_iter is not None:
                        os.chdir(case.run_dir)
                        # Create a control_file in run folder
                        if not os.path.exists('control_file'):
                            control_file = open('control_file','w')
                            control_file.write("time_step_limit " + str(self.__n_iter) + "\n")
                            # Flush to ensure that control_file content is seen
                            # when control_file is copied to the run directory on all systems
                            control_file.flush()
                            control_file.close

                    self.reporting('    - running %s ...' % case.title,
                                   stdout=True, report=False, status=True)

                    error = case.run(resource_name = self.__resource_name)
                    if case.is_time:
                        is_time = "%s s" % case.is_time
                    else:
                        is_time = "existed already"

                    if not error:
                        self.reporting('    - run %s --> OK (%s)' \
                                       % (case.title, is_time))
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
                        self.reporting('    - run %s --> FAILED (%s)' \
                                       % (case.title, is_time))
                        self.reporting('      * see run_case.log in dest/' + \
                                       case.title)

                    self.__log_file.flush()

        os.chdir(home)

        self.reporting('')

    #---------------------------------------------------------------------------

    def check_compare(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to make a comparison with a file which does not exist.
        """
        for case in self.graph.graph_dict:
            check_msg = "  o Check compare of case: " + case.title
            self.report_action_location(check_msg, destination)

            # reference directory passed in studymanager command line overwrites
            # destination in all cases (even if compare is defined by a compare
            # markup with a non empty destination)

            ref = None
            if self.__ref:
                ref = os.path.join(self.__ref, case.study)
            cases_to_disable = []

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
                msg = case.disable()
                self.reporting(msg)

        self.reporting('')

    #---------------------------------------------------------------------------

    def compare_case_and_report(self, case, repo, dest, threshold, args,
                                reference=None):
        """
        Compare the results for one computation and report
        """
        case.is_compare = "done"

        if not case.disabled:
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

        if case.disabled:
            self.reporting('    - compare %s (%s) --> DISABLED'
                           %(case.title, s_args))
        elif not m_size_eq:
            self.reporting('    - compare %s (%s) --> DIFFERENT MESH SIZES FOUND'
                           %(case.title, s_args))
        elif diff_value:
            self.reporting('    - compare %s (%s) --> DIFFERENCES FOUND'
                           %(case.title, s_args))
        else:
            self.reporting('    - compare %s (%s) --> NO DIFFERENCES FOUND'
                           %(case.title, s_args))


    #---------------------------------------------------------------------------

    def compare(self):
        """
        Compare the results of the new computations with those from the
        Repository.
        """
        if self.__compare:
            for case in self.graph.graph_dict:
                self.reporting('  o Compare case: ' + case.title)
                # reference directory passed in studymanager command line
                # overwrites destination in all cases (even if compare is
                # defined by a compare markup with a non empty destination)
                ref = None
                if self.__ref:
                    ref = os.path.join(self.__ref, case.study)
                if case.compare == 'on' and case.is_run != "KO":
                    is_compare, nodes, repo, dest, t, args = \
                                         self.__parser.getCompare(case.node)
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
                msg = case.disable()
                self.reporting(msg)

        if scripts_checked:
            self.reporting('')

        return scripts_checked

    #---------------------------------------------------------------------------

    def scripts(self):
        """
        Launch external additional scripts with arguments.
        """
        # create smgr_post_pro.log in dest/STUDY

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
                                r = os.path.join(self.__repo,  l, case.label,
                                                 "RESU", repo[i])
                                cmd += " -r " + r
                            if dest[i]:
                                d = os.path.join(self.__dest, l, case.label,
                                                 "RESU", dest[i])
                                cmd += " -d " + d

                            retcode, t = run_studymanager_command(cmd,
                                                                  self.__log_post_file)
                            stat = "FAILED" if retcode != 0 else "OK"

                            self.reporting('    - script %s --> %s (%s s)'
                                           %(stat, sc_name, t),
                                           stdout=True, report=False)

                            self.reporting('    - script %s --> %s (%s s)'
                                           %(stat, cmd, t),
                                           stdout=False, report=True)
                        else:
                            self.reporting('    - script %s not found' % cmd)

        self.reporting('')

    #---------------------------------------------------------------------------

    def postpro(self):
        """
        Launch external additional scripts with arguments.
        """
        for l, s in self.studies:
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

                        retcode, t = run_studymanager_command(cmd, self.__log_post_file)
                        stat = "FAILED" if retcode != 0 else "OK"

                        self.reporting('    - postpro %s --> %s (%s s)' \
                                       % (stat, sc_name, t),
                                       stdout=True, report=False)

                        self.reporting('    - postpro %s --> %s (%s s)' \
                                       % (stat, cmd, t),
                                       stdout=False, report=True)
                    else:
                        self.reporting('    - postpro %s not found' % cmd)

        # erase empty log file
        self.__log_post_file.close()
        log_post_file = open(self.__log_post_name, mode='r')
        content = log_post_file.read()
        log_post_file.close()
        if not content:
            os.system("rm -f " + self.__log_post_name)
        else:
            self.reporting(" /!\ ERROR during post. See %s\n" \
                           %self.__log_post_name, stdout=True, report=True)

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
                msg = case.disable()
                self.reporting(msg)

        self.reporting('')

    #---------------------------------------------------------------------------

    def plot(self):
        """
        Plot data.
        """
        if self.__plotter:
            for l, s in self.studies:
                if s.cases:
                    self.reporting('  o Plot study: ' + l)
                    self.__plotter.plot_study(l, s,
                                              self.__dis_tex,
                                              self.__default_fmt)

        self.reporting('')

    #---------------------------------------------------------------------------

    def report_input(self, doc, i_nodes, s_label, c_label=None):
        """
        Add input to report detailed.
        """
        for i_node in i_nodes:
            f, dest, repo, tex = self.__parser.getInput(i_node)
            doc.appendLine("\\subsubsection{%s}" % f)

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
                fd = os.path.join(dd, s_label, c_label, 'RESU')
            else:
                fd = os.path.join(dd, s_label, 'POST')

            if d == '':
                ff = os.path.join(fd, f)
                if not os.path.isfile(ff):
                    l = os.listdir(fd)
                    if len(l) == 1:
                        ff = os.path.join(fd, l[0], f)
            else:
                ff = os.path.join(fd, d, f)

            if not os.path.isfile(ff):
                self.reporting("\n Warning: this file does not exist: %s\n" %ff)
            elif ff[-4:] in ('.png', '.jpg', '.pdf') or ff[-5:] == '.jpeg':
                doc.addFigure(ff)
            elif tex == 'on':
                doc.addTexInput(ff)
            else:
                doc.addInput(ff)

    #---------------------------------------------------------------------------

    def build_reports(self, report_fig):
        """
        @type report_fig: C{String}
        @param report_fig: name of the figures report.
        @rtype: C{List} of C{String}
        @return: list of file to be attached to the report.
        """
        attached_files = []

        # figures report
        doc = Report(self.__dest,
                     report_fig,
                     self.__pdflatex)

        for l, s in self.studies:
            if not s.needs_report_detailed(self.__postpro):
                continue

            doc.appendLine("\\section{%s}" % l)

            if s.matplotlib_figures or s.input_figures:
                doc.appendLine("\\subsection{Graphical results}")
                for g in s.matplotlib_figures:
                    doc.addFigure(g)
                for g in s.input_figures:
                    doc.addFigure(g)

            for case in s.cases:

                # handle the input nodes that are inside case nodes
                if case.plot == "on" and case.is_run != "KO":
                    nodes = self.__parser.getChildren(case.node, "input")
                    if nodes:
                        doc.appendLine("\\subsection{Results for "
                                       "case %s}" % case.label)
                        self.report_input(doc, nodes, l, case.label)

            # handle the input nodes that are inside postpro nodes

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
                doc.appendLine("\\subsection{Results for "
                               "post-processing cases}")
                for i in range(len(label)):
                    if script[i]:
                        input_nodes = \
                            self.__parser.getChildren(nodes[i], "input")
                        if input_nodes:
                            self.report_input(doc, input_nodes, l)

        attached_files.append(doc.close())

        return attached_files

    #---------------------------------------------------------------------------

    def getlabel(self):
        return self.labels

#-------------------------------------------------------------------------------
# class dependency_graph
#-------------------------------------------------------------------------------

class dependency_graph(object):

    def __init__(self):
        """ Initializes a dependency graph object to an empty dictionary
        """
        self.graph_dict = {}

    def add_dependency(self, dependency):
        """ Defines dependency between two cases as an edge in the graph
        """
        (node1, node2) = dependency
        # TODO: Add error message as only one dependency is possible
        if node1 in self.graph_dict:
            self.graph_dict[node1].append(node2)
        else:
            self.graph_dict[node1] = [node2]

    def add_node(self, case):
        """ Add a case in the graph if not already there.
            Add a dependency when depends parameters is defined
        """
        if case not in self.graph_dict:
            self.graph_dict[case] = []

            if case.depends:
                for neighbor in self.graph_dict:
                    neighbor_name = neighbor.study + '/' + neighbor.label + '/' \
                                  + neighbor.run_id
                    if neighbor_name == case.depends:
                        # cases with dependency are level > 0 and connected to the dependency
                        self.add_dependency((case, neighbor))
                        case.level = neighbor.level + 1
                        break
                if case.level is None:
                    msg = "Problem in graph construction : dependency " \
                          + case.depends + " is not found.\n"
                    sys.exit(msg)
            else:
                # cases with no dependency are level 0
                case.level = 0

    def nodes(self):
        """ returns the cases of the dependency graph """
        return list(self.graph_dict.keys())

    def dependencies(self):
        """ returns the dependencies between the cases of the graph """
        dependencies = []
        for node in self.graph_dict:
            for neighbor in self.graph_dict[node]:
                if (neighbor, node) not in dependencies:
                    dependencies.append((node, neighbor))
        return dependencies

    def extract_sub_graph(self, filter_level, filter_n_procs):
        """ extracts a sub_graph based on level and n_procs criteria"""
        sub_graph = dependency_graph()
        for node in self.graph_dict:
            keep_node = True
            if filter_level is not None:
                keep_node = node.level is int(filter_level)
            if filter_n_procs is not None:
                keep_node = node.n_procs is int(filter_n_procs)
            if keep_node:
                sub_graph.add_node(node)
        return sub_graph

    def __str__(self):
        res = "\nList of cases: "
        for node in self.nodes():
            res += str(node) + " "
        res += "\nList of dependencies: "
        for dependency in self.dependencies():
            (node1, node2) = dependency
            res += '\n' + str(node1.name) + ' depends on ' + str(node2.name)
        return res

#-------------------------------------------------------------------------------
