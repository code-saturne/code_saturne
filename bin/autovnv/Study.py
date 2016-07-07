# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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

from cs_exec_environment import get_shell_type, enquote_arg
from autovnv.Parser import Parser
from autovnv.TexMaker import Report1, Report2
try:
    from autovnv.Drawing import Plotter
except Exception:
    pass
try:
    from autovnv.PlotVTK import PlotVTK
except Exception:
    pass
from autovnv.Command import run_autovnv_command
import cs_runcase

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
        self.compare    = data['compare']

        self.is_compil  = "not done"
        self.is_run     = "not done"
        self.is_time    = "not done"
        self.is_plot    = "not done"
        self.is_compare = "not done"
        self.threshold  = "default"
        self.diff_value = []
        self.subdomains = None
        self.run_dir    = ""

        self.resu = 'RESU'

        # Specific case for coupling

        coupling = os.path.join(self.__repo, self.label, "coupling_parameters.py")
        if os.path.isfile(coupling):
            import cs_case_coupling
            try:
                exec(compile(open(coupling).read(), '<string>', 'exec'))
            except Exception:
                execfile(coupling)
            run_ref = os.path.join(self.__repo, self.label, "runcase")
            self.exe, self.pkg = self.__get_exe(pkg, run_ref)
            self.subdomains = []
            for d in locals()['domains']:
                if d['solver'] == self.pkg.code_name:
                    self.subdomains.append(d['domain'])
            self.resu = 'RESU_COUPLING'

        else:
            run_ref = os.path.join(self.__repo, self.label, "SCRIPTS", "runcase")
            self.exe, self.pkg = self.__get_exe(pkg, run_ref)


    def __get_exe(self, old_pkg, run_ref):
        """
        Return the name of the exe of the case, in order to mix
        Code_Saturne and NEPTUNE_CFD test cases in the same study.
        """

        # Read the runcase script from the Repository

        runcase = cs_runcase.runcase(run_ref)

        if runcase.cmd_name == "code_saturne":
            from cs_package import package
            pkg = package(old_pkg.scriptdir)
        elif runcase.cmd_name == "neptune_cfd":
            from nc_package import package
            pkg = package(old_pkg.scriptdir)

        return runcase.cmd_name, pkg


    def __update_domain(self, subdir, xmlonly=False):
        """
        Update path for the script in the Repository.
        """
        # 1) Load the xml file of parameters in order to update it
        #    with the __backwardCompatibility method.

        from Base.XMLengine import Case
        from Pages.ScriptRunningModel import ScriptRunningModel
        from cs_exec_environment import \
            separate_args, get_command_single_value, update_command_single_value, assemble_args

        if self.exe == "code_saturne":
            from Base.XMLinitialize import XMLinit
        elif self.exe == "neptune_cfd":
            from core.XMLinitialize import XMLinit

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
                        print("This file is not in accordance with XML specifications.")
                        sys.exit(1)

                    case['xmlfile'] = fp
                    case.xmlCleanAllBlank(case.xmlRootNode())
                    XMLinit(case).initialize()
                    case.xmlSaveDocument()

        # 2) Create RESU directory if needed
        if not xmlonly:
            r = os.path.join(self.__repo, subdir, "RESU")
            if not os.path.isdir(r):
                os.makedirs(r)

        # 3) Update the GUI script from the Repository
        ref = os.path.join(self.__repo, subdir, "DATA", self.pkg.guiname)

        have_gui = 1
        try:
            f = open(ref, mode = 'r')
        except IOError:
           print("Warning SaturneGUI does not exist: %s\n" % ref)
           have_gui = 0

        if have_gui:
           lines = f.readlines()
           f.close()

           for i in range(len(lines)):
               if re.search(r'^prefix=', lines[i]):
                   if xmlonly:
                       lines[i] = "prefix=\n"
                   else:
                       lines[i] = "prefix=" + self.pkg.get_dir('prefix') + "\n"
               if re.search(r'^export PATH=', lines[i]):
                   if xmlonly:
                       lines[i] = 'export PATH="":$PATH\n'
                   else:
                       lines[i] = 'export PATH="' + self.pkg.get_dir('bindir') +'":$PATH\n'
           f = open(ref, mode = 'w')
           f.writelines(lines)
           f.close()

        # 4) Update the runcase script from the Repository
        self.update_runcase_path(os.path.join(self.__repo, subdir, "SCRIPTS"), None, xmlonly)


    def update_runcase_path(self, subdir, destdir=None, xmlonly=False):
        """
        Update path for the script in the Repository.
        """
        if not destdir:
            ref = os.path.join(self.__repo, subdir, "runcase")
        else:
            ref = os.path.join(destdir, subdir, "runcase")

        try:
            f = open(ref, mode = 'r')
        except IOError:
            print("Error: can not open %s\n" % ref)
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
                j = lines[i].find(self.pkg.name)

        f = open(ref, mode = 'w')
        f.writelines(lines)
        f.close()


    def update(self, xmlonly=False):
        """
        Update path for the script in the Repository.
        """
        # 1) Load the xml file of parameters in order to update it
        #    with the __backwardCompatibility method.

        from Pages.ScriptRunningModel import ScriptRunningModel
        from cs_exec_environment import \
            separate_args, get_command_single_value, update_command_single_value, assemble_args, enquote_arg

        if self.exe == "code_saturne":
            from Base.XMLinitialize import XMLinit
        elif self.exe == "neptune_cfd":
            from core.XMLinitialize import XMLinit

        if self.subdomains:
            cdirs = []
            for d in self.subdomains:
                cdirs.append(os.path.join(self.label, d))
        else:
            cdirs = (self.label,)

        for d in cdirs:
            self.__update_domain(d, xmlonly)

        # Update the runcase script from the Repository in case of coupling

        if self.subdomains:
            ref = os.path.join(self.__repo, self.label, "runcase")
            try:
                f = open(ref, mode = 'r')
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
                f = open(ref, mode = 'w')
                f.writelines(lines)
                f.close()
            except IOError:
                pass


    def compile(self, d):
        """
        Just compile user sources if exist.
        @rtype: C{String}
        @return: the status of the succes of the compilation.
        """

        if self.subdomains:
            sdirs = []
            for sd in self.subdomains:
                sdirs.append(os.path.join(d, self.label, sd, 'SRC'))
        else:
            sdirs = (os.path.join(d, self.label, 'SRC'),)

        home = os.getcwd()

        self.is_compil = "OK"

        for s in sdirs:

            os.chdir(s)
            cmd = os.path.join(self.pkg.get_dir('bindir'), self.exe) + " compile -t"
            p = subprocess.Popen(cmd,
                                 shell=True,
                                 executable=get_shell_type(),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
            output = p.communicate()
            o = output[1]
            if o.find('erreur') != -1 or o.find('error') != -1:
                self.is_compil = "KO"

        os.chdir(home)

        return self.is_compil


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


    def __updateRuncase(self, run_id):
        """
        Update the command line in the launcher C{runcase}.
        """
        from cs_exec_environment import separate_args, \
                                 get_command_single_value

        if self.subdomains:
            run_ref = os.path.join(self.__repo, self.label, "runcase")
            run_new = os.path.join(self.__dest, self.label, "runcase")
            run_ref_win = os.path.join(self.__repo, self.label, "runcase.bat")
        else:
            run_ref = os.path.join(self.__repo, self.label, "SCRIPTS", "runcase")
            run_new = os.path.join(self.__dest, self.label, "SCRIPTS", "runcase")
            run_ref_win = os.path.join(self.__repo, self.label, "SCRIPTS", "runcase.bat")

        if sys.platform.startswith('win'):
            if self.subdomains:
                run_new = os.path.join(self.__dest, self.label, "runcase.bat")
            else:
                run_new = os.path.join(self.__dest, self.label, "SCRIPTS", "runcase.bat")

        # Read the runcase from the Repository

        try:
            f = open(run_ref, mode = 'r')
        except:
            try:
                f = open(run_ref_win, mode = 'r')
            except IOError:
                print("Error: can not open %s\n" % run_ref)
                sys.exit(1)

        for line in f.readlines():
            if sys.platform.startswith('win'):
                if re.search(self.exe, line):
                    #we suppress "\" if needed
                    content = line.split()
                    if re.search(r'^\\' + self.exe, line):
                        content[0] = self.exe
                    run_cmd = " ".join(content)
            else:
                if re.search(r'^\\' + self.exe, line):
                    run_cmd = " ".join(line.split())
        f.close()

        # Write the new runcase

        try:
            f = open(run_new, mode = 'r')
        except IOError:
            print("Error: can not open %s\n" % run_new)
            sys.exit(1)

        lines = f.readlines()
        f.close()

        for i in range(len(lines)):
            l = lines[i].strip()
            if len(l) == 0:
                continue
            if l[0] == '#':
                continue
            if sys.platform.startswith('win'):
                if re.search(self.exe, lines[i]):
                    lines[i] = run_cmd + " --id=" + run_id
                elif re.search(r' cd ', lines[i]): # do not switch to batch submit directory
                    lines[i] = '# ' + lines[i]
            else:
                if re.search(r'^\\' + self.exe, lines[i]):
                    lines[i] = run_cmd + " --id=" + run_id
                elif re.search(r' cd ', lines[i]): # do not switch to batch submit directory
                    lines[i] = '# ' + lines[i]

        f = open(run_new, mode = 'w')
        f.writelines(lines)
        f.close()


    def run(self):
        """
        Run the case a thread.
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
                print("Warning: the directory %s already exists in the destination." % run_dir)

                if os.path.isfile(os.path.join(run_dir, "error")):
                    self.is_run = "KO"
                    self.is_time = 0.
                    error = 1
                else:
                    self.is_run = "OK"
                    self.is_time = 0.
                    error = 0
                os.chdir(home)

                return error

        else:
            run_id, run_dir = self.__suggest_run_id()

            while os.path.isdir(run_dir):
                time.sleep(5)
                run_id, run_dir = self.__suggest_run_id()

        self.run_id  = run_id
        self.run_dir = run_dir

        self.__updateRuncase(run_id)

        if sys.platform.startswith('win'):
            error, self.is_time = run_autovnv_command("runcase.bat", self.__log)
        else:
            error, self.is_time = run_autovnv_command("./runcase", self.__log)

        if not error:
            self.is_run = "OK"
        else:
            self.is_run = "KO"

        os.chdir(home)

        return error


    def runCompare(self, studies, r, d, threshold, args, reference=None):
        home = os.getcwd()

        node = None

        if reference:
            result = os.path.join(reference, self.label, self.resu)
        else:
            result = os.path.join(self.__repo, self.label, self.resu)
        repo = self.check_dir(studies, node, result, r, "repo")
        repo = os.path.join(result, repo, 'checkpoint', 'main')

        result = os.path.join(self.__dest, self.label, self.resu)
        dest = self.check_dir(studies, node, result, d, "dest")
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

        tab = []
        for i in range(len(lines)):
            if lines[i].find("Support") != -1:
                line = [x.replace("\""," ").strip() for x in lines[i].split(";")]
                name = line[0]
                info = [x.split(":") for x in line[1:]]
                info = [[x[0].strip(),x[1].strip()] for x in info]
                if info[1][1] not in ['i4', 'u4', 'c']:
                    line = [x.strip() for x in lines[i+1].split(";")]
                    vals = [x.split(":") for x in line]
                    vals = [[x[0].strip(),x[1].strip()] for x in vals]
                    tab.append([name.replace("_", "\_"),
                                vals[1][1],
                                vals[2][1],
                                self.threshold])

        os.chdir(home)

        return tab


    def check_dir(self, studies, node, result, rep, attr):
        """
        Check coherency between xml file of parameters and repository or destination.
        """
        # 1. Check if the given result directory exists.
        if rep != "":
            rep_f = os.path.join(result, rep)
            rep_e = os.path.join(result, rep, 'error')

            if not os.path.isdir(rep_f):
                msg = "Study %s case %s:\nthe directory %s does not exist.\nStop.\n" % \
                    (self.__study, self.label, rep_f)
                studies.reporting(msg)
                sys.exit(1)

            if os.path.isfile(rep_e):
                msg = "Study %s case %s:\nthe directory %s contains an error file.\nStop.\n" % \
                    (self.__study, self.label, rep_f)
                studies.reporting(msg)
                sys.exit(1)

            return rep

        # 2. The result directory must be read automatically;
        #    check if there is a single result directory.
        elif rep == "":
            if len(list(filter(nodot, os.listdir(result)))) == 0:
                msg = "Study %s case %s:\nthere is no result directory in %s\nStop.\n" % \
                    (self.__study, self.label, result)
                studies.reporting(msg)
                sys.exit(1)

            # if no run_id is specified in the xml file
            # only one result directory allowed in RESU
            if len(list(filter(nodot, os.listdir(result)))) > 1 and self.run_id == "":
                msg = "Study %s case %s:\nthere are several directories in %s\nStop.\n" % \
                    (self.__study, self.label, result)
                studies.reporting(msg)
                sys.exit(1)

            rep = self.run_id
            # if no run_id is specified in the xml file
            # the only result directory present in RESU is taken
            if rep == "":
                rep = list(filter(nodot, os.listdir(result)))[0]

            rep_f = os.path.join(result, rep)
            if not os.path.isdir(rep_f):
                msg = "Study %s case %s:\nthe directory %s does not exist.\nStop.\n" % \
                    (self.__study, self.label, rep_f)
                studies.reporting(msg)
                sys.exit(1)

            rep_e = os.path.join(result, rep, 'error')
            if os.path.isfile(rep_e):
                msg = "Study %s case %s:\nthe directory %s contains an error file.\nStop.\n" % \
                    (self.__study, self.label, rep_f)
                studies.reporting(msg)
                sys.exit(1)

            # 3. Update the file of parameters with the name of the result directory
            if node:
                self.__parser.setAttribute(node, attr, rep)
            return rep

        else:
            studies.reporting('Error: check compare/script/plot/probe/resu/input failed.')
            sys.exit(1)


    def check_dirs(self, studies, node, repo, dest, reference=None):
        """
        Check coherency between xml file of parameters and repository and destination.
        """
        if repo != None:
            # build path to RESU directory with path to study and case label in the repo
            if reference:
                result = os.path.join(reference, self.label, self.resu)
            else:
                result = os.path.join(self.__repo, self.label, self.resu)
            self.check_dir(studies, node, result, repo, "repo")

        if dest != None:
            # build path to RESU directory with path to study and case label in the dest
            result = os.path.join(self.__dest, self.label, self.resu)
            self.check_dir(studies, node, result, dest, "dest")

#-------------------------------------------------------------------------------

class Study(object):
    """
    Create, run and compare all cases for a given study.
    """
    def __init__(self, pkg, parser, study, exe, dif, rlog):
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
        """
        # Initialize attributes
        self.__parser   = parser
        self.__main_exe = exe
        self.__diff     = dif
        self.__log      = rlog

        self.__repo = os.path.join(self.__parser.getRepository(),  study)
        self.__dest = os.path.join(self.__parser.getDestination(), study)

        if not os.path.isdir(self.__repo):
            print("Error: the directory %s does not exist" % self.__repo)
            sys.exit(1)

        self.label = study

        self.Cases = []
        self.matplotlib_figures = []
        self.vtk_figures = []

        # build the list of the cases
        self.cases = parser.getCasesLabel(study)
        if not self.cases:
            print("\n\n\nWarning: no case defined in %s study\n\n\n" % study)
        else:
            for data in self.__parser.getCasesKeywords(self.label):
                c = Case(pkg,
                         self.__log,
                         self.__diff,
                         self.__parser,
                         self.label,
                         data,
                         self.__repo,
                         self.__dest)
                self.Cases.append(c)


    def getCasesLabel(self):
        """
        Return the list of the labels of the cases of the current study.
        @rtype: C{List}
        @return: labels of the cases
        """
        return self.cases


    def createCases(self):
        """
        Create a single study with its all cases.
        """
        repbase = os.getcwd()

        # Create study if necessary
        if not os.path.isdir(self.__dest):
            cmd = self.__main_exe + " create --quiet --study " + self.__dest
            retval, t = run_autovnv_command(cmd, self.__log)
            shutil.rmtree(os.path.join(self.__dest, "CASE1"))

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

        # Create cases
        os.chdir(self.__dest)

        log_lines = []
        for c in self.Cases:
            if not os.path.isdir(c.label):
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
                            node_retval, t = run_autovnv_command(cmd, self.__log)
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
                    retval, t = run_autovnv_command(cmd, self.__log)
                if retval == 0:
                    log_lines += ['    - create case: ' + c.label]
                else:
                    log_lines += ['    - create case: %s --> FAILED' % c.label]
            else:
                print("Warning: the case %s already exists in the destination." % c.label)

        os.chdir(repbase)

        return log_lines


    def getRunDirectories(self):
        list_cases = []
        list_dir   = []

        for case in self.Cases:
            if case.is_run != "KO":
                list_cases.append(case.label)
                list_dir.append(case.run_dir)

        return " ".join(list_cases), " ".join(list_dir)

#-------------------------------------------------------------------------------

class Studies(object):
    """
    Manage all Studies and all Cases described in the files of parameters.
    """
    def __init__(self, pkg, f, v, x, r, n, c, d, p, exe, dif, log, dis_tex):
        """
        Constructor.
          1. create if necessary the destination directory,
          2. initialize the parser and the plotter,
          3. build the list of the studies,
          4. start the report.
        @type f: C{String}
        @param f: xml file of parameters.
        @type v: C{True} or C{False}
        @param v: verbose mode.
        @type x: C{True} or C{False}
        @param x: update xml in repo only.
        @type r: C{True} or C{False}
        @param r: have to run the case.
        @type c: C{True} or C{False}
        @param c: have to do a comparison.
        @type exe: C{String}
        @param exe: name of the solver executable: C{code_saturne} or C{neptune_cfd}.
        @type dif: C{String}
        @param dif: name of the diff executable: C{cs_io_dump -d}.
        """

        # create a first xml parser only for
        #   the repository verification and
        #   the destination creation

        self.__parser = Parser(f)

        # Test if the repository exists

        self.repo = self.getRepository()

        # create if necessary the destination directory

        self.__xmlupdate = x
        if not self.__xmlupdate:
            self.dest = self.getDestination()
        else:
            self.dest = self.getRepository()
        if not os.path.isdir(self.dest):
            os.makedirs(self.dest)

        # copy the xml file of parameters for update and restart

        file = os.path.join(self.dest, os.path.basename(f))
        try:
            shutil.copyfile(f, file)
        except:
            pass

        # create a new parser, which is definitive and the plotter

        self.__parser  = Parser(file)
        try:
            self.__plotter = Plotter(self.__parser)
        except Exception:
            self.__plotter = None
        try:
            self.__plotvtk = PlotVTK(self.__parser)
        except Exception:
            self.__plotvtk = None

        # build the list of the studies

        doc = os.path.join(self.dest, log)
        self.__log = open(doc, "w")
        self.labels  = self.__parser.getStudiesLabel()
        self.studies = []
        for l in self.labels:
            self.studies.append( [l, Study(pkg, self.__parser, l, \
                                  exe, dif, self.__log)] )

        # start the report
        self.report = os.path.join(self.dest, "report.txt")
        self.reportFile = open(self.report, mode='w')
        self.reportFile.write('\n')

        # attributes

        self.__verbose = v
        self.__running = r
        self.__n_iter  = n
        self.__compare = c
        self.__ref     = d
        self.__postpro = p
        self.__dis_tex = dis_tex

        # in case of restart
        iok = 0
        for l, s in self.studies:
            for case in s.Cases:
                if case.compute == 'on':
                   iok+=1
        if not iok:
            self.__running = False

        if self.__xmlupdate:
            os.remove(file)
            os.remove(doc)

    def reporting(self, msg, screen_only=False):
        """
        Write the report.
        http://www.developpez.net/forums/d651999/autres-langages/python-zope/general-python/console-ecrire-ligne-precedente/
        Carriage return "chr(13)" renvoi en debut de ligne et line feed "chr(10)" effectue un saut de ligne.
        @type l: C{String}
        @param l: the sentence to be writing.
        """
        if not screen_only:
            s = chr(10)
        else:
            s  = ""

        if not self.__verbose:
            sys.stdout.write (msg + chr(13) + s)
            sys.stdout.flush()

        if not screen_only:
            self.reportFile.write(msg + '\n')
            self.reportFile.flush()


    def getRepository(self):
        """
        Return the directory of the repository.
        @rtype: C{String}
        @return: the directory of the repository of all studies.
        """
        r = self.__parser.getRepository()

        if not os.path.isdir(r):
            self.reporting('Error: the repository does not exist.')
            sys.exit(1)

        return r


    def getDestination(self):
        """
        Return the directory of the destination.
        @rtype: C{String}
        @return: the directory of the destination.
        """
        return self.__parser.getDestination()


    def updateRepository(self, xml_only=False):
        """
        Create all studies and all cases.
        """
        iok = 0
        for l, s in self.studies:
            self.reporting('  o Update repository: ' + l)
            for case in s.Cases:
                self.reporting('    - update  %s' % case.label)
                case.update(xml_only)
                if case.compile(os.path.join(self.repo, l)) == "OK":
                    self.reporting('    - compile %s --> OK' % case.label)
                else:
                    self.reporting('    - compile %s --> FAILED' % case.label)
                    iok+=1
        if iok:
            self.reporting('Error: compilation failed. Number of failed cases: %s' % iok)
            sys.exit(1)


    def createStudies(self):
        """
        Create all studies and all cases.
        """
        for l, s in self.studies:
            self.reporting('  o Create study: ' + l)
            log_lines = s.createCases()
            for line in log_lines:
                self.reporting(line)


    def compilation(self):
        """
        Compile sources of all cases.
        """
        iok = 0
        for l, s in self.studies:
            self.reporting('  o Compile study: ' + l)
            for case in s.Cases:
                if case.compute == 'on':
                    if case.compile(os.path.join(self.dest, l)) == "OK":
                        self.reporting('    - compile %s --> OK' % case.label)
                    else:
                        self.reporting('    - compile %s --> FAILED' % case.label)
                        iok+=1
        if iok:
            self.reporting('Error: compilation failed. Number of failed cases: %s' % iok)
            sys.exit(1)


    def prepro(self, l, s, case):
        """
        Launch external additional scripts with arguments.
        """
        pre, label, nodes, args = self.__parser.getPrepro(case.node)
        for i in range(len(label)):
            if pre[i]:
                # search if the script is the directoty MESH
                # if not, the script is searched in the directories of the current case
                cmd = os.path.join(self.dest, l, "MESH", label[i])
                if not os.path.isfile(cmd):
                    filePath = ""
                    for root, dirs, fs in os.walk(os.path.join(self.dest, l, case.label)):
                        if label[i] in fs:
                            filePath = root
                            break

                    cmd = os.path.join(filePath, label[i])

                if os.path.isfile(cmd):
                    cmd += " " + args[i]
                    cmd += " -c " + os.path.join(self.dest, l, case.label)
                    repbase = os.getcwd()
                    os.chdir(os.path.join(self.dest, l, "MESH"))
                    # Prepro external script might need pythondir and pkgpythondir
                    pdir = case.pkg.get_dir('pythondir')
                    pdir = pdir + ":" + case.pkg.get_dir('pkgpythondir')
                    # Pythondirs might be different for NEPTUNE_CFD
                    # and have to be added as well
                    if case.pkg.code_name == "NEPTUNE_CFD":
                        sdir = case.pkg.get_cs_dir('pythondir')
                        sdir = sdir + ":" + case.pkg.get_cs_dir('pkgpythondir')
                        pdir = pdir + ":" + sdir
                    retcode, t = run_autovnv_command(cmd, self.__log, pythondir = pdir)
                    os.chdir(repbase)
                    self.reporting('    - script %s --> OK (%s s)' % (cmd, t))
                else:
                    self.reporting('    - script %s not found' % cmd)


    def run(self):
        """
        Update and run all cases.
        Warning, if the makup of the case is repeated in the xml file of parameters,
        the run of the case is also repeated.
        """
        for l, s in self.studies:
            self.reporting('  o Script prepro and run for study: ' + l)
            for case in s.Cases:
                self.prepro(l, s, case)
                if self.__running:
                    if case.compute == 'on' and case.is_compil != "KO":

                        if self.__n_iter:
                            if case.subdomains:
                                case_dir = os.path.join(self.dest, s.label, case.label,
                                                        case.subdomains[0], "DATA")
                            else:
                                case_dir = os.path.join(self.dest, s.label, case.label, "DATA")
                            os.chdir(case_dir)
                            # Create a control_file in each case DATA
                            if not os.path.exists('control_file'):
                                control_file = open('control_file','w')
                                control_file.write("time_step_limit " + str(self.__n_iter) + "\n")
                                # Flush to ensure that control_file content is seen
                                # when control_file is copied to the run directory on all systems
                                control_file.flush()
                                control_file.close

                        self.reporting('    - running %s ...' % case.label, True)
                        error = case.run()
                        if not error:
                            if not case.run_id:
                                self.reporting('    - run %s --> Warning suffix is not read' % case.label)

                            self.reporting('    - run %s --> OK (%s s) in %s' % (case.label, case.is_time, case.run_id))
                            self.__parser.setAttribute(case.node, "compute", "off")

                            # update dest="" attribute
                            n1 = self.__parser.getChilds(case.node, "compare")
                            n2 = self.__parser.getChilds(case.node, "script")
                            n3 = self.__parser.getChilds(case.node, "data")
                            n4 = self.__parser.getChilds(case.node, "probe")
                            n5 = self.__parser.getChilds(case.node, "resu")
                            n6 = self.__parser.getChilds(case.node, "input")
                            for n in n1 + n2 + n3 + n4 + n5 + n6:
                                if self.__parser.getAttribute(n, "dest") == "":
                                    self.__parser.setAttribute(n, "dest", case.run_id)
                        else:
                            self.reporting('    - run %s --> FAILED' % case.label)

                        self.__log.flush()


    def check_compare(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to make a comparison with a file which does not exist.
        """
        for l, s in self.studies:
            # reference directory passed in autovnv command line overwrites
            # destination in all cases (even if compare is defined by a compare
            # markup with a non empty destination)
            ref = None
            if self.__ref:
                ref = os.path.join(self.__ref, s.label)
            for case in s.Cases:
                if case.compare == 'on' and case.is_run != "KO":
                    compare, nodes, repo, dest, threshold, args = self.__parser.getCompare(case.node)
                    if compare:
                        is_checked = False
                        for i in range(len(nodes)):
                            if compare[i]:
                                is_checked = True
                                if destination == False:
                                    dest[i]= None
                                case.check_dirs(self, nodes[i], repo[i], dest[i], reference=ref)
                    if not compare or not is_checked:
                        node = None
                        repo = ""
                        dest = ""
                        if destination == False:
                            dest = None
                        case.check_dirs(self, node, repo, dest, reference=ref)


    def compare(self):
        """
        Compare the results of the new computations with those from the Repository.
        """
        if self.__compare:
            for l, s in self.studies:
                self.reporting('  o Compare study: ' + l)
                # reference directory passed in autovnv command line overwrites
                # destination in all cases (even if compare is defined by a
                # compare markup with a non empty destination)
                ref = None
                if self.__ref:
                    ref = os.path.join(self.__ref, s.label)
                for case in s.Cases:
                    if case.compare == 'on' and case.is_run != "KO":
                        is_compare, nodes, repo, dest, t, args = self.__parser.getCompare(case.node)
                        if is_compare:
                            for i in range(len(nodes)):
                                if is_compare[i]:
                                    case.is_compare = "done"
                                    diff_value = case.runCompare(self, repo[i], dest[i], t[i], args[i], reference=ref)
                                    case.diff_value += diff_value
                                    if diff_value:
                                        self.reporting('    - compare %s (with args: %s) --> DIFFERENCES FOUND' % (case.label, args[i]))
                                    else:
                                        self.reporting('    - compare %s (with args: %s) --> NO DIFFERENCES FOUND' % (case.label, args[i]))
                        if not is_compare or case.is_compare != "done":
                            repo = ""
                            dest = ""
                            t    = None
                            args = None
                            case.is_compare = "done"
                            case.diff_value += case.runCompare(self, repo, dest, t, args, reference=ref)
                            if case.diff_value:
                                self.reporting('    - compare %s (default mode) --> DIFFERENCES FOUND' % (case.label))
                            else:
                                self.reporting('    - compare %s (default mode) --> NO DIFFERENCES FOUND' % (case.label))



    def check_script(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to run a script with a file which does not exist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                script, label, nodes, args, repo, dest = self.__parser.getScript(case.node)
                for i in range(len(nodes)):
                    if script[i] and case.is_run != "KO":
                        if destination == False:
                            dest[i] = None
                        case.check_dirs(self, nodes[i], repo[i], dest[i])


    def scripts(self):
        """
        Launch external additional scripts with arguments.
        """
        for l, s in self.studies:
            self.reporting('  o Script postpro study: ' + l)
            for case in s.Cases:
                script, label, nodes, args, repo, dest = self.__parser.getScript(case.node)
                for i in range(len(label)):
                    if script[i] and case.is_run != "KO":
                        cmd = os.path.join(self.dest, l, "POST", label[i])
                        if os.path.isfile(cmd):
                            cmd += " " + args[i]
                            if repo[i]:
                                r = os.path.join(self.repo,  l, case.label, "RESU", repo[i])
                                cmd += " -r " + r
                            if dest[i]:
                                d = os.path.join(self.dest, l, case.label, "RESU", dest[i])
                                cmd += " -d " + d
                            retcode, t = run_autovnv_command(cmd, self.__log)
                            self.reporting('    - script %s --> OK (%s s)' % (cmd, t))
                        else:
                            self.reporting('    - script %s not found' % cmd)


    def postpro(self):
        """
        Launch external additional scripts with arguments.
        """
        for l, s in self.studies:
            # fill results directories and ids for the cases of the current study
            # that were not run by the current autovnv command
            for case in s.Cases:
                if case.is_run != "KO":
                    if case.run_dir == "":
                        resu = os.path.join(self.dest, l, case.label, case.resu)
                        rep = case.check_dir(self, None, resu, "", "dest")
                        case.run_id = rep
                        case.run_dir = os.path.join(resu, rep)

            self.reporting('  o Postprocessing cases of study: ' + l)
            script, label, nodes, args = self.__parser.getPostPro(l)
            for i in range(len(label)):
                if script[i]:
                    cmd = os.path.join(self.dest, l, "POST", label[i])
                    if os.path.isfile(cmd):
                        list_cases, list_dir = s.getRunDirectories()
                        cmd += ' ' + args[i] + ' -c "' + list_cases + '" -d "' + list_dir + '" -s ' + l
                        retcode, t = run_autovnv_command(cmd, self.__log)
                        self.reporting('    - postpro %s --> OK (%s s)' % (cmd, t))
                    else:
                        self.reporting('    - postpro %s not found' % cmd)


    def check_plot(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to make a plot of a file which does not exist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                if case.plot == "on" and case.is_run != "KO":
                    for node in self.__parser.getChilds(case.node, "data"):
                        plots, file, dest, repo = self.__parser.getResult(node)
                        if destination == False:
                            dest = None
                        case.check_dirs(self, node, repo, dest)

                    for node in self.__parser.getChilds(case.node, "probes"):
                        file, dest, fig = self.__parser.getProbes(node)
                        if destination == False:
                            dest = None
                        repo = None
                        case.check_dirs(self, node, repo, dest)

                    for node in self.__parser.getChilds(case.node, "resu"):
                        plots, file, dest, repo = self.__parser.getResult(node)
                        if destination == False:
                            dest = None
                        case.check_dirs(self, node, repo, dest)

                    for node in self.__parser.getChilds(case.node, "input"):
                        file, dest, repo, tex = self.__parser.getInput(node)
                        if destination == False:
                            dest = None
                        case.check_dirs(self, node, repo, dest)


    def plot(self):
        """
        Plot data.
        """
        if self.__plotter:
            for l, s in self.studies:
                self.reporting('  o Plot study: ' + l)
                self.__plotter.plot_study(l, s, self.__dis_tex)

        if self.__plotvtk:
            for l, s in self.studies:
                self.__plotvtk.plot_study(l, s)


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

        # Fisrt global report
        doc1 = Report1(self.dest,
                       report1,
                       self.__log,
                       self.report,
                       self.__parser.write())

        for l, s in self.studies:
            for case in s.Cases:
                if case.diff_value:
                    is_nodiff = "KO"
                else:
                    is_nodiff = "OK"

                doc1.add_row(l,
                             case.label,
                             case.is_compil,
                             case.is_run,
                             case.is_time,
                             case.is_compare,
                             is_nodiff)

        attached_files.append(doc1.close())

        # Second detailed report
        if self.__compare or self.__postpro:
            doc2 = Report2(self.dest, report2, self.__log)

            for l, s in self.studies:
                doc2.appendLine("\\section{%s}" % l)

                if s.matplotlib_figures or s.vtk_figures:
                    doc2.appendLine("\\subsection{Graphical results}")
                    for g in s.matplotlib_figures:
                        doc2.addFigure(g)
                    for g in s.vtk_figures:
                        doc2.addFigure(g)

                for case in s.Cases:
                    if case.is_compare == "done":
                        run_id = None
                        if case.run_id != "":
                            run_id = case.run_id
                        doc2.appendLine("\\subsection{Comparison for case %s (run_id: %s)}" % (case.label, run_id))
                        if case.diff_value:
                            doc2.add_row(case.diff_value, l, case.label)
                        elif self.__compare:
                            doc2.appendLine("No difference between the repository and the destination.")

                    # handle the input nodes that are inside case nodes
                    if case.plot == "on" and case.is_run != "KO":
                        nodes = self.__parser.getChilds(case.node, "input")
                        if nodes:
                            doc2.appendLine("\\subsection{Results for case %s}" % case.label)
                            for node in nodes:
                                f, dest, repo, tex = self.__parser.getInput(node)
                                doc2.appendLine("\\subsubsection{%s}" % f)

                                if dest:
                                    d = dest
                                    dd = self.dest
                                elif repo:
                                    d = repo
                                    dd = self.repo
                                else:
                                    d = ""
                                    dd = ""

                                ff = os.path.join(dd, l, case.label, "RESU", d, f)

                                if not os.path.isfile(ff):
                                    print("\n\nWarning: this file does not exist: %s\n\n" % ff)
                                elif tex == 'on':
                                    doc2.addTexInput(ff)
                                else:
                                    doc2.addInput(ff)

                # handle the input nodes that are inside postpro nodes
                if self.__postpro:
                    script, label, nodes, args = self.__parser.getPostPro(l)
                    doc2.appendLine("\\subsection{Results for post-processing cases}")
                    for i in range(len(label)):
                        if script[i]:
                            input_nodes = self.__parser.getChilds(nodes[i], "input")
                            if input_nodes:
                                for node in input_nodes:
                                    f, dest, repo, tex = self.__parser.getInput(node)
                                    doc2.appendLine("\\subsubsection{%s}" % f)

                                    if dest:
                                        d = dest
                                        dd = self.dest
                                    elif repo:
                                        d = repo
                                        dd = self.repo
                                    else:
                                        d = ""
                                        dd = ""

                                    ff = os.path.join(dd, l, "POST", d, f)

                                    if not os.path.isfile(ff):
                                        print("\n\nWarning: this file does not exist: %s\n\n" % ff)
                                    elif tex == 'on':
                                        doc2.addTexInput(ff)
                                    else:
                                        doc2.addInput(ff)

            attached_files.append(doc2.close())

        return attached_files


    def getlabel(self):
        return self.labels


    def logs(self):
        try:
            self.reportFile.close()
        except:
            pass
        self.reportFile = open(self.report, mode='r')
        return self.reportFile.read()


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
