# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

from autovnv.Parser import Parser
from autovnv.TexMaker import Report1, Report2
from autovnv.Drawing import Plotter
from autovnv.PlotVTK import PlotVTK
from autovnv.Command import run_command

#-------------------------------------------------------------------------------
# log config.
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
#log.setLevel(logging.DEBUG)
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------

class RunThread(threading.Thread):
    """
    Run an external script as a simple thread.
    """

    def __init__(self, cmd, rlog):

        self.__cmd = cmd
        self.__log = rlog
        self.__retcode = -999
        self.__time = 0.0
        threading.Thread.__init__ (self, target=self.run)
        self.__stopevent = threading.Event()

    def run(self):
        try:
            if not self.__stopevent.isSet():
                self.__retcode, self.__time = self.__runcommand()
        except:
            if not self.__stopevent.is_set():
                self.__retcode, self.__time = self.__runcommand()

    def stop(self):
        self.__stopevent.set()
        return self.__retcode, self.__time

    def __runcommand(self):
        retcode, t = run_command(self.__cmd, self.__log)
        self.__stopevent.wait(1.0)
        return retcode, t

#-------------------------------------------------------------------------------

class Case(object):
    def __init__(self, pkg, rlog, diff, study, data, repo, dest):
        """
        @type data: C{Dictionary}
        @param data: contains all keyword and value read in the parameters file
        """
        self.__log      = rlog
        self.__diff     = diff
        self.__study    = study
        self.__data     = data
        self.__repo     = repo
        self.__dest     = dest

        self.node       = data['node']
        self.label      = data['label']
        self.compute    = data['compute']
        self.plot       = data['post']

        self.is_compil  = "not done"
        self.is_run     = "not done"
        self.is_time    = "not done"
        self.is_plot    = "not done"
        self.is_compare = "not done"
        self.threshold  = "default"
        self.diff_value = []

        self.exe, self.pkg = self.__get_exe()


    def __get_exe(self):
        """
        Return the name of the exe of the case, in order to mix
        Code_Saturne and NEPTUNE_CFD test cases in the same study.
        """
        run_ref = os.path.join(self.__repo, self.label, "SCRIPTS", "runcase")

        # Read the runcase script from the Repository

        try:
            f = file(run_ref, mode = 'r')
        except IOError:
            print "Error: can not opening %s\n" % run_ref
            sys.exit(1)

        lines = f.readlines()
        f.close()

        exe    = ""
        for name in ('code_saturne', 'neptune_cfd'):
            for line in lines:
                if re.search(r'^\\' + name, line):
                    exe = name

        if not exe:
            print "Error: name of the executable for the case %s not found. Stop." % self.label
            sys.exit(1)

        if exe == "code_saturne":
            from Base.XMLinitialize import XMLinit
            from cs_package import package
            pkg = package()
        elif exe == "neptune_cfd":
            from core.XMLinitialize import XMLinit
            from nc_package import package
            pkg = package()

        return exe, pkg


    def update(self):
        """
        Update path for the script in the Repository.
        """
        # 1) Load the xml file of parameters in order to update it
        #    with the __backwardCompatibility method.

        from Base.XMLengine import Case

        if self.exe == "code_saturne":
            from Base.XMLinitialize import XMLinit
        elif self.exe == "neptune_cfd":
            from core.XMLinitialize import XMLinit

        try:
            fn = os.path.join(self.__repo, self.label, "DATA", "liu.xml")
            case = Case(package = self.pkg, file_name = fn)
        except:
            print "File of parameters reading error.\n"
            print "This file is not in accordance with XML specifications."
            sys.exit(1)

        case['xmlfile'] = fn
        case.xmlCleanAllBlank(case.xmlRootNode())
        XMLinit(case).initialize()
        case.xmlSaveDocument()

        # 2) Create RESU directory if needed
        r = os.path.join(self.__repo, self.label, "RESU")
        if not os.path.isdir(r):
            os.makedirs(r)

        # 3) Update the GUI script from the Repository
        ref = os.path.join(self.__repo, self.label, "DATA", self.pkg.guiname)

        try:
            f = file(ref, mode = 'r')
        except IOError:
            print "Error: can not opening %s\n" % ref
            sys.exit(1)

        lines = f.readlines()
        f.close()

        for i in range(len(lines)):
            if re.search(r'^prefix=', lines[i]):
                 lines[i] = "prefix=" + self.pkg.prefix + "\n"

        f = file(ref, mode = 'w')
        f.writelines(lines)
        f.close()

        # 4) Update the runcase script from the Repository
        ref = os.path.join(self.__repo, self.label, "SCRIPTS", "runcase")

        try:
            f = file(ref, mode = 'r')
        except IOError:
            print "Error: can not opening %s\n" % ref
            sys.exit(1)

        lines = f.readlines()
        f.close()

        for i in range(len(lines)):
            if re.search(r'^export PATH=', lines[i]):
                 lines[i] = 'export PATH="' + self.pkg.bindir +'":$PATH\n'

        f = file(ref, mode = 'w')
        f.writelines(lines)
        f.close()


    def compile(self, d):
        """
        Just compile user sources if exist.
        @rtype: C{String}
        @return: the status of the succes of the compilation.
        """
        home = os.getcwd()
        os.chdir(os.path.join(d, self.label, 'SRC'))
        cmd = os.path.join(self.pkg.bindir, self.exe) + " compile -t"

        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output = p.communicate()
        o = output[1]
        if o.find('erreur') != -1 or o.find('error') != -1:
            self.is_compil = "KO"
        else:
            self.is_compil = "OK"

        os.chdir(home)

        return self.is_compil


    def __suggest_run_id(self):

        cmd = os.path.join(self.pkg.bindir, self.exe) + " run --suggest-id"
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        i = p.communicate()[0]
        run_id = string.join(i.split())

        return run_id, os.path.join(self.__dest, self.label, "RESU", run_id)


    def __updateRuncase(self, run_id):
        """
        Update the command line in the launcher C{runcase}.
        """
        run_ref = os.path.join(self.__repo, self.label, "SCRIPTS", "runcase")
        run_new = os.path.join(self.__dest, self.label, "SCRIPTS", "runcase")

        # Read the runcase from the Repository

        try:
            f = file(run_ref, mode = 'r')
        except IOError:
            print "Error: can not opening %s\n" % run_ref
            sys.exit(1)

        for line in f.readlines():
            if re.search(r'^\\' + self.exe, line):
                run_cmd = string.join(line.split())
        f.close()

        # Write the new runcase

        try:
            f = file(run_new, mode = 'r')
        except IOError:
            print "Error: can not opening %s\n" % run_new
            sys.exit(1)

        lines = f.readlines()
        f.close()

        for i in range(len(lines)):
            if re.search(r'^\\' + self.exe, lines[i]):
                lines[i] = run_cmd + " --id=" + run_id

        f = file(run_new, mode = 'w')
        f.writelines(lines)
        f.close()


    def run(self):
        """
        Run the case a thread.
        """
        home = os.getcwd()
        os.chdir(os.path.join(self.__dest, self.label, 'SCRIPTS'))

        run_id, run_dir = self.__suggest_run_id()

        while os.path.isdir(run_dir):
            time.sleep(5)
            run_id, run_dir = self.__suggest_run_id()

        self.__updateRuncase(run_id)

        error, self.is_time = run_command("./runcase", self.__log)
        #t1 = RunThread("runcase", self.__log)
        #t1.start()
        #t1.join()
        #error, self.is_time = t1.stop()

        if not error:
            self.is_run = "OK"
        else:
            self.is_run = "KO"

        os.chdir(home)

        return error, run_id


    def compare(self, r, d, threshold, args):
        home = os.getcwd()

        repo = os.path.join(self.__repo, self.label,
                            'RESU', r, 'checkpoint', 'main')
        if not os.path.isfile(repo):
            msg = "In %s \nthe directory %s does not exist.\n" % (repo, r)
            raise ValueError, msg

        dest = os.path.join(self.__dest, self.label,
                            'RESU', d, 'checkpoint', 'main')
        if not os.path.isfile(dest):
            msg = "In %s \nthe directory %s does not exist.\n" % (dest, d)
            raise ValueError, msg

        cmd = self.__diff + ' ' + repo + ' ' + dest

        self.threshold = "default"
        if threshold != None:
            cmd += ' --threshold ' + threshold
            self.threshold = threshold

        if args != None:
            cmd += (" " + args)
            l = string.split(args)
            try:
                i = l.index('--threshold')
                self.threshold = l[i+1]
            except:
                pass

        l = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
        info = l.read().replace("\"", " ").replace(";", " ").replace(":", " ").split()

        tab = []

        for i in range(len(info)):
            if info[i][:4] == 'Diff':
                if info[i-3] not in ['i4', 'u4']:
                    tab.append( [info[i-7].replace("_", "\_"),
                                 info[i+3],
                                 info[i+5],
                                 self.threshold] )

        os.chdir(home)

        return tab

#-------------------------------------------------------------------------------

class Study(object):
    """
    Create, run and compare all cases fir a given study.
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
        self.__study    = study
        self.__main_exe = exe
        self.__diff     = dif
        self.__log      = rlog

        self.__repo = os.path.join(self.__parser.getRepository(),  study)
        self.__dest = os.path.join(self.__parser.getDestination(), study)

        if not os.path.isdir(self.__repo):
            print "Error: the directory %s does not exist" % self.__repo
            sys.exit(1)

        # build the list of the cases
        self.cases = parser.getCasesLabel(study)
        if not self.cases:
            print "Error: no case defined in %s study" % study
            sys.exit(1)

        self.Cases = []
        for data in self.__parser.getCasesKeywords(self.__study):
            c = Case(pkg,
                     self.__log,
                     self.__diff,
                     self.__study,
                     data,
                     self.__repo,
                     self.__dest)
            self.Cases.append(c)

            self.matplotlib_figures = []
            self.vtk_figures = []


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
            retval, t = run_command(cmd, self.__log)
            shutil.rmtree(os.path.join(self.__dest, "CASE1"))

            # Link meshes and copy other files
            ref = os.path.join(self.__repo, "MESH")
            if os.path.isdir(ref):
                l = os.listdir(ref)
                meshes = []
                for cpr in ["", ".gz"]:
                    for fmt in ["unv",
                                "med",
                                "ngeom",
                                "ccm",
                                "cgns",
                                "neu",
                                "msh",
                                "des"]:
                        meshes += fnmatch.filter(l, "*." + fmt + cpr)
                des = os.path.join(self.__dest, "MESH")
                for m in l:
                    if m in meshes:
                        os.symlink(os.path.join(ref, m), os.path.join(des, m))
                    elif m != ".svn":
                        shutil.copy2(os.path.join(ref, m), des)

            # Copy external scripts for post-processing
            ref = os.path.join(self.__repo, "POST")
            if os.path.isdir(ref):
                des = os.path.join(self.__dest, "POST")
                shutil.rmtree(des)
                shutil.copytree(ref, des, symlinks=True)

        # Create cases
        os.chdir(self.__dest)

        for c in self.Cases:
            if not os.path.isdir(c.label):
                e = os.path.join(c.pkg.bindir, c.exe)
                cmd = e + " create --case " + c.label  \
                      + " --quiet --noref --copy-from "    \
                      + os.path.join(self.__repo, c.label)
                retval, t = run_command(cmd, self.__log)
            else:
                print "Warning: the case %s already exists in the destination." % c.label

        os.chdir(repbase)

#-------------------------------------------------------------------------------

class Studies(object):
    """
    Manage all Studies and all Cases described in the files of parameters.
    """
    def __init__(self, pkg, f, v, r, c, p, exe, dif):
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

        self.dest = self.getDestination()
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
        self.__plotter = Plotter(self.__parser)
        self.__plotvtk = PlotVTK(self.__parser)

        # build the list of the studies

        doc = os.path.join(self.dest, "auto_vnv.log")
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
        self.__compare = c
        self.__postpro = p

        # in case of restart
        iok = 0
        for l, s in self.studies:
            for case in s.Cases:
                if case.compute == 'on':
                   iok+=1
        if not iok:
            self.__running = False

        iok = 0
        for l, s in self.studies:
            for case in s.Cases:
                if case.plot == 'on':
                   iok+=1
        if not iok:
            self.__postpro = False


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


    def updateRepository(self):
        """
        Create all studies and all cases.
        """
        iok = 0
        for l, s in self.studies:
            self.reporting('  o Update repository: ' + l)
            for case in s.Cases:
                self.reporting('    - update  %s' % case.label)
                case.update()
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
            s.createCases()
            for c in s.getCasesLabel():
                self.reporting('    - create case: ' + c)


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


    def prepro(self):
        """
        Launch external additional scripts with arguments.
        """
        for l, s in self.studies:
            self.reporting('  o Script prepro study: ' + l)
            for case in s.Cases:
                pre, label, nodes, args = self.__parser.getPrepro(case.node)
                for i in range(len(label)):
                    if pre[i]:
                        cmd = os.path.join(self.dest, l, "MESH", label[i])
                        if os.path.isfile(cmd):
                            cmd += " " + args[i]
                            repbase = os.getcwd()
                            os.chdir(os.path.join(self.dest, l, "MESH"))
                            retcode, t = run_command(cmd, self.__log)
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
        if self.__running:
            for l, s in self.studies:
                self.reporting('  o Run study: ' + l)
                for case in s.Cases:
                    if case.compute == 'on' and case.is_compil != "KO":
                        self.reporting('    - running %s ...' % case.label, True)
                        error, run_id = case.run()
                        if not error:
                            if not run_id:
                                self.reporting('    - run %s --> Warning suffixe is not read' % case.label)

                            self.reporting('    - run %s --> OK (%s s) in %s' % (case.label, case.is_time, run_id))
                            self.__parser.setAttribute(case.node, "compute", "off")

                            # update dest="" attribute
                            n1 = self.__parser.getChilds(case.node, "compare")
                            n2 = self.__parser.getChilds(case.node, "script")
                            n3 = self.__parser.getChilds(case.node, "data")
                            n4 = self.__parser.getChilds(case.node, "probe")
                            n5 = self.__parser.getChilds(case.node, "resu")
                            for n in n1 + n2 + n3 + n4 + n5:
                                if self.__parser.getAttribute(n, "dest") == "":
                                    self.__parser.setAttribute(n, "dest", run_id)
                        else:
                            self.reporting('    - run %s --> FAILED' % case.label)

                        self.__log.flush()


    def __check_dir(self, study_label, case_label, node, result, rep, attr):
        """
        Check coherency between xml file of parameters and repository or destination.
        """
        # 1. Check if the given result directory exists.
        if rep != "":
            rep_f = os.path.join(result, rep, 'checkpoint', 'main')
            if not os.path.isfile(rep_f):
                msg = "Study %s case %s:\nthe directory %s does not exist.\nStop.\n" % \
                    (study_label, case_label, rep_f)
                self.reporting(msg)
                sys.exit(1)

        # 2. The result directory must be read automatically;
        #    check if there is a single result directory.
        elif rep == "":
            if not (len(os.listdir(result)) == 1):
                msg = "Study %s case %s:\nthere is not a single result directory in %s\nStop.\n" % \
                    (study_label, case_label, result)
                self.reporting(msg)
                sys.exit(1)

        # 3. Update the file of parameters with the name of the result directory
            self.__parser.setAttribute(node, attr, os.listdir(result)[0])
        else:
            self.reporting('Error: check compare/script/plot/probe/resu failed.')
            sys.exit(1)


    def __check_dirs(self, study_label, case_label, node, repo, dest):
        """
        Check coherency between xml file of parameters and repository and destination.
        """
        if repo != None:
            result = os.path.join(self.repo, study_label, case_label, 'RESU')
            self.__check_dir(study_label, case_label, node, result, repo, "repo")

        if dest != None:
            result = os.path.join(self.dest, study_label, case_label, 'RESU')
            self.__check_dir(study_label, case_label, node, result, dest, "dest")


    def check_compare(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to make a comparison with a file which does not exsist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                compare, nodes, repo, dest, threshold, args = self.__parser.getCompare(case.node)
                for i in range(len(nodes)):
                    if compare[i] and case.is_run != "KO":
                        if destination == False:
                            dest[i]= None
                        self.__check_dirs(l, case.label, nodes[i], repo[i], dest[i])


    def compare(self):
        """
        Compare the results of the new computations with thoses from the Repository.
        """
        if self.__compare:
            for l, s in self.studies:
                self.reporting('  o Compare study: ' + l)
                for case in s.Cases:
                    is_compare, nodes, repo, dest, t, args = self.__parser.getCompare(case.node)
                    for i in range(len(nodes)):
                        if is_compare[i] and case.is_run != "KO":
                            case.is_compare = "done"
                            self.reporting('    - compare %s (%s)' % (case.label, args[i]))
                            case.diff_value += case.compare(repo[i], dest[i], t[i], args[i])


    def check_script(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to run a script with a file which does not exsist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                script, label, nodes, args, repo, dest = self.__parser.getScript(case.node)
                for i in range(len(nodes)):
                    if script[i] and case.is_run != "KO":
                        if destination == False:
                            dest[i] = None
                        self.__check_dirs(l, case.label, nodes[i], repo[i], dest[i])


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
                            retcode, t = run_command(cmd, self.__log)
                            self.reporting('    - script %s --> OK (%s s)' % (cmd, t))
                        else:
                            self.reporting('    - script %s not found' % cmd)


    def check_plot(self, destination=True):
        """
        Check coherency between xml file of parameters and repository.
        Stop if you try to make a plot of a file which does not exsist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                if case.plot == "on" and case.is_run != "KO":
                    for node in self.__parser.getChilds(case.node, "data"):
                        plots, file, dest, repo = self.__parser.getResult(node)
                        if destination == False:
                            dest = None
                        self.__check_dirs(l, case.label, node, repo, dest)

                    for node in self.__parser.getChilds(case.node, "probes"):
                        file, dest, fig = self.__parser.getProbes(node)
                        if destination == False:
                            dest = None
                        repo = None
                        self.__check_dirs(l, case.label, node, repo, dest)

                    for node in self.__parser.getChilds(case.node, "resu"):
                        plots, file, dest, repo = self.__parser.getResult(node)
                        if destination == False:
                            dest = None
                        self.__check_dirs(l, case.label, node, repo, dest)


    def plot(self):
        """
        Plot data.
        """
        for l, s in self.studies:
            self.reporting('  o Plot study: ' + l)
            self.__plotter.plot_study(l, s)

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
        dest = self.getDestination()

        # Fisrt global report
        doc1 = Report1(dest,
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
            doc2 = Report2(dest, report2, self.__log)

            for l, s in self.studies:
                doc2.appendLine("\\section{%s}" % l)

                if s.matplotlib_figures or s.vtk_figures:
                    doc2.appendLine("\\subsection{Graphical results}")
                    for png in s.matplotlib_figures:
                        doc2.addFigure(png)
                    for png in s.vtk_figures:
                        doc2.addFigure(png)

                for case in s.Cases:
                    if case.is_compare == "done":
                        doc2.appendLine("\\subsection{Comparison for case %s}" % case.label)
                        if case.diff_value:
                            doc2.add_row(case.diff_value, l, case.label)
                        elif self.__compare:
                            doc2.appendLine("No difference between the repository and the destination.")

            attached_files.append(doc2.close())

        return attached_files


    def logs(self):
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
