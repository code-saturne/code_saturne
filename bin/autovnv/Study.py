# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne Scripts, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2011 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
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

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from autovnv.Parser import Parser
from autovnv.TexMaker import Report1, Report2
from autovnv.Drawing import Plotter
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
    def __init__(self, rlog, exe, diff, study, data, repo, dest):
        """
        @type data: C{Dictionary}
        @param data: contains all keyword and value read in the parameters file
        """
        self.__log      = rlog
        self.__exe      = exe
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

        self.figures    = []
        self.verbatim   = []


    def compile(self):
        """
        Just compile user sources if exist.
        @rtype: C{String}
        @return: the status of the succes of the compilation.
        """
        home = os.getcwd()
        os.chdir(os.path.join(self.__dest, self.label, 'SRC'))
        cmd = self.__exe + " compile -t"

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

        cmd = self.__exe + " run --suggest-id"
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
            run_ref_f = file(run_ref, mode='r')
        except IOError:
            print "Error: can not opening %s\n" % run_ref
            sys.exit(1)

        patterns = r'^\\code_saturne', r'^\\neptune_cfd'

        for line in run_ref_f.readlines():
            for p in patterns:
                if re.search(p, line):
                    run_cmd = string.join(line.split())
        run_ref_f.close()

        # Write the new runcase

        try:
            run_new_f = file(run_new, mode='r')
        except IOError:
            print "Error: can not opening %s\n" % run_new
            sys.exit(1)

        lines = run_new_f.readlines()
        run_new_f.close()

        for i in range(len(lines)):
            for p in patterns:
                if re.search(p, lines[i]):
                    lines[i] = run_cmd + " --id=" + run_id

        run_new_f = file(run_new, mode='w')
        run_new_f.writelines(lines)
        run_new_f.close()


    def run(self):
        """
        Run the case a thread.
        """
        run_id, run_dir = self.__suggest_run_id()

        while os.path.isdir(run_dir):
            time.sleep(5)
            run_id, run_dir = self.__suggest_run_id()

        self.__updateRuncase(run_id)

        home = os.getcwd()
        os.chdir(os.path.join(self.__dest, self.label, 'SCRIPTS'))

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


    def compare(self, r, d, thresold):
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

        if thresold != None:
            cmd += ' --threshold ' + thresold

        l = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
        info = l.read().replace("\"", " ").replace(";", " ").replace(":", " ").split()

        tab = []

        for i in range(len(info)):
            if info[i][:4] == 'Diff':
                if info[i-3] not in ['i4', 'u4']:
                    tab.append([info[i-7].replace("_", "\_"), info[i+3], info[i+5]])

        os.chdir(home)

        return tab

#-------------------------------------------------------------------------------

class Study(object):
    """
    Create, run and compare all cases fir a given study.
    """
    def __init__(self, parser, study, exe, dif, rlog):
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
        self.__parser  = parser
        self.__study   = study
        self.__exe     = exe
        self.__diff    = dif
        self.__log     = rlog

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
            c = Case(self.__log,
                     self.__exe,
                     self.__diff,
                     self.__study,
                     data,
                     self.__repo,
                     self.__dest)
            self.Cases.append(c)

        self.matplotlib_figures = []


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
            cmd = self.__exe + " create --quiet --study " + self.__dest
            retval, t = run_command(cmd, self.__log)
            shutil.rmtree(os.path.join(self.__dest, "CASE1"))

            # Link meshes
            ref = os.path.join(self.__repo, "MESH")
            if os.path.isdir(ref):
                des = os.path.join(self.__dest, "MESH")
                for m in os.listdir(ref):
                    os.symlink(os.path.join(ref, m), os.path.join(des, m))

            # Copy external scripts for post-processing
            ref = os.path.join(self.__repo, "POST")
            if os.path.isdir(ref):
                des = os.path.join(self.__dest, "POST")
                shutil.rmtree(des)
                shutil.copytree(ref, des, symlinks=True)

        # Create cases
        os.chdir(self.__dest)

        for c in self.cases:
            if not os.path.isdir(c):
                cmd = self.__exe + " create --case " + c  \
                      + " --quiet --noref --copy-from " \
                      + os.path.join(self.__repo, c)
                retval, t = run_command(cmd, self.__log)
            else:
                print "Warning: the case %s already exists in the destination." % c

        os.chdir(repbase)

#-------------------------------------------------------------------------------

class Studies(object):
    """
    Manage all Studies and all Cases described in the files of parameters.
    """
    def __init__(self, f, v, r, c, p, exe, dif):
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

        # Create the xml parser

        self.__parser = Parser(f)

        # Test if the repository exists

        self.getRepository()

        # create if necessary the destination directory

        try:
            os.makedirs(self.getDestination())
        except:
            pass

        # initialize the parser and the plotter

        file = os.path.join(self.getDestination(), f)
        try:
            shutil.copyfile(f, file)
        except:
            pass
        self.__parser  = Parser(file)
        self.__plotter = Plotter(self.__parser)

        # build the list of the studies

        doc = os.path.join(self.getDestination(), "auto_vnv.log")
        self.__log = open(doc, "w")
        self.labels  = self.__parser.getStudiesLabel()
        self.studies = []
        for l in self.labels:
            self.studies.append( [l, Study(self.__parser, l, exe, dif, self.__log)] )

        # start the report
        self.report = os.path.join(self.getDestination(), "report.txt")
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
                    if case.compile() == "OK":
                        self.reporting('    - compile %s --> OK' % case.label)
                    else:
                        self.reporting('    - compile %s --> FAILED' % case.label)
                        iok+=1
        if iok:
            self.reporting('Error: compilation failed. Number of failed cases: %s' % iok)
            sys.exit(1)


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
                            for n in n1 + n2 + n3:
                                if self.__parser.getAttribute(n, "dest", False) == "":
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
                msg = "Study %s case %s:\nthe directory %s does not exist.\n" % \
                    (study_label, case_label, rep_f)
                raise ValueError, msg

        # 2. The result directory must be read automatically;
        #    check if there is a single result directory.
        elif rep == "":
            if not (len(os.listdir(result)) == 1):
                msg = "Study %s case %s:\nthere is not a single result directory in %s\n" % \
                    (study_label, case_label, result)
                raise ValueError, msg

        # 3. Update the file of parameters with the name of the result directory
            self.__parser.setAttribute(node, attr, os.listdir(result)[0])
        else:
            self.reporting('Error: check compare/script/plot failed.')
            sys.exit(1)


    def __check_dirs(self, study_label, case_label, node, repo, dest):
        """
        Check coherency between xml file of parameters and repository and destination.
        """
        if repo != None:
            result = os.path.join(self.getRepository(), study_label, case_label, 'RESU')
            self.__check_dir(study_label, case_label, node, result, repo, "repo")

        if dest != None:
            result = os.path.join(self.getDestination(), study_label, case_label, 'RESU')
            self.__check_dir(study_label, case_label, node, result, dest, "dest")


    def check_compare(self):
        """
        Check coherency between xml file of parameters and repository.
        Stops if one try to make a comparison with a file which does not exsist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                bool, repo, dest, threshold = self.__parser.getCompare(case.node)
                if bool and case.is_run != "KO":
                    node = self.__parser.getChild(case.node, "compare")
                    self.__check_dirs(l, case.label, node, repo, dest)


    def compare(self):
        """
        Compare the results of the new computations with thoses from the Repository.
        """
        if self.__compare:
            for l, s in self.studies:
                self.reporting('  o Compare study: ' + l)
                for case in s.Cases:
                    case.is_compare, repo, dest, case.threshold = self.__parser.getCompare(case.node)
                    if case.is_compare and case.is_run != "KO":
                        self.reporting('    - compare %s' % case.label)
                        case.diff_value = case.compare(repo, dest, case.threshold)


    def check_script(self):
        """
        Check coherency between xml file of parameters and repository.
        Stops if one try to run a script with a file which does not exsist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                bool, label, nodes, args, repo, dest = self.__parser.getScript(case.node)
                if bool and case.is_run != "KO":
                    for i in range(len(label)):
                        self.__check_dirs(l, case.label, nodes[i], repo[i], dest[i])


    def scripts(self):
        """
        Launch external additional scripts with arguments.
        """
        for l, s in self.studies:
            self.reporting('  o Script study: ' + l)
            for case in s.Cases:
                bool, label, nodes, args, repo, dest = self.__parser.getScript(case.node)
                if bool and case.is_run != "KO":
                    for i in range(len(label)):
                        cmd = os.path.join(self.getDestination(), l, "POST", label[i])
                        if os.path.isfile(cmd):
                            cmd += " " + args[i]
                            if repo[i]:
                                r = os.path.join(self.getRepository(),  l, case.label, "RESU", repo[i])
                                cmd += " -r " + r
                            if dest[i]:
                                d = os.path.join(self.getDestination(), l, case.label, "RESU", dest[i])
                                cmd += " -d " + d
                            retcode, t = run_command(cmd, self.__log)
                            self.reporting('    - script %s --> OK (%s s)' % (cmd, t))
                        else:
                            self.reporting('    - script %s not found' % cmd)


    def check_plot(self):
        """
        Check coherency between xml file of parameters and repository.
        Stops if one try to make a plot of a file which does not exsist.
        """
        for l, s in self.studies:
            for case in s.Cases:
                if case.plot == "on" and case.is_run != "KO":
                    for node in self.__parser.getChilds(case.node, "data"):
                        plots, file, dest, repo = self.__parser.getResult(node)
                        self.__check_dirs(l, case.label, node, repo, dest)

    def plot(self):
        """
        Draw plots.
        """
        for l, s in self.studies:
            self.reporting('  o Plot study: ' + l)
            self.__plotter.plot_study(l, s)


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
                             is_nodiff,
                             case.threshold)

        attached_files.append(doc1.close())

        # Second detailed report
        if self.__compare or self.__postpro:
            doc2 = Report2(dest, report2, self.__log)

            for l, s in self.studies:
                doc2.appendLine("\\section{%s}" % l)

                if s.matplotlib_figures:
                    doc2.appendLine("\\subsection{Results}")
                    for png in s.matplotlib_figures:
                        doc2.addFigure(png)

                for case in s.Cases:
                    if case.is_compare:
                        doc2.appendLine("\\subsection{%s}" % case.label)
                        if case.diff_value:
                            doc2.add_row(case.diff_value,
                                         l,
                                         case.label,
                                         case.threshold)
                        else:
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
