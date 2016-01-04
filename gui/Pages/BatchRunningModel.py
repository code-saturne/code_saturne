# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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

"""
This module modify the batch file
- BatchRunningModel
"""
#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, os.path, shutil, sys, string, types, re

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import Base.Toolbox as Tool
from Pages.SolutionDomainModel import MeshModel, SolutionDomainModel
from Pages.CoalCombustionModel import CoalCombustionModel
from Pages.AtmosphericFlowsModel import AtmosphericFlowsModel
from Base.XMLvariables import Variables, Model

#-------------------------------------------------------------------------------
# Class BatchRunningModel
#-------------------------------------------------------------------------------

class BatchRunningModel(Model):
    """
    This class modifies the batch file (runcase_batch)
    """
    def __init__(self, parent, case):
        """
        Constructor.
        """
        self.parent = parent
        self.case = case

        if not self.case['batch']:
            self.case['batch'] = ""

        self.dictValues = {}

        self.dictValues['job_name'] = None
        self.dictValues['job_nodes'] = None
        self.dictValues['job_ppn'] = None
        self.dictValues['job_procs'] = None
        self.dictValues['job_walltime'] = None
        self.dictValues['job_class'] = None
        self.dictValues['job_account'] = None
        self.dictValues['job_wckey'] = None

        # Is a batch file present ?

        self.batch1 = os.path.join(self.case['scripts_path'], self.case['batch'])

        if os.path.isfile(self.batch1):
            if self.parent.batch_lines:
                self.parseBatchFile()
            else:
                self.readBatchFile()


    def preParse(self, s):
        """
        Pre-parse batch file lines
        """
        r = ' '
        i = s.find('#')
        if i > -1:
            s = s[:i]
        s = r.join(s.split())

        return s


    def parseBatchCCC(self):
        """
        Parse Platform LSF batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#MSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                if len(tok) < 2:
                    continue
                kw = tok[0]
                val = tok[1].split(',')[0].strip()
                if kw == '-r':
                    self.dictValues['job_name'] = val
                elif kw == '-n':
                    self.dictValues['job_procs'] = int(val)
                elif kw == '-N':
                    self.dictValues['job_nodes'] = int(val)
                elif kw == '-T':
                    self.dictValues['job_walltime'] = int(val)
                elif kw == '-q':
                        self.dictValues['job_class'] = val


    def updateBatchCCC(self):
        """
        Update the Platform LSF batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#MSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                if len(tok) < 2:
                    continue
                kw = tok[0]
                if kw == '-r':
                    val = str(self.dictValues['job_name'])
                elif kw == '-n':
                    val = str(self.dictValues['job_procs'])
                elif kw == '-N':
                    val = str(self.dictValues['job_nodes'])
                elif kw == '-T':
                    val = str(self.dictValues['job_walltime'])
                elif kw == '-q':
                    val = self.dictValues['job_class']
                else:
                    continue
                batch_lines[i] = '#MSUB ' + kw + ' ' + str(val) + '\n'


    def parseBatchLOADL(self):
        """
        Parse LoadLeveler batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0] == '#':
                batch_args = self.preParse(batch_lines[i][1:])
                try:
                    if batch_args[0] == '@':
                        kw, val = batch_args[1:].split('=')
                        kw = kw.strip()
                        val = val.split(',')[0].strip()
                        if kw == 'job_name':
                            self.dictValues['job_name'] = val
                        elif kw == 'node':
                            self.dictValues['job_nodes'] = val
                        elif kw == 'tasks_per_node':
                            self.dictValues['job_ppn'] = val
                        elif kw == 'total_tasks':
                            self.dictValues['job_procs'] = val
                        elif kw == 'wall_clock_limit':
                            wt = (val.split(',')[0].rstrip()).split(':')
                            if len(wt) == 3:
                                self.dictValues['job_walltime'] \
                                    = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                            elif len(wt) == 2:
                                self.dictValues['job_walltime'] \
                                    = int(wt[0])*60 + int(wt[1])
                            elif len(wt) == 1:
                                self.dictValues['job_walltime'] = int(wt[0])
                        elif kw == 'class':
                            self.dictValues['job_class'] = val
                        elif kw == 'group':
                            self.dictValues['job_account'] = val
                except Exception:
                    pass


    def updateBatchLOADL(self):
        """
        Update the LoadLeveler batch file from dictionary self.dictValues.
        """

        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0] == '#':
                batch_args = self.preParse(batch_lines[i][1:])
                try:
                    if batch_args[0] == '@':
                        kw, val = batch_args[1:].split('=')
                        kw = kw.strip()
                        val = val.split(',')[0].strip()
                        if kw == 'job_name':
                            val = self.dictValues['job_name']
                        elif kw == 'node':
                            val = self.dictValues['job_nodes']
                        elif kw == 'tasks_per_node':
                            val = self.dictValues['job_ppn']
                        elif kw == 'total_tasks':
                            val = self.dictValues['job_procs']
                        elif kw == 'wall_clock_limit':
                            wt = self.dictValues['job_walltime']
                            val = '%d:%02d:%02d' % (wt/3600,
                                                    (wt%3600)/60,
                                                    wt%60)
                        elif kw == 'class':
                            val = self.dictValues['job_class']
                        elif kw == 'group':
                            val = self.dictValues['job_account']
                        else:
                            continue
                        batch_lines[i] = '# @ ' + kw + ' = ' + str(val) + '\n'
                except Exception:
                    pass


    def parseBatchLSF(self):
        """
        Parse Platform LSF batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#BSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                kw = tok[0]
                val = tok[1].split(',')[0].strip()
                if kw == '-J':
                    self.dictValues['job_name'] = val
                elif kw == '-n':
                    self.dictValues['job_procs'] = int(val)
                elif kw == '-W' or kw == '-wt' or kw == '-We':
                    wt = val.split(':')
                    if len(wt) == 1:
                        self.dictValues['job_walltime'] = int(wt[0])*60
                    elif len(wt) == 2:
                        self.dictValues['job_walltime'] \
                            = int(wt[0])*3600 + int(wt[1])*60
                elif kw == '-q':
                        self.dictValues['job_class'] = val


    def updateBatchLSF(self):
        """
        Update the Platform LSF batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#BSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                kw = tok[0]
                if kw == '-J':
                    val = str(self.dictValues['job_name'])
                elif kw == '-n':
                    val = str(self.dictValues['job_procs'])
                elif kw == '-W' or kw == '-wt' or kw == '-We':
                    wt = self.dictValues['job_walltime']
                    val = '%d:%02d' % (wt/3600, (wt%3600)/60)
                elif kw == '-q':
                    val = self.dictValues['job_class']
                else:
                    continue
                batch_lines[i] = '#BSUB ' + kw + ' ' + str(val) + '\n'


    def parseBatchPBS(self):
        """
        Parse PBS batch file lines
        """
        batch_lines = self.parent.batch_lines

        # TODO: specialize for PBS Professional and TORQUE (OpenPBS has not been
        # maintained since 2004, so we do not support it).
        # We use the "-l nodes=N:ppn=P" syntax here, which is common to all PBS
        # variants, but PBS Pro considers the syntax depecated, and prefers its
        # own "-l select=N:ncpus=P:mpiprocs=P" syntax.
        # We do not have access to a PBS Professional system, but according to
        # its documentation, it has commands such as "pbs-report" or "pbs_probe"
        # which are not part of TORQUE, while the latter has "pbsnodelist" or
        # #pbs-config". The presence of either could help determine which
        # system is available.

        for i in range(len(batch_lines)):
            if batch_lines[i][0:4] == '#PBS':
                batch_args = ' ' + self.preParse(batch_lines[i][4:])
                index = string.rfind(batch_args, ' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        self.dictValues['job_name'] = arg.split()[1]
                    elif arg[0:9] == '-l nodes=':
                        arg_tmp = arg[9:].split(':')
                        self.dictValues['job_nodes'] = arg_tmp[0]
                        self.dictValues['job_ppn'] = arg_tmp[1].split('=')[1]
                    elif arg[0:12] == '-l walltime=':
                        wt = (arg.split('=')[1]).split(':')
                        if len(wt) == 3:
                            self.dictValues['job_walltime'] \
                                = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                        elif len(wt) == 2:
                            self.dictValues['job_walltime'] \
                                = int(wt[0])*60 + int(wt[1])
                        elif len(wt) == 1:
                            self.dictValues['job_walltime'] \
                                = int(wt[0])
                    elif arg[0:2] == '-q':
                            self.dictValues['job_class'] = arg.split()[1]
                    index = string.rfind(batch_args, ' -')


    def updateBatchPBS(self):
        """
        Update the PBS batch file from dictionary self.dictValues.
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:4] == '#PBS':
                ch = '\n'
                batch_args = ' ' + self.preParse(batch_lines[i][4:])
                index = string.rfind(batch_args, ' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        ch = ' -N ' + self.dictValues['job_name'] + ch
                    elif arg[0:9] == '-l nodes=':
                        ch = ' -l nodes=' + self.dictValues['job_nodes'] \
                            +  ':ppn=' + self.dictValues['job_ppn'] + ch
                    elif arg[0:12] == '-l walltime=':
                        wt = self.dictValues['job_walltime']
                        s_wt = '%d:%02d:%02d' % (wt/3600,
                                                 (wt%3600)/60,
                                                 wt%60)
                        ch = ' -l walltime=' + s_wt + ch
                    elif arg[0:2] == '-q':
                        ch = ' -q ' + self.dictValues['job_class'] + ch
                    else:
                        ch = ' ' + arg + ch
                    index = string.rfind(batch_args, ' -')
                ch = '#PBS' + ch
                batch_lines[i] = ch


    def parseBatchSGE(self):
        """
        Parse Sun Grid Engine batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:2] == '#$':
                batch_args = ' ' + self.preParse(batch_lines[i][2:])
                index = string.rfind(batch_args, ' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        self.dictValues['job_name'] = arg.split()[1]
                    elif arg[0:3] == '-pe':
                        try:
                            arg_tmp = arg[3:].split(' ')
                            self.dictValues['job_procs'] = arg_tmp[2]
                        except Exception:
                            pass
                    elif arg[0:8] == '-l h_rt=':
                        wt = (arg.split('=')[1]).split(':')
                        if len(wt) == 3:
                            self.dictValues['job_walltime'] \
                                = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                        elif len(wt) == 2:
                            self.dictValues['job_walltime'] \
                                = int(wt[0])*60 + int(wt[1])
                        elif len(wt) == 1:
                            self.dictValues['job_walltime'] = int(wt[0])
                    elif arg[0:2] == '-q':
                        self.dictValues['job_class'] = arg.split()[1]
                    index = string.rfind(batch_args, ' -')


    def updateBatchSGE(self):
        """
        Update the Sun Grid Engine batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:2] == '#$':
                ch = '\n'
                batch_args = ' ' + self.preParse(batch_lines[i][2:])
                index = string.rfind(batch_args, ' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        ch = ' -N ' + self.dictValues['job_name'] + ch
                    elif arg[0:3] == '-pe':
                        try:
                            arg_tmp = arg[3:].split(' ')
                            ch = ' -pe ' + arg_tmp[1] + ' ' \
                                + str(self.dictValues['job_procs']) + ch
                        except Exception:
                            pass
                    elif arg[0:8] == '-l h_rt=':
                        wt = self.dictValues['job_walltime']
                        s_wt = '%d:%02d:%02d' % (wt/3600,
                                                 (wt%3600)/60,
                                                 wt%60)
                        ch = ' -l h_rt=' + s_wt + ch
                    elif arg[0:2] == '-q':
                        ch = ' -q ' + self.dictValues['job_class'] + ch
                    else:
                        ch = ' ' + arg + ch
                    index = string.rfind(batch_args, ' -')
                    ch = '#$' + ch
                    batch_lines[i] = ch


    def parseBatchSLURM(self):
        """
        Parse SLURM batch file lines
        """
        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:7] == '#SBATCH':
                batch_args = self.preParse(batch_lines[i][7:])
                if batch_args[0:2] == '--':
                    tok = batch_args.split('=')
                    if len(tok) < 2:
                        continue
                    kw = tok[0] + '='
                    val = tok[1].split(',')[0].strip()
                elif batch_args[0] == '-':
                    kw = batch_args[0:2]
                    val = batch_args[2:].split(',')[0].strip()
                else:
                    continue
                if kw == '--job-name=' or kw == '-J':
                    self.dictValues['job_name'] = val
                elif kw == '--ntasks=' or kw == '-n':
                    self.dictValues['job_procs'] = val
                elif kw == '--nodes=' or kw == '-N':
                    self.dictValues['job_nodes'] = val
                elif kw == '--ntasks-per-node=':
                    self.dictValues['job_ppn'] = val
                elif kw == '--time=' or kw == '-t':
                    wt0 = val.split('-')
                    if len(wt0) == 2:
                        th = int(wt0[0])*3600*24
                        wt = wt0[1].split(':')
                    else:
                        th = 0
                        wt = wt0[0].split(':')
                    if len(wt) == 3:
                        self.dictValues['job_walltime'] \
                            = th + int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                    elif len(wt) == 2:
                        if len(wt0) == 2:
                            self.dictValues['job_walltime'] \
                                = th + int(wt[0])*3600 + int(wt[1])*60
                        else:
                            self.dictValues['job_walltime'] \
                                = th + int(wt[0])*60 + int(wt[1])
                    elif len(wt) == 1:
                        if len(wt0) == 2:
                            self.dictValues['job_walltime'] \
                                = th + int(wt[0])*3600
                        else:
                            self.dictValues['job_walltime'] \
                                = th + int(wt[0])*60
                elif kw == '--partition=' or kw == '-p':
                    self.dictValues['job_class'] = val
                elif kw == '--account=' or kw == '-A':
                    self.dictValues['job_account'] = val
                elif kw == '--wckey=':
                    self.dictValues['job_wckey'] = val


    def updateBatchSLURM(self):
        """
        Update the SLURM batch file from dictionary self.dictValues.
        """
        batch_lines = self.parent.batch_lines

        batch_lines = self.parent.batch_lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:7] == '#SBATCH':
                batch_args = self.preParse(batch_lines[i][7:])
                if batch_args[0:2] == '--':
                    tok = batch_args.split('=')
                    if len(tok) < 2:
                        continue
                    kw = tok[0] + '='
                    val = tok[1].split(',')[0].strip()
                elif batch_args[0] == '-':
                    kw = batch_args[0:2]
                    val = batch_args[2:].split(',')[0].strip()
                if kw == '--job-name=' or kw == '-J':
                    val = str(self.dictValues['job_name'])
                elif kw == '--ntasks=' or kw == '-n':
                    val = str(self.dictValues['job_procs'])
                elif kw == '--nodes=' or kw == '-N':
                    val = str(self.dictValues['job_nodes'])
                elif kw == '--ntasks-per-node=':
                    val = self.dictValues['job_ppn']
                elif kw == '--time=' or kw == '-t':
                    wt = self.dictValues['job_walltime']
                    if wt > 86400: # 3600*24
                        val = '%d-%d:%02d:%02d' % (wt/86400,
                                                   (wt%86400)/3600,
                                                   (wt%3600)/60,
                                                   wt%60)
                    else:
                        val = '%d:%02d:%02d' % (wt/3600,
                                                (wt%3600)/60,
                                                wt%60)
                elif kw == '--partition=' or kw == '-p':
                    val = self.dictValues['job_class']
                elif kw == '--account=' or kw == '-A':
                    val = self.dictValues['job_account']
                elif kw == '--wckey=':
                    val = self.dictValues['job_wckey']
                else:
                    continue
                batch_lines[i] = '#SBATCH ' + kw + str(val) + '\n'


    def readBatchFile(self):
        """
        Fill self.dictValues reading the batch file.
        """

        if not self.case['batch']:
            return

        self.batch1 = os.path.join(self.case['scripts_path'], self.case['batch'])
        if not os.path.isfile(self.batch1):
            return

        # Read the batch file line by line.
        # All lines are stored in a list called "self.batch_lines".

        f = open(self.batch1, 'r')
        batch_lines = f.readlines()
        f.close()

        for i in range(len(batch_lines)):
            batch_lines[i] = batch_lines[i].rstrip(' \t')

        self.parent.batch_lines = batch_lines


        # Parse lines depending on batch type

        self.parseBatchFile()


    def parseBatchFile(self):
        """
        Fill self.dictValues reading the batch file.
        """

        # Parse lines depending on batch type

        if self.case['batch_type'] == None:
            return

        elif self.case['batch_type'][0:3] == 'CCC':
            self.parseBatchCCC()
        elif self.case['batch_type'][0:5] == 'LOADL':
            self.parseBatchLOADL()
        elif self.case['batch_type'][0:3] == 'LSF':
            self.parseBatchLSF()
        elif self.case['batch_type'][0:3] == 'PBS':
            self.parseBatchPBS()
        elif self.case['batch_type'][0:3] == 'SGE':
            self.parseBatchSGE()
        elif self.case['batch_type'][0:5] == 'SLURM':
            self.parseBatchSLURM()

    def updateBatchFile(self, keyword=None):
        """
        Update the batch file from reading dictionary self.dictValues.
        If keyword == None, all keywords are updated
        If keyword == key, only key is updated.
        """
        l = list(self.dictValues.keys())
        l.append(None) # Add 'None' when no keyword is specified in argument.
        for k in list(self.dictValues.keys()):
            if self.dictValues[k] == 'None':
                self.dictValues[k] = None
        self.isInList(keyword, l)

        batch_type = self.case['batch_type']
        if batch_type:
            if batch_type[0:3] == 'CCC':
                self.updateBatchCCC()
            elif batch_type[0:5] == 'LOADL':
                self.updateBatchLOADL()
            elif batch_type[0:3] == 'LSF':
                self.updateBatchLSF()
            elif batch_type[0:3] == 'PBS':
                self.updateBatchPBS()
            elif batch_type[0:3] == 'SGE':
                self.updateBatchSGE()
            elif batch_type[0:5] == 'SLURM':
                self.updateBatchSLURM()


    def _getRegex(self, word):
        """
        Get regular expression to extract line without comment
        """
        # handles both comments and "$word":
        regex = re.compile(r"""(^(?#)^\s*(?<!$)""" + word + r""".*$)""")

        return regex


#-------------------------------------------------------------------------------
# BatchRunningModel test class
#-------------------------------------------------------------------------------

class BatchRunningModelTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all 'check' methods.
        """
        from Base.XMLengine import Case
        from Base.XMLinitialize import XMLinit
        from Base.Toolbox import GuiParam
        GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case).initialize()

        self.case['batch_type'] = None
        self.case['scripts_path'] = os.getcwd()
        self.case['batch'] = 'runcase_batch'

        lance_PBS = '# test \n'\
        '#\n'\
        '#                  CARTES BATCH POUR CLUSTERS sous PBS\n'\
        '#\n'\
        '#PBS -l nodes=16:ppn=1,walltime=34:77:22\n'\
        '#PBS -j eo -N super_toto\n'

        lance_LSF = '# test \n'\
        '#\n'\
        '#        CARTES BATCH POUR LE CCRT (Platine sous LSF)\n'\
        '#\n'\
        '#BSUB -n 2\n'\
        '#BSUB -c 00:05\n'\
        '#BSUB -o super_tataco.%J\n'\
        '#BSUB -e super_tatace.%J\n'\
        '#BSUB -J super_truc\n'

        self.f = open('lance_PBS','w')
        self.f.write(lance_PBS)
        self.f.close()
        self.f = open('lance_LSF','w')
        self.f.write(lance_LSF)
        self.f.close()


    def tearDown(self):
        """
        This method is executed after all 'check' methods.
        """
        f = self.case['batch']
        if os.path.isfile(f): os.remove(f)
        if os.path.isfile(f+"~"): os.remove(f+"~")


    def checkReadBatchPBS(self):
        """ Check whether the BatchRunningModel class could be read file"""
        self.case['batch_type'] = 'PBS'
        mdl = BatchRunningModel(self.case)
        mdl.readBatchFile()

        dico_PBS = {\
        'job_nodes': '16',
        'job_name': 'super_toto',
        'job_ppn': '1',
        'job_walltime': '34:77:22'}

        for k in list(dico_PBS.keys()):
            if mdl.dictValues[k] != dico_PBS[k] :
                print("\nwarning for key: ", k)
                print("  read value in the batch description:", mdl.dictValues[k])
                print("  reference value:", dico_PBS[k])
            assert  mdl.dictValues[k] == dico_PBS[k], 'could not read the batch file'


    def checkUpdateBatchFile(self):
        """ Check whether the BatchRunningModel class could update file"""
        mdl = BatchRunningModel(self.case)
        mdl.readBatchFile()
        mdl.dictValues['N_PROCS']=48
        dico_updated = mdl.dictValues
        mdl.updateBatchFile()
        mdl.readBatchFile()
        dico_read = mdl.dictValues

        assert dico_updated == dico_read, 'error on updating batch script file'


    def checkUpdateBatchPBS(self):
        """ Check whether the BatchRunningModel class could update file"""
        mdl = BatchRunningModel(self.case)
        mdl.readBatchFile()
        mdl.dictValues['job_walltime']='12:42:52'
        dicojob_updated = mdl.dictValues
        mdl.updateBatchFile()
        mdl.readBatchFile()
        dicojob_read = mdl.dictValues

        assert dicojob_updated == dicojob_read, 'error on updating PBS batch script file'


def suite():
    testSuite = unittest.makeSuite(BatchRunningModelTestCase, "check")
    return testSuite


def runTest():
    print("BatchRunningModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End of BatchRunningModel
#-------------------------------------------------------------------------------
