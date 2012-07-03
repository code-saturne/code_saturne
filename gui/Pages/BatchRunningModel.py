# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2007 EDF S.A., France
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

"""
This module modify the "runcase" script file
- RuncaseModel
- BatchRunningModel
"""
#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, sys, string, types, re

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import Base.Toolbox as Tool
from SolutionDomainModel import SolutionDomainModel
from CoalCombustionModel import CoalCombustionModel
from AtmosphericFlowsModel import AtmosphericFlowsModel
from ProfilesModel import ProfilesModel
from Base.XMLvariables import Variables, Model

import cs_config

#-------------------------------------------------------------------------------
# Class BatchRunningModel
#-------------------------------------------------------------------------------

class RuncaseModel:
    """
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        if not self.case['batchScript']:
            self.case['batchScript'] = ""

        if not self.case['backupBatchScript']:
            self.case['backupBatchScript'] = "no"

        # get type of batch file

        self.case['batch_type'] = ''
        batch_template = cs_config.batch_template.template
        if batch_template and batch_template != 'no':
            if os.path.isabs(batch_template):
                self.case['batch_type'] = os.path.basename(batch_template)
            else:
                self.case['batch_type'] = batch_template


class BatchRunningModel(Model):
    """
    This class modify saturne running file (runcase)
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        RuncaseModel(self.case)

        self.script1 = self.case['scripts_path'] + "/" + self.case['batchScript']

        # Save a backup file before any modification.
        # There is only one backup for the entire session.

        script2 = self.script1 + "~"

        if self.case['backupBatchScript'] == "no" \
               or not os.path.isfile(script2):
            os.popen('cp ' + self.script1 + " " +script2)
            self.case['backupBatchScript'] = "yes"

        # Read the batch script file line by line.
        # All lines are stored in a list called "self.lines".

        f = open(self.script1, 'rw')
        self.lines = f.readlines()
        f.close()

        for i in range(len(self.lines)):
            self.lines[i] = self.lines[i].rstrip(' \t')

        # DicoValues's initialisation

        self.dicoValues = {}

        self.dicoValues['job_name'] = None
        self.dicoValues['job_nodes'] = None
        self.dicoValues['job_ppn'] = None
        self.dicoValues['job_procs'] = None
        self.dicoValues['job_walltime'] = None
        self.dicoValues['job_class'] = None
        self.dicoValues['job_group'] = None

        self.dicoValues['SOLCOM'] = '0'
        self.dicoValues['PARAM'] = ""
        self.dicoValues['NUMBER_OF_PROCESSORS'] = ""
        self.dicoValues['PROCESSOR_LIST'] = ""
        self.dicoValues['PARTITION_LIST'] = ""
        self.dicoValues['USER_INPUT_FILES'] = ""
        self.dicoValues['USER_OUTPUT_FILES'] = ""
        self.dicoValues['CS_TMP_PREFIX'] = ""
        self.dicoValues['MESH'] = ""
        self.dicoValues['COMMAND_REORIENT'] = ""
        self.dicoValues['COMMAND_JOIN'] = ""
        self.dicoValues['COMMAND_CWF'] = ""
        self.dicoValues['COMMAND_PERIO'] = ""
        self.dicoValues['EXEC_PREPROCESS'] = "yes"
        self.dicoValues['EXEC_PARTITION'] = "yes"
        self.dicoValues['EXEC_KERNEL'] = "yes"
        self.dicoValues['CS_LIB_ADD'] = ""
        self.dicoValues['VALGRIND'] = ""
        self.dicoValues['ARG_CS_OUTPUT'] = ""
        self.dicoValues['ARG_CS_VERIF'] = ""
        self.dicoValues['THERMOCHEMISTRY_DATA'] = ""
        self.dicoValues['METEO_DATA'] = ""

        if self.case['salome']:
            self.dicoValues['ARG_CS_OUTPUT'] = "--log 0"


    def _getRegex(self, word):
        """
        Get regular expression to extract line without comment
        """
##  fonctionne mais incomplet:      regex = re.compile(r"""(^\s*""" + word + r""".*$)""")
##  fonctionne en tenant compte des lignes commencant par # :
##      regex = re.compile(r"""(^(?#)^\s*""" + word + r""".*$)""")
        #tient compte a la fois des commentaires et des "$word":
        regex = re.compile(r"""(^(?#)^\s*(?<!$)""" + word + r""".*$)""")

        return regex


    def _getLineToModify(self, regex, txt):
        """
        Search word in txt if and only if it's not a comment
        """
        pattern = None
        if regex != None :
            pattern = regex.search(txt)
        return pattern


    def _substituteLine(self, pattern, newword, txt):
        """
        Substitute pattern by newword
        """
        new_pattern = pattern.re.sub(newword, txt, re.VERBOSE)
        return new_pattern


    def _getValueInPattern(self, pattern, word, dico):
        """
        Return value of pattern
        """
        self.dicoValues[word] = dico
        resu = pattern.group().split('=')
        L = resu[1]
        for i in range(2,len(resu)): L = L  + "="+ resu[i]
        resu = L

        if resu:
            if resu.split(' ')[0] == '': resu = resu.split(' ')[1]
            if resu == '""': resu =''
        else:
            resu =''

        self.dicoValues[word] = resu
        return self.dicoValues[word]


    def _removeQuotes(self, line, index=0, word=""):
        """
        1) Delete quotes and return caracters,
        2) Return  the associated value of the word
        """
        if not line: return ""

        ch = line[index+len(word):]
        if not ch: return ""
        ch = string.join(string.split(ch))

        try:
            if ch[-1:] == '\n': ch = ch[:-1]
        except IndexError:
            pass

        try:
            if ch[-1:] == '"': ch = ch[:-1]
        except IndexError:
            pass

        try:
            if ch[0] == '"': ch = ch[1:]
        except IndexError:
            pass

        ch = string.join(string.split(ch))

        return ch


    def _addQuotes(self, ch):
        """
        Add quotes in front of and behind string c.
        """
        ch = string.join(string.split(ch))
        if string.rfind(ch, " ") != -1:
            ch = '"' + ch + '"'
        return ch


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
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#MSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                if len(tok) < 2:
                    continue
                kw = tok[0]
                val = tok[1].split(',')[0].strip()
                if kw == '-r':
                    self.dicoValues['job_name'] = val
                elif kw == '-n':
                    self.dicoValues['job_procs'] = int(val)
                elif kw == '-N':
                    self.dicoValues['job_nodes'] = int(val)
                elif kw == '-T':
                    self.dicoValues['job_walltime'] = int(val)
                elif kw == '-q':
                        self.dicoValues['job_class'] = val


    def updateBatchCCC(self):
        """
        Update the Platform LSF batch file lines
        """
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#MSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                if len(tok) < 2:
                    continue
                kw = tok[0]
                if kw == '-r':
                    val = str(self.dicoValues['job_name'])
                elif kw == '-n':
                    val = str(self.dicoValues['job_procs'])
                elif kw == '-N':
                    val = str(self.dicoValues['job_nodes'])
                elif kw == '-T':
                    val = str(self.dicoValues['job_walltime'])
                elif kw == '-q':
                    val = self.dicoValues['job_class']
                else:
                    continue
                batch_lines[i] = '#MSUB ' + kw + ' ' + str(val) + '\n'


    def parseBatchLOADL(self):
        """
        Parse LoadLeveler batch file lines
        """
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0] == '#':
                batch_args = self.preParse(batch_lines[i][1:])
                try:
                    if batch_args[0] == '@':
                        kw, val = batch_args[1:].split('=')
                        kw = kw.strip()
                        val = val.split(',')[0].strip()
                        if kw == 'job_name':
                            self.dicoValues['job_name'] = val
                        elif kw == 'node':
                            self.dicoValues['job_nodes'] = val
                        elif kw == 'tasks_per_node':
                            self.dicoValues['job_ppn'] = val
                        elif kw == 'total_tasks':
                            self.dicoValues['job_procs'] = val
                        elif kw == 'wall_clock_limit':
                            wt = (val.split(',')[0].rstrip()).split(':')
                            if len(wt) == 3:
                                self.dicoValues['job_walltime'] \
                                    = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                            elif len(wt) == 2:
                                self.dicoValues['job_walltime'] \
                                    = int(wt[0])*60 + int(wt[1])
                            elif len(wt) == 1:
                                self.dicoValues['job_walltime'] = int(wt[0])
                        elif kw == 'class':
                            self.dicoValues['job_class'] = val
                        elif kw == 'group':
                            self.dicoValues['job_group'] = val
                except Exception:
                    pass


    def updateBatchLOADL(self):
        """
        Update the LoadLeveler batch file from dictionary self.dicoValues.
        """

        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0] == '#':
                batch_args = self.preParse(batch_lines[i][1:])
                try:
                    if batch_args[0] == '@':
                        kw, val = batch_args[1:].split('=')
                        kw = kw.strip()
                        val = val.split(',')[0].strip()
                        if kw == 'job_name':
                            val = self.dicoValues['job_name']
                        elif kw == 'node':
                            val = self.dicoValues['job_nodes']
                        elif kw == 'tasks_per_node':
                            val = self.dicoValues['job_ppn']
                        elif kw == 'total_tasks':
                            val = self.dicoValues['job_procs']
                        elif kw == 'wall_clock_limit':
                            wt = self.dicoValues['job_walltime']
                            val = '%d:%02d:%02d' % (wt/3600,
                                                    (wt%3600)/60,
                                                    wt%60)
                        elif kw == 'class':
                            val = self.dicoValues['job_class']
                        elif kw == 'group':
                            val = self.dicoValues['job_group']
                        else:
                            continue
                        batch_lines[i] = '# @ ' + kw + ' = ' + str(val) + '\n'
                except Exception:
                    pass


    def parseBatchLSF(self):
        """
        Parse Platform LSF batch file lines
        """
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#BSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                kw = tok[0]
                val = tok[1].split(',')[0].strip()
                if kw == '-J':
                    self.dicoValues['job_name'] = val
                elif kw == '-n':
                    self.dicoValues['job_procs'] = int(val)
                elif kw == '-W' or kw == '-wt' or kw == '-We':
                    wt = val.split(':')
                    if len(wt) == 1:
                        self.dicoValues['job_walltime'] = int(wt[0])*60
                    elif len(wt) == 2:
                        self.dicoValues['job_walltime'] \
                            = int(wt[0])*3600 + int(wt[1])*60
                elif kw == '-q':
                        self.dicoValues['job_class'] = val


    def updateBatchLSF(self):
        """
        Update the Platform LSF batch file lines
        """
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#BSUB':
                batch_args = self.preParse(batch_lines[i][5:])
                tok = batch_args.split()
                kw = tok[0]
                if kw == '-J':
                    val = str(self.dicoValues['job_name'])
                elif kw == '-n':
                    val = str(self.dicoValues['job_procs'])
                elif kw == '-W' or kw == '-wt' or kw == '-We':
                    wt = self.dicoValues['job_walltime']
                    val = '%d:%02d' % (wt/3600, (wt%3600)/60)
                elif kw == '-q':
                    val = self.dicoValues['job_class']
                else:
                    continue
                batch_lines[i] = '#BSUB ' + kw + ' ' + str(val) + '\n'


    def parseBatchPBS(self):
        """
        Parse PBS batch file lines
        """
        batch_lines = self.lines

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
                        self.dicoValues['job_name'] = arg.split()[1]
                    elif arg[0:9] == '-l nodes=':
                        arg_tmp = arg[9:].split(':')
                        self.dicoValues['job_nodes'] = arg_tmp[0]
                        self.dicoValues['job_ppn'] = arg_tmp[1].split('=')[1]
                    elif arg[0:12] == '-l walltime=':
                        wt = (arg.split('=')[1]).split(':')
                        if len(wt) == 3:
                            self.dicoValues['job_walltime'] \
                                = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                        elif len(wt) == 2:
                            self.dicoValues['job_walltime'] \
                                = int(wt[0])*60 + int(wt[1])
                        elif len(wt) == 1:
                            self.dicoValues['job_walltime'] \
                                = int(wt[0])
                    elif arg[0:2] == '-q':
                            self.dicoValues['job_class'] = arg.split()[1]
                    index = string.rfind(batch_args, ' -')


    def updateBatchPBS(self):
        """
        Update the PBS batch file from dictionary self.dicoValues.
        """
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:4] == '#PBS':
                ch = '\n'
                batch_args = ' ' + self.preParse(batch_lines[i][4:])
                index = string.rfind(batch_args, ' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        ch = ' -N ' + self.dicoValues['job_name'] + ch
                    elif arg[0:9] == '-l nodes=':
                        ch = ' -l nodes=' + self.dicoValues['job_nodes'] \
                            +  ':ppn=' + self.dicoValues['job_ppn'] + ch
                    elif arg[0:12] == '-l walltime=':
                        wt = self.dicoValues['job_walltime']
                        s_wt = '%d:%02d:%02d' % (wt/3600,
                                                 (wt%3600)/60,
                                                 wt%60)
                        ch = ' -l walltime=' + s_wt + ch
                    elif arg[0:2] == '-q':
                        ch = ' -q ' + self.dicoValues['job_class'] + ch
                    else:
                        ch = ' ' + arg + ch
                    index = string.rfind(batch_args, ' -')
                ch = '#PBS' + ch
                batch_lines[i] = ch


    def parseBatchSGE(self):
        """
        Parse Sun Grid Engine batch file lines
        """
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:2] == '#$':
                batch_args = ' ' + self.preParse(batch_lines[i][2:])
                index = string.rfind(batch_args, ' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        self.dicoValues['job_name'] = arg.split()[1]
                    elif arg[0:3] == '-pe':
                        try:
                            arg_tmp = arg[3:].split(' ')
                            self.dicoValues['job_procs'] = arg_tmp[2]
                        except Exception:
                            pass
                    elif arg[0:8] == '-l h_rt=':
                        wt = (arg.split('=')[1]).split(':')
                        if len(wt) == 3:
                            self.dicoValues['job_walltime'] \
                                = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                        elif len(wt) == 2:
                            self.dicoValues['job_walltime'] \
                                = int(wt[0])*60 + int(wt[1])
                        elif len(wt) == 1:
                            self.dicoValues['job_walltime'] = int(wt[0])
                    elif arg[0:2] == '-q':
                        self.dicoValues['job_class'] = arg.split()[1]
                    index = string.rfind(batch_args, ' -')


    def updateBatchSGE(self):
        """
        Update the Sun Grid Engine batch file lines
        """
        batch_lines = self.lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:2] == '#$':
                ch = '\n'
                batch_args = ' ' + self.preParse(batch_lines[i][2:])
                index = string.rfind(batch_args, ' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        ch = ' -N ' + self.dicoValues['job_name'] + ch
                    elif arg[0:3] == '-pe':
                        try:
                            arg_tmp = arg[3:].split(' ')
                            ch = ' -pe ' + arg_tmp[1] + ' ' \
                                + str(self.dicoValues['job_procs']) + ch
                        except Exception:
                            pass
                    elif arg[0:8] == '-l h_rt=':
                        wt = self.dicoValues['job_walltime']
                        s_wt = '%d:%02d:%02d' % (wt/3600,
                                                 (wt%3600)/60,
                                                 wt%60)
                        ch = ' -l h_rt=' + s_wt + ch
                    elif arg[0:2] == '-q':
                        ch = ' -q ' + self.dicoValues['job_class'] + ch
                    else:
                        ch = ' ' + arg + ch
                    index = string.rfind(batch_args, ' -')
                    ch = '#$' + ch
                    batch_lines[i] = ch


    def parseBatchSLURM(self):
        """
        Parse SLURM batch file lines
        """
        batch_lines = self.lines

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
                    self.dicoValues['job_name'] = val
                elif kw == '--ntasks=' or kw == '-n':
                    self.dicoValues['job_procs'] = val
                elif kw == '--nodes=' or kw == '-N':
                    self.dicoValues['job_nodes'] = val
                elif kw == '--ntasks-per-node=':
                    self.dicoValues['job_ppn'] = val
                elif kw == '--time=' or kw == '-t':
                    wt0 = val.split('-')
                    if len(wt0) == 2:
                        th = int(wt0[0])*3600*24
                        wt = wt0[1].split(':')
                    else:
                        th = 0
                        wt = wt0[0].split(':')
                    if len(wt) == 3:
                        self.dicoValues['job_walltime'] \
                            = th + int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                    elif len(wt) == 2:
                        if len(wt0) == 2:
                            self.dicoValues['job_walltime'] \
                                = th + int(wt[0])*3600 + int(wt[1])*60
                        else:
                            self.dicoValues['job_walltime'] \
                                = th + int(wt[0])*60 + int(wt[1])
                    elif len(wt) == 1:
                        if len(wt0) == 2:
                            self.dicoValues['job_walltime'] \
                                = th + int(wt[0])*3600
                        else:
                            self.dicoValues['job_walltime'] \
                                = th + int(wt[0])*60
                elif kw == '--partition=' or kw == '-p':
                    self.dicoValues['job_class'] = val


    def updateBatchSLURM(self):
        """
        Update the SLURM batch file from dictionary self.dicoValues.
        """
        batch_lines = self.lines

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
                    val = str(self.dicoValues['job_name'])
                elif kw == '--ntasks=' or kw == '-n':
                    val = str(self.dicoValues['job_procs'])
                elif kw == '--nodes=' or kw == '-N':
                    val = str(self.dicoValues['job_nodes'])
                elif kw == '--ntasks-per-node=':
                    val = self.dicoValues['job_ppn']
                elif kw == '--time=' or kw == '-t':
                    wt = self.dicoValues['job_walltime']
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
                    val = self.dicoValues['job_class']
                else:
                    continue
                batch_lines[i] = '#SBATCH ' + kw + str(val) + '\n'


    def readBatchScriptFile(self):
        """
        Fill self.dicoValues reading the backup file.
        """

        if self.case['batch_type'] != '':
            if self.case['batch_type'][0:3] == 'CCC':
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

        lines = self.lines

        # Other keywords
        for k in self.dicoValues.keys():
            if k not in ('THERMOCHEMISTRY_DATA', 'METEO_DATA') \
                    and k[0:4] != 'job_':
                nbkey = 0
                for i in range(len(lines)):
                    reg = self._getRegex(k)
                    if reg != None:
                        pat = self._getLineToModify(reg, lines[i])
                        if pat != None:
                            nbkey = nbkey + 1
                            if nbkey == 1:
                                ch = self._getValueInPattern(pat, k, self.dicoValues)
                                ch = self._removeQuotes(str(ch))
                                self.dicoValues[k] = ch
                            else:
                                # If there are more than one occurence of the keyword in the
                                # batch script file, only the first one is modified
                                #
                                pass

        self.initializeBatchScriptFile()


    def initializeBatchScriptFile(self):
        """
        Initialize the backup file from reading dictionary self.dicoValues the first time.
        """
        # Basic verification
        #
        #node_ecs = self.case.xmlGetNode('solution_domain')
        #if not Tool.GuiParam.matisse:
            #if not node_ecs:
                #raise ValueError, "No preprocessor heading!"

        sdm = SolutionDomainModel(self.case)
        prm = ProfilesModel(self.case)

        # MESH
        meshes = string.join(sdm.getMeshList())
        if meshes:
            self.dicoValues['MESH'] = meshes

        self.dicoValues['COMMAND_REORIENT'] = sdm.getReorientCommand()
        self.dicoValues['COMMAND_JOIN'] = sdm.getJoinCommand()
        self.dicoValues['COMMAND_CWF'] = sdm.getCutCommand()
        self.dicoValues['COMMAND_PERIO'] = sdm.getPerioCommand()
        self.dicoValues['PARAM'] = os.path.basename(self.case['xmlfile'])

        # User 1D profiles are loaded as user result files

        if prm.getProfilesLabelsList():
            if self.dicoValues['USER_OUTPUT_FILES']:
                v = string.split(self.dicoValues['USER_OUTPUT_FILES'])
            else:
                v = []
            vlist = prm.updateOutputFiles(v)
            self.dicoValues['USER_OUTPUT_FILES'] = string.join(vlist, " ")

        # Specific data file for specific physics

        model = CoalCombustionModel(self.case).getCoalCombustionModel()
        if model == 'coal_homo' or model == 'coal_homo2':
            self.dicoValues['THERMOCHEMISTRY_DATA'] = 'dp_FCP'

        atmo = AtmosphericFlowsModel(self.case)
        if atmo.getAtmosphericFlowsModel() != 'off':
            self.dicoValues['METEO_DATA'] = atmo.getMeteoDataFileName()


    def updateBatchScriptFile(self, keyword=None):
        """
        Update the backup file from reading dictionary self.dicoValues.
        If keyword == None, all keywords are updated
        If keyword == key, only key is updated.
        """
        # update the name of the param, useful when the xml file is new
        # and was never saved before
        self.dicoValues['PARAM'] = os.path.basename(self.case['xmlfile'])

        l = self.dicoValues.keys()
        l.append(None) # Add 'None' when no keyword is specified in argument.
        for k in self.dicoValues.keys():
            if self.dicoValues[k] == 'None':
                self.dicoValues[k] = None

        lines = self.lines

        if keyword == None or keyword[0:4] == 'job_':
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

        else:
            self.isInList(keyword, l)

        if self.case['batch_type'] != "" \
          or str(self.dicoValues['NUMBER_OF_PROCESSORS']) == "1":
            self.dicoValues['NUMBER_OF_PROCESSORS'] = ""

        for k in self.dicoValues.keys():
            if k[0:4] == 'job_':
                continue
            nbkey = 0
            if keyword: k = keyword
            for i in range(len(lines)):
                if self._getRegex(k) != None:
                    pat = self._getLineToModify(self._getRegex(k),lines[i])
                    if pat != None:
                        nbkey = nbkey + 1
                        if nbkey == 1:
                            ch = self._addQuotes(str(self.dicoValues[k]))
                            new = k + "=" + ch
                            lines[i] = self._substituteLine(pat, new, lines[i])
            if keyword: break


        f = open(self.script1, 'w')
        f.writelines(lines)
        f.close()
        os.system('chmod +x ' + self.script1)


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
        XMLinit(self.case)

        domain = SolutionDomainModel(self.case)
        domain.addMesh('mail1.des', 'des')
        domain.addMesh('mail2.des', 'des')
        domain.addMesh('mail3.des', 'des')
        domain.setCutStatus('on')
        domain.setCutAngle(0.0321)
        domain.setOrientation('on')

        self.case['xmlfile'] = 'NEW.xml'
        self.case['scripts_path'] = os.getcwd()
        self.case['backupBatchScript'] = ''
        lance_test = '# test \n'\
        '#SOLCOM=6\n'\
        'SOLCOM=999\n'\
        'PARAM=NEW.xml\n'\
        'NUMBER_OF_PROCESSORS=2\n'\
        'PROCESSOR_LIST=\n'\
        'PARTITION_LIST=\n'\
        'USER_INPUT_FILES=data\n'\
        'USER_OUTPUT_FILES=titi\n'\
        'CS_TMP_PREFIX=/home/toto\n'\
        'MESH=\n'\
        'COMMAND_REORIENT=\n'\
        'COMMAND_JOIN=" --j  --color 98 99  --fraction 0.1  --plane 0.8"\n'\
        'COMMAND_CWF="--cwf 0.001"\n'\
        'COMMAND_PERIO=\n'\
        'EXEC_PREPROCESS=yes\n'\
        'EXEC_PARTITION=yes\n'\
        'EXEC_KERNEL=yes\n'\
        'CS_LIB_ADD=''\n'\
        'VALGRIND=''\n'\
        'ARG_CS_OUTPUT=''\n'\
        'ARG_CS_VERIF=''\n'

        self.f = open('lance_test','w')
        self.f.write(lance_test)
        self.f.close()


    def tearDown(self):
        """
        This method is executed after all 'check' methods.
        """
        f = self.case['batchScript']
        if os.path.isfile(f): os.remove(f)
        if os.path.isfile(f+"~"): os.remove(f+"~")


    def checkGetRegexAndGetLineToModify(self):
        """ Check whether the BatchRunningModel class could be get line"""
        mdl = BatchRunningModel(self.case)
        txt1 = '# fic = 1  '
        txt2 = '      fic=2'
        txt3 = 'fic=33'
        txt4 = '   fic =55'
        txt5 = '   fic = 55'
        txt6 = '   fic = " fic jklm    " '
        regex1 = mdl._getRegex('fic')
        regex2 = mdl._getRegex('fic')
        regex3 = mdl._getRegex('fic')
        regex4 = mdl._getRegex('fic')
        regex5 = mdl._getRegex('fic')
        regex6 = mdl._getRegex('fic')

        pat1 = mdl._getLineToModify(regex1,txt1)
        pat2 = mdl._getLineToModify(regex2,txt2)
        pat3 = mdl._getLineToModify(regex3,txt3)
        pat4 = mdl._getLineToModify(regex4,txt4)
        pat5 = mdl._getLineToModify(regex5,txt5)
        pat6 = mdl._getLineToModify(regex6,txt6)

        assert pat1 == None, 'Could not get pat1 to modify text'
        assert pat2.group() == '      fic=2', 'Could not get pat2 to modify text'
        assert pat3.group() == 'fic=33', 'Could not get pat3 to modify text'
        assert pat4.group() == '   fic =55', 'Could not get pat4 to modify text'
        assert pat5.group() == '   fic = 55', 'Could not get pat5 to modify text'
        assert pat6.group() == '   fic = " fic jklm    " ', 'Could not get pat6 to modify text'


    def checkGetValueInPattern(self):
        """ Check whether the class could be get value from regular expression"""
        mdl = BatchRunningModel(self.case)
        dico = {}
        txt = 'fic=33'
        txt1 = '# fic = 1  '
        txt2 = '      fic=2'
        txt5 = '   fic = 55'
        regex = mdl._getRegex('fic')
        pat = mdl._getLineToModify(regex,txt)
        value = mdl._getValueInPattern(pat, 'fic', dico)
        regex1 = mdl._getRegex('fic')
        pat1 = mdl._getLineToModify(regex1,txt1)
        regex2 = mdl._getRegex('fic')
        pat2 = mdl._getLineToModify(regex2,txt2)
        value2 = mdl._getValueInPattern(pat2, 'fic', dico)
        regex5 = mdl._getRegex('fic')
        pat5 = mdl._getLineToModify(regex5,txt5)
        value5 = mdl._getValueInPattern(pat5, 'fic', dico)

        assert value == '33','could not get value from regular expression'
        assert pat1 == None,'could not get value1 from regular expression'
        assert value2 == '2','could not get value2 from regular expression'
        assert value5 == '55','could not get value5 from regular expression'


    def checkSubstituteLine(self):
        """ Check whether the BatchRunningModel class could be substitute line"""
        mdl = BatchRunningModel(self.case)
        txt1 = '      fic='
        txt2 = ' fic= rien'
        pat1 = mdl._getLineToModify(mdl._getRegex('fic'),txt1)
        new1 = mdl._substituteLine(pat1,'vacances',txt1)
        pat2 = mdl._getLineToModify(mdl._getRegex('fic'),txt2)
        new_pat2 = 'fic=' + 'vacances'
        new2 = mdl._substituteLine(pat2,new_pat2,txt2)

        assert new1 == 'vacances','could not substitute line from regular expression'
        assert new2 == 'fic=vacances','could not substitute line from regular expression'


    def checkReadBatchScriptFile(self):
        """ Check whether the BatchRunningModel class could be read file"""
        mdl = BatchRunningModel(self.case)
        mdl.readBatchScriptFile()

        # The following keywords from the batch script file
        # are cancelled by the informations from the case !
        #   MESH
        #   COMMAND_JOIN
        #   COMMAND_CWF
        #   COMMAND_PERIO
        #   COMMAND_REORIENT
        #
        dico = {\
        'PROCESSOR_LIST': '',
        'PARTITION_LIST': '',
        'MESH': 'mail1.des mail2.des mail3.des',
        'COMMAND_CWF': ' --cwf 0.0321',
        'SOLCOM': '999',
        'USER_OUTPUT_FILES': 'titi',
        'PARAM': 'NEW.xml',
        'NUMBER_OF_PROCESSORS': '2',
        'USER_INPUT_FILES': 'data',
        'COMMAND_JOIN': '',
        'COMMAND_REORIENT': ' --reorient ',
        'CS_TMP_PREFIX': '/home/toto',
        'COMMAND_PERIO': '',
        'EXEC_PREPROCESS':'yes',
        'EXEC_PARTITION':'yes',
        'EXEC_KERNEL':'yes',
        'CS_LIB_ADD':'',
        'VALGRIND':'',
        'ARG_CS_OUTPUT':'',
        'ARG_CS_VERIF':'',
        'THERMOCHEMISTRY_DATA':'',
        'METEO_DATA':''}
        for k in mdl.dicoValues.keys():
            if mdl.dicoValues[k] != dico[k]:
                print "\nwarning for key: ", k
                print "  read value in the script:", mdl.dicoValues[k]
                print "  reference value:", dico[k]
            assert  mdl.dicoValues[k] == dico[k], 'could not read the batch script file'
        assert  mdl.dicoValues == dico, 'could not read batch script file'


    def checkUpdateBatchScriptFile(self):
        """ Check whether the BatchRunningModel class could update file"""
        mdl = BatchRunningModel(self.case)
        mdl.readBatchScriptFile()
        mdl.dicoValues['NUMBER_OF_PROCESSORS']='48'
        dico_updated = mdl.dicoValues
        mdl.updateBatchScriptFile()
        mdl.readBatchScriptFile()
        dico_read = mdl.dicoValues

        assert dico_updated == dico_read, 'error on updating batch script file'


def suite():
    testSuite = unittest.makeSuite(BatchRunningModelTestCase, "check")
    return testSuite


def runTest():
    print "BatchRunningModelTestCase"
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End of BacthRunningModel
#-------------------------------------------------------------------------------
