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

"""
This module modify the batch file
- BatchRunningModel
"""
#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

from __future__ import print_function

import sys, unittest
import os, os.path, shutil, sys, types, re

import cs_batch

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Variables, Model

import cs_exec_environment

#-------------------------------------------------------------------------------
# Class BatchRunningModel
#-------------------------------------------------------------------------------

class BatchRunningModel(Model):
    """
    This class modifies the batch file (runcase)
    """
    def __init__(self, parent, case):
        """
        Constructor.
        """
        self.parent = parent
        self.case = case

        if not self.case['runcase']:
            self.case['runcase'] = None

        self.batch = cs_batch.batch(self.case['package'])

        self.dictValues = {}

        # Do we force a number of MPI ranks ?

        self.dictValues['run_nprocs'] = None
        self.dictValues['run_nthreads'] = None
        self.dictValues['run_id'] = None
        self.dictValues['run_build'] = None
        self.dictValues['run_stage_init'] = None

        # Is a batch file present ?

        if self.case['runcase']:
            self.parseBatchFile()


    def parseBatchRunOptions(self):
        """
        Get info from the run command
        """

        self.dictValues['run_nprocs'] = self.case['runcase'].get_nprocs()
        self.dictValues['run_nthreads'] = self.case['runcase'].get_nthreads()
        self.dictValues['run_id'] = self.case['runcase'].get_run_id()[0]
        self.dictValues['run_build'] = self.case['runcase'].get_compute_build()
        self.dictValues['run_stage_init'] = self.case['runcase'].get_run_stage('initialize')


    def updateBatchRunOptions(self, keyword=None):
        """
        Update the run command
        """

        if (keyword == 'run_nprocs' or not keyword) and self.case['runcase']:
            self.case['runcase'].set_nprocs(self.dictValues['run_nprocs'])
        if (keyword == 'run_nthreads' or not keyword) and self.case['runcase']:
            self.case['runcase'].set_nthreads(self.dictValues['run_nthreads'])
        if (keyword == 'run_id' or not keyword) and self.case['runcase']:
            self.case['runcase'].set_run_id(run_id = self.dictValues['run_id'])
        if (keyword == 'run_build' or not keyword) and self.case['runcase']:
            self.case['runcase'].set_compute_build(self.dictValues['run_build'])
        if (keyword == 'run_stage_init' or not keyword) and self.case['runcase']:
            self.case['runcase'].set_run_stage('initialize',
                                               self.dictValues['run_stage_init'])


    def parseBatchFile(self):
        """
        Fill self.dictValues reading the batch file.
        """

        # Parse lines depending on batch type

        self.parseBatchRunOptions()

        if self.batch.rm_type == None:
            return

        self.batch.parse_lines(self.case['runcase'].lines)


    def updateBatchFile(self, keyword=None):
        """
        Update the batch file from reading dictionary self.dictValues.
        If keyword == None, all keywords are updated
        If keyword == key, only key is updated.
        """

        k = str(keyword)[0:3]

        if k[0:3] == 'run':
            l = list(self.dictValues.keys())
            l.append(None) # Add 'None' when no keyword is specified in argument.
            for k in list(self.dictValues.keys()):
                if self.dictValues[k] == 'None':
                    self.dictValues[k] = None
            self.isInList(keyword, l)

        if k[0:3] != 'job':
            self.updateBatchRunOptions()

        if k[0:3] != 'run':
            self.batch.update_lines(self.case['runcase'].lines, keyword)


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
        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit
        from code_saturne.model.Common import GuiParam
        GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case).initialize()

        self.case['scripts_path'] = os.getcwd()
        self.case['runcase'] = cs_runcase.runcase('runcase')

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
        f = self.case['runcase'].path
        if os.path.isfile(f): os.remove(f)


    def checkReadBatchPBS(self):
        """ Check whether the BatchRunningModel class could be read file"""
        # TODO update: self.case['batch_type'] = 'PBS'
        mdl = BatchRunningModel(self.case)

        dico_PBS = {\
        'job_nodes': '16',
        'job_name': 'super_toto',
        'job_ppn': '1',
        'job_walltime': '34:77:22'}

        for k in list(dico_PBS.keys()):
            if mdl.batch.params[k] != dico_PBS[k] :
                print("\nwarning for key: ", k)
                print("  read value in the batch description: ", mdl.batch.params[k])
                print("  reference value: ", dico_PBS[k])
            assert  mdl.batch.params[k] == dico_PBS[k], 'could not read the batch file'


    def checkUpdateBatchFile(self):
        """ Check whether the BatchRunningModel class could update file"""
        mdl = BatchRunningModel(self.case)
        mdl.batch.params['job_procs']=48
        dico_updated = mdl.batch.params
        mdl.updateBatchFile()
        dico_read = mdl.batch.params

        assert dico_updated == dico_read, 'error on updating batch script file'


    def checkUpdateBatchPBS(self):
        """ Check whether the BatchRunningModel class could update file"""
        mdl = BatchRunningModel(self.case)
        mdl.batch.params['job_walltime']='12:42:52'
        dicojob_updated = mdl.batch.params
        mdl.updateBatchFile()
        dicojob_read = mdl.batch.params

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
