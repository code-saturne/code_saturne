# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
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
            self.case['batchScript'] = {'station': "",
                                        'pbs': "",
                                        'lsf': "",
                                        'sge': ""}

        if not self.case['backupBatchScript']:
            self.case['backupBatchScript'] = {'station': "no",
                                              'pbs': "no",
                                              'lsf': "no",
                                              'sge': "no"}


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

        # we get up batch script file
        key = self.case['computer']

        self.script1 = self.case['scripts_path'] + "/" + self.case['batchScript'][key]

        # Save a backup file before any modification.
        # There is only one backup for the entire session.

        script2 = self.script1 + "~"

        if self.case['backupBatchScript'][key] == "no" \
               or not os.path.isfile(script2):
            os.popen('cp ' + self.script1 + " " +script2)
            self.case['backupBatchScript'][key] = "yes"

        # Read the batch script file line by line.
        # All lines are stored in a list called "self.lines".

        f = open(self.script1, 'rw')
        self.lines = f.readlines()
        f.close()

        # DicoValues's initialisation

        self.dicoValues = {}
        self.dicoValues['MESHES'] = None
        self.dicoValues['REORIENT'] = False
        self.dicoValues['JOIN'] = None
        self.dicoValues['PERIODICITY'] = None
        self.dicoValues['THERMOCHEMISTRY_DATA'] = None
        self.dicoValues['METEO_DATA'] = None
        self.dicoValues['USER_INPUT_FILES'] = None
        self.dicoValues['USER_OUTPUT_FILES'] = None
        self.dicoValues['N_PROCS'] = None
        self.dicoValues['PARTITION_LIST'] = None
        self.dicoValues['PARAMETERS'] = None
        self.dicoValues['SOLCOM'] = False
        self.dicoValues['CUT_WARPED_FACES'] = None
        self.dicoValues['CHECK_ARGS'] = None
        self.dicoValues['OUTPUT_ARGS'] = None
        self.dicoValues['HOSTS_LIST'] = None
        self.dicoValues['EXEC_PREPROCESS'] = True
        self.dicoValues['EXEC_PARTITION'] = True
        self.dicoValues['EXEC_SOLVER'] = True

        self.dicoValues['CS_TMP_PREFIX'] = None
        self.dicoValues['PBS_JOB_NAME'] = ""
        self.dicoValues['PBS_nodes'] = '1'
        self.dicoValues['PBS_ppn'] = '2'
        self.dicoValues['PBS_walltime'] = '1:00:00'
        self.dicoValues['PBS_mem'] = '320'
        self.dicoValues['VALGRIND'] = None

        if self.case['salome']:
            self.dicoValues['OUTPUT_ARGS'] = "--log 0"


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

        try:
            self.dicoValues[word] = eval(resu)
        except Exception:
            self.dicoValues[word] = None

        return self.dicoValues[word]


    def _addQuotes(self, ch):
        """
        Add quotes in front of and behind string c.
        """

        if ch[0] != "'" and ch[0] != '"' and ch != 'True' and ch != 'False':
            ch = "'" + ch + "'"

        return ch


    def readBatchScriptFile(self):
        """
        Fill self.dicoValues reading the backup file.
        """
        lines = self.lines

        list = ['PBS_JOB_NAME','PBS_nodes','PBS_ppn','PBS_walltime','PBS_mem']

        for k in self.dicoValues.keys():
            if k not in list and k not in ('THERMOCHEMISTRY_DATA', 'METEO_DATA'):
                nbkey = 0
                for i in range(len(lines)):
                    reg = self._getRegex(k)
                    if reg != None:
                        pat = self._getLineToModify(reg, lines[i])
                        if pat != None:
                            nbkey = nbkey + 1
                            if nbkey == 1:
                                ch = self._getValueInPattern(pat, k, self.dicoValues)
                            else:
                                # If there are more than one occurence of the keyword in the
                                # batch script file, only the first one is modified
                                pass

        if self.case['computer'] == "pbs":
            for (word, ind, lab) in [('#PBS -j eo -N ', 0, 'PBS_JOB_NAME')]:
                for i in range(len(lines)):
                    index = string.rfind(lines[i], word)
                    if index == ind:
                        self.dicoValues[lab] = self._removeQuotes(lines[i], ind, word)

            for (word, next, lab) in [('#PBS -l nodes', ':ppn'    , 'PBS_nodes'),
                                      (':ppn'         ,',walltime', 'PBS_ppn'),
                                      (',walltime'    ,',mem'    , 'PBS_walltime'),
                                      (',mem'        , 'mb'     , 'PBS_mem')]:
                word = word + "="
                for i in range(len(lines)):
                    ind1 = string.rfind(lines[i], word)
                    ind2 = string.rfind(lines[i], next)
                    if ind1 != -1 and ind2 != -1:
                        ch = lines[i][ind1+len(word):ind2]
                        self.dicoValues[lab] = ch

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
        self.dicoValues['MESH'] = sdm.getMeshList()

        self.dicoValues['REORIENT'] = sdm.getReorientCommand()
        self.dicoValues['JOIN'] = sdm.getJoinCommand()
        self.dicoValues['CUT_WARPED_FACES'] = sdm.getCutCommand()
        self.dicoValues['PERIODICITY'] = sdm.getPerioCommand()
        self.dicoValues['PARAMETERS'] = os.path.basename(self.case['xmlfile'])

        # User 1D profiles are loaded as user result files

        if prm.getProfilesLabelsList():
            if self.dicoValues['USER_OUTPUT_FILES']:
                v = self.dicoValues['USER_OUTPUT_FILES']
                vlist = prm.updateOutputFiles(v)
            else:
                vlist = None
            self.dicoValues['USER_OUTPUT_FILES'] = vlist

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
        self.dicoValues['PARAMETERS'] = os.path.basename(self.case['xmlfile'])

        l = self.dicoValues.keys()
        l.append(None) # Add 'None' when no keyword is specified in argument.
        for k in self.dicoValues.keys():
            if self.dicoValues[k] == 'None':
                self.dicoValues[k] = None
        self.isInList(keyword, l)
        lines = self.lines

        if self.case['computer'] == "pbs" \
          or self.dicoValues['N_PROCS'] == 1:
            self.dicoValues['N_PROCS'] = None

        for k in self.dicoValues.keys():
            nbkey = 0
            if keyword: k = keyword
            for i in range(len(lines)):
                if self._getRegex(k) != None:
                    pat = self._getLineToModify(self._getRegex(k),lines[i])
                    if pat != None:
                        nbkey = nbkey + 1
                        if nbkey == 1:
                            if isinstance(self.dicoValues[k], str):
                                ch = self.dicoValues[k]
                                if len(ch) > 0:
                                    ch = self._addQuotes(self.dicoValues[k])
                                else:
                                    ch = 'None'
                            else:
                                ch = str(self.dicoValues[k])
                                if ch == '[]':
                                    ch = 'None'
                            new = k + " = " + ch
                            lines[i] = self._substituteLine(pat, new, lines[i])
            if keyword: break

        #  keywords only for the PBS Cluster
        if self.case['computer'] == "pbs":
            for (word, ind, var) in [('#PBS -j eo -N ', 0, self.dicoValues['PBS_JOB_NAME'])]:
                for i in range(len(lines)):
                    index = string.rfind(lines[i], word)
                    if index == ind:
                        if type(var) != types.StringType :
                            var = str(var)
                            if var == "0" : var=""
                        lines[i] = word + var + '\n'

            for (word, ind, next) in [('#PBS -l nodes', 0, ':ppn')]:
                for i in range(len(lines)):
                    ind1 = string.rfind(lines[i], word)
                    ind2 = string.rfind(lines[i], next)
                    if ind1 == ind and ind2 != -1:
                        ch = ""
                        for (w, dic) in [('#PBS -l nodes=', self.dicoValues['PBS_nodes']),
                                         (':ppn=',          self.dicoValues['PBS_ppn']),
                                         (',walltime=',     self.dicoValues['PBS_walltime']),
                                         (',mem=',          self.dicoValues['PBS_mem'])]:
                            ch = ch + w + str(dic)
                        lines[i] = ch + "mb\n"

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
        self.case['computer'] = 'station'
        self.case['scripts_path'] = os.getcwd()
        self.case['batchScript'] = {'pbs': 'lance_PBS', 'lsf': 'lance_LSF', 'station': 'lance_test'}
        self.case['backupBatchScript'] = {'pbs': 'yes', 'lsf': 'yes', 'station': 'yes'}
        lance_test = "# test \n"\
        "#SOLCOM=6\n"\
        "SOLCOM=False\n"\
        "PARAMETERS='NEW.xml'\n"\
        "N_PROCS=2\n"\
        "HOSTS_LIST=None\n"\
        "PARTITION_LIST=None\n"\
        "USER_INPUT_FILES=['data']\n"\
        "USER_OUTPUT_FILES=['titi']\n"\
        "CS_TMP_PREFIX='/home/toto'\n"\
        "MESHES=None\n"\
        "REORIENT=\n"\
        "JOIN='--j  --color 98 99  --fraction 0.1  --plane 0.8'\n"\
        "CUT_WARPED_FACES='--cwf 0.001'\n"\
        "PERIODICITY=\n"\
        "EXEC_PREPROCESS=True\n"\
        "EXEC_PARTITION=True\n"\
        "EXEC_SOLVER=True\n"\
        "VALGRIND=None\n"\
        "OUTPUT_ARGS=None\n"\
        "CHECK_ARGS=None\n"

        lance_PBS = '# test \n'\
        '#\n'\
        '#                  CARTES BATCH POUR CLUSTERS sous PBS\n'\
        '#\n'\
        '#PBS -l nodes=16:ppn=1,walltime=34:77:22,mem=832mb\n'\
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

        self.f = open('lance_test','w')
        self.f.write(lance_test)
        self.f.close()
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
        for plateform in ('station', 'pbs','lsf'):
            f = self.case['batchScript'][plateform]
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
        self.case['computer'] = 'station'
        mdl = BatchRunningModel(self.case)
        mdl.readBatchScriptFile()

        # The following keywords from the batch script file
        # are cancelled by the informations from the case !
        #   MESHES
        #   JOIN
        #   CUT_WARPED_FACES
        #   PERIODICITY
        #   REORIENT
        #
        dico = {\
        'HOSTS_LIST': '',
        'PARTITION_LIST': None,
        'PBS_nodes': '1',
        'MESHES': ['mail1.des', 'mail2.des', 'mail3.des'],
        'PBS_JOB_NAME': '',
        'CUT_WARPED_FACES': ' --cwf 0.0321',
        'SOLCOM': False,
        'USER_OUTPUT_FILES': ['titi'],
        'PARAMETERS': 'NEW.xml',
        'N_PROCS': 2,
        'USER_INPUT_FILES': ['data'],
        'JOIN': None,
        'REORIENT': ' --reorient ',
        'CS_TMP_PREFIX': '/home/toto',
        'PBS_ppn': '2',
        'PBS_walltime': '1:00:00',
        'PBS_mem': '320',
        'PERIODICITY': None,
        'EXEC_PREPROCESS':True,
        'EXEC_PARTITION':True,
        'EXEC_SOLVER':True,
        'VALGRIND':None,
        'OUTPUT_ARGS':None,
        'CHECK_ARGS':None,
        'THERMOCHEMISTRY_DATA':None,
        'METEO_DATA':None}
        for k in mdl.dicoValues.keys():
            if mdl.dicoValues[k] != dico[k]:
                print "\nwarning for key: ", k
                print "  read value in the script:", mdl.dicoValues[k]
                print "  reference value:", dico[k]
            assert  mdl.dicoValues[k] == dico[k], 'could not read the batch script file'
        assert  mdl.dicoValues == dico, 'could not read batch script file'


    def checkReadBatchScriptPBS(self):
        """ Check whether the BatchRunningModel class could be read file"""
        self.case['computer'] = 'pbs'
        mdl = BatchRunningModel(self.case)
        mdl.readBatchScriptFile()

        # The following keywords from the batch script file
        # are cancelled by the informations from the case !
        #   MESHES
        #   JOIN
        #   CUT_WARPED_FACES
        #   PERIODICITY
        #
        #
        dico_PBS = {\
        'PBS_nodes': '16',
        'MESHES': ['mail1.des', 'mail2.des', 'mail3.des'],
        'PBS_JOB_NAME': 'super_toto',
        'CUT_WARPED_FACES': ' --cwf 0.0321',
        'SOLCOM': False,
        'PARAMETERS': 'NEW.xml',
        'N_PROCS': None,
        'USER_INPUT_FILES': None,
        'JOIN': '',
        'REORIENT': ' --reorient ',
        'CS_TMP_PREFIX': '',
        'PBS_ppn': '1',
        'PBS_walltime': '34:77:22',
        'PBS_mem': '832',
        'PERIODICITY': ''}

        for k in dico_PBS.keys():
            if mdl.dicoValues[k] != dico_PBS[k] :
                print "\nwarning for key: ", k
                print "  read value in the script:", mdl.dicoValues[k]
                print "  reference value:", dico_PBS[k]
            assert  mdl.dicoValues[k] == dico_PBS[k], 'could not read the batch script file'


    def checkUpdateBatchScriptFile(self):
        """ Check whether the BatchRunningModel class could update file"""
        mdl = BatchRunningModel(self.case)
        mdl.readBatchScriptFile()
        mdl.dicoValues['N_PROCS']=48
        dico_updated = mdl.dicoValues
        mdl.updateBatchScriptFile()
        mdl.readBatchScriptFile()
        dico_read = mdl.dicoValues

        assert dico_updated == dico_read, 'error on updating batch script file'


    def checkUpdateBatchScriptPBS(self):
        """ Check whether the BatchRunningModel class could update file"""
        mdl = BatchRunningModel(self.case)
        mdl.readBatchScriptFile()
        mdl.dicoValues['PBS_mem']='512'
        mdl.dicoValues['PBS_walltime']='12:42:52'
        dicoPBS_updated = mdl.dicoValues
        mdl.updateBatchScriptFile()
        mdl.readBatchScriptFile()
        dicoPBS_read = mdl.dicoValues

        assert dicoPBS_updated == dicoPBS_read, 'error on updating PBS batch script file'


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
