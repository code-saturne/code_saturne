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
        self.dicoValues['SOLCOM'] = '0'
        self.dicoValues['PARAM'] = ""
        self.dicoValues['VERSION'] = ""
        self.dicoValues['NUMBER_OF_PROCESSORS'] = '1'
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
        self.dicoValues['COMMAND_SYRTHES'] = ""
        self.dicoValues['PBS_JOB_NAME'] = ""
        self.dicoValues['PBS_nodes'] = '1'
        self.dicoValues['PBS_ppn'] = '2'
        self.dicoValues['PBS_walltime'] = '1:00:00'
        self.dicoValues['PBS_mem'] = '320'
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
            self.mdl.dicoValues['ARG_CS_OUTPUT'] = "--log 0"


    def _getRegex(self, word):
        """
        Get regular expression to extract line without comment
        """
##  fonctionne mais incomplet:      regex = re.compile(r"""(^\s*""" + word + r""".*$)""")
##  fonctionne en tenant compte des lignes commencant par # :
##      regex = re.compile(r"""(^(?#)^\s*""" + word + r""".*$)""")
        #tient compte a  la fois des commentaires et des "$word":
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
                                ch = self._removeQuotes(str(ch))
                                self.dicoValues[k] = ch
                            else:
                                # If there are more than one occurence of the keyword in the
                                # batch script file, only the first one is modified
                                #
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
        meshes = string.join(sdm.getMeshList())
        if meshes:
            self.dicoValues['MESH'] = meshes

        self.dicoValues['COMMAND_REORIENT'] = sdm.getReorientCommand()
        self.dicoValues['COMMAND_JOIN'] = sdm.getJoinCommand()
        self.dicoValues['COMMAND_CWF'] = sdm.getCutCommand()
        self.dicoValues['COMMAND_PERIO'] = sdm.getPerioCommand()
        self.dicoValues['COMMAND_SYRTHES'] = sdm.getSyrthesCommand()
        self.dicoValues['PARAM'] = os.path.basename(self.case['xmlfile'])

        list = prm.getProfilesLabelsList()
        if list:
            if self.dicoValues['USER_OUTPUT_FILES']:
                vlist = string.split(self.dicoValues['USER_OUTPUT_FILES'])
            else:
                vlist = []
            for file in list:
                if file not in vlist:
                    vlist.append(file)
            self.dicoValues['USER_OUTPUT_FILES'] = string.join(vlist, " ")

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
        l = self.dicoValues.keys()
        l.append(None)
        self.isInList(keyword, l)
        lines = self.lines

        if self.case['computer'] == "pbs":
            self.dicoValues['NUMBER_OF_PROCESSORS'] = ""

        for k in self.dicoValues.keys():
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
        lance_test = '# test \n'\
        '#SOLCOM=6\n'\
        'SOLCOM=999\n'\
        'PARAM=NEW.xml\n'\
        'VERSION=tutu\n'\
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
        'COMMAND_SYRTHES=\n'\
        'EXEC_PREPROCESS=yes\n'\
        'EXEC_PARTITION=yes\n'\
        'EXEC_KERNEL=yes\n'\
        'CS_LIB_ADD=''\n'\
        'VALGRIND=''\n'\
        'ARG_CS_OUTPUT=''\n'\
        'ARG_CS_VERIF=''\n'

        lance_PBS = '# test \n'\
        '#\n'\
        '#                  CARTES BATCH POUR CLUSTERS sous PBS\n'\
        '#\n'\
        '#PBS -l nodes=16:ppn=1,walltime=34:77:22,mem=832mb\n'\
        '#PBS -j eo -N super_toto\n'

        lance_LSF = '# test \n'\
        '#\n'\
        '#        CARTES BATCH POUR LE CCRT (Nickel/Chrome/Tantale sous LSF)\n'\
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
        #   MESH
        #   COMMAND_JOIN
        #   COMMAND_CWF
        #   COMMAND_SYRTHES
        #   COMMAND_PERIO
        #   COMMAND_REORIENT
        #
        dico = {\
        'PROCESSOR_LIST': '',
        'PARTITION_LIST': '',
        'PBS_nodes': '1',
        'MESH': 'mail1.des mail2.des mail3.des',
        'PBS_JOB_NAME': '',
        'COMMAND_CWF': ' --cwf 0.0321',
        'SOLCOM': '999',
        'USER_OUTPUT_FILES': 'titi',
        'PARAM': 'NEW.xml',
        'NUMBER_OF_PROCESSORS': '2',
        'USER_INPUT_FILES': 'data',
        'COMMAND_JOIN': '',
        'COMMAND_REORIENT': ' --reorient ',
        'VERSION': 'tutu',
        'CS_TMP_PREFIX': '/home/toto',
        'COMMAND_SYRTHES': '',
        'PBS_ppn': '2',
        'PBS_walltime': '1:00:00',
        'PBS_mem': '320',
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


    def checkReadBatchScriptPBS(self):
        """ Check whether the BatchRunningModel class could be read file"""
        self.case['computer'] = 'pbs'
        mdl = BatchRunningModel(self.case)
        mdl.readBatchScriptFile()

        # The following keywords from the batch script file
        # are cancelled by the informations from the case !
        #   MESH
        #   COMMAND_JOIN
        #   COMMAND_CWF
        #   COMMAND_SYRTHES
        #   COMMAND_PERIO
        #
        #
        dico_PBS = {\
        'PBS_nodes': '16',
        'MESH': 'mail1.des mail2.des mail3.des',
        'PBS_JOB_NAME': 'super_toto',
        'COMMAND_CWF': ' --cwf 0.0321',
        'SOLCOM': '0',
        'PARAM': 'NEW.xml',
        'NUMBER_OF_PROCESSORS': '1',
        'USER_INPUT_FILES': '',
        'COMMAND_JOIN': '',
        'COMMAND_REORIENT': ' --reorient ',
        'CS_TMP_PREFIX': '',
        'COMMAND_SYRTHES': '',
        'PBS_ppn': '1',
        'PBS_walltime': '34:77:22',
        'PBS_mem': '832',
        'COMMAND_PERIO': ''}

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
        mdl.dicoValues['NUMBER_OF_PROCESSORS']='48'
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
