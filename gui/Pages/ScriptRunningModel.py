# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
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

"""
This module modifies the "runcase" script file
- ScriptModel
- ScriptRunningModel
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
# Class ScriptRunningModel
#-------------------------------------------------------------------------------

class ScriptModel:
    """
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        if not self.case['script']:
            self.case['script'] = ""

        if not self.case['backupScript']:
            self.case['backupScript'] = False


class ScriptRunningModel(Model):
    """
    This class modifies the script file (runcase)
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        ScriptModel(self.case)

        self.script1 = self.case['scripts_path'] + "/" + self.case['script']

        # Save a backup file before any modification.
        # There is only one backup for the entire session.

        script2 = self.script1 + "~"

        if self.case['backupScript'] == False:
            shutil.copy2(self.script1, script2)
            self.case['backupScript'] = True

        # Read the script file line by line.
        # All lines are stored in a list called "self.script_lines".

        f = open(self.script1, 'rw')
        self.script_lines = f.readlines()
        f.close()

        # DicoValues's initialisation

        self.dictValues = {}
        self.dictValues['MESHES'] = None
        self.dictValues['REORIENT'] = False
        self.dictValues['THERMOCHEMISTRY_DATA'] = None
        self.dictValues['METEO_DATA'] = None
        self.dictValues['USER_INPUT_FILES'] = None
        self.dictValues['N_PROCS'] = None
        self.dictValues['PARTITION_LIST'] = None
        self.dictValues['PARAMETERS'] = None
        self.dictValues['CHECK_ARGS'] = None
        self.dictValues['OUTPUT_ARGS'] = None
        self.dictValues['HOSTS_LIST'] = None
        self.dictValues['EXEC_PREPROCESS'] = True
        self.dictValues['EXEC_PARTITION'] = True
        self.dictValues['EXEC_SOLVER'] = True

        self.dictValues['CS_TMP_PREFIX'] = None
        self.dictValues['VALGRIND'] = None

        if self.case['salome']:
            self.dictValues['OUTPUT_ARGS'] = "--log 0"


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
        self.dictValues[word] = dico
        resu = pattern.group().split('=')
        L = resu[1]
        for i in range(2,len(resu)): L = L  + "="+ resu[i]
        resu = L

        try:
            self.dictValues[word] = eval(resu)
        except Exception:
            self.dictValues[word] = None

        return self.dictValues[word]


    def _addQuotes(self, ch):
        """
        Add quotes in front of and behind string c.
        """

        if ch[0] != "'" and ch[0] != '"' and ch != 'True' and ch != 'False':
            ch = "'" + ch + "'"

        return ch


    def readScriptFile(self):
        """
        Fill self.dictValues reading the backup file.
        """
        script_lines = self.script_lines

        for k in self.dictValues.keys():
            if k not in ('THERMOCHEMISTRY_DATA', 'METEO_DATA'):
                nbkey = 0
                for i in range(len(script_lines)):
                    reg = self._getRegex(k)
                    if reg != None:
                        pat = self._getLineToModify(reg, script_lines[i])
                        if pat != None:
                            nbkey = nbkey + 1
                            if nbkey == 1:
                                ch = self._getValueInPattern(pat, k, self.dictValues)
                            else:
                                # If there are more than one occurence of the
                                # keyword, only the first one is modified
                                pass

        self.initializeScriptFile()


    def initializeScriptFile(self):
        """
        Initialize the backup file from reading dictionary self.dictValues the first time.
        """

        sdm = SolutionDomainModel(self.case)
        mm = MeshModel()

        # MESHES
        self.dictValues['MESHES'] = []
        for m in sdm.getMeshList():
            l_args = []
            extension = mm.getMeshExtension(m)
            if not extension in mm.getExtensionFileList():
                l_args.append('--format ' + sdm.getMeshFormat(m))
            n = sdm.getMeshNumber(m)
            if n != None:
                l_args.append('--num ' + str(n))
            gc = sdm.getMeshGroupCells(m)
            if gc != 'off':
                l_args.append('--grp-cel ' + gc)
            gf = sdm.getMeshGroupFaces(m)
            if gf != 'off':
                l_args.append('--grp-fac ' + gf)
            if len(l_args) >  0:
                l_args.insert(0, m)
                self.dictValues['MESHES'].append(tuple(l_args))
            else:
                self.dictValues['MESHES'].append(m)

        self.dictValues['REORIENT'] = sdm.getReorientCommand()
        self.dictValues['PARAMETERS'] = os.path.basename(self.case['xmlfile'])

        # Specific data file for specific physics

        model = CoalCombustionModel(self.case).getCoalCombustionModel()
        if model == 'coal_homo' or model == 'coal_homo2':
            self.dictValues['THERMOCHEMISTRY_DATA'] = 'dp_FCP'

        atmo = AtmosphericFlowsModel(self.case)
        if atmo.getAtmosphericFlowsModel() != 'off':
            self.dictValues['METEO_DATA'] = atmo.getMeteoDataFileName()


    def updateScriptFile(self, keyword=None):
        """
        Update the backup file from reading dictionary self.dictValues.
        If keyword == None, all keywords are updated
        If keyword == key, only key is updated.
        """
        # update the name of the param, useful when the xml file is new
        # and was never saved before
        self.dictValues['PARAMETERS'] = os.path.basename(self.case['xmlfile'])

        l = self.dictValues.keys()
        l.append(None) # Add 'None' when no keyword is specified in argument.
        for k in self.dictValues.keys():
            if self.dictValues[k] == 'None':
                self.dictValues[k] = None
        self.isInList(keyword, l)
        script_lines = self.script_lines

        if self.dictValues['N_PROCS'] == 1:
            self.dictValues['N_PROCS'] = None

        for k in self.dictValues.keys():
            nbkey = 0
            if keyword: k = keyword
            for i in range(len(script_lines)):
                if self._getRegex(k) != None:
                    pat = self._getLineToModify(self._getRegex(k),script_lines[i])
                    if pat != None:
                        nbkey = nbkey + 1
                        if nbkey == 1:
                            if isinstance(self.dictValues[k], str):
                                ch = self.dictValues[k]
                                if len(ch) > 0:
                                    ch = self._addQuotes(self.dictValues[k])
                                else:
                                    ch = 'None'
                            else:
                                ch = str(self.dictValues[k])
                                if ch == '[]':
                                    ch = 'None'
                            new = k + " = " + ch
                            script_lines[i] = self._substituteLine(pat, new, script_lines[i])
            if keyword: break

        f = open(self.script1, 'w')
        f.writelines(script_lines)
        f.close()
        os.system('chmod +x ' + self.script1)


#-------------------------------------------------------------------------------
# ScriptRunningModel test class
#-------------------------------------------------------------------------------

class ScriptRunningModelTestCase(unittest.TestCase):
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
        domain.setOrientation('on')

        self.case['xmlfile'] = 'NEW.xml'
        self.case['scripts_path'] = os.getcwd()
        self.case['script'] = 'runcase_test'
        self.case['backupScript'] = True
        self.case['script'] = 'runcase'
        runcase_test = "# test \n"\
        "PARAMETERS='NEW.xml'\n"\
        "N_PROCS=2\n"\
        "HOSTS_LIST=None\n"\
        "PARTITION_LIST=None\n"\
        "USER_INPUT_FILES=['data']\n"\
        "CS_TMP_PREFIX='/home/toto'\n"\
        "MESHES=None\n"\
        "REORIENT=False\n"\
        "EXEC_PREPROCESS=True\n"\
        "EXEC_PARTITION=True\n"\
        "EXEC_SOLVER=True\n"\
        "VALGRIND=None\n"\
        "OUTPUT_ARGS=None\n"\
        "CHECK_ARGS=None\n"

        lance_PBS = '# test \n'\
        '#\n'\
        '#                  CARTES SCRIPT POUR CLUSTERS sous PBS\n'\
        '#\n'\
        '#PBS -l nodes=16:ppn=1,walltime=34:77:22\n'\
        '#PBS -j eo -N super_toto\n'

        lance_LSF = '# test \n'\
        '#\n'\
        '#        CARTES SCRIPT POUR LE CCRT (Platine sous LSF)\n'\
        '#\n'\
        '#BSUB -n 2\n'\
        '#BSUB -c 00:05\n'\
        '#BSUB -o super_tataco.%J\n'\
        '#BSUB -e super_tatace.%J\n'\
        '#BSUB -J super_truc\n'

        self.f = open('runcase_test','w')
        self.f.write(runcase_test)
        self.f.close()


    def checkGetRegexAndGetLineToModify(self):
        """ Check whether the ScriptRunningModel class could be get line"""
        mdl = ScriptRunningModel(self.case)
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
        mdl = ScriptRunningModel(self.case)
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
        """ Check whether the ScriptRunningModel class could be substitute line"""
        mdl = ScriptRunningModel(self.case)
        txt1 = '      fic='
        txt2 = ' fic= rien'
        pat1 = mdl._getLineToModify(mdl._getRegex('fic'),txt1)
        new1 = mdl._substituteLine(pat1,'vacances',txt1)
        pat2 = mdl._getLineToModify(mdl._getRegex('fic'),txt2)
        new_pat2 = 'fic=' + 'vacances'
        new2 = mdl._substituteLine(pat2,new_pat2,txt2)

        assert new1 == 'vacances','could not substitute line from regular expression'
        assert new2 == 'fic=vacances','could not substitute line from regular expression'


    def checkReadScriptScriptFile(self):
        """ Check whether the ScriptRunningModel class could be read file"""
        mdl = ScriptRunningModel(self.case)
        mdl.readScriptScriptFile()

        # The following keywords from the script script file
        # are cancelled by the informations from the case !
        #   MESHES
        #   REORIENT
        #
        dico = {\
        'HOSTS_LIST': '',
        'PARTITION_LIST': None,
        'MESHES': ['mail1.des', 'mail2.des', 'mail3.des'],
        'PARAMETERS': 'NEW.xml',
        'N_PROCS': 2,
        'USER_INPUT_FILES': ['data'],
        'REORIENT': True,
        'CS_TMP_PREFIX': '/home/toto',
        'EXEC_PREPROCESS':True,
        'EXEC_PARTITION':True,
        'EXEC_SOLVER':True,
        'VALGRIND':None,
        'OUTPUT_ARGS':None,
        'CHECK_ARGS':None,
        'THERMOCHEMISTRY_DATA':None,
        'METEO_DATA':None}
        for k in mdl.dictValues.keys():
            if mdl.dictValues[k] != dico[k]:
                print("\nwarning for key: ", k)
                print("  read value in the script:", mdl.dictValues[k])
                print("  reference value:", dico[k])
            assert  mdl.dictValues[k] == dico[k], 'could not read the script script file'
        assert  mdl.dictValues == dico, 'could not read script script file'


    def checkUpdateScriptScriptFile(self):
        """ Check whether the ScriptRunningModel class could update file"""
        mdl = ScriptRunningModel(self.case)
        mdl.readScriptScriptFile()
        mdl.dictValues['N_PROCS']=48
        dico_updated = mdl.dictValues
        mdl.updateScriptScriptFile()
        mdl.readScriptScriptFile()
        dico_read = mdl.dictValues

        assert dico_updated == dico_read, 'error on updating script script file'


def suite():
    testSuite = unittest.makeSuite(ScriptRunningModelTestCase, "check")
    return testSuite


def runTest():
    print("ScriptRunningModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End of ScriptRunningModel
#-------------------------------------------------------------------------------
