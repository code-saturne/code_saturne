#============================================================================
#
#     This file is part of the Code_Saturne Kernel, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2008 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne Kernel is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne Kernel is distributed in the hope that it will be
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
#============================================================================
#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, os, shutil
import re, time

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

import Autovalidation.Parser as Parser
import Autovalidation.Common as Common
import Autovalidation.Listing as Listing
import Autovalidation.Chrono as Chrono

class Case :
    """
    Describe a Saturne case
    """
    
    def __init__(self, parser, studyLabel, caseLabel, variables, reportListing, reportChrono):
        """
        constructor
        """
        self.studyLabel = studyLabel
        self.caseLabel = caseLabel
        self.variables = variables
        self.parser = parser
        self.studyPath = Common.localDirectory+"/"+studyLabel.upper()
        self.casePath = self.studyPath+"/"+caseLabel.upper()
        self.reportListing = reportListing
        self.reportChrono = reportChrono
        #
        # if case doesn't exist, create "case" with cree_sat
        #
        if not os.path.isdir(self.studyPath+"/"+caseLabel.upper()):
            os.chdir(self.studyPath)
            proc = os.popen("cree_sat -noihm -cas "+caseLabel)
            proc.close()
            os.chdir(Common.localDirectory)
            #
            # Si le cas est deja cree l'indiquer dans le fichier report
            # Si il y a deja des resultats dans RESU, ils sont supprimes
            # Ajouter en parametre le fichier report
            #
            # ln or cp :
            #    - users
            #    - data.xml
            refFortPath = Common.referencePath+"/"+studyLabel.upper()+"/"+caseLabel.upper()+"/FORT"
            try:
                fortransList = os.listdir(refFortPath)
                for fortran in fortransList :
                    if fortran != "USERS" :
                        shutil.copyfile(refFortPath+"/"+fortran,self.studyPath+"/"+caseLabel.upper()+"/FORT/"+fortran)
            except:
                pass

            refDataPath = Common.referencePath+"/"+studyLabel.upper()+"/"+caseLabel.upper()+"/DATA"
            try:
                datasList = os.listdir(refDataPath)
                for data in datasList :
                    if data != "THCH" and data != "SaturneGUI" :
                        shutil.copyfile(refDataPath+"/"+data,self.studyPath+"/"+caseLabel.upper()+"/DATA/"+data)
            except:
                pass
            #
            # mise a jour du lance
            refScriptsPath = Common.referencePath+"/"+studyLabel.upper()+"/"+caseLabel.upper()+"/SCRIPTS"
            try:
                lanceFile = file(refScriptsPath+'/lance', mode='r')
            except IOError:
                print "Error : opening "+refScriptsPath+'/lance'
                print sys.exit(1)           
            
            keywordsLance = ["SOLCOM","LONGIA","LONGRA","PARAM","MAILLAGE",
                             "COMMANDE_RC","COMMANDE_DF","COMMANDE_PERIO",
                             "COMMANDE_SYRTHES","DONNEES_THERMOCHIMIE",
                             "NOMBRE_DE_PROCESSEURS","LISTE_PROCESSEURS",
                             "FICHIERS_DONNEES_UTILISATEURS",
                             "FICHIERS_RESULTATS_UTILISATEUR",
                             "OPTIMISATION","LISTE_LIB_SAT","OPTION_LIB_EXT",
                             "VALGRIND","ARG_ECS_VERIF","ARG_CS_VERIF","ECHOCOMM",
                             "PILOTAGE_ADAPTATION"]
            
            keywordsBsub = ["#BSUB -n","#BSUB -c"]
            
            keywordsValues = {}

            keywordsValBsub = {}
            
            while 1:
                line = lanceFile.readline()
                if (line == ""):
                    break
            
                for keywordLance in keywordsLance :
                    kw = re.compile("^"+keywordLance+"=")

                    lineVar = kw.match(line)
                    if lineVar :
                        tmp = line.split("=")

                        try:
                            keywordsValues[keywordLance] = str(tmp[1])
                        except:
                            keywordsValues[keywordLance] = None

                for keywordBsub in keywordsBsub :
                    kw = re.compile("^"+keywordBsub+" ")

                    lineVar = kw.match(line)
                    if lineVar :
                        tmp = line.split(" ")

                        try:
                            keywordsValBsub[keywordBsub] = str(tmp[2])
                        except:
                            keywordsValBsub[keywordBsub] = None
        
            lanceFile.close()

            testScriptsPath = self.studyPath+"/"+caseLabel.upper()+"/SCRIPTS"
            try:
                lanceTestFile = file(testScriptsPath+'/lance', mode='r')
            except IOError:
                print "Error : opening "+testScriptsPath+'/lance'
                print sys.exit(1)

            lanceTmpFile = file(testScriptsPath+'/lance.tmp', mode='w')

            while 1:
                line = lanceTestFile.readline()
                if (line == ""):
                    break

                indic = False
                for keywordLance in keywordsLance :
                    kw = re.compile("^"+keywordLance+"=")

                    lineVar = kw.match(line)
                    if lineVar :
                        lanceTmpFile.write(keywordLance+"="+keywordsValues[keywordLance]+'\n')
                        indic = True

                for keywordBsub in keywordsBsub :
                    kw = re.compile("^"+keywordBsub+" ")

                    lineVar = kw.match(line)
                    if lineVar :
                        lanceTmpFile.write(keywordBsub+" "+keywordsValBsub[keywordBsub]+'\n')
                        indic = True

                if Common.tmpDirectory != 'Default directory' :
                    kw = re.compile("^CS_TMP_PREFIX=")

                    lineVar = kw.match(line)
                    if lineVar :
                        lanceTmpFile.write("CS_TMP_PREFIX="+Common.tmpDirectory+'\n')
                        indic = True

                if not indic:
                    lanceTmpFile.write(line)

            lanceTmpFile.close()
            os.rename(testScriptsPath+'/lance.tmp',testScriptsPath+'/lance')
            

    def run(self):
        """
        run case
        """
        #
        # verifier que le cas ne devant pas etre recalcule a
        # bien des resultats sinon envoie d'un message a l'utilisateur
        # et relance du calcul
        
        run = True
        
        if not self.parser.getIfComputeCase(self.studyLabel,self.caseLabel):
            run = False
            
            if not self.isResult(self.caseLabel) :
                print "  -- No results for "+self.studyLabel+"/"+self.caseLabel
                run = True

        if run :
            resuList = os.listdir(self.casePath+"/RESU")
            for resu in resuList :
                try :
                    shutil.rmtree(self.casePath+"/RESU/"+resu)
                except :
                    os.remove(self.casePath+"/RESU/"+resu)

            print "        "+self.studyLabel+"/"+self.caseLabel+ " is running"
            # Ajouter la commande de lancement avec os.system
            testRunPath = self.studyPath+"/"+self.caseLabel.upper()+"/SCRIPTS"
            os.chdir(testRunPath)
            os.chmod(testRunPath+"/lance",0777)
            arch = os.uname()
            if (arch[0]=='OSF1' or (arch[1].find("tantal") >=0)):
                proc = os.popen("bsub < lance")
                proc.close()
            elif (arch[0]=='Linux_Ch'):
                proc = os.popen("qsub < lance")
                proc.close()
            else:
                proc = os.popen("nohup ./lance > list&")
                proc.close()
            end = False
            while 1:
                resuList = os.listdir(self.casePath+"/RESU")
                for resuFile in resuList:
                    if resuFile.find("error") >= 0:
                        return "Execution error"
                    elif resuFile.find("resume") >= 0:
                        end = True
                    elif resuFile.find("compil.log") >= 0:
                        compilFile = file(self.casePath+'/RESU/'+resuFile, mode='r')
                        while 1:
                            line = compilFile.readline()
                            if (line == ""):
                                break
                            error = "ERREUR DE COMPILATION OU D'EDITION DE LIENS"
                            err = re.compile(error)
                            lineVar = err.match(line)
                            if lineVar :
                                return "Compilation error"
                if end:
                    break
                time.sleep(10.)
                    
            os.chdir(Common.localDirectory)

        return "OK"

    def listingCompare(self):
        """
        compare listing with reference listing
        """
        refCasePath = Common.referencePath+"/"+self.studyLabel.upper()+"/"+self.caseLabel.upper()
        
        refListing = Listing.Listing(refCasePath)
        testListing = Listing.Listing(self.casePath)

        varMinRef, varMaxRef, clipMinRef, clipMaxRef = refListing.getMinMaxVariables(self.variables)
        varMinTest, varMaxTest, clipMinTest, clipMaxTest = testListing.getMinMaxVariables(self.variables)

        #
        # Creation d'un fichier cas_listing.report
        #
        fa = open(self.reportListing,'a')
        fa.write('\n')
        fa.write('Case: '+self.caseLabel.upper())
        fa.write('\n')
        fa.write('-------------------------------------------------------------------------------------------------\n')
        fa.write('    Variable          Ref Value        Test Value            Norm       ClipRef   ClipTest  <tol\n')
        fa.write('-------------------------------------------------------------------------------------------------\n')
        fa.close()
        print "        Listing files analysing"

        #
        # Gestion des erreurs
        #
        if len(varMinRef)==0:
            # Verification que la variable a bien ete trouvee
	    print '%12.12s' %self.variables, ' X -> Variable not find in listing'
	    return
        
        if len(varMinRef)!=len(varMinTest):
            # Verification que l'on effectue une comparaison avec des "fichiers comparables"
	    print '%12.12s' %self.variables, ' -> Different listing files ! case 1:',len(varMinRef), ' values - case 2:',len(varMinTest), ' values'		
            fa = open(self.reportListing,'a')
	    fa.write('%12.12s' %self.variables)
            fa.write(" X -> different listing files!")
	    fa.write("  case 1: ")
	    fa.write(str(len(varMinRef)))
	    fa.write(" values")
	    fa.write(" - case 2: ")
	    fa.write(str(len(varMinTest)))
	    fa.write(" values\n")
	    fa.close()
	    indic=1
	    return
        #
        # Calcul des normes
        #
        epsilon = 1E-12
        for i in self.variables :
            try:
                norme1 = abs(varMinRef[i]-varMinTest[i])/abs(varMaxRef[i]-varMinRef[i] + epsilon)
                norme2 = abs(varMaxRef[i]-varMaxTest[i])/abs(varMaxRef[i]-varMinRef[i] + epsilon)
            except:
                norme1 = None
                norme2 = None
                
        for i in self.variables :
            tolerance = self.variables[i]
            try:
                if norme1 > tolerance :
                    test1 ='NOK'
                else :
                    test1 ='OK'
            except:
                test1 = None

            try:
                if norme2 > tolerance :
                    test2 ='NOK'
                else :
                    test2 ='OK'
            except:
                test2 = None

            fa = open(self.reportListing,'a')
            fa.write('%12s' % i)
            fa.write('  min  ')
            try:
                fa.write('%12.5e' % varMinRef[i])
            except:
                fa.write('%12s' % varMinRef[i])
            fa.write('      ')
            try:
                fa.write('%12.5e' % varMinTest[i])
            except:
                fa.write('%12s' % varMinTest[i]) 
            fa.write('      ')
            try:
                fa.write('%12.5e' % norme1)
            except:
                fa.write('%12s' % norme1)
            fa.write('      ')
            try:
                fa.write('%5i' % clipMinRef[i])
            except:
                fa.write('%5s' % clipMinRef[i])
            fa.write('      ')
            try:
                fa.write('%5i' % clipMinTest[i])
            except:
                fa.write('%5s' % clipMinTest[i])
            fa.write('      ')
            fa.write(test1)	
            fa.write('\n')
            fa.write('            ')
            fa.write('  max  ')
            try:
                fa.write('%12.5e' % varMaxRef[i])
            except:
                fa.write('%12s' % varMaxRef[i])
            fa.write('      ')
            try:
                fa.write('%12.5e' % varMaxTest[i])
            except:
                fa.write('%12s' % varMaxTest[i])
            fa.write('      ')
            try:
                fa.write('%12.5e' % norme2)
            except:
                fa.write('%12s' % norme2)
            fa.write('      ')
            try:
                fa.write('%5i' % clipMaxRef[i])
            except:
                fa.write('%5s' % clipMaxRef[i])
            fa.write('      ')
            try:
                fa.write('%5i' % clipMaxTest[i])
            except:
                fa.write('%5s' % clipMaxTest[i])
            fa.write('      ')
            fa.write(test2)	
            fa.write('\n')
            fa.close()

        return test1,test2
            

    def chrCompare(self):
        """
        compare chrono file with reference file
        """
        refCasePath = Common.referencePath+"/"+self.studyLabel.upper()+"/"+self.caseLabel.upper()
        refResuPath = refCasePath+"/RESU"
        testResuPath = self.casePath+"/RESU"

        refFilesNames = os.listdir(refResuPath)
        testFilesNames = os.listdir(testResuPath)

        for fileName in refFilesNames :
            if fileName.find('CHR.ENSIGHT') >= 0 :
                refChrName = fileName

        for fileName in testFilesNames :
            if fileName.find('CHR.ENSIGHT') >= 0 :
                testChrName = fileName

        refChrono = Chrono.Chrono(refCasePath,self.variables)
        testChrono = Chrono.Chrono(self.casePath,self.variables)

        timeRef = refChrono.getTime()
        timeTest = testChrono.getTime()

        #
        # Creation d'un fichier cas_chrono.report
        #
        fa = open(self.reportChrono,'a')
        fa.write('Case: '+self.caseLabel.upper())
        fa.write('\n')
        fa.close()
        print "        Chrono files analysing \n"

        #
        # Gestion des erreurs
        #
        
        # Verification que les instants sont identiques a epsilon pres
        epsilon=1.0 # en secondes
 
        if abs(timeRef-timeTest)>epsilon :
            print " -> ERROR: chrono files incompatibles: different times !"
            fa = open(self.reportChrono,'a')
            fa.write(' X -> ERROR: chrono files incompatibles: different times)')
            fa.write('\n')    
            fa.close()
            return

        fa = open(self.reportChrono,'a')
        fa.write('--------------------------------------------------------------------\n')
        fa.write('Time = ')
        fa.write(str(timeRef))
        fa.write(' s\n')
        fa.write('--------------------------------------------------------------------\n')
        fa.write('    Variable     Delta_max     Delta_moy.      Norm     Norm<tol.  \n')
        fa.write('--------------------------------------------------------------------\n')
        fa.close()

        for variable in self.variables:
            tablRef = refChrono.getValues(variable)
            tablTest = testChrono.getValues(variable)

            delta = []

            # Calcul de delta(variable)
            for i in range (len(tablRef)):
                try:
                    delta.append(abs(tablRef[i] - tablTest[i]))
                except:
                    print "Warning : incompatibility between the files"
                    return "No possible comparison"
                
            vminRef = min(tablRef)
            vmaxRef = max(tablRef)
            
            # Calcul de l'ecart max
            deltaMax = max(delta)
            
            somVar = 0.0
              
            for i in range (len(delta)):
                somVar += delta[i]

            # Calcul de l'ecart moyen
            deltaMoy = somVar/len(delta)

            # Calcul de la norme
            epsilon= 1e-12
            norme = deltaMax/(vmaxRef-vminRef+epsilon)

            tolerance = self.variables[variable]
            if norme > tolerance:
                test='NOK'
                if self.parser.getIfComputeCase(self.studyLabel,self.caseLabel):
                    typ = testChrono.getTyp(variable)
                    testChrono.addVariable(delta, typ, "delta_"+variable, variable)
            else:
                test='OK'
                	
            #
            # Ecriture du fichier cas_chrono.report
            #
            fa = open(self.reportChrono,'a')
            fa.write('%12s' % variable)
            fa.write('  ')
            try:
                fa.write('%12.5e' % deltaMax)
            except:
                fa.write('%12s' % deltaMax)
            fa.write('  ')
            try:
                fa.write('%12.5e' % deltaMoy)
            except:
                fa.write('%12s' % deltaMoy)
            fa.write('  ')
            try:
                fa.write('%12.5e' % norme)
            except:
                fa.write('%12s' % norme)
            fa.write('      ')
            fa.write(test)	
            fa.write('\n')
            fa.close()    

        return test

        
    def getCaseLabel(self):
        return self.caseLabel
    

    def isResult(self, caseLabel):
        """
        find a result file in RESU
        """
        testResuPath = self.studyPath+"/"+caseLabel.upper()+"/RESU"
        
        try:
            resusList = os.listdir(testResuPath)
            if resusList != None and resusList != [] :
                return True
            else :
                return False
        except:
            return False
