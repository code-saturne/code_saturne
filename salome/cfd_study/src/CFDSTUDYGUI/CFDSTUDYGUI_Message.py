#Qt4/Qt5
try :
  from PyQt4.QtGui import QDialog, QMessageBox
except :
  from PyQt5.QtWidgets import QDialog, QMessageBox

# --- Get SALOME PyQt interface
import SalomePyQt
sgPyQt = SalomePyQt.SalomePyQt()

class cfdStudyMessageImpl(QDialog):
    """
    """
    MESSAGE="CFDSTUDY Module Message: "
    def __init__(self):
        QDialog.__init__(self)
        self.msg=QMessageBox()

    def aboutMessage(self,message):
        messageInformation = cfdStudyMessageImpl.MESSAGE+ "Information"
        return self.msg.about(sgPyQt.getDesktop(),messageInformation,message)

    def criticalMessage(self,message):
        criticalMess = cfdStudyMessageImpl.MESSAGE + "Error"
        return self.msg.critical(sgPyQt.getDesktop(), criticalMess, message, QMessageBox.Ok, QMessageBox.No)

    def warningMessage(self,message):
        warningMess = cfdStudyMessageImpl.MESSAGE + "Warning"
        return self.msg.warning(sgPyQt.getDesktop(), warningMess, message, QMessageBox.Ok, QMessageBox.No)

    def trMessage(self,mess,listeMots):
        for i in listeMots :
            try:#Qt4
                mess = mess.replace(mess.indexOf("%s"),2,i)
            except:#Qt5
                mess = mess.replace("%s",i,1)
        return mess

cfdstudyMess = cfdStudyMessageImpl()
