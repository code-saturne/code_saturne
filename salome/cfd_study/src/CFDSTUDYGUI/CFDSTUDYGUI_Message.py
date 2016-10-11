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
    MESSAGE="CFDSTUDY Module Message"
    def __init__(self,):
        QDialog.__init__(self)
        print cfdStudyMessageImpl.MESSAGE

    def aboutMessage(self,message):
        msg=QMessageBox()
        texte=message
        msg.about(sgPyQt.getDesktop(),cfdStudyMessageImpl.MESSAGE,texte)

cfdstudyMess = cfdStudyMessageImpl()
