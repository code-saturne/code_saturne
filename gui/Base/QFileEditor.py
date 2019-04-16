import sys
from PyQt4 import QtGui, QtCore

class QFileEditor(QtGui.QMainWindow):

    def __init__(self, parent=None):
        super(QFileEditor, self).__init__(parent)
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("Text editor")
#        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))

        openFile = QtGui.QAction("Open", self)
        openFile.setShortcut("Ctrl+O")
        openFile.setStatusTip('Open File')
        openFile.triggered.connect(self.openFile)

        newFile = QtGui.QAction("New", self)
        newFile.setShortcut("Ctrl+E")
        newFile.setStatusTip('Create new file')
        newFile.triggered.connect(self.newFile)

        saveFile = QtGui.QAction("Save", self)
        saveFile.setShortcut("Ctrl+S")
        saveFile.setStatusTip('Save file')
        saveFile.triggered.connect(self.saveFile)

        quitAction = QtGui.QAction("Quit", self)
        quitAction.setShortcut("Ctrl+Q")
        quitAction.setStatusTip('Quit the editor')
        quitAction.triggered.connect(self.closeApplication)

        self.statusBar()

        mainMenu = self.menuBar()

        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(newFile)
        fileMenu.addAction(openFile)
        fileMenu.addAction(saveFile)
        fileMenu.addAction(quitAction)

    def openFile(self):
        name = QtGui.QFileDialog.getOpenFileName(self, 'Open File')

        if name != None and name != '':
            file = open(name,'r')

            self.editor()
            with file:
                text = file.read()
                self.textEdit.setText(text)

    def newFile(self):
        self.textEdit = QtGui.QTextEdit()
        self.setCentralWidget(self.textEdit)


    def saveFile(self):
        name = QtGui.QFileDialog.getSaveFileName(self, 'Save File')

        if name != None and name != '':
            file = open(name,'w')
            text = self.textEdit.toPlainText()
            file.write(text)
            file.close()


    def closeApplication(self):
        choice = QtGui.QMessageBox.question(self, 'Built-in editor',
                                            "Exit text editor?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass


