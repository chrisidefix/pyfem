import sys
from PyQt4 import QtCore, QtGui
 
 
class MainWindow(QtGui.QMainWindow):
 
    def __init__(self, *args):
         
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle("Adding and deleting widgets")
        self.resize(800, 600)
         
        self.mainWidget = QtGui.QWidget(self)
        self.setCentralWidget(self.mainWidget)
         
        self.buttonLayout = QtGui.QGridLayout(self.mainWidget)
         
        self.bnAdd = QtGui.QPushButton(self.mainWidget)
        self.bnAdd.setText("Ok")
        self.buttonLayout.addWidget(self.bnAdd,1,0,1,1)
 
        self.tree = QtGui.QTreeWidget(self.mainWidget)
        self.root = QtGui.QTreeWidgetItem(self.tree)
        self.tree.topLevelItem(0).setText(0, "Moteur")
         
        for i in range(10) :
            tt = QtGui.QTreeWidgetItem(self.root)
            self.tree.topLevelItem(0).child(i).setText(0, "%d" % i)
            self.tree.topLevelItem(0).child(i).setCheckState(0, QtCore.Qt.Unchecked)

        ss=QtGui.QTreeWidgetItem(tt)
        ss.setText(0,"A")
        ss.setCheckState(0, QtCore.Qt.Unchecked)
        self.tree.setColumnWidth(0,250)
        header=QtGui.QTreeWidgetItem(["Tree","First","secondo"])
        self.tree.setHeaderItem(header)
             
     
        self.buttonLayout.addWidget(self.tree,0,0,1,2)
        QtCore.QObject.connect(self.bnAdd, QtCore.SIGNAL("clicked()"), self.killChildreen)
     
    def killChildreen(self):
     
        tmp = []
        for i in range(self.tree.topLevelItem(0).childCount()) :
            if self.tree.topLevelItem(0).child(i).checkState(0) :
                tmp.append(self.tree.topLevelItem(0).takeChild(i))
                 
def main(args):
    app=QtGui.QApplication(args)
    win=MainWindow()
    win.show()
    app.connect(app, QtCore.SIGNAL("lastWindowClosed()")
                                 , app
                                 , QtCore.SLOT("quit()")
                                 )
    sys.exit(app.exec_())
 
if __name__=="__main__":
    main(sys.argv)

