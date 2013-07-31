import wx

class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=(500,400))
        panel = wx.Panel(self)
        self.tree  = wx.TreeCtrl(panel, wx.ID_ANY, None, size=(200,300))
        root = self.tree.AddRoot('PyFEM Analysis')
        self.tree.SetItemHasChildren(root)
        ch_mesh = self.tree.AppendItem(root, 'Mesh')
        ch_anys  = self.tree.AppendItem(root, 'Analysis')

        #self.control = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self.CreateStatusBar() # A Statusbar in the bottom of the window

        # Setting up the menu.
        filemenu= wx.Menu()

        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        filemenu.Append(wx.ID_ABOUT, "&About"," Information about this program")
        filemenu.AppendSeparator()
        filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")

        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        self.Show(True)

app = wx.App(False)
frame = MainWindow(None, "Sample editor")
app.MainLoop()

