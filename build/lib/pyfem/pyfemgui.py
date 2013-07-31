import pyfem
from pyfem.gui import *

app = wx.PySimpleApp(0)
wx.InitAllImageHandlers()
frame = MyFrame(None, -1, "")
app.SetTopWindow(frame)
frame.Show()
app.MainLoop()

