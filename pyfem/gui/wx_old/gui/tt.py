import wx
import wx.lib.agw.customtreectrl as CT

class MyFrame(wx.Frame):

    def __init__(self, parent):

        wx.Frame.__init__(self, parent, -1, "CustomTreeCtrl Demo")

        # Create a CustomTreeCtrl instance
        custom_tree = CT.CustomTreeCtrl(self, agwStyle=wx.TR_DEFAULT_STYLE)

        # Add a root node to it
        root = custom_tree.AddRoot("The Root Item")

        # Create an image list to add icons next to an item
        il = wx.ImageList(16, 16)
        fldridx     = il.Add(wx.ArtProvider_GetBitmap(wx.ART_FOLDER,      wx.ART_OTHER, 16))
        fldropenidx = il.Add(wx.ArtProvider_GetBitmap(wx.ART_FILE_OPEN,   wx.ART_OTHER, 16))
        fileidx     = il.Add(wx.ArtProvider_GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, 16))

        custom_tree.SetImageList(il)

        custom_tree.SetItemImage(root, fldridx, wx.TreeItemIcon_Normal)
        custom_tree.SetItemImage(root, fldropenidx, wx.TreeItemIcon_Expanded)

        for x in range(15):
            child = custom_tree.AppendItem(root, "Item %d" % x)
            custom_tree.SetItemImage(child, fldridx, wx.TreeItemIcon_Normal)
            custom_tree.SetItemImage(child, fldropenidx, wx.TreeItemIcon_Expanded)

            for y in range(5):
                last = custom_tree.AppendItem(child, "item %d-%s" % (x, chr(ord("a")+y)))
                custom_tree.SetItemImage(last, fldridx, wx.TreeItemIcon_Normal)
                custom_tree.SetItemImage(last, fldropenidx, wx.TreeItemIcon_Expanded)

                for z in range(5):
                    item = custom_tree.AppendItem(last,  "item %d-%s-%d" % (x, chr(ord("a")+y), z))
                    custom_tree.SetItemImage(item, fileidx, wx.TreeItemIcon_Normal)
                    custom_tree.SetItemImage(item, smileidx, wx.TreeItemIcon_Selected)

        custom_tree.Expand(self.root)


# our normal wxApp-derived class, as usual

app = wx.PySimpleApp()

frame = MyFrame(None)
app.SetTopWindow(frame)
frame.Show()

app.MainLoop()

