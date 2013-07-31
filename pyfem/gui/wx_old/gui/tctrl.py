#!/usr/bin/env python

import wx
import os.path, dircache

class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # this is just setup boilerplate
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)

        # our tree object, self.tree
        self.tree = wx.TreeCtrl(self, -1, style=wx.TR_HAS_BUTTONS|wx.TR_DEFAULT_STYLE|wx.SUNKEN_BORDER)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.tree, 1, wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer)
        sizer.Fit(self)
        sizer.SetSizeHints(self)
        self.Layout()

        # register the self.onExpand function to be called
        wx.EVT_TREE_ITEM_EXPANDING(self.tree, self.tree.GetId(), self.onExpand)
        # initialize the tree
        self.buildTree('C:\\Dropbox')

    def onExpand(self, event):
        '''onExpand is called when the user expands a node on the tree
        object. It checks whether the node has been previously expanded. If
        not, the extendTree function is called to build out the node, which
        is then marked as expanded.'''

        # get the wxID of the entry to expand and check it's validity
        itemID = event.GetItem()
        if not itemID.IsOk():
            itemID = self.tree.GetSelection()

        # only build that tree if not previously expanded
        old_pydata = self.tree.GetPyData(itemID)
        if old_pydata[1] == False:
            # clean the subtree and rebuild it
            self.tree.DeleteChildren(itemID)
            self.extendTree(itemID)
            self.tree.SetPyData(itemID,(old_pydata[0], True))

    def buildTree(self, rootdir):
        '''Add a new root element and then its children'''
        self.rootID = self.tree.AddRoot(rootdir)
        self.tree.SetPyData(self.rootID, (rootdir,1))
        self.extendTree(self.rootID)
        self.tree.Expand(self.rootID)

    def extendTree(self, parentID):
        '''extendTree is a semi-lazy directory tree builder. It takes
        the ID of a tree entry and fills in the tree with its child
        subdirectories and their children - updating 2 layers of the
        tree. This function is called by buildTree and onExpand methods'''

        # This is something to work around, because Windows will list
        # this directory but throw a WindowsError exception if you
        # try to use the listdir() command on it. I need a better workaround
        # for this...this is a temporary kludge.
        excludeDirs=["c:\\System Volume Information","/System Volume Information"]

        # retrieve the associated absolute path of the parent
        parentDir = self.tree.GetPyData(parentID)[0]


        subdirs = dircache.listdir(parentDir)
        subdirs.sort()
        for child in subdirs:
            child_path = os.path.join(parentDir,child)
            if os.path.isdir(child_path) and not os.path.islink(child):
                if child_path in excludeDirs:
                    continue
                # add the child to the parent
                childID = self.tree.AppendItem(parentID, child)
                # associate the full child path with its tree entry
                self.tree.SetPyData(childID, (child_path, False))

                # Now the child entry will show up, but it current has no
                # known children of its own and will not have a '+' showing
                # that it can be expanded to step further down the tree.
                # Solution is to go ahead and register the child's children,
                # meaning the grandchildren of the original parent
                newParent = child
                newParentID = childID
                newParentPath = child_path
                newsubdirs = dircache.listdir(newParentPath)
                newsubdirs.sort()
                for grandchild in newsubdirs:
                    grandchild_path = os.path.join(newParentPath,grandchild)
                    if os.path.isdir(grandchild_path) and not os.path.islink(grandchild_path):
                        grandchildID = self.tree.AppendItem(newParentID, grandchild)
                        self.tree.SetPyData(grandchildID, (grandchild_path,False))

if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    frame = MyFrame(None, -1, "")
    app.SetTopWindow(frame)
    frame.Show()
    app.MainLoop()

