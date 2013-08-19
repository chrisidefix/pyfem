#!/usr/bin/env python

import wx
import wx.gizmos
import os.path, dircache
from inspect import isfunction, isclass

from pyfem.mesh import *

class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):

        # Base class constructor
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        kwds["size"]  = (800,600)
        wx.Frame.__init__(self, *args, **kwds)


        # Setting up menu bar
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
        menubar.Append(fileMenu, '&File')
        self.SetMenuBar(menubar)
        
        #self.Bind(wx.EVT_MENU, self.OnQuit, fitem)
        #self.SetMenuBar(menubar)

        # Setting up toolbar
        #self.toolbar = self.CreateToolBar()
        #self.toolbar.AddLabelTool(1, '', wx.Bitmap('icons//back.png'))
        #self.toolbar.AddLabelTool(1, '', wx.Bitmap('icons//save.png'))
        #self.toolbar.AddLabelTool(1, '', wx.Bitmap('icons//delete.png'))
        #self.toolbar.AddLabelTool(1, '', wx.Bitmap('icons//settings.png'))
        #self.toolbar.AddLabelTool(1, '', wx.Bitmap('icons//play.png'))
        #self.toolbar.Realize()

        #self.Bind(wx.EVT_TOOL, self.OnQuit, qtool)
        
        self.panel = wx.Panel(self, size=(-1,-1))

        # Setting up the Tree control
        self.tree = wx.gizmos.TreeListCtrl(self.panel, -1,
                size=(900,500),style=wx.TR_HAS_BUTTONS|wx.TR_DEFAULT_STYLE|wx.SUNKEN_BORDER|wx.TR_EDIT_LABELS)
                #pos=(0,100),
                #size=(800,500),style=wx.TR_EDIT_LABELS)



        #sizer = wx.BoxSizer(wx.VERTICAL)


        self.txt = wx.TextCtrl(self.panel, -1, pos=(0,0), size=(800,100) )
        #self.txt = wx.TextCtrl(self.panel, -1, style=wx.TE_MULTILINE, pos=(0,0), size=(800,100) )
        #sizer.Add(self.txt, 1, wx.EXPAND, 0)


        #sizer.Add(self.tree, 1, wx.EXPAND, 0)
        #self.SetAutoLayout(True)
        #self.SetSizer(sizer)
        #sizer.Fit(self)
        #sizer.SetSizeHints(self)


        #sizer.Add(panel,0, wx.EXPAND, 0)
        #panel.SetSizer(sizer)
        #panel.Layout()

        #self.Layout()

        # register the self.onExpand function to be called
        wx.EVT_TREE_ITEM_EXPANDING(self.tree, self.tree.GetId(), self.onExpand)
        # initialize the tree
        self.buildTree()

        print self.panel.GetSize()
        wx.EVT_TREE_ITEM_RIGHT_CLICK(self.tree, self.tree.GetId(), self.onRightClick) 
        wx.EVT_TREE_SEL_CHANGED(self.tree, self.tree.GetId(), self.onSelectionChange) 

        wx.EVT_SIZE(self, self.onResize)

    def onResize(self, event):
        print "onResize"
        h, v = self.GetClientSize()
        v_tbox = 40
        v_tree = v - v_tbox
        self.panel.SetDimensions(0,0,h,v)
        self.txt.SetDimensions(0,0, h, v_tbox)
        self.tree.SetDimensions(0,v_tbox,h,v_tree)
        self.Refresh()


    def onExpand(self, event):
        '''onExpand is called when the user expands a node on the tree
        object. It checks whether the node has been previously expanded. If
        not, the extendTree function is called to build out the node, which
        is then marked as expanded.'''

        ## get the wxID of the entry to expand and check it's validity
        #itemID = event.GetItem()
        #if not itemID.IsOk():
        #    itemID = self.tree.GetSelection()

        ## only build that tree if not previously expanded
        #old_pydata = self.tree.GetPyData(itemID)
        #if old_pydata[1] == False:
        #    # clean the subtree and rebuild it
        #    self.tree.DeleteChildren(itemID)
        #    self.extendTree(itemID)
        #    self.tree.SetPyData(itemID,(old_pydata[0], True))


    def buildTree(self):
        '''Add a new root element and then its children'''
        self.tree.AddColumn("Class")
        self.tree.AddColumn("Value")
        self.tree.SetColumnWidth(0,200)
        self.tree.SetColumnWidth(1,300)
        self.rootID = self.tree.AddRoot("Analysis")
        #self.tree.AppendItem(self.rootID, "Mesh")
        self.tree.AppendItem(self.rootID, "Domain")
        self.tree.AppendItem(self.rootID, "Solver")
        self.tree.Expand(self.rootID)
        #self.rootID = self.tree.AddRoot(rootdir)
        #self.tree.SetPyData(self.rootID, (rootdir,1))
        #self.extendTree(self.rootID)
        #self.tree.Expand(self.rootID)


        # Load the mesh branch
        #self.buildBranch(Mesh, self.rootID)
        #mesh = Mesh()
        #data = mesh.gui_get_default_data()
        
        #mesh_item = self.tree.AppendItem(self.rootID, Mesh.__name__)
        #self.buildBranch(Mesh, mesh_item)

        self.buildBranch(Mesh, self.rootID)

        #for key, feature in data.iteritems():
        #    display = feature.get('display')
        #    value   = str(feature.get('value',''))
        #    expand  = feature.get('expand', False)
        #    print key, feature
        #    item = self.tree.AppendItem(mesh_item, display)
        #    if not expand:
        #        self.tree.SetItemText(item, value, 1)
        #    else:
        #        # call recursively for value
        #        edit = False
            
    def buildBranch(self, obj_type, root_item):
        print obj_type
        obj = obj_type()
        data = obj.gui_get_default_data()
        #branch_item = self.tree.AppendItem(root_item, obj_type.__name__)
        root_item = self.tree.AppendItem(root_item, obj_type.__name__)
        for key, feature in data.iteritems():
            display = feature.get('display')
            value   = str(feature.get('value',''))
            expand  = feature.get('expand', False)
            #print key, feature
            if not expand:
                item = self.tree.AppendItem(root_item, display)
                self.tree.SetItemText(item, value, 1)
            else:
                # call recursively for value
                edit = False
                print ">",value
                typ  = feature['type']
                
                item = self.tree.AppendItem(root_item, typ.__name__)
                self.buildBranch(typ, item)

                #item = self.buildBranch(typ, root_item)  #<<<<<<<<<<<<<<
                
                pydata = {}
                self.tree.SetPyData(item, pydata)
                pydata['context'] = typ().gui_context()
                menu = wx.Menu()
                mitem = menu.Append(-1, "Pop item")
                self.Bind(wx.EVT_MENU, self.onPopupSel, mitem)
                self.tree.Bind(wx.EVT_CONTEXT_MENU, self.onShowPopup)

    #def onTreePopupItemSel(self, event):
    #    item = self.tree.GetSelection()
    #    if isclass(self.popup_command):
    #        bind_func = self.buildBranch(self.popup_command, item)

    #    pass


    def onRightClick(self, event):
        print "onRightClick"
        itemID = event.GetItem()
        self.tree.SelectItem(itemID)
        data = self.tree.GetPyData(itemID)

        if data:
            context = data.get('context', None)
            if context:
                menu = wx.Menu()
                for key, menu_item in data['context'].iteritems():
                    mitem = menu.Append(-1, menu_item['display'])
                    command = menu_item['command']
                    if isclass(command):
                        handle = lambda event: self.buildBranch(command, itemID)
                    else:
                        handle = command

                    self.Bind(wx.EVT_MENU, handle, mitem)


                pos = self.GetPosition()
                pos = event.GetPoint()
                pos = pos[0]+20, pos[1]+20
                self.tree.PopupMenu(menu, pos)
        pass


    def extendTree(self, parentID):
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


    def onPopupSel(self):
        print "onPopupSel"
        pass

    def onShowPopup(self):
        print "onShowPopup"
        pass

    def onSelectionChange(self, event):
        print "onSelectionChange"
        itemID = event.GetItem()
        #self.tree.SelectItem(itemID)
        self.txt.SetValue(self.tree.GetItemText(itemID, 1))
        pass

def on_test(event):
    print "on_test"


def launch_gui():
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    frame = MyFrame(None, -1, "")
    app.SetTopWindow(frame)
    frame.Show()
    app.MainLoop()

if __name__ == "__main__":
    launch_gui()

