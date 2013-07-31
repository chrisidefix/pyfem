import wx, random


class LazyTree(wx.TreeCtrl):
    ''' LazyTree is a simple "Lazy Evaluation" tree, that is, it only adds 
        items to the tree view when they are needed.'''

    def __init__(self, *args, **kwargs):
        super(LazyTree, self).__init__(*args, **kwargs)
        self.Bind(wx.EVT_TREE_ITEM_EXPANDING, self.OnExpandItem)
        self.Bind(wx.EVT_TREE_ITEM_COLLAPSING, self.OnCollapseItem)
        self.__collapsing = False
        root = self.AddRoot('root')
        self.SetItemHasChildren(root)

    def OnExpandItem(self, event):
        # Add a random number of children and randomly decide which 
        # children have children of their own.
        nrChildren = random.randint(1, 6)
        for childIndex in range(nrChildren):
            child = self.AppendItem(event.GetItem(), 'child %d'%childIndex)
            self.SetPyData(child, ["asljfdlsdfjl"])
            self.SetItemHasChildren(child, random.choice([True, False]))

    def OnCollapseItem(self, event):
        # Be prepared, self.CollapseAndReset below may cause
        # another wx.EVT_TREE_ITEM_COLLAPSING event being triggered.
        if self.__collapsing:
            event.Veto()
        else:
            self.__collapsing = True
            item = event.GetItem()
            self.CollapseAndReset(item)
            self.SetItemHasChildren(item)
            self.__collapsing = False


class LazyTreeFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super(LazyTreeFrame, self).__init__(*args, **kwargs)
        self.__tree = LazyTree(self)


if __name__ == "__main__":
    app = wx.App(False)
    frame = LazyTreeFrame(None)
    frame.Show()
    app.MainLoop() 

