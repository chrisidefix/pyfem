class LNode:
    def __init__(self, data):
        self.data = data
        self.prev = None
        self.next = None

class LList:
    def __init__(self, alist=None):
        self.begin  = None
        self.end    = None
        self.lenght = 0

        if alist:
            for data in alist:
                self.add_node(data)

    def add_node(self, data):
        new_node = LNode(data)

        if self.begin==None:
            self.begin = new_node

        if self.end:
            new_node.prev = self.end
            self.end.next = new_node

        self.end = new_node
        self.lenght += 1

    def insert_node(self, node, data):
        if not self.begin: raise Exception("LList::__iter__: Empty linked list.")

        new_node = LNode(data)

        if node==self.end:
            new_node.prev = self.end
            new_node.next = None
            self.end.next = new_node
            self.end      = new_node
        else:
            if node:
                new_node.next   = node.next
                new_node.prev   = node
                node.next.prev  = new_node
                node.next       = new_node
            else:   # node==None: insert at begin
                new_node.prev   = None
                new_node.next   = self.begin
                self.begin.prev = new_node
                self.begin      = new_node

        self.lenght += 1

    def remove_node(self, node):
        if not self.begin: raise Exception("LList::__iter__: Empty linked list.")

        if node==self.begin:
            self.begin      = node.next
            self.begin.prev = None
        else:
            node.prev.next  = node.next

        if node==self.end:
            self.end        = node.prev
            self.end.next   = None
        else:
            node.next.prev = node.prev

        node.prev = None
        node.next = None

        self.lenght -= 1


    def __iter__(self):
        curr_node = self.begin
        while curr_node:
            yield curr_node
            curr_node = curr_node.next

    def __repr__(self):
        return "LList("+str([n.data for n in self])+")"



if __name__=="__main__":


    def print_LL(L):
        for n in L:
            print "node:", n.data,
            if n.prev:
                print "prev:", n.prev.data,
            if n.next:
                print "next", n.next.data,
            print
        print "---"

    # Test
    L = LList()
    L.add_node(1)
    L.add_node(2)
    L.add_node(3)

    print_LL(L)

    L.remove_node(L.begin.next)
    print_LL(L)

    L.insert_node(L.begin, 2)
    print_LL(L)

    L.insert_node(L.end, 4)
    print_LL(L)

    L.insert_node(L.begin.prev, 0)
    print_LL(L)

    print_LL(LList([5,6,7,8,9]) )

    print L

