class A():
    def aa(self,i):
        return i
    def a1(self,i):
        return i
    def a2(self,i):
        return i
    def a3(self,i):
        return i
    def a4(self,i):
        return i
    def a5(self,i):
        return i
    def a6(self,i):
        return i

B = A()

for i in range(10000):
    #print getattr(B,"aa")(i), "\r",
    print B.aa(i), "\r",



