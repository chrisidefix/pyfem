#!/usr/bin/python
import os, string
import sys
import traceback

try:
    import pyfem
except:
    import os,sys; sys.path.insert(0, os.getcwd()+"/..")


count=0
setpath = "./"
for root, dirs, files in os.walk(setpath, topdown=False):
      for fname in files:
          s = fname.find('.py')
          if fname[:5]=='test_' and fname[-3:]=='.py':
              fpath = os.path.join(root,fname)
              fsize = os.path.getsize(fpath)
              try:
                  print
                  print '\033[96mRunning test in:', fpath, '...\033[0m'
                  execfile(fpath)
              except:
                  print '\a', # Beep
                  print '\033[91mError during test: \033[0m'
                  print '-'*60
                  traceback.print_exc(file=sys.stdout)
                  print '-'*60
              else:
                  print "Test", fname, "finalized."

