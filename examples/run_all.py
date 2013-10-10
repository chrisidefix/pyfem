#!/usr/bin/python
import os, string
import sys
import traceback

count=0
rootpath = os.path.dirname(os.path.abspath(__file__))

setpath  = "./"
for root, dirs, files in os.walk(setpath, topdown=False):
      os.chdir(rootpath)
      print "root", root
      print "files", files

      for fname in files:
          if not fname.endswith(".py"): continue
          os.chdir( os.path.dirname( os.path.abspath(root) ) )

          if fname[:4]=='run_': continue
          #if fname[:5]=='test_' and fname[-3:]=='.py':

          fpath = os.path.join(root,fname)
          fsize = os.path.getsize(fpath)
          try:
              print
              print '\033[96mRunning test in:', fpath, '...\033[0m'
              execfile(fpath)
          except KeyboardInterrupt:
              exit()
          except:
              print '\a', # Beep
              print '\033[91mError during test: \033[0m'
              print '-'*60
              traceback.print_exc(file=sys.stdout)
              print '-'*60
          else:
              print "Test", fname, "finalized."

os.chdir(rootpath)
