from distutils.core import setup

# Try to remove old installation
import os, shutil, site
from distutils.sysconfig import get_python_lib

sitepackages_dir = site.getsitepackages()[0]
osname = os.name

if osname=='posix':
    pyfem_path = os.path.join(sitepackages_dir, 'pyfem')

if osname=='nt':
    pyfem_path = os.path.join(sitepackages_dir, 'Lib', 'site-packages', 'pyfem')

if os.path.isdir(pyfem_path):
    try:
        shutil.rmtree(pyfem_path)
        print "Deleting current installation path : ", pyfem_path
    except:
        print "Permission denied."
        pass

