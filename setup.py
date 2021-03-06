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
    print "Current installation path : ", pyfem_path
    shutil.rmtree(pyfem_path, ignore_errors=True)

# Installation setup
setup(name='pyfem',
      version="1.0",
      #package_dir={'pyfem':'../'},
      packages=['pyfem', 'pyfem.mesh', 'pyfem.tools', 'pyfem.equilib', 'pyfem.seepage', 'pyfem.hydromec']
      )
