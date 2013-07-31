from distutils.core import setup
setup(name='pyfem',
      version='1.0',
      #package_dir={'pyfem':'../'},
      packages=['pyfem', 'pyfem.mesh', 'pyfem.tools', 'pyfem.equilib', 'pyfem.seepage'],
      #packages=['pyfem']
      #package_dir={'': 'src'},
      )
