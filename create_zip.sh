# Create a zip file containing
# pyfem directory
# test  directory
# setup.py    file 
# install.bat file 

FILENAME=pyfem
FILENAME="$FILENAME-`date +%Y-%m-%d`.zip"

find . -name "*.pyc" -delete -print

zip -r $FILENAME pyfem test setup.py install.bat install.sh

mkdir -p dist
mv -f $FILENAME ./dist/

