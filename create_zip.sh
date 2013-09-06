# Create a zip file containing
# pyfem directory
# test  directory
# setup.py    file 
# install.bat file 

FILENAME=pyfem
FILENAME="$FILENAME-`date +%Y-%m-%d`.zip"

zip -r $FILENAME pyfem test setup.py install.bat install.sh

mkdir -p dist
mv $FILENAME ./dist/

