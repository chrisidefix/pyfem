cd doc

touch index.rst
make html

mkdir -p ~/Dropbox/Public/pyfem-doc
unalias cp
cp -r _build/html ~/Dropbox/Public/pyfem-doc/ > /dev/null

cd ..
