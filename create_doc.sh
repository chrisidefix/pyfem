cd doc

touch index.rst
mkdir -p _build
mkdir -p _build/html

make html

mkdir -p ~/Dropbox/Public/pyfem-doc
unalias cp
cp -r _build/html ~/Dropbox/Public/pyfem-doc/ 
#> /dev/null

cd ..
