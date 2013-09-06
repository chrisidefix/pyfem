#!/bin/bash

#find . -type f -name "*.pyc" -exec rm -f  '{}' \;
find . -iname "*.pyc" -delete -print
rm -rf build/
rm -rf dist/

