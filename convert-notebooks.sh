#!/bin/sh

mv temp/src/*.py src

mv temp/requirements.txt .
mv temp/README.md .

rm -d -r scripts
mv temp scripts
