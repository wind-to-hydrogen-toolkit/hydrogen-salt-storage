#!/bin/sh

mv temp/src/*.py src
mv temp/tests/*.py tests

mv temp/requirements.txt .
mv temp/README.md .

rm -d -r scripts
mv temp scripts
