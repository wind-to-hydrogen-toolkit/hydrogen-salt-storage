#!/bin/sh

mv temp/src/*.py src

mv temp/requirements.txt .

rm -d -r scripts
mv temp scripts
