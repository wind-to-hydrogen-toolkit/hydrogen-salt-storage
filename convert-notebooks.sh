#!/bin/sh

# convert Jupyter Notebooks to Python scripts
jupyter nbconvert --to script docs/notebooks/*.ipynb

# clear scripts/ directory
rm -rf scripts/*

# move scripts to scripts/ directory
mv docs/notebooks/*.py scripts

# remove "# In []" and multiple blank lines in converted scripts
for f in scripts/*.py;
do sed -i -e '/^# In\[/d' $f
cat -s $f > $f.txt
mv $f.txt $f
done
