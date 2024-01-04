#!/bin/sh

# format notebooks
black -l 79 */*.ipynb

# sort imports in notebooks
isort */*.ipynb

# convert Jupyter Notebooks to Python scripts
jupyter nbconvert --to script */*.ipynb

# move files to temporary directory
mkdir -p temp
mv notebooks/*.py temp

# remove "# In []" and multiple blank lines
for f in temp/*.py;
do sed -i -e '/^# In\[/d' $f
cat -s $f > $f.txt
mv $f.txt $f
done

# format scripts
black -l 79 **/*.py

# sort imports
isort **/*.py

# copy other stuff to main
cp -r src temp
cp -r tests temp
cp requirements.txt temp
cp README.md temp
