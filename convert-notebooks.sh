#!/bin/sh

# activate Python environment
source .venv/bin/activate

# convert Jupyter Notebooks to Python scripts
jupyter nbconvert --to script */*.ipynb

# remove "# In []" and multiple blank lines in converted scripts
for f in temp/*.py;
do sed -i -e '/^# In\[/d' $f
cat -s $f > $f.txt
mv $f.txt $f
done
