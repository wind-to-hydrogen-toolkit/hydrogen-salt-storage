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
black -l 79 temp/*.py

# sort imports
isort temp/*.py
