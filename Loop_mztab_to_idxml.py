name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
/Users/jensvandeperre/opt/anaconda3/bin/python3 /Users/jensvandeperre/Master-Thesis-2022/mztab_to_idxml.py /Users/jensvandeperre/Desktop/Inputs/ALL_mzTab/${array[0]}.mztab /Users/jensvandeperre/Desktop/Inputs/idxml/ALL_idxml/${array[0]}.idxml
done
