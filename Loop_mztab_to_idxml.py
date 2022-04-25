#DO NOT RUN as PYTHON
#Run in Terminal

#loop mztab to idxml
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time /Users/jensvandeperre/opt/anaconda3/bin/python3 /Users/jensvandeperre/Master-Thesis-2022/mztab_to_idxml.py /Users/jensvandeperre/Desktop/Inputs/ALL_mzTab/${array[0]}.mztab /Users/jensvandeperre/Desktop/Inputs/idxml/ALL_idxml/${array[0]}.idxml
done


#loop peptideindex n idxmls
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/idxml_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time peptideindexer -in /Users/jensvandeperre/Desktop/Inputs/idxml_test/${array[0]}.idxml -fasta /Users/jensvandeperre/Desktop/Inputs/Seq_database_fasta/Human_Proteome.fasta -missing_decoy_action "silent" -IL_equivalent -enzyme:name 'Trypsin' -enzyme:specificity "semi" -out /Users/jensvandeperre/Desktop/Outputs/peptideindexer/${array[0]}.idxml
done
