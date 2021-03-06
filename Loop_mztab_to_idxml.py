#DO NOT RUN as PYTHON
#Run in Terminal

#loop mztab to idxml
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time /Users/jensvandeperre/opt/anaconda3/bin/python3 /Users/jensvandeperre/Master-Thesis-2022/mztab_to_idxml.py /Users/jensvandeperre/Desktop/Inputs/ALL_mzTab_pure_seq/${array[0]}.mztab /Users/jensvandeperre/Desktop/Inputs/idxml/${array[0]}.idxml
done

/Users/jensvandeperre/opt/anaconda3/bin/python3 /Users/jensvandeperre/Master-Thesis-2022/mztab_to_idxml.py /Users/jensvandeperre/Desktop/Inputs/ALL_mzTab_pure_seq/01CPTAC_COprospective_W_PNNL_20170123_B1S1_f02.mztab /Users/jensvandeperre/Desktop/Inputs/idxml/TEST_WORKS.idxml



#loop peptideindex on idxmls
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/idxml_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time peptideindexer -in /Users/jensvandeperre/Desktop/Inputs/idxml/${array[0]}.idxml -fasta /Users/jensvandeperre/Desktop/Inputs/Seq_database_fasta/Human_Proteome.fasta -missing_decoy_action "silent" -IL_equivalent -enzyme:name 'Trypsin' -enzyme:specificity "semi" -out /Users/jensvandeperre/Desktop/Inputs/peptideindexer/${array[0]}.idxml
done

peptideindexer -in /Users/jensvandeperre/Desktop/Inputs/idxml/WORKS.idxml -fasta /Users/jensvandeperre/Desktop/Inputs/Seq_database_fasta/Human_Proteome.fasta -missing_decoy_action "silent" -IL_equivalent -enzyme:name 'Trypsin' -enzyme:specificity "semi" -out /Users/jensvandeperre/Desktop/Outputs/peptideindexer/WORKS.idxml


#loop PIA_compile
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/idxml_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time java -jar /Applications/pia-1/pia.jar /Users/jensvandeperre/Desktop/Inputs/peptideindexer/${array[0]}.idxml -o /Users/jensvandeperre/Desktop/Inputs/PIA_compile/${array[0]}.xml --compile 
done

java -jar /Applications/pia-1/pia.jar /Users/jensvandeperre/Desktop/Outputs/peptideindexer/WORKS.idxml -o /Users/jensvandeperre/Desktop/Inputs/PIA_compile/WORKS.xml --compile 


#loop PIA_analysis
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/idxml_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time java -jar /Applications/pia-1/pia.jar /Users/jensvandeperre/Desktop/Inputs/PIA_parameter_files/${array[0]}.json /Users/jensvandeperre/Desktop/Inputs/PIA_compile/${array[0]}.xml
done

java -jar /Applications/pia-1/pia.jar /Users/jensvandeperre/Desktop/Inputs/PIA_parameter_files/01CPTAC_COprospective_W_PNNL_20170123_B1S1_f01.json /Users/jensvandeperre/Desktop/Inputs/PIA_compile/01CPTAC_COprospective_W_PNNL_20170123_B1S1_f01.xml


#loop PTM identification mass tolerance 10
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time /Users/jensvandeperre/opt/anaconda3/bin/python3 /Users/jensvandeperre/Master-Thesis-2022/PTM_identification.py /Users/jensvandeperre/Desktop/Inputs/PTM_mass_differences/Mass_tolerance_10/${array[0]}.csv /Users/jensvandeperre/Desktop/Inputs/Unimod_library/unimod_py.csv /Users/jensvandeperre/Desktop/Outputs/PTM_identification_tol_10/Mod_${array[0]}.csv
done

#loop PTM identification mass tolerance 20
name_file=/Users/jensvandeperre/Desktop/Inputs/file_names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
IFS=';' read -r -a array <<< "$line"
time /Users/jensvandeperre/opt/anaconda3/bin/python3 /Users/jensvandeperre/Master-Thesis-2022/PTM_identification.py /Users/jensvandeperre/Desktop/Inputs/PTM_mass_differences/Mass_tolerance_20/${array[0]}.csv
/Users/jensvandeperre/Desktop/Inputs/Unimod_library/
/Users/jensvandeperre/Desktop/Outputs/PTM_identification/Mass_tolerance_20/Mod_${array[0]}.csv
done



