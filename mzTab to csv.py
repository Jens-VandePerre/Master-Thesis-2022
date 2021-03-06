import sys
#print(sys.executable)  #sys.executable contains full path of the currently running Python interpreter
import pandas
import sys
import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('mztab', type=pathlib.Path)
parser.add_argument('bait', type=str)
args = parser.parse_args()

path = '~/Desktop/mzTab/Stored files/Test 6 mzTabs'
output_path = path + 'csv/psm/' + args.bait + '.csv'
output_path = args.bait + '.csv'

with open (args.mztab,'r') as f:
    outF = open(output_path, 'w')
    for line in f:
        if line.startswith('MTD'):
            continue
        if line.startswith('PRT'):
            continue
        if line.startswith('PRH'):
            continue
        else:
            outF.write(line)
    outF.close()
    
#mztabList.append(line)