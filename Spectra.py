import sys

input_file = '' #input file uses 17mix_test2.mzxml.gz in the repository, which must be unzipped first.
input_scan_number = 1301
peptide = 'TYDSYLGDDYVR'

#Add Sys command line arguments
#input_file = sys.argv[1]
#input_scan_number = sys.argv[2]
#peptide = sys.argv[3]

#Establish molecular weights
MW = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07, 'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08, 'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10, 'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

#compute b-ion and y-ion m/z values
ions = {}
for i in range(1, len(peptide)):
    b_ion = peptide[0:i]
    b_ion_mass = sum([MW[aa] for aa in b_ion]) + 1
    ions['b'+str(i)] = b_ion_mass
    #Repeat for y_ion
    y_ion = peptide[i:len(peptide)]
    y_ion_mass = sum([MW[aa] for aa in y_ion]) + 19
    ions['y' + str(len(peptide)-i)] = y_ion_mass

#parse the mzXML file to get the correct scan and its peaks element
import xml.etree.ElementTree as ET
ns = r'{http://sashimi.sourceforge.net/schema/}'
for event, ele in ET.iterparse(input_file):
    if ele.tag == ns+'scan':
        number = ele.attrib.get('num')
        if int(number) == int(input_scan_number):
            peakselt = ele.find(ns+'peaks')
            break

#Extract the binary spectra data from the peaks element
from array import array
from base64 import b64decode
peaks = array('f',b64decode(peakselt.text)) # peakselt is the XML element corresponding to the peaks list
if sys.byteorder != 'big':
   peaks.byteswap()
mzs = peaks[::2]
ints = peaks[1::2]
ints = [x/max(ints)*100 for x in ints] #Converting absolute intensities to relative abundance (%)

#Group the mz, int pairs by their matches to b-ions, y-ions, or non-matches
b_mzs = []
b_labels = []
b_ints = []
y_mzs = []
y_ints = []
y_labels = []
non_mzs = []
non_ints = []
non_labels = []
for mz, int in zip(mzs, ints):
    for key, mass in ions.items():
        if abs(mz - mass) < 0.3 and int > max(ints)*0.05:  #cutoff values
            if key[0] == 'b':
                b_mzs.append(mz)
                b_ints.append(int)
                b_labels.append(key)
            elif key[0] == 'y':
                y_mzs.append(mz)
                y_ints.append(int)
                y_labels.append(key)
            #break #Allows only one ion match for each mz
    if mz not in b_mzs and mz not in y_mzs:
        non_mzs.append(mz)
        non_ints.append(int)
        non_labels.append(key)

#Plot the peaks
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
plt.xlabel('m/z')
plt.ylabel('Relative abundance (%)')
plt.title('scan #' + input_scan_number + ', Peptide: ' + peptide + ', length = ' + str(len(peptide)))
if b_mzs:
    ax.stem(b_mzs, b_ints, basefmt=' ', markerfmt='', linefmt='blue')
if y_mzs:
    ax.stem(y_mzs, y_ints, basefmt=' ', markerfmt='', linefmt='red')
if non_mzs:
    ax.stem(non_mzs, non_ints, basefmt=' ', markerfmt='', linefmt='grey')

#label b-ion matches:
for i in range(len(b_mzs)):
    ax.text(b_mzs[i], b_ints[i] + 5, b_labels[i], color='blue')
#label y-ion matches:
for i in range(len(y_mzs)):
    ax.text(y_mzs[i], y_ints[i], y_labels[i], color='red')

plt.show()




