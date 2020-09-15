{\rtf1\ansi\ansicpg1252\cocoartf2513
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fswiss\fcharset0 Helvetica-Bold;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww18220\viewh11320\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\tx13105\pardirnatural\partightenfactor0

\f0\fs24 \cf0 Code to generate clones, integrate, and classify TNOs.\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 Instructions:\
0. Download all of the files in the repository.\
1. Create a folder for each TNO to be classified. This folder should contain an aei file and abg file for the object. \
2. Create an empty folder called \'91res\'92 in each TNO folder.\
3. Create a csv file with names of all folders.\
4. Generate clones for each TNO by running 
\f1\b make_clones_aei.py 
\f0\b0 (adjust path names in this file). This script outputs a small.in file for each object.\
5. Run integrations. The folder \'93sample_integration_params\'94 has the *.in files that were used in the classification paper.\
6. Separately generate a table of barycentric aei elements for the best fit of each TNO orbit.\
7. Classify TNOs by running 
\f1\b main.py
\f0\b0  (adjust path names in this file). The script main.py calls 
\f1\b classify.py
\f0\b0  (generates plots of the integrations and does the straightforward classifications) and  
\f1\b find_resonance.py 
\f0\b0 (looks for resonances).\
\
I am providing an example by using resonant object \'91s12_good_3\'92. The folder \'91s12_good_3\'92 has the files needed for step 1 above. The folder \'91s12_good_3_output\'92 is the end result of the full classification process (what you should end up with).\
\
All of the scripts should be fairly straightforward except for find_resonance.py. This script has quite a few parameters that might need to be tuned depending on the data. For example, since the code relies on identifying dense regions of the resonance angle plots, the simulation output frequency is important. I would recommend using the *.in files provided in the sample_integration_params folder for this reason.\
\
Additional details about outputs.\
-the folder \'91res\'92 will hold plots of each instance of a resonance found for a TNO clone, a summary plot, and a pandas table of the resonance information.\
-elements.txt stores the best fit orbit\
-scattering.txt stores number of clones that are scattering\
-classification.txt stores classification of TNO IF resonances are NOT taken into account\
-ejected.txt stores number of clones that were ejected during the integration\
-resonance.txt stores less detailed resonance information than the resonance_table in \'91res\'92 \
\
}