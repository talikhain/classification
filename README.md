Code to generate clones, integrate, and classify TNOs.
Instructions:
0. Download all of the files in the repository.
1. Create a folder for each TNO to be classified. This folder should contain an aei file and abg file for the object. 
2. Create an empty folder called ‘res’ in each TNO folder.
3. Create a csv file with names of all folders.
4. Generate clones for each TNO by running make_clones_aei.py (adjust path names in this file). This script outputs a small.in file for each object.
5. Run integrations. The folder “sample_integration_params” has the *.in files that were used in the classification paper.
6. Separately generate a table of barycentric aei elements for the best fit of each TNO orbit.
7. Classify TNOs by running main.py (adjust path names in this file). The script main.py calls classify.py (generates plots of the integrations and does the straightforward classifications) and  find_resonance.py (looks for resonances).

I am providing an example by using resonant object ‘s12_good_3’. The folder ‘s12_good_3’ has the files needed for step 1 above. The folder ‘s12_good_3_output’ is the end result of the full classification process (what you should end up with).

All of the scripts should be fairly straightforward except for find_resonance.py. This script has quite a few parameters that might need to be tuned depending on the data. For example, since the code relies on identifying dense regions of the resonance angle plots, the simulation output frequency is important. I would recommend using the *.in files provided in the sample_integration_params folder for this reason.

Additional details about outputs.
-the folder ‘res’ will hold plots of each instance of a resonance found for a TNO clone, a summary plot, and a pandas table of the resonance information.
-elements.txt stores the best fit orbit
-scattering.txt stores number of clones that are scattering
-classification.txt stores classification of TNO IF resonances are NOT taken into account
-ejected.txt stores number of clones that were ejected during the integration
-resonance.txt stores less detailed resonance information than the resonance_table in ‘res’ 

