# gamma_dist
Fitting Whitson's Gamma Distribution to C10+ PVT sample composition

The algorithm is based on Curtis Whitson's work primarily extracted
from Phase Behavior SPE Monograph (Volume 20) by Whitson and Brule.
Additional information sources are numerous papers on Gamma distribution
which can be found at Whitson's website https://whitson.com/publications/
Particularly useful was the 2019 paper by Bilal Younus, Curtis Whitson et al
"Field-Wide Equation of State Model Development" that is also available from
the download section of Whitson's website.
Gamma distribution fitting assumes C10 as the starting component.
The script take a .csv file with the following columns with names (in square brakets):
  1st column [SCN]: SCN identifiers (e.g., C10, C11, C12 etc.);
  2nd column [mfi_lab]: mole fraction of component as per full composition;
  3rd column [wfi_lab]: weight fraction of component as per full composition.
  Average sample molecular weight is entered in the main section of the code.
