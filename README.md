# Computational-Chemistry-Research-Project
Final year undergraduate university research project in computational chemistry searching for novel organic molecules that will display room temperature phosphorescence, (RTP), to be used in photoluminescence applications.

This python files contain the code written for the data analysis/computational chemistry part this project, the final result being a ~6,900 word dissertation, the abstract of the project is as follows: 

Organic materials pose an environmentally friendly, low cost and low toxicity option for photoluminescence applications; however, room temperature phosphorescence, (RTP), is rare in organic materials in solution. This paper explores the various phenomena that enhance a molecules photoluminescent properties, for example thermally activated delayed fluorescence, (TADF) along with the characteristics that are suggestive of photoluminescence in a molecule, such as the energy of the first excited triplet state, E(T1). A database consisting of 48,182 organic semiconductors has been analysed using computational methods and suggestions for molecules to be investigated further will be made in this thesis. The heterofission mechanism was considered for each molecule where carbazole, (Cz), was the host in the system and the molecules being tested were potential guests. A total of 214 Cz derivatives were found in the database, the geometries of which were optimised at the T1 minima to increase the accuracy of the heterofission mechanism rule. Suggestions were then made based on this mechanism, providing valuable suggestions for the guest in a Cz host-guest system that can be synthesised for potential photoluminescence applications.

The database was sourced from the University of Liverpool and can be found at the following link in the file 'CSD_EES_DB.csv': https://datacat.liverpool.ac.uk/1472/ 

The conclusion of this project was as follows: 

This research demonstrates that carbazole, (Cz), derivatives present a suitable option for the dopant in dopant-host photoluminescent systems, where Cz is the host molecule. 214 Cz derivatives were found in the database of 48182 organic semiconductors, where the data was originally obtained from the Cambridge Structural Database. In addition the research identifies molecules worth further investigation to be the dopant in the host-dopant system. The basis of these predictions came from first optimising the geometries of the Cz derivatives at the T1 minima, before calculating their deviation from the heterofission mechanism and comparing it to a system proven to emit RTP, Bd/Cz, which had a deviation of 0.34 eV. It is important to note that out of the 214 Cz derivatives, 76 of the geometry optimisation calculations were run successfully at the TDA-DFT/M06-2X/def2-SVP basis set and level of theory. This method was chosen as the same level of theory was used in the database, where the geometries of the molecules were optimised at the S0 minima. The T1 minima was chosen as this is one of the parameters in the heterofission mechanism, allowing for the dopant T1 energy level to be compared to the Cz host T1 energy level at a greater level of accuracy. The geometry optimisation calculations were run on Guassian16 and submitted to the Queen Mary HPC cluster on 05/04/2023, the remaining calculations are in the queue to be completed, unfortunately meaning the result of these could not be included in the final analysis (09/05/2023). Therefore, the time required to run large numbers of Guassian16 calculations should be a key consideration for future work.

Upon geometry optimisation at the T1 minima, five Cz derivatives were found that showed a deviation smaller than 0.34 eV, all of which possessed triplet values significantly lower than Cz, the host, which has also been proven to be an important factor in RTP as it reduces triplet migration. Based off previous research, the concentration of the dopant in the host-dopant system can be kept as low as 0.1 mol% as a Bd/Lab-Cz host-dopant system has shown to emit efficient RTP with this ratio. Therefore, when doping Cz with the five proposed molecules, a series of concentrations can be tested, starting with 0.1 mol%, in order to find the most effective systems.

Thus, this thesis has been successful in analysing a database of organic semiconductors, containing compounds that are stable in solid state whose electronic properties were not recorded for a specific purpose, making it unbiased with respect to application.24 For this reason, the five Cz derivatives that have been proposed as a suitable dopant in the Cz host- dopant system provide a novel suggestion in the RTP field as they have not been suggested for this purpose previously. Furthermore, an innovative insight into organic compounds for RTP photoluminescence applications has been achieved, and we hope to continue making break throughs in further explorations.

References (Royal Society of Chemistry format):

F. J. Hernández and R. Crespo-Otero, J. Mater. Chem. C, 2021, 9, 11882–11892
O. Omar, T. Nematiaram, A. Troisi and D. Padula, Organic Materials Repurposing, https://datacat.liverpool.ac.uk/1472/, (accessed May 2023)
Q. Sun, J. Ren, Q. Peng, Z Shuai, 2021, Research Square: https://assets.researchsquare.com/files/rs2520159/v1/d16e6bd7d2c450c8db398934.pdf?c=1676443353 (accessed May 2023)
T. King, S. Butcher and L. Zalewski, DOI:10.5281/ZENODO.438045
Q. Wang and H. Aziz, Appl. Phys. Lett., 2014, 105, 053304
C. Chen, Z. Chi, K. C. Chong, A. S. Batsanov, Z. Yang, Z. Mao, Z. Yang and B. Liu, Nat. Mater., 2021, 20, 175–180









