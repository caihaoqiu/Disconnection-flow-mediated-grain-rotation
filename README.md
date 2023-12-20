# Disconnection flow-mediated grain rotation
Python codes to model the grain rotation whithin a multi-disconnection-mode framework published in PNAS (DOI needed to be added).

Here, we provide example python codes considering 2 reference interfaces (for the results in Fig.2) and 4 reference interfaces (examples of Σ17 and Σ29 shown in Fig.3, whose rotation directions are clockwise and counter-clockwise, respectively).

In directory "phase diagram", you could find the code and data files to plot the figures in Fig.4. Please be minded that we adopted the opposite conventions to define the rotation directions in Figs.2-3 and the modified Cahn-Taylor model (fig.4). Therefore, we provide two columns of rotation rates obtained by our continuum model in each data file corresponding to the ordinary convention used in Figs.2-3 and the opposite convention (to be consistent with the convention in modified Cahn-Taylor model). We plot Fig.4b, c with the opposite convention in order to directly compare with the mCT model.  

Please be reminded that all the parameters used in the codes are dimensionless and the definition of the reduced units can be found in the main text and supporting information. Due to the existence of sigularity points, the current codes might be instable. If one would like to change the parameters in the codes, please use smaller timestep and avoid low temperature and large Burgers vectors (which will lead to instability of the code).

We will continuously improve the readability and stability of the codes. Please pay attention to our future updated versions.
