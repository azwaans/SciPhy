SciPhy
-------

SciPhy (Sequential Cas-9 Insertion based Phylogenetics) is a [BEAST 2](http://www.beast2.org/) package to estimate time-scaled lineage trees and cell population dynamics from ordered Cas9-mediated insertion data. This package is designed to handle any alignment of such order-aware insertion/prime-editing based lineage tracing constructs, such as: 
- DNA Typewriter lineage tracing constructs, as described [here](https://doi.org/10.1038/s41586-022-04922-8).
- peCHYRON, as described [here](https://doi.org/10.1101/2021.11.05.467507).

The preprint describing the model behind SciPhy and its applications can be found [here](https://doi.org/10.1101/2024.10.01.615771).

Installation
-------
SciPhy requires an installation of BEAST 2.7 which can be obtained from https://www.beast2.org/. SciPhy can then be installed using BEAUti with the following instructions: 

- Open BEAUti.
- In the `File` menu, select `Manage Packages`.
- Click the `Package repositories` button at the bottom of the dialog box.
- Click on`Add URL` and enter the following URL:
https://raw.githubusercontent.com/azwaans/SciPhy/refs/heads/master/package.xml. Close this window.
- `sciphy` should now appear as a package in the list of available packages. Select it and click the `Install/Upgrade` button. 
- SciPhy should now be available for use in BEAST 2 and Beauti. Close and restart Beauti to start generating an xml. 

Running Sciphy on your data
-------

Start from a csv file, where the rows are cells, the columns are the consecutive sites in the barcode, and the entries denote the edits as integers. Then use the `write_nexus.R` script under `scripts/` to generate a `[your_filename].sciphy` file that can be read by BEAUTi. 

Open BEAUTi, click on File -> Template -> Sciphy to open the Sciphy template. Then, click on the "plus" button in the bottom left corner and read your `[your_filename].sciphy` into BEAUTi. Now, you can specify your preferred Site Model, Clock Model etc. For general information on setting up a Bayesian phylogenetic analysis, refer to the [BEAST tutorials](https://taming-the-beast.org/).

You can find an example xml file in `examples/example.xml` file, that is based on the `scripts/example_data.csv`.

To run the analysis from the command line, run the following command:
`beast -seed 1 -overwrite example.xml`

You can also run the analysis using the GUI. To do that, open BEAST. When asked for the input file, choose your xml and hit the "Run" button.

License
-------

SciPhy is free software.  It is distributed under the terms of version 3
of the GNU General Public License.  A copy of this license should
be found in the file COPYING located in the root directory of this repository.
If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.
