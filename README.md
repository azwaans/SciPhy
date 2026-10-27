# SciPhy


SciPhy (Sequential Cas9 Insertion-based Phylogenetics) is a BEAST 2 package designed to estimate time-scaled lineage trees and cell population dynamics from lineage-tracing data characterized by the sequential accumulation of mutations within lineage-tracing constructs.

While SciPhy provides a framework built to generically handle any type of mutation events (e.g., insertions, deletions, substitutions, or complex edits) that accumulate in an ordered sequence, it was natively developed for modalities relying on Cas9 prime-editing-mediated insertions, such as:

The manuscript describing the base model behind SciPhy and its applications can be found [here](https://www.nature.com/articles/s41467-026-73377-6). 
The base model was subsequently extended to handle incomplete alignments inputs, as described [here](https://www.biorxiv.org/content/10.64898/2026.10.07.757283v1).

---

## System requirements

SciPhy runs on any hardware and operating system supported by BEAST 2 (Windows, macOS, Linux). BEAST 2.8 or newer is required. No non-standard hardware is needed. SciPhy has been tested on BEAST 2.6, 2.7 and 2.8.

---

## Installation

Install BEAST 2.8 or newer from https://www.beast2.org/. Then install SciPhy using BEAUti:

1. Open BEAUti.
2. In the `File` menu, select `Manage Packages`.
3. Click the `Package repositories` button at the bottom of the dialog box.
4. Click on`Add URL` and enter the following URL:
https://raw.githubusercontent.com/azwaans/SciPhy/refs/heads/master/package.xml. Close this window.
5. `sciphy` should now appear as a package in the list of available packages. Select it and click the `Install/Upgrade` button. 
6. Close and restart BEAUti.

**Typical installation time:** ~1 minute if BEAST 2 is already installed; ~5 minutes including installation of BEAST 2.

---

## Demo


An example XML file is provided at `examples/example.xml`, based on `scripts/example_data.csv`. To run the demo from the command line:

```bash
beast -seed 1 -overwrite examples/example.xml
```

**Expected run time:** ~5 minutes on a standard desktop computer.

---


## Running SciPhy on your data


Start from a csv file, where rows are cells, columns are the consecutive sites in the barcode, and the entries denote the edits as integers. Then use the `write_nexus.R` script under `scripts/` to generate a `[your_filename].sciphy` file that can be read by BEAUTi.

1. Open BEAUTi, click on File -> Template -> Sciphy to open the Sciphy template.
2. Click on the "+" button in the bottom left corner and load your `[your_filename].sciphy` file.
3. Specify your preferred Site Model, Clock Model etc. For general information on setting up a Bayesian phylogenetic analysis, refer to the [BEAST tutorials](https://taming-the-beast.org/).
4. Export your XML file


To run your analysis from the command line:

```bash
beast -seed 1 -overwrite your_analysis.xml`
```

Or use the BEAST GUI: open BEAST, select your XML file, and click "Run".


A detailed description of the model and algorithm is provided in the Methods section of the accompanying manuscript ([link](https://www.nature.com/articles/s41467-026-73377-6).

---

## License

SciPhy is free software.  It is distributed under the terms of version 3
of the GNU General Public License.  A copy of this license should
be found in the file COPYING located in the root directory of this repository.
If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.
