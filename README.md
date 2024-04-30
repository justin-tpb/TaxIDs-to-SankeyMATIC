# TaxIDs to SankeyMATIC
This script fetches taxonomic information from the Entrez database using a list of TaxIDs and outputs user-selected taxonomic ranks to a file.
A ready-to-use code of the taxonomic distribution is then generated for visualization with [SankeyMATIC](https://sankeymatic.com/).


## Features
* Fetch taxonomic information from Entrez for specified TaxIDs.
* Output the fetched data with user-selectable taxonomic ranks.
* Generate code compatible with SankeyMATIC for visualizing taxonomic distribution, sorted by hierarchy and count.
* Optional grouping of taxa with low counts to reduce clutter in the Sankey diagram.


## Prerequisites
Before running this script, it is necessary to have the following installed:
* Python 3.6 or later
* `pandas` and `biopython` libraries. To install these libraries, use:
```bash
pip install pandas biopython
```


## Usage
```bash
python taxids2sankey.py <input_file> [options]
```


### Arguments
* `<input_file>`: File containing a column of TaxIDs, one per line.
  * The first line must contain column headers.
  * The delimiter will be automatically detected.


### Options
* `-c`, `--header <name>`: Header of the column containing the TaxIDs. Default is `#Taxid` from the BLAST text output.
* `-t`, `--tax_ranks <taxonomic_ranks>`: Space-separated list of taxonomic ranks. Default is `class order genus`.
* `-e`, `--email <address>`: Email for identification by Entrez. Will be saved to `entrez_config.ini` for future use.
Entrez will show a warning without an email and might block access in case of excessive usage.
* `-g`, `--group <threshold>`: Group together ranks which are below this threshold for less cluttered SankeyMATIC diagrams,
which will then be named `<parent_rank> (grouped)`. Default is no grouping.
* `-s`, `--skip`: Skip the generation of SankeyMATIC compatible code (only output taxonomic information).
* `-h`, `--help`: Show this help message.


### Examples
Running the script with the default parameters:
```bash
python taxids2sankey.py example.csv -e email@example.com
```

Running the script after an email was saved to `entrez_config.ini`:
```bash
python taxids2sankey.py example.csv
```

Specifying a custom column header for the TaxIDs and skipping SankeyMATIC code generation:
```bash
python taxids2sankey.py example.csv -c TaxIDs -s
```

Specifying custom taxonomic ranks and enabling grouping for SankeyMATIC:
```bash
python taxids2sankey.py example.csv -t "phylum class order" -g 10
```


### Note
Entrez will sometimes cause an `HTTP Error 400 Bad Request` error. If this happens, just try again after a few seconds.


## Output
* A `<input_filename>.taxonomy.csv` file containing the input TaxIDs along with the fetched taxonomic information.
* If not skipped, a `<input_filename>.sankey.txt` file with SankeyMATIC code based on the taxonomic distribution.
  * The code will be sorted by taxonomic hierarchy and count.
  * If grouping is enabled, the grouping threshold will be added to the filename.
  * The file content can be copied to the input field of [SankeyMATIC](https://sankeymatic.com/build/) to generate a Sankey diagram.


## Author
Justin Teixeira Pereira Bassiaridis


## License
Distributed under the MIT License. See `LICENSE` for more information.
