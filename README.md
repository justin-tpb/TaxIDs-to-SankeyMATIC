# TaxIDs to SankeyMATIC
This script retrieves taxonomic information from Entrez for a list of TaxIDs, outputs user-selected taxonomic ranks to a file, and generates ready-to-use code for visualizing the taxonomic distribution with [SankeyMATIC](https://sankeymatic.com/).


## Features
* Fetch taxonomic information from the Entrez database for up to 10,000 specified TaxIDs.
* Export the retrieved data to a CSV file, allowing users to select specific taxonomic ranks.
* Generate SankeyMATIC-compatible code for visualizing the taxonomic distribution of the selected ranks, sorted by hierarchy and count.
* Group taxa with low counts to declutter the Sankey diagram using a user-defined threshold.
* Specify an email address, as required by NCBI.


## Prerequisites
Before running this script, it is necessary to have the following installed:
* Python 3.6 or later
* `pandas` and `biopython` libraries. To install these libraries, execute the following command:
```
pip install pandas biopython
```


## Usage
```
python taxids2sankey.py <input_file> [options]
```
or
```
./taxids2sankey.py <input_file> [options]
```

### Arguments
* `<input_file>`: File containing a column of TaxIDs, one per line.
  * A maximum of 10,000 TaxIDs can be processed.
  * The first line must contain column headers.
  * The delimiter will be automatically detected.


### Options
* `-c`, `--header <name>`: Header of the column containing the TaxIDs.
  * Default is `#Taxid` from the BLAST text output.
* `-t`, `--tax_ranks <taxonomic_ranks>`: Comma-separated list of taxonomic ranks.
  * Default is `class,order,genus`.
* `-e`, `--email <address>`: Email address for identification by NCBI.
  * Will be saved to `entrez_config.ini` for future use.
  * E-utilities will show a warning without an email address.
* `-g`, `--group <threshold>`: Group ranks below this threshold to declutter Sankey diagrams.
  * Groups will be named `<parent_rank> (grouped)`.
  * Default is no grouping.
* `-s`, `--skip`: Skip the generation of SankeyMATIC compatible code (only output taxonomic information).
* `-h`, `--help`: Show this help message.


### Examples
Running the script with the default parameters:
```
python taxids2sankey.py example.csv -e email@example.com
```

Running the script after an email was saved to `entrez_config.ini`:
```
python taxids2sankey.py example.csv
```

Specifying a custom column header for the TaxIDs and skipping SankeyMATIC code generation:
```
python taxids2sankey.py example.csv -c TaxIDs -s
```

Specifying custom taxonomic ranks and enabling grouping for SankeyMATIC:
```
python taxids2sankey.py example.csv -t phylum,class,order -g 10
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
Distributed under the MIT License. See [`LICENSE`](https://github.com/justin-tpb/TaxIDs-to-SankeyMATIC/blob/main/LICENSE) for more information.
