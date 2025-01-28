#!/usr/bin/env python3

import argparse
import collections
import configparser
import csv
import textwrap
import os
import sys
import pandas as pd
from Bio import Entrez

def parse_args():
    """ Set up command-line argument parsing """
    description = textwrap.dedent("""
    This script fetches taxonomic information from the Entrez database using a list of TaxIDs and outputs selected taxonomic ranks to a file.
    A ready-to-use code of the taxonomic distribution is then generated for visualization with SankeyMATIC (https://sankeymatic.com/).
    Taxa with low counts can be grouped together for less clutter. The Sankey diagram is sorted by taxonomic hierarchy and count.
    Two non-default python modules are required, pandas and biopython. Install with: pip install pandas biopython
    """)
    parser = argparse.ArgumentParser(description=description, epilog="Made by Justin Teixeira Pereira Bassiaridis from the Hess Lab at TU Darmstadt.",
                                     add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input_file", metavar="<input_file>", type=str, help="File containing a column of TaxIDs, one per line. First line must contain column headers.")
    parser.add_argument("-c", "--header", metavar="<name>", type=str, help="Header of the column containing the TaxIDs (default: '#Taxid' from BLAST text output).", default="#Taxid")
    parser.add_argument("-t", "--tax_ranks", metavar="<taxonomic_ranks>", type=str, help="Comma-separated list of taxonomic ranks in hierarchical order (default: class,order,genus).\n"
                        "Commonly used ranks: superkingdom,kingdom,phylum,class,order,family,genus,species\n"
                        "All available ranks: clade,domain,superkingdom,kingdom,subkingdom,infrakingdom,superphylum,phylum,subphylum,infraphylum,superclass,class,subclass,infraclass,\n"
                        "                     superorder,order,suborder,infraorder,parvorder,superfamily,family,subfamily,tribe,subtribe,genus,subgenus,species,subspecies,variety,form\n"
                        "Note: The rank 'clade' can only be used to fetch taxonomic information, but not to generate SankeyMATIC code.\n"
                        "      This is because 'clade' does not follow any hierarchical structure and can also have multiple entries.\n"
                        "      When 'clade' is selected, '-s' is applied automatically.", default="class,order,genus")
    parser.add_argument("-e", "--email", metavar="<address>", type=str, help="Email for identification by Entrez. Will be saved to 'entrez_config.ini' for future use.")
    parser.add_argument("-g", "--group", metavar="<threshold>", type=int, help="Group together ranks which are below this threshold for SankeyMATIC (default: no grouping).", default=0)
    parser.add_argument("-s", "--skip", help="Skip the generation of SankeyMATIC compatible code (only output taxonomic information).", action="store_true")
    parser.add_argument("-h", "--help", action="help", help="Show this help message.")
    args = parser.parse_args()

    # Check if the input file exists
    if not os.path.isfile(args.input_file):
        print(f"\nError: The input file '{args.input_file}' does not exist.\n")
        sys.exit(1)

    # Automatically apply -s/--skip if "clade" is in taxonomic ranks
    if "clade" in args.tax_ranks.split(","):
        args.skip = True
        print("\nSkipping SankeyMATIC code generation because the 'clade' rank was selected.")

    return args

def setup_config(email):
    """ Write or read the configuration file """
    config = configparser.ConfigParser()
    config_file = "entrez_config.ini"
    if os.path.isfile(config_file):
        config.read(config_file)

    if email:  # Write email to config file
        Entrez.email = email
        if not config.has_section("Entrez"):
            config.add_section("Entrez")
        config.set("Entrez", "email", email)
        with open(config_file, "w") as file:
            config.write(file)
        print(f"\nEmail saved to '{config_file}'.\n")
    elif os.path.isfile(config_file):  # Use existing email if available
        if config.has_section("Entrez"):
            Entrez.email = config["Entrez"].get("email", None)
            if Entrez.email:
                print(f"\nEmail read from '{config_file}'.\n")
            else:
                print(f"\nNo email found in '{config_file}'. Please provide an email using the -e option to avoid issues with Entrez.\n")
        else:
            print(f"\nNo email found in '{config_file}'. Please provide an email using the -e option to avoid issues with Entrez.\n")
    else:
        print("\nPlease provide an email using the -e option to avoid issues with Entrez.\n")

def detect_delimiter(input_file):
    """ Automatically detect the input table delimiter """
    with open(input_file, "r", newline="") as file:
        try:
            dialect = csv.Sniffer().sniff(file.read(8192))  # Read a small portion of the file
            return dialect.delimiter
        except csv.Error:
            print("Input table delimiter could not be detected. Exiting...\n")
            sys.exit(1)

def fetch_taxonomies(taxids):
    """ Fetch taxonomic information for all TaxIDs in one request from the Entrez taxonomy database """
    print("Fetching taxonomic information from the Entrez database...")
    handle = Entrez.efetch(db="taxonomy", id=taxids, retmax=10000)  # Limited to 10000 TaxIDs
    records = Entrez.read(handle)
    handle.close()
    return records

def append_taxonomies(input_file, header, tax_ranks, delimiter):
    """ Append the taxonomic information from Entrez to the TaxIDs from the input file and output the result in a new file """
    # Load DataFrame from input file, filter entries and fetch taxonomies
    tax_df = pd.read_csv(input_file, sep=delimiter)
    tax_df.columns = tax_df.columns.str.strip()
    header = header.strip()
    if header not in tax_df.columns:
        print(f"Header '{header}' not found in input file. Exiting...\n")
        sys.exit(1)
    header_prefix = "#" if header.startswith("#") else ""
    tax_df[header] = tax_df[header].astype(str).str.strip()
    valid_taxids = tax_df[tax_df[header].apply(lambda x: x.isdigit())][header].tolist()
    invalid_taxids = tax_df[~tax_df[header].apply(lambda x: x.isdigit())].copy()
    if not valid_taxids:
        print("No valid TaxIDs found in the input file. Please check your data and try again.\n")
        sys.exit(1)
    print(f"Valid TaxIDs: {len(valid_taxids)}")
    invalid_taxids[""] = invalid_taxids[header].apply(lambda x: "(empty)" if x.strip() == "" else "(nondigit)").str.ljust(10)
    print(f"Invalid TaxIDs: {len(invalid_taxids)}\n{invalid_taxids}\n")
    tax_records = fetch_taxonomies(valid_taxids)

    # Build taxid_to_record mapping, considering "AkaTaxIds"
    taxid_to_record = {}
    for record in tax_records:
        record_taxid = record["TaxId"]
        aka_taxids = record.get("AkaTaxIds", [])
        if not isinstance(aka_taxids, list):
            aka_taxids = [aka_taxids]
        all_taxids = [record_taxid] + aka_taxids
        for taxid in all_taxids:
            taxid_to_record[taxid] = record

    # Process records to extract taxonomic hierarchy and append to the DataFrame
    tax_ranks = tax_ranks.split(",")
    max_clade_count = 0
    for taxid in valid_taxids:
        record = taxid_to_record.get(taxid)
        if not record:
            print(f"TaxID {taxid} does not exist.")
            continue
        tax_data = {rank: [] for rank in tax_ranks}  # Default to empty list for all ranks

        # Extract scientific names from "LineageEx"
        for taxon in record["LineageEx"]:
            if taxon["Rank"] in tax_ranks:
                tax_data[taxon["Rank"]].append(taxon["ScientificName"])

        # Fetch the rank of the current TaxID
        taxid_rank = record.get("Rank")
        if taxid_rank in tax_ranks:
            tax_data[taxid_rank].append(record.get("ScientificName", f"Unknown_{taxid_rank}"))

        # Set placeholders if rank is unknown
        for i, rank in enumerate(tax_ranks):
            if not tax_data[rank]:
                if i == 0:  # Default for most basal rank
                    tax_data[rank].append(f"Unknown_{rank}")
                else:  # Default for other ranks
                    prev_rank = tax_ranks[i - 1]
                    if tax_data[prev_rank]:
                        prev_taxon_name = tax_data[prev_rank][0].split("_")[0]
                        if prev_taxon_name == "Unknown":
                            tax_data[rank].append(f"{prev_taxon_name}_{rank}")
                        else:
                            tax_data[rank].append(f"{prev_taxon_name}_unclassified_{rank}")

        # Update the DataFrame with taxonomy information
        idx = tax_df[tax_df[header] == str(taxid)].index  # Index of rows with matching TaxIDs
        for rank in tax_ranks:
            if rank == "clade":
                for i, value in enumerate(tax_data[rank]):
                    column_name = f"{header_prefix}{rank.capitalize()}_{i + 1}"
                    tax_df.loc[idx, column_name] = value
                max_clade_count = max(max_clade_count, len(tax_data[rank]))
            else:
                if tax_data[rank]:
                    column_name = f"{header_prefix}{rank.capitalize()}"
                    tax_df.loc[idx, column_name] = tax_data[rank][0]

    # Ensure that all clade columns up to the maximum count are present in the DataFrame in the correct order
    for i in range(1, max_clade_count + 1):
        clade_column = f"{header_prefix}Clade_{i}"
        if clade_column not in tax_df.columns:
            tax_df[clade_column] = pd.NA

    # Reorder columns to match the taxonomic rank order specified by -t/--tax_ranks
    ordered_columns = []
    for rank in tax_ranks:
        if rank == "clade":
            for i in range(1, max_clade_count + 1):
                ordered_columns.append(f"{header_prefix}Clade_{i}")
        else:
            ordered_columns.append(f"{header_prefix}{rank.capitalize()}")

    # Include original columns that are not taxonomic ranks
    non_tax_columns = [col for col in tax_df.columns if col not in ordered_columns]
    tax_df = tax_df[non_tax_columns + ordered_columns]
    tax_df.fillna("-", inplace=True)

    # Remove placeholder underscores from the taxonomic columns before saving
    for column in tax_df.columns:
        if header_prefix in column:
            pattern = f"_{column.strip('#').lower()}"
            replacement = f" {column.strip('#').lower()}"
            tax_df[column] = tax_df[column].astype(str).replace(pattern, replacement, regex=True)
            tax_df[column] = tax_df[column].astype(str).replace("_unclassified", " unclassified", regex=True)

    # Save the enhanced DataFrame
    output_filename = f"{os.path.splitext(input_file)[0]}.taxonomy.csv"
    tax_df.to_csv(output_filename, sep=delimiter, index=False)
    print(f"\nTaxonomic information has been saved to '{output_filename}'.\n")

    # Replace square brackets in taxon names to avoid errors in SankeyMATIC code generation
    tax_df = tax_df.map(lambda x: x.replace("[", "(").replace("]", ")") if isinstance(x, str) else x)

    return tax_df, header_prefix

def generate_sankeymatic_code(tax_df, header_prefix, tax_ranks, group_threshold):
    """ Generate SankeyMATIC compatible code for the taxonomic distribution with optional grouping """
    sankeymatic_code = []
    node_presence = set()  # To avoid orphaned entries

    # Setup for linking "All" to the highest rank taxa
    tax_ranks = tax_ranks.split(",")
    highest_rank = f"{header_prefix}{tax_ranks[0].capitalize()}"
    filtered_df = tax_df[tax_df[highest_rank] != "-"]  # Exclude entries with no taxonomic information
    total_counts = filtered_df[highest_rank].value_counts()
    normal_entries = {}
    other_entries = {}

    # Process totals for the highest rank
    for rank, count in total_counts.items():
        if count >= group_threshold:  # Handle as normal entry
            normal_entries[rank] = count
            node_presence.add(rank)
        else:  # Handle as entry to be grouped
            other_entries[rank] = count

    # Link "All" to the highest rank taxa
    for rank, count in normal_entries.items():
        sankeymatic_code.append(f"All [{count}] {rank}")
        node_presence.add(rank)

    # If only one entry is below the threshold, do not group to avoid unnecessary loss of taxonomic information
    if len(other_entries) == 1 and sum(other_entries.values()) < group_threshold:
        for rank, count in other_entries.items():
            sankeymatic_code.append(f"All [{count}] {rank}")
            node_presence.add(rank)
    elif other_entries:  # Group multiple low count entries
        total_other = sum(other_entries.values())
        sankeymatic_code.append(f"All [{total_other}] Other {tax_ranks[0].lower()} (grouped)")

    # Process each rank, handling linkage based on node presence
    for i in range(len(tax_ranks) - 1):
        current_rank = f"{header_prefix}{tax_ranks[i].capitalize()}"
        next_rank = f"{header_prefix}{tax_ranks[i + 1].capitalize()}"

        # Exclude entries with no taxonomic information in the current or next rank
        filtered_df = tax_df[(tax_df[current_rank] != "-") & (tax_df[next_rank] != "-")]

        group_data = filtered_df.groupby([current_rank, next_rank]).size()
        grouped_entries = {}

        # Determine entries to be grouped
        for (higher, lower), count in group_data.items():
            if higher in node_presence:  # Only process if the parent node is present to avoid orphaned entries
                if count >= group_threshold:
                    sankeymatic_code.append(f"{higher} [{count}] {lower}")
                    node_presence.add(lower)
                else:
                    if higher not in grouped_entries:
                        grouped_entries[higher] = {}
                    grouped_entries[higher][lower] = count

        # Group low count entries under existing parents
        for higher, lower_counts in grouped_entries.items():
            if higher in node_presence:  # Only process if the parent node is present to avoid orphaned entries
                # If only one entry below threshold, do not group to avoid unnecessary loss of taxonomic information
                if len(lower_counts) == 1 and sum(lower_counts.values()) < group_threshold:
                    for lower, count in lower_counts.items():
                        sankeymatic_code.append(f"{higher} [{count}] {lower}")
                        node_presence.add(lower)
                else:  # Group multiple low count entries
                    total_count = sum(lower_counts.values())
                    sankeymatic_code.append(f"{higher} [{total_count}] {higher} (grouped)")

    return sankeymatic_code

def sort_and_output_sankeymatic_code(input_file, group_threshold, sankeymatic_code):
    """ Sort the SankeyMATIC compatible code by taxonomic hierarchy and count """

    def parse_line(line):
        """ Extract the count and ranks from a line of the SankeyMATIC code """
        count = int(line.split("[")[1].split("]")[0])
        parts = line.split("] ")
        left_part = parts[0].split("[")[0].strip()
        right_part = parts[1].strip()
        return left_part, right_part, count

    def build_hierarchy(data):
        """ Build a hierarchical tree from the data """
        tree = collections.defaultdict(list)

        # Populate the tree with data
        for left_part, right_part, count in data:
            tree[left_part].append((right_part, count))

        # Sort children of each node by their count in descending order
        for key in tree:
            tree[key].sort(key=lambda x: x[1], reverse=True)

        return tree

    def sort_hierarchy(tree):
        """ Sort the hierarchy, remove cycles, and generate the output iteratively for a Sankey-compatible structure """
        sorted_output = []
        visited_nodes = set()
        stack = [("All", 0, None, None)]

        # Stack entries are tuples of (node, depth, parent, count)
        while stack:
            current_node, depth, parent, count = stack.pop()

            # Detect and skip cycles
            if current_node in visited_nodes:
                continue  # Skip nodes already processed to avoid cycles
            visited_nodes.add(current_node)

            # Add the current node's relationship to the sorted output if it has a parent
            if parent is not None and count is not None:
                sorted_output.append(f"{'  ' * int(depth - 1)}{parent} [{count}] {current_node}")

            # Process children if they exist in the hierarchy
            if current_node in tree:
                # Sort children by count in descending order
                children = sorted(tree[current_node], key=lambda x: x[1], reverse=True)
                # Push children onto the stack in reverse order to process in correct order
                for child, child_count in reversed(children):
                    stack.append((child, depth + 1, current_node, child_count))

        return sorted_output

    # Process each line to parse and collect data
    sankey_data = [parse_line(line.strip()) for line in sankeymatic_code if line.strip()]
    hierarchy = build_hierarchy(sankey_data)  # Build hierarchy from data
    sorted_sankeymatic_code = sort_hierarchy(hierarchy)  # Sort and prepare the final output

    # Write the sorted SankeyMATIC code to a file
    if group_threshold == 0:
        sankey_filename = f"{os.path.splitext(input_file)[0]}.sankey.txt"        
    else:
        sankey_filename = f"{os.path.splitext(input_file)[0]}.sankey-group{group_threshold}.txt"
    with open(sankey_filename, "w") as file:
        for line in sorted_sankeymatic_code:
            file.write(line + "\n")
    print(f"SankeyMATIC code has been saved to '{sankey_filename}'.\n")

def main():
    # Parse command-line arguments
    args = parse_args()

    # Write or read configuration file
    setup_config(args.email)

    # Automatically detect input file delimiter
    delimiter = detect_delimiter(args.input_file)

    # Fetch and append taxonomic information from Entrez database to TaxIDs
    tax_df, header_prefix = append_taxonomies(args.input_file, args.header, args.tax_ranks, delimiter)

    # Run SankeyMATIC code generation if not skipped
    if not args.skip:
        sankeymatic_code = generate_sankeymatic_code(tax_df, header_prefix, args.tax_ranks, args.group)
        sort_and_output_sankeymatic_code(args.input_file, args.group, sankeymatic_code)

if __name__ == "__main__":
    main()
