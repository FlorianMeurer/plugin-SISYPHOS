import csv
import os

import os
import csv


def parse_block(block):
    """
    Parse a single block of text into a dictionary.
    The block is expected to contain lines with either a simple key:value pair or a composite field
    where the value contains semicolon-separated sub key:value pairs.
    """
    record = {}

    for line in block.splitlines():
        line = line.strip()
        if not line or ":" not in line:
            continue

        key, rest = line.split(":", 1)
        key = key.strip()
        rest = rest.strip()

        # ------------------------------------------------------------------
        # Special handling: anomalous dispersion
        # ------------------------------------------------------------------
        if key.startswith("Refined Disps from"):
            # Decide which column suffix to use
            if "cif" in key.lower():
                # Refined Disps from cif:
                target = "_cif"
            else:
                # Refined Disps from xray structure (or others):
                target = ""

            # There may be multiple entries separated by semicolons
            for sub in rest.split(";"):
                sub = sub.strip()
                if not sub or ":" not in sub:
                    continue

                site_label, sub_val = sub.split(":", 1)
                site_label = site_label.strip()
                sub_val = sub_val.strip()

                # Use FULL site label including numbering (e.g. Mo1, Mo2, Mo1A)
                col_name = f"{site_label}"
                record[col_name] = sub_val

            # Done with this special line
            continue

        # ------------------------------------------------------------------
        # Generic handling
        # ------------------------------------------------------------------
        if ";" in rest:
            # composite field: value contains semicolon-separated sub key:val pairs
            if rest == "":
                record[key] = ""
            else:
                sub_pairs = rest.split(";")
                for sub in sub_pairs:
                    sub = sub.strip()
                    if not sub or ":" not in sub:
                        continue

                    sub_key, sub_val = sub.split(":", 1)
                    sub_key = sub_key.strip()
                    sub_val = sub_val.strip()

                    # Make sure the *category* (key) is never at the beginning:
                    # use subkey_category instead of category_subkey.
                    combined_key = f"{sub_key}_{key}"

                    # Keep values as text (tuples or scalars) for CSV display
                    record[combined_key] = sub_val
        else:
            # Simple key:value (no semicolons)
            record[key] = rest

    return record


def main(inputloc, outputloc):
    """
    Generate a .csv file from a SYSout.txt output file from a SISYPHOS run.
    """
    with open(inputloc, "r", encoding="utf-8") as infile:
        content = infile.read()

    blocks = content.split("+++++++++++++++++++")
    records = []
    for block in blocks:
        block = block.strip()
        if not block:
            continue
        record = parse_block(block)
        records.append(record)

    # Collect all keys used in any record
    all_keys = set()
    for rec in records:
        all_keys.update(rec.keys())
    all_keys = sorted(all_keys)

    out = os.path.join(outputloc, "SYSoutput.csv")
    delimiter = ";"

    with open(out, "w", newline="", encoding="utf-8") as csvfile:
        # Excel-compatible separator declaration
        csvfile.write(f"sep={delimiter}\n")
        writer = csv.DictWriter(csvfile, fieldnames=all_keys, delimiter=delimiter)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


if __name__ == "__main__":
    # Example:
    main("SYSout.txt", ".")
    #main()

    
