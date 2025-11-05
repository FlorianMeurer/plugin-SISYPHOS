import csv
import os

def parse_block(block):
    """
    Parse a single block of text into a dictionary.
    The block is expected to contain lines with either a simple key:value pair or a composite field
    where the value contains semicolon-separated sub key:value pairs.
    """
    record = {}
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue
        if ":" not in line:
            continue
        key, rest = line.split(":", 1)
        key = key.strip()
        rest = rest.strip()
        if ";" in rest:
            # composite field: value contains semicolon-separated sub key:val pairs
            if rest == "":
                record[key] = ""
            else:
                sub_pairs = rest.split(";")
                for sub in sub_pairs:
                    sub = sub.strip()
                    if not sub:
                        continue
                    if ":" in sub:
                        sub_key, sub_val = sub.split(":", 1)
                        # include parent key to avoid collisions between similarly-named subkeys
                        combined_key = f"{sub_key.strip()}_{key}"
                        record[combined_key] = sub_val.strip()
                    else:
                        continue
        else:
            record[key] = rest
    return record

def main(inputloc, outputloc):
    """
    Generate a .csv file from a SYSout.txt output file from a SISYPHOS run.
    """
    with open(inputloc, "r") as infile:
        content = infile.read()

    blocks = content.split("+++++++++++++++++++")
    records = []
    for block in blocks:
        block = block.strip()
        if not block:
            continue
        record = parse_block(block)
        records.append(record)

    all_keys = set()
    for rec in records:
        all_keys.update(rec.keys())
    all_keys = sorted(all_keys)  
    out = os.path.join(outputloc, 'SYSoutput.csv')
    # choose delimiter used in CSV (semicolon by default)
    delimiter = ';'
    with open(out, "w", newline="", encoding='utf-8') as csvfile:
        # Write Excel-compatible separator declaration on the first line so
        # applications (Excel in many locales) can detect the delimiter.
        # Example: "sep=;" tells Excel that fields are separated by semicolons.
        csvfile.write(f"sep={delimiter}\n")
        writer = csv.DictWriter(csvfile, fieldnames=all_keys, delimiter=delimiter)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)

if __name__ == '__main__':
    main()
    
