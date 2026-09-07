#!/usr/bin/env python3

import re
import sys

if len(sys.argv) != 4:
    print("Usage: append_go_names.py go-basic.obo enrichment.txt output.txt")
    sys.exit(1)

obo_file = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

# -----------------------------
# Read GO ontology
# -----------------------------
go_info = {}

with open(obo_file) as f:
    go_id = None
    name = None
    namespace = None
    definition = None

    for line in f:
        line = line.rstrip()

        if line == "[Term]":
            if go_id:
                go_info[go_id] = {
                    "name": name or "",
                    "namespace": namespace or "",
                    "definition": definition or ""
                }

            go_id = None
            name = None
            namespace = None
            definition = None

            continue

        if line.startswith("id: GO:"):
            go_id = line.split("id: ")[1]

        elif line.startswith("name: "):
            name = line.split("name: ", 1)[1]

        elif line.startswith("namespace: "):
            namespace = line.split("namespace: ", 1)[1]

        elif line.startswith("def: "):
            m = re.search(r'"([^"]+)"', line)
            if m:
                definition = m.group(1)

# Save final term
if go_id:
    go_info[go_id] = {
        "name": name or "",
        "namespace": namespace or "",
        "definition": definition or ""
    }

# -----------------------------
# Append information
# -----------------------------
with open(infile) as fin, open(outfile, "w") as fout:

    fout.write(
        "GO_ID\tName\tNamespace\tDefinition\tResults\n"
    )

    for line in fin:

        line = line.rstrip()

        if not line:
            continue

        m = re.match(r'(GO:\d+)', line)

        if not m:
            continue

        go = m.group(1)

        info = go_info.get(
            go,
            {
                "name": "",
                "namespace": "",
                "definition": ""
            }
        )

        fout.write(
            f"{go}\t"
            f"{info['name']}\t"
            f"{info['namespace']}\t"
            f"{info['definition']}\t"
            f"{line}\n"
        )
