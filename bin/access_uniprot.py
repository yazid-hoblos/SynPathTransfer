import os
import requests
import csv
import time
import sys

if len(sys.argv) != 3:
    print("Usage: python access_uniprot.py <ec_file> <output_file>")
    sys.exit(1)

ec_file = sys.argv[1]
output_file = sys.argv[2]

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/"

with open(ec_file, 'r') as f:
    ec_numbers = [line.strip() for line in f if line.strip()]

with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["EC_number", "UniProt_ID", "Protein_Name", "Pfam_ID", "Pfam_Description"])

    for ec in ec_numbers:
        print(f"\n🔍 EC {ec}")
        params = {
            "query": f"ec:{ec} AND reviewed:true",
            "format": "json",
            "fields": "accession,protein_name",
            "size": "1"
        }

        try:
            r = requests.get(UNIPROT_SEARCH, params=params)
            r.raise_for_status()
            results = r.json().get("results", [])

            if not results:
                print(f"❌ No reviewed entry found.")
                continue

            entry = results[0]
            accession = entry["primaryAccession"]
            protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")
            print(f"✅ {accession} - {protein_name}")

            # Get full UniProt record to extract domains
            full_entry = requests.get(f"{UNIPROT_ENTRY}{accession}.json").json()

            features = full_entry.get("uniProtKBCrossReferences", [])
            pfam_domains = []
            for feature in features:
                if feature.get("database") == "Pfam":
                    pfam_id = feature.get("id")
                    pfam_domains.append(pfam_id)

                    # Write to CSV
                    writer.writerow([ec, accession, protein_name, pfam_id])
                    
            print("Pfam domains found:", pfam_domains)
        except Exception as e:
            print(f"🚨 Error with EC {ec}: {e}")

        time.sleep(1)  # gentle on the server
