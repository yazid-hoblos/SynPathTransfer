hits_file = 'hits_table.txt'
threshold = 1e-5
hit_ids = []

with open(hits_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        # Format fields from hmmsearch tblout:
        # target name, accession, query name, accession, E-value, score, etc.
        target_name = parts[0]
        evalue = float(parts[4])
        if evalue < threshold:
            hit_ids.append(target_name)

print(f"Found {len(hit_ids)} hits with E-value < {threshold}")
print(hit_ids)

unique_ids = sorted(set(hit_ids))

with open('unique_hits_ids.txt', 'w') as f:
    for uid in unique_ids:
        f.write(uid + '\n')

print(f"Saved {len(unique_ids)} unique IDs to unique_ids.txt")
