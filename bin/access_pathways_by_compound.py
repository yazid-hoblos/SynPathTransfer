import requests

def find_kegg_pathways_by_compound_name(compound_name):
    # Step 1: Search compounds
    search_url = f"http://rest.kegg.jp/find/compound/{compound_name}"
    response = requests.get(search_url)
    if response.status_code != 200 or not response.text.strip():
        return []
    
    compounds = [line.split('\t')[0] for line in response.text.strip().split('\n')]
    
    # Step 2: Get pathways for each compound
    pathways = set()
    for cid in compounds:
        link_url = f"http://rest.kegg.jp/link/pathway/{cid}"
        link_response = requests.get(link_url)
        if link_response.status_code != 200:
            continue
        for line in link_response.text.strip().split('\n'):
            _, pid = line.split('\t')
            pathways.add(pid)
    
    return list(pathways)

compound_name = input("Enter compound or abbreviation: ") 
pathways = find_kegg_pathways_by_compound_name(compound_name)
print("Associated KEGG pathways:", pathways)
