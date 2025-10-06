import find_best_pw_for_compound as bpw

if __name__ == "__main__":
    map = "map00720" # Other carbon fixation pathways
    compound = "C00022" # Pyruvate
    maps = bpw.find_maps_from_compound(compound)
    #for element in maps:
         #bpw.visualize_best_subpathway(map, compound) 
    bpw.visualize_best_subpathway(map, compound) 
    # Very slow but actually displays the highlighted subpathway on the KEGG website

