import requests
import webbrowser
import cost

def find_reactions_from_compound(map: str, compound: str) -> dict:
    """
        Takes map and compounds as parameters.
        Returns the map with every reaction using the given compound as start product
    """
    # Getting every reactions involving the given compound
    reactions_c = "https://rest.kegg.jp/link/reaction/{}".format(compound)
    # Getting every reactions of the given pathway
    reactions = "https://rest.kegg.jp/link/reaction/{}".format(map)

    # Transforming these texts into lists so we can use them later
    rc = requests.get(reactions_c, headers=None)
    all_reactions_c = rc.text.replace(f"cpd:{compound}", "").replace("\t", "").replace("rn:", "").split("\n")
    rs = requests.get(reactions, headers=None)
    all_reactions = rs.text.replace(f"path:{map}", "").replace("\t", "").replace("rn:", "").split("\n")

    # Looking for every reaction of the pathway involving the compound
    map_reactions_c = []
    for element in all_reactions_c:
        if element in all_reactions:
            map_reactions_c.append(element)
    if "" in map_reactions_c:
        map_reactions_c.remove("")

    # Checking if the compound is actually a substrate or not
    for reaction in map_reactions_c:
        check_start = "https://rest.kegg.jp/get/reaction:{}".format(reaction)
        raw = requests.get(check_start).text
        for info in raw.splitlines():
            if info.startswith("EQUATION"):
                equation = info.replace("EQUATION", "").strip()
                if "=>" in equation:
                    products_part = equation.split("=>")  # I'm expecting a non reversible reaction, even though I didn't see one yet
                else:
                    products_part = equation.split("<=>")
        if compound in products_part[1]:
            map_reactions_c.remove(reaction)

    return {
        'compound': compound,
        'reactions': map_reactions_c,
        'map_reactions': all_reactions
    }


def clean_equation(eq: str):
    cofactors = ["ATP", "ADP", "AMP", "Pi", "PPi", "NAD+", "NADH", "NADP+", "NADPH", "H2O", "H+", "CO2", "CoA", "NH"]
    for c in cofactors:
        eq = eq.replace(c, "")
    return eq.replace("  ", " ").strip()


def parse_equation(equation: str):
    if "=>" in equation:
        sub, prod = equation.split("=>")
        sign = -1
    elif "<=>" in equation:
        sub, prod = equation.split("<=>")
        sign = 1
    else:
        return set(), set(), 0
    sub = {x.strip() for x in sub.split("+") if x.strip()}
    prod = {x.strip() for x in prod.split("+") if x.strip()}
    return sub, prod, sign


def load_equations(all_reactions):
    """
        Loads all reactions at once returns dict {reaction: (subs, prods, sign)}
        Without this, the amount of time taken to analyze is longer
    """
    equations = {}
    for r in all_reactions:
        try:
            raw = requests.get(f"https://rest.kegg.jp/list/reaction:{r}").text
            if ";" not in raw:
                equations[r] = (set(), set(), 0)
                continue
            eq = raw.split(";")[1].strip()
            clean = clean_equation(eq)
            equations[r] = parse_equation(clean)
        except Exception:
            equations[r] = (set(), set(), 0)
    return equations


def find_subpathways_from_reac(start_reac, all_reactions, equations):
    """
        Returns every reaction linked to start_reac
    """
    visited = set()
    stack = [start_reac]
    pathway = []

    while stack:
        reac = stack.pop()
        if reac in visited:
            continue
        visited.add(reac)
        sub, prod, sign = equations[reac]
        pathway.append((reac, sign))

        for r in list(all_reactions):
            if r in visited:
                continue
            sub_next, prod_next, _ = equations[r]
            if prod & sub_next or prod_next & sub:
                stack.append(r)

    return pathway


def find_all_subpathways(pathway, compound):
    """
        Explore a KEGG pathway and return all possible subpathways
            starting from each reaction involving the given compound
    """
    f = find_reactions_from_compound(pathway, compound)
    all_reac = f["map_reactions"]
    equations = load_equations(all_reac)

    all_subpaths = []
    for reac in all_reac:
        subpaths = find_subpathways_from_reac(reac, all_reac, equations)
        if len(subpaths) > 1:
            all_subpaths.append(subpaths)
    return all_subpaths


def best_pw(subpaths):
    """
        Explore a list of subpathways and returns the best one (when it comes to cost)
    """
    best_pw = None
    min_cost = 0
    for element in subpaths:
        total, step = cost.subpathway_cost_relative(element)
        if total < min_cost:
            min_cost = total
            best_pw = element
    
    reactions = ""
    for reaction in best_pw:
        reactions += f"{reaction}/"
    return reactions


def visualize_best_subpathway(pathway, compound):
    """
        Finds all subpathways for a given KEGG pathway and compound,
            selects the best one, and opens it in the KEGG Website
    """
    all_subpaths = find_all_subpathways(pathway, compound)
    reactions = best_pw(all_subpaths)
    url = f"https://www.kegg.jp/kegg-bin/show_pathway?{pathway}/{reactions}/{compound}%20%23ff0000"
    webbrowser.open(url)


# Example usage
map = "map00720"
compound = "C00022"
visualize_best_subpathway(map, compound)
