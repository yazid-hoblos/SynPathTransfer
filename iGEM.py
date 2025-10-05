import requests

#INITIALIZATION
ec_num=input("What is the EC number of the reaction you want to implement ? Example: 1.1.1.1 \n ec:")


#org_prompt=input("""
                 #Give me the list of organisms you want this reaction to be compatible with.
                 #Copy-paste the text from the GOLD data base as follows:
                 #"Go0678575	Polycladomyces zharkentensis ZKZ2	BACILLOTA	Aerobe	55
                #Go0672758	Metallosphaera javensis J1	THERMOPROTEOTA	Aerobe	55
                #Go0357320	Melghirimyces algeriensis DSM 45474	BACILLOTA	Aerobe	55
                #Go0357291	Thermoflavifilum aggregans DSM 27268	BACTEROIDOTA	Aerobe	60"
                 #""")

org_prompt="""Go0678575	Polycladomyces zharkentensis ZKZ2	BACILLOTA	Aerobe	55
Go0672758	Metallosphaera javensis J1	THERMOPROTEOTA	Aerobe	55
Go0615179	Thermorudis pharmacophila Wkt50.2	THERMOMICROBIOTA	Aerobe	50
Go0579817	Xylanibacillus composti K13	BACILLOTA	Aerobe	45
Go0543959	Thermaurantimonas aggregans LA	BACTEROIDOTA	Aerobe	50
Go0511066	Microvirga flocculans DSM 15743	PSEUDOMONADOTA	Aerobe	43
Go0511007	Anoxybacillus voinovskiensis DSM 17075	BACILLOTA	Aerobe	45
Go0511006	Amphiplicatus metriothermophilus DSM 105738	PSEUDOMONADOTA	Aerobe	50
Go0456155	Aeribacillus composti KCTC 33824	BACILLOTA	Obligate aerobe	50
Go0357320	Melghirimyces algeriensis DSM 45474	BACILLOTA	Aerobe	55
Go0357291	Thermoflavifilum aggregans DSM 27268	BACTEROIDOTA	Aerobe	60
Go0143341	Halobiforma haloterrestris DSM 13078	EURYARCHAEOTA	Aerobe	42
Go0120274	Alicyclobacillus montanus USBA-503	BACILLOTA	Aerobe	45
Go0103850	Thermoactinomyces sp. Gus2-1	BACILLOTA	Aerobe	60
Go0091374	Limisphaera ngatamarikiensis NGM72.4	VERRUCOMICROBIOTA	Aerobe	60
Go0088733	Thermosporothrix hazakensis ATCC BAA-1881	CHLOROFLEXOTA	Aerobe	50
Go0050713	Thermoflavifilum aggregans P373	BACTEROIDOTA	Obligate aerobe	60
Go0046751	Aigarchaeota archaeon OS1	NITROSOSPHAEROTA	Aerobe	80
Go0036963	Natronorubrum sulfidifaciens JCM 14089	EURYARCHAEOTA	Aerobe	44
Go0036947	Natrinema limicola JCM 13563	EURYARCHAEOTA	Aerobe	45
Go0036933	Natrialba chahannaoensis JCM 10990	EURYARCHAEOTA	Aerobe	50
Go0032038	Sulfolobus acidocaldarius 98-3	THERMOPROTEOTA	Obligate aerobe	70
Go0029898	Rhodothermus marinus JCM 9785	RHODOTHERMOTA	Aerobe	80
Go0025906	Natronorubrum tibetense GA33	EURYARCHAEOTA	Aerobe	45
Go0024539	Methylothermalis aethiopiae gen. nov. sp. nov. LS7-MT	PSEUDOMONADOTA	Aerobe	50
Go0020890	Thermomyces lanuginosus ATCC 200065	FUNGI-ASCOMYCOTA	Aerobe	50
Go0019017	Thermus igniterrae ATCC 700962	DEINOCOCCOTA	Aerobe	65
Go0018176	Thermomicrobium sp. HL1	THERMOMICROBIOTA	Aerobe	70
Go0015989	Thermoflexus hugenholtzii JAD2	CHLOROFLEXOTA	Obligate aerobe	70
Go0014275	Alicyclobacillus acidocaldarius HSU9	BACILLOTA	Aerobe	50
Go0014145	Thermus oshimai DSM 12092	DEINOCOCCOTA	Aerobe	70
Go0014143	Thermus igniterrae DSM 12459	DEINOCOCCOTA	Aerobe	65
Go0014119	Hydrogenobaculum acidophilum 3H-1, DSM 11251	AQUIFICOTA	Aerobe	65
Go0014019	Alicyclobacillus herbarius DSM 13609	BACILLOTA	Aerobe	55
Go0014018	Alicyclobacillus contaminans DSM 17975	BACILLOTA	Aerobe	50
Go0013801	Actinopolyspora mortivallis HS-1, DSM 44261	ACTINOMYCETOTA	Aerobe	45
Go0013583	Caldibacillus debilis DSM 16016	BACILLOTA	Aerobe	70
Go0013491	Ferrithrix thermotolerans DSM 19514	ACTINOMYCETOTA	Aerobe	43
Go0013444	Thermus sp. CCB_US3_UF1	DEINOCOCCOTA	Aerobe	92.4
Go0013280	Silanimonas lenta DSM 16282	PSEUDOMONADOTA	Aerobe	47
Go0013233	Rubritepida flocculans DSM 14296	PSEUDOMONADOTA	Aerobe	50
Go0013230	Rubellimicrobium thermophilum DSM 16684	PSEUDOMONADOTA	Aerobe	45
Go0013206	Pseudoxanthomonas taiwanensis DSM 22914	PSEUDOMONADOTA	Aerobe	50
Go0013194	Pseudomonas thermotolerans DSM 14292	PSEUDOMONADOTA	Aerobe	47
Go0013146	Erythrobacter cryptus DSM 12079	PSEUDOMONADOTA	Aerobe	50
Go0013138	Picrophilus oshimae DSM 9789	EURYARCHAEOTA	Aerobe	60
Go0013035	Calidithermus timidus DSM 17022	DEINOCOCCOTA	Aerobe	55
Go0013034	Meiothermus taiwanensis DSM 14542	DEINOCOCCOTA	Aerobe	55
Go0013033	Meiothermus rufus DSM 22234	DEINOCOCCOTA	Aerobe	55
Go0013015	Elioraea tepidiphila DSM 17972	PSEUDOMONADOTA	Aerobe	45
Go0012925	Alicyclobacillus pomorum DSM 14955	BACILLOTA	Aerobe	45
Go0012543	Alicyclobacillus acidocaldarius acidocaldarius Tc-4-1	BACILLOTA	Aerobe	70
Go0011220	Thermocrispum agreste DSM 44070	ACTINOMYCETOTA	Aerobe	45
Go0011201	Natronorubrum tibetense DSM 13204	EURYARCHAEOTA	Aerobe	45
Go0011166	Thermocrispum municipale DSM 44069	ACTINOMYCETOTA	Aerobe	45
Go0011135	Brevibacillus laterosporus NRS 682, LMG 15441	BACILLOTA	Aerobe	51
Go0011067	Hydrogenobaculum sp. SO	AQUIFICOTA	Aerobe	58
Go0010487	Bacillus methanolicus PB1	BACILLOTA	Aerobe	45
Go0010387	Caldalkalibacillus thermarum TA2.A1	BACILLOTA	Obligate aerobe	70
Go0008893	Hydrogenobaculum sp. SHO	AQUIFICOTA	Aerobe	58
Go0008665	Halobacterium sp. DL1	EURYARCHAEOTA	Aerobe	42
Go0008470	Hydrogenobaculum sp. 3684	AQUIFICOTA	Aerobe	58
Go0007990	Thermus thermophilus SG0.5JP17-16	DEINOCOCCOTA	Obligate aerobe	70
Go0007944	Saccharomonospora glauca K62, DSM 43769	ACTINOMYCETOTA	Aerobe	45
Go0007397	Thermus thermophilus JL-18	DEINOCOCCOTA	Obligate aerobe	70
Go0007395	Microvirga lotononidis WSM3557	PSEUDOMONADOTA	Aerobe	41
Go0007234	Bacillus methanolicus MGA3	BACILLOTA	Aerobe	55
Go0007190	Halorubrum arcis AJ201	EURYARCHAEOTA	Aerobe	42
Go0007187	Haloferax sp. D1227, ATCC 51408	EURYARCHAEOTA	Aerobe	45
Go0007164	Natrinema saccharevitans AB14	EURYARCHAEOTA	Aerobe	42
Go0007163	Haloterrigena longa JCM 13562	EURYARCHAEOTA	Aerobe	45
Go0007140	Natronorubrum aibiense JCM 13488	EURYARCHAEOTA	Aerobe	45
Go0007113	Halobiforma haloterrestris JCM 11627	EURYARCHAEOTA	Aerobe	42
Go0007112	Natrialba hulunbeirensis JCM 10989	EURYARCHAEOTA	Aerobe	50
Go0006841	Natronorubrum tibetense JCM 10636	EURYARCHAEOTA	Aerobe	45
Go0006582	Hydrogenobaculum sp. HO	AQUIFICOTA	Aerobe	58
Go0006309	Marinithermus hydrothermalis T1, DSM 14884	DEINOCOCCOTA	Obligate aerobe	67.5
Go0006278	Brevibacillus sp. JXL	BACILLOTA	Aerobe	57
Go0005858	Sphaerobacter thermophilus 4ac11, DSM 20745	THERMOMICROBIOTA	Obligate aerobe	55
Go0005831	Alicyclobacillus sp. A4	BACILLOTA	Aerobe	57
Go0005502	Haloferax volcanii DS2	EURYARCHAEOTA	Aerobe	42
Go0005483	Hydrogenobaculum sp. SN	AQUIFICOTA	Aerobe	58
Go0005159	Thermaerobacter subterraneus C21, DSM 13965	BACILLOTA	Obligate aerobe	70
Go0005139	Thermobacillus composti KWC4, DSM 18247	BACILLOTA	Aerobe	50
Go0005109	Thermus brockianus	DEINOCOCCOTA	Obligate aerobe	70
Go0005097	Candidatus Nitrososphaera gargensis Ga9-2	NITROSOSPHAEROTA	Aerobe	50
Go0004668	Chthonomonas calidirosea T49, DSM 23976	ARMATIMONADOTA	Aerobe	68
Go0004365	Thermaerobacter marianensis 7p75a, DSM 12885	BACILLOTA	Obligate aerobe	75
Go0004361	Rubrobacter radiotolerans DSM 5868	ACTINOMYCETOTA	Aerobe	60
Go0004276	Saccharolobus solfataricus 98/2	THERMOPROTEOTA	Obligate aerobe	85
Go0003696	Ruegeria lacuscaerulensis ITI-1157	PSEUDOMONADOTA	Aerobe	45
Go0003690	Truepera radiovictrix RQ-24, DSM 17093	DEINOCOCCOTA	Obligate aerobe	50
Go0003649	Pyrolobus fumarii 1A, DSM 11204	THERMOPROTEOTA	Aerobe	106
Go0003483	Saccharomonospora viridis P101, DSM 43017	ACTINOMYCETOTA	Aerobe	55
Go0003405	Methylacidiphilum infernorum V4	VERRUCOMICROBIOTA	Aerobe	60
Go0003221	Hydrogenobacter thermophilus TK-6, DSM 6534	AQUIFICOTA	Aerobe	70
Go0002339	Sulfolobus islandicus L.D.8.5	THERMOPROTEOTA	Aerobe	75
Go0002338	Sulfolobus islandicus L.S.2.15	THERMOPROTEOTA	Aerobe	75
Go0002336	Sulfolobus islandicus M.16.4	THERMOPROTEOTA	Aerobe	75
Go0002335	Sulfolobus islandicus M.16.27	THERMOPROTEOTA	Aerobe	75
Go0002334	Sulfolobus islandicus M.14.25	THERMOPROTEOTA	Aerobe	75
Go0002300	Halobacterium salinarum R1	EURYARCHAEOTA	Obligate aerobe	50
Go0002120	Thermus thermophilus	DEINOCOCCOTA	Aerobe	85
Go0002108	Thermomicrobium roseum DSM 5159	THERMOMICROBIOTA	Aerobe	70
Go0002007	Alicyclobacillus acidocaldarius LAA1	BACILLOTA	Aerobe	57
Go0001773	Hydrogenobaculum sp. Y04AAS1	AQUIFICOTA	Aerobe	58
Go0001225	Halorhabdus utahensis AX-2, DSM 12940	EURYARCHAEOTA	Aerobe	50
Go0001214	Thermomonospora curvata DSM 43183	ACTINOMYCETOTA	Aerobe	65
Go0001213	Thermobispora bispora R51, DSM 43833	ACTINOMYCETOTA	Aerobe	55
Go0001192	Kyrpidia tusciae T2, DSM 2912	BACILLOTA	Aerobe	50
Go0001191	Acidimicrobium ferrooxidans ICP, DSM 10331	ACTINOMYCETOTA	Aerobe	48
Go0001190	Alicyclobacillus acidocaldarius acidocaldarius 104-IA, DSM 446	BACILLOTA	Aerobe	60
Go0001179	Allomeiothermus silvanus DSM 9946	DEINOCOCCOTA	Aerobe	55
Go0001173	Rhodothermus marinus R-10, DSM 4252	RHODOTHERMOTA	Aerobe	65
Go0001170	Meiothermus ruber 21, DSM 1279	DEINOCOCCOTA	Obligate aerobe	50
Go0001165	Thermobaculum terrenum YNP1, ATCC BAA-798	CHLOROFLEXOTA	Obligate aerobe	65
Go0001163	Isosphaera pallida IS1B, ATCC 43644	PLANCTOMYCETOTA	Obligate aerobe	45
Go0000767	Aeropyrum pernix K1	THERMOPROTEOTA	Aerobe	70
Go0000749	Halobacterium salinarum NRC-1	EURYARCHAEOTA	Aerobe	42
Go0000728	Saccharolobus solfataricus P2	THERMOPROTEOTA	Obligate aerobe	50
Go0000723	Sulfurisphaera tokodaii 7	THERMOPROTEOTA	Aerobe	80
Go0000624	Chloracidobacterium thermophilum B	ACIDOBACTERIOTA	Aerobe	53
Go0000583	Thermus thermophilus HB27	DEINOCOCCOTA	Aerobe	68
Go0000574	Picrophilus torridus DSM 9790	EURYARCHAEOTA	Aerobe	60
Go0000544	Methylococcus capsulatus Bath	PSEUDOMONADOTA	Aerobe	45
Go0000523	Thermus thermophilus HB8	DEINOCOCCOTA	Aerobe	85
Go0000519	Geobacillus kaustophilus HTA426	BACILLOTA	Aerobe	60
Go0000473	Sulfolobus acidocaldarius 98-3	THERMOPROTEOTA	Obligate aerobe	70
Go0000457	Thermobifida fusca YX	ACTINOMYCETOTA	Aerobe	50
Go0000350	Deinococcus geothermalis DSM 11300	DEINOCOCCOTA	Aerobe	47
Go0000338	Rubrobacter xylanophilus DSM 9941	ACTINOMYCETOTA	Aerobe	60
Go0000258	Acidothermus cellulolyticus 11B	ACTINOMYCETOTA	Aerobe	58
Go0000163	Metallosphaera sedula DSM 5348	THERMOPROTEOTA	Aerobe	70
Go0000055	Candidatus Nitrosocaldus yellowstonensis HL72	NITROSOSPHAEROTA	Aerobe	72
Go0000025	Caldivirga maquilingensis IC-167	THERMOPROTEOTA	Aerobe	85
Go0000004	Thermus aquaticus Y51MC23	DEINOCOCCOTA	Obligate aerobe	70"""

selected_genes={}

#TRANSFORM PASTED TEXT INTO LIST OF DESIRED ORGANISMS
def list_organisms(prompt):
    desired_org={}
    prompt_lines=str(prompt).split("\n")
    for line in prompt_lines:

        line=line[line.index("\t")+1:]
        line=line[:line.index("\t")]
        desired_org[line]=True

    return desired_org

desired_org=list_organisms(org_prompt)

#GET KO NUMBERS CORRESPONDING TO THE EC NUMBER
ko_url = "https://rest.kegg.jp/link/ko/ec:"+ec_num
ko_resp= requests.get(ko_url).text
ko_lines=ko_resp.split("\n")
ko_list=[]

for line in ko_lines[:-1]:

    i=line.index("ko:")
    ko_list.append(line[i:i+9])

print("KO numbers associated with this ec number:",ko_list)

#FETCH ORGANISM NAME CORRESPONDING TO GENE PREFIX 
prefix_to_org={}

org_url="https://rest.kegg.jp/list/organism"
org_text=requests.get(org_url).text
org_list=org_text.split("\n")

for line in org_list:
    split_line=line.split("\t")
    if len(split_line)>2:
        prefix,org_name=split_line[1],split_line[2]
        prefix_to_org[prefix]=org_name

#GET GENES CORRESPONDING TO EACH KO NUMBER FOUND

for ko in ko_list:

    genes_url = "https://rest.kegg.jp/link/genes/"+ko
    gene_resp=requests.get(genes_url).text
    gene_text = gene_resp.split("\n")

    for gene_line in gene_text[:-1]:

        i=gene_line.index("\t")+1
        new_gene=gene_line[i:]
        #print("new_gene:"+new_gene)

        #CHECK WHETHER THIS NEW GENE COMES FROM A DESIRED ORGANISM
        new_gene_prefix=new_gene[:new_gene.index(":")]
        new_gene_org=prefix_to_org[new_gene_prefix]
        
        if new_gene_org in desired_org:
            selected_genes[new_gene]=new_gene_org

        #print(selected_genes)

print("Genes that fulfill your requirements:\n",list(selected_genes.keys()))
#print(desired_org)