import re, random
import pandas as pd, numpy as np

def get_which_kins(files):
    unq_kins = set()
    for f in files:
        for uk in pd.read_csv(f)['lab_name'].to_list():
            unq_kins.add(uk)
    return unq_kins

def get_motifs_and_matrices(which_kins):
    if not isinstance(which_kins, list):
        which_kins = list(which_kins)
    uniprot_to_fam = pd.read_csv("uniprot_to_family.csv").set_index('Uniprot').to_dict()['Family']
    kin_to_uniprot = pd.read_csv("kin_to_uniprot.csv").set_index('GENE').to_dict()['KIN_ACC_ID']
    for k in kin_to_uniprot:
        kin_to_uniprot[k] = re.sub("(.*)-[0-9]+", r"\1", kin_to_uniprot[k])

    uniprot_to_kin = {y:x for x, y in kin_to_uniprot.items()}

    kin_to_fam = {}
    not_found = []
    for wk in which_kins:
        uniprot = kin_to_uniprot[wk]
        if uniprot not in uniprot_to_fam:
            print("Not found --", wk, uniprot)
            not_found.append((wk, uniprot))
        else:
            kin_to_fam[uniprot_to_kin[uniprot]] = uniprot_to_fam[uniprot].upper()

    def one_hot_it(arr, pos): 
        if not isinstance(pos, int) and isinstance(pos, int):
            pos = list(pos)
        arr[pos] = 1
        return arr

    fams = list(set(kin_to_fam.values()))
    with open("data/large_fams.csv", "w") as f:
        f.write("\n".join(fams))
    fam_to_vec = {f: str(one_hot_it(np.zeros(len(fams), dtype = int), [i]).tolist()) for i, f in enumerate(fams)}
    
    mapping_df = pd.DataFrame({"Kinase": which_kins, "Family": [kin_to_fam[wk] if wk in kin_to_fam else None for wk in which_kins], "Vector": [fam_to_vec[kin_to_fam[wk]] if wk in kin_to_fam else None for wk in which_kins]}).sort_values(by = ["Vector", "Kinase"], ascending=[False, True])
    mapping_df.to_csv("kin_to_vec.csv", index=False)

def process_from_formatted(infile):
    kin_to_vec = pd.read_csv("kin_to_vec.csv")[['Kinase', 'Vector']]
    kin_to_vec_dict = kin_to_vec.set_index('Kinase').to_dict()['Vector']
    indata_full = pd.read_csv(infile)
    indata = indata_full['orig_lab_name'].tolist()
    outdata_vec = [eval(kin_to_vec_dict[x]) if not (isinstance(kin_to_vec_dict[x], float) and np.isnan(kin_to_vec_dict[x])) else [float("nan")] for x in indata]
    outdata_vec = [",".join(map(str, x)) for x in outdata_vec]
    outdata_seq = indata_full['seq'].tolist()

    while "nan" in outdata_vec:
        idx = outdata_vec.index('nan')
        del outdata_vec[idx]
        del outdata_seq[idx]

    assert len(outdata_seq) == len(outdata_vec), "Seqs and Vecs are not the same length."

    shuffle_inds = list(range(len(outdata_seq)))
    random.seed(0)
    random.shuffle(shuffle_inds)

    outdata_seq = [outdata_seq[i] for i in shuffle_inds]
    outdata_vec = [outdata_vec[i] for i in shuffle_inds]

    with open(f"data/{len(outdata_seq)}_motifs.csv", "w") as f:
        f.write("\n".join(outdata_seq))
    with open(f"data/{len(outdata_vec)}_motifxFamMatrix.csv", "w") as f:
        f.write("\n".join(outdata_vec))

if __name__ == "__main__":
    infiles = ["raw_data_1647_no_overlaps_formatted_95.csv", "raw_data_7530_no_overlaps_formatted_50.csv"]
    which = get_which_kins(infiles)
    get_motifs_and_matrices(which)
    for infile in infiles:
        process_from_formatted(infile)