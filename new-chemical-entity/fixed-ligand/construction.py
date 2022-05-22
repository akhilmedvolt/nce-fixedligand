import pandas as pd

import random
import itertools

import rdkit.Chem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import BRICS
import rdkit.Chem.QED


def fragment_construction():
    results_df = pd.read_csv('results.txt', sep="\t", names=['Frag', 'Smiles', 'Score'])
    smiles = results_df['Smiles'].tolist()
    frag = results_df['Frag'].tolist()
    score = results_df['Score'].tolist()

    new_smiles = []
    for x in smiles:
        x = x.replace("()", "")
        new_smiles.append(x)
    smiles = new_smiles

    rings = []
    wgt = []
    smiles_n = []
    frag_n = []
    score_n = []

    k = 0
    while k < len(smiles):
        if '+' not in str(smiles[k]):
            if "[O-]" not in str(smiles[k]):
                try:
                    mol = Chem.MolFromSmiles(smiles[k])

                    ri = mol.GetRingInfo()
                    nu = ri.NumRings()
                    rings.append(nu)

                    wgt.append(Descriptors.MolWt(mol))
                    smiles_n.append(smiles[k])
                    frag_n.append(frag[k])
                    score_n.append(score[k])

                    print(smiles[k])
                except:
                    pass
        k += 1

        new_df = pd.DataFrame(list(zip(smiles_n, rings, wgt, frag_n, score_n)),
                              columns=['smiles', 'Number of Rings', 'Mol.Weight', 'Frag', 'Score'])
        new_df = new_df.sort_values('Score', ascending=True)

        # user parameters
        df1 = new_df[new_df['Number of Rings'] <= 2]
        # df1 = df1[df1['Score'] >= -5]
        df1 = df1[df1['Mol.Weight'] <= 180]
        fraglist = df1['Frag'].tolist()

        g_dict = {}
        for x in range(1, 26):
            g_dict[str(x)] = []
        for g in fraglist:
            if "[1*]" in g:
                g_dict["1"].append(g)

            if "[2*]" in g:
                g_dict["2"].append(g)

            if "[3*]" in g:
                g_dict["3"].append(g)

            if "[4*]" in g:
                g_dict["4"].append(g)

            if "[5*]" in g:
                g_dict["5"].append(g)

            if "[6*]" in g:
                g_dict["6"].append(g)

            if "[7*]" in g:
                g_dict["7"].append(g)

            if "[8*]" in g:
                g_dict["8"].append(g)

            if "[9*]" in g:
                g_dict["9"].append(g)

            if "[10*]" in g:
                g_dict["10"].append(g)

            if "[11*]" in g:
                g_dict["11"].append(g)

            if "[12*]" in g:
                g_dict["12"].append(g)

            if "[13*]" in g:
                g_dict["13"].append(g)

            if "[14*]" in g:
                g_dict["14"].append(g)

            if "[15*]" in g:
                g_dict["15"].append(g)

            if "[16*]" in g:
                g_dict["16"].append(g)

            if "[17*]" in g:
                g_dict["17"].append(g)

            if "[18*]" in g:
                g_dict["18"].append(g)

            if "[19*]" in g:
                g_dict["19"].append(g)

            if "[20*]" in g:
                g_dict["20"].append(g)

            if "[21*]" in g:
                g_dict["21"].append(g)

            if "[22*]" in g:
                g_dict["22"].append(g)

            if "[23*]" in g:
                g_dict["23"].append(g)

            if "[24*]" in g:
                g_dict["24"].append(g)

            if "[25*]" in g:
                g_dict["25"].append(g)

        n_emp = []
        for x in range(1, 26):
            if len(g_dict[str(x)]) != 0:
                n_emp.append(str(x))
        if len(n_emp) > 16:
            n_emp = n_emp[:16]

        new_comb = list(itertools.combinations(n_emp, 4))
        new_comb = [list(x) for x in new_comb]

        for x in new_comb:
            new_frag = []
            new_frag.extend(random.sample(g_dict[x[0]], 3))
            new_frag.extend(random.sample(g_dict[x[1]], 3))
            new_frag.extend(random.sample(g_dict[x[2]], 3))
            new_frag.extend(random.sample(g_dict[x[3]], 3))

            fragms = [Chem.MolFromSmiles(x) for x in new_frag]
            ms = BRICS.BRICSBuild(fragms)
            prods = list(ms)
            new_smi = []

            for prod in prods:
                prod.UpdatePropertyCache(strict=False)

            for prod in prods:
                smi = Chem.MolToSmiles(prod, isomericSmiles=False, canonical=True)
                mol = Chem.MolFromSmiles(smi)
                if smi not in new_smi:
                    pieces = BRICS.BRICSDecompose(mol)
                    res2 = list(pieces)
                    res1 = []
                    for p in res2:
                        if p not in res1:
                            res1.append(p)
                    mol = Chem.MolFromSmiles(smi)
                    ri = mol.GetRingInfo()
                    nu = ri.NumRings()
                    qed = rdkit.Chem.QED.qed(mol)
                    wgt = (Descriptors.MolWt(mol))
                    if (nu == 3) or (nu == 4):
                        bv = 0
                        if wgt < 500:
                            if len(res1) == len(res2):
                                if qed >= 0.5:
                                    new_smi.append(smi)
                                    hi5 = open("molecules.csv", "a+")
                                    hi5.write(smi + '\n')
                                    hi5.close()


