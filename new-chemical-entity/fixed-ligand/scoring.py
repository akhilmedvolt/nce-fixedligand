from biopandas.pdb import PandasPdb
from bs4 import BeautifulSoup as bs

import copy
import os

from rdkit import Chem


def fragment_scoring():
    vg = open("zinc_frag_newtry_unduplicated_new.txt", "r")
    fragments = vg.readlines()

    smiles = fragments
    smiles_act = []
    smiles_brics = []

    for smii in smiles:
        mol = Chem.MolFromSmiles(smii)
        if mol is not None:
            smiles_act.append(Chem.MolToSmiles(mol, isomericSmiles=False))
            smiles_brics.append(smii)

    new_smiles = []
    for s in smiles_act:
        a1 = s.replace('[*]', '').replace("\n", "").replace("*", "").replace('()', '')
        new_smiles.append(a1)

    smiles_act = new_smiles

    # User inputs
    PDB = str('7RFW')  # PDB ID
    LIG_ID = str('4WI')  # User input of pocket information required for docking in the form of reference co-crystalized ligand
    smiles = str('C1=CC2=C(C(=C(C=C2Cl)I)O)N=C1')  # dummy smiles required for grid creation and does not come from user

    ppdb = PandasPdb().fetch_pdb(PDB)
    ppdb_orig = copy.deepcopy(ppdb)
    p_chain = ppdb.df['ATOM']['chain_id'].tolist()
    p_chain_id = str(p_chain[0])

    # Saving as protein.pdb
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == p_chain_id]
    ppdbpath = './protein.pdb'
    ppdb.to_pdb(path=ppdbpath,
                records=['ATOM'],
                gz=False,
                append_newline=True)
    protein_pdb = copy.deepcopy(ppdb)
    os.mkdir("receptor")

    prepprot = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -A bonds_hydrogens -U nphs_lps_waters_nonstdres_deleteAltB -r protein.pdb -o finalmodel.pdbqt'
    os.system(prepprot)

    ldb = copy.deepcopy(ppdb)
    ldb.df['HETATM'] = ldb.df['HETATM'][ldb.df['HETATM']['chain_id'] == p_chain_id]
    ldb.df['HETATM'] = ldb.df['HETATM'][ldb.df['HETATM']['residue_name'] == LIG_ID]
    minx = ldb.df['HETATM']['x_coord'].min()
    maxx = ldb.df['HETATM']['x_coord'].max()
    cent_x = round((maxx + minx) / 2, 2)
    size_x = round(abs(maxx - minx) + 12, 2)
    miny = ldb.df['HETATM']['y_coord'].min()
    maxy = ldb.df['HETATM']['y_coord'].max()
    cent_y = round((maxy + miny) / 2, 2)
    size_y = round(abs(maxy - miny) + 12, 2)
    minz = ldb.df['HETATM']['z_coord'].min()
    maxz = ldb.df['HETATM']['z_coord'].max()
    cent_z = round((maxz + minz) / 2, 2)
    size_z = round(abs(maxz - minz) + 12, 2)

    # Preparing for the gpf files
    s1 = 'obabel -:"' + str(smiles) + '" -opdb -O lig.pdb --gen3d'
    os.system(s1)
    os.system('obminimize -ff GAFF lig.pdb > lig_opt.pdb')
    os.system(
        '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -A '
        'bonds_hydrogens -U nphs_lps -l lig_opt.pdb -o lig1.pdbqt')
    gpf2 = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -l ' \
           'lig1.pdbqt -r finalmodel.pdbqt -p npts=75,75,75 -p ligand_types=A,C,HD,F,I,NA,Cl,OA,N,P,S,Br,' \
           'SA -p gridcenter=' + str(int(cent_x)) + ',' + str(int(cent_y)) + ',' + str(int(cent_z)) + ' -y TRUE'
    os.system(gpf2)

    hi1 = open("finalmodel.gpf", "r")
    Lines = hi1.readlines()
    Lines[6] = 'gridcenter ' + str(cent_x) + ' ' + str(cent_y) + ' ' + str(cent_z)
    os.remove('finalmodel.gpf')

    vg1 = open("finalmodel.gpf", "a+")
    for l in Lines:
        vg1.write(l + '\n')
    vg1.close()
    os.system('autogrid4 -p finalmodel.gpf')

    a = 1
    g = 0
    print("Scoring - Docking Begining")
    while g < len(smiles_act):
        try:
            print(g)
            print()
            print(smiles_act[g])
            s1 = 'obabel -:"' + str(smiles_act[g]) + '" -opdb -O ligact.pdb --gen3d'

            os.system(s1)

            os.system('obminimize -ff GAFF ligact.pdb > lig_optact.pdb')

            os.system(
                '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24'
                '/prepare_ligand4.py -A bonds_hydrogens -U nphs_lps -l lig_optact.pdb -o ligact.pdbqt')

            os.system('autodock_gpu_64wi -ffile finalmodel.maps.fld -lfile ligact.pdbqt -nrun 100')

            with open("docking.xml", "r") as file:
                content = file.readlines()
                content = "".join(content)
                bs_content = bs(content, "html.parser")
            results = bs_content.findAll("cluster", {"cluster_rank": "1"})
            score = 0
            for r in results:
                a = str(r).split()
                b = str(a[2]).replace("lowest_binding_energy=", "")
                score = b.replace('"', "")

            print("Score of " + str(g + 1) + " is " + str(score))
            asd = str(smiles_brics[g]).replace('\n', '')

            results_file = open("results.txt", "a+")
            results_file.write(str(asd) + '\t' + str(smiles_act[g]) + '\t' + str(score) + '\n')
            results_file.close()
        except Exception as e:
            print("Error Occurred")
            print(e)
            pass
        g += 1
