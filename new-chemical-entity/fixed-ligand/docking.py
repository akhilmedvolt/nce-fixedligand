import pandas as pd
from biopandas.pdb import PandasPdb

import os
import time

from bs4 import BeautifulSoup as bs


def fragment_docking():
    smile_df = pd.read_csv("molecules.csv", names=["smiles"])
    smile_df.drop_duplicates(['smiles'], inplace=True)

    smiles_act = smile_df["smiles"].tolist()

    os.mkdir("complex")
    os.mkdir("VG")

    g = 0
    print("Main Docking Starts")
    while g < len(smiles_act):
        try:
            print(g)
            print()
            print(smiles_act[g])

            # Conversion of smiles
            s1 = 'obabel -:"' + str(smiles_act[g]) + '" -opdb -O ligact.pdb --gen3d'
            os.system(s1)

            os.system('obminimize -ff GAFF ligact.pdb > lig_optact.pdb')
            os.system(
                '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24'
                '/prepare_ligand4.py -A bonds_hydrogens -U nphs_lps -l lig_optact.pdb -o ligact.pdbqt')

            # Docking Starts
            os.system('autodock_gpu_64wi -ffile finalmodel.maps.fld -lfile ligact.pdbqt -nrun 100')

            # writing the output of the best docked pose of ligand

            os.system(
                '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24'
                '/write_lowest_energy_ligand.py -f docking.dlg -o lig_out2.pdbqt')

            os.system('obabel lig_out2.pdbqt -opdb -O lig_out2.pdb')

            # creating the complex

            ligo = PandasPdb()
            ligo.read_pdb('lig_out2.pdb')

            ldb_2 = PandasPdb()
            ldb_2.read_pdb('receptor/protein.pdb')

            # Creating the Complex Files
            ldb_2.df['ATOM'] = ldb_2.df['ATOM'].append(ligo.df['HETATM'])
            ldb_2.to_pdb(path='complex/complex_{}.pdb'.format(smiles_act[g]),
                         records=['ATOM', 'HETATM'], gz=False, append_newline=True)

            # creating files needed for VG's pIC50 prediction (merged pocket and ligand file in pdb format)
            ldb2 = PandasPdb().read_pdb('complex/complex_{}.pdb'.format(smiles_act[g]))

            maxx = ldb2.df['HETATM']['x_coord'].max() + 3.5
            minx = ldb2.df['HETATM']['x_coord'].min() + -3.5
            maxy = ldb2.df['HETATM']['y_coord'].max() + 3.5
            miny = ldb2.df['HETATM']['y_coord'].min() + -3.5
            maxz = ldb2.df['HETATM']['z_coord'].max() + 3.5
            minz = ldb2.df['HETATM']['z_coord'].min() + -3.5

            ldb = PandasPdb().read_pdb('complex/complex_{}.pdb'.format(smiles_act[g]))
            ldb.df['ATOM'] = ldb.df['ATOM'][(ldb.df['ATOM']['x_coord'] < maxx) & (ldb.df['ATOM']['x_coord'] > minx) & (
                    ldb.df['ATOM']['y_coord'] < maxy) & (ldb.df['ATOM']['y_coord'] > miny) & (
                                                    ldb.df['ATOM']['z_coord'] < maxz) & (
                                                    ldb.df['ATOM']['z_coord'] > minz)]
            hi = ldb.df['ATOM']['residue_number'].tolist()

            res = []
            for i in hi:
                if i not in res:
                    res.append(i)

            ldb2.to_pdb(path='tempdb11.pdb', records=['HETATM'], gz=False, append_newline=True)
            ldb2.read_pdb('tempdb11.pdb')

            e1 = 0

            while e1 < len(res):
                ldb = PandasPdb().read_pdb('complex/complex_{}.pdb'.format(smiles_act[g]))
                ldb2.df['ATOM'] = ldb2.df['ATOM'].append(ldb.df['ATOM'][ldb.df['ATOM']['residue_number'] == res[e1]])
                e1 += 1

            ldb2.to_pdb('VG/pocketlig_{}.pdb'.format(smiles_act[g]), records=['ATOM', 'HETATM'], gz=False,
                        append_newline=True)

            content = []
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

            res_file = open("final_results.txt", "a+")
            res_file.write(str(smiles_act[g]) + '\t' + str(score) + '\n')
            res_file.close()

        except Exception as e:
            print(e)
            pass
        g += 1
    os.mkdir("final_result")
    results_df = pd.read_csv('final_results.txt', sep="\t", names=['Smile', 'Score'])
    results_df.to_csv("final_result/final_results.csv")

    print("FINAL Docking Completed")