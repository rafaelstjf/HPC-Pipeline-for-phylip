import parsl
import os
import random
import json
import argparse
import logging
from parsl.channels import LocalChannel
from parsl.launchers import SrunLauncher, SingleNodeLauncher
from parsl.addresses import address_by_interface, address_by_hostname
from parsl.executors import HighThroughputExecutor, ThreadPoolExecutor
from parsl.providers import LocalProvider, SlurmProvider
from apps import *
from config import GenerateSettings




def main(params):
    parsl.load(GenerateSettings(params))
    if params["align"] != False:
        r_mafft = Mafft(params)
    else:
        r_mafft = None
    r_folders = CreateFolders(params, inputs=[r_mafft])
    r_trees = list()
    for i in range(1, params["bootstrap"] + 1):
        path = os.path.join(params["dir"], os.path.join("bootstrap",  str(i)))
        seed = 0
        while True:
            seed = random.randint(0, 30000)
            if seed % 2 == 1:
                break
        r_boot = SeqBoot(params, path, seed, inputs=[r_folders])
        r_ren0 = Rename(params, path, "outfile",
                        "alignment.phy", inputs=[r_boot])
        if params["method"] == "NJ":
            r_dnadist = DnaDist(params, path, inputs=[r_ren0])
            r_ren1 = Rename(params, path, "outfile",
                        "dna_dist.data", inputs=[r_dnadist])
            r_tree = CreateTree(params, path, inputs=[r_ren1])
        elif params["method"] == "MP":
            r_tree = CreateTree(params, path, inputs=[r_ren0])
        r_ren2 = Rename(params, path, "outtree",
                        "bstree.newick", inputs=[r_tree])
        r_ren3 = Rename(params, path, "outfile", "tree.log", inputs=[r_ren2])
        r_trees.append(r_ren3)
    r_conc = ConcatenateBSTrees(params, inputs=r_trees)
    r_consense = ConsenseTree(params, inputs=[r_conc])
    r_ren4 = Rename(params, params['dir'], "outfile", "con_tree.log", inputs=[r_consense])
    r_ren5 = Rename(params, params['dir'], "outtree", "con_tree.newick", inputs=[r_ren4])
    r_root = AssignBSvalues(params, inputs = [r_ren5])
    return r_root.result()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='HPC-Phylip-Pipeline',
        description='Create a phylogenetic tree using phylip from a alignment file',
    )
    parser.add_argument(
        'json_file', help="Json file containing all the settings and parameters")
    args = parser.parse_args()
    json_data = dict()
    with open(args.json_file, 'r') as json_file:
        json_data = json.load(json_file)
    try:
        main(json_data)
    except Exception as e:
        pass
    finally:
        parsl.dfk().cleanup()