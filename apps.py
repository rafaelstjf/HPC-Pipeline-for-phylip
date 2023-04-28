import parsl
from Bio import Phylo
from parsl.app.app import python_app, bash_app

@bash_app
def Mafft(params, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os
    seq_file = os.path.join(params["dir"], params["input"])
    align_file = os.path.join(params["dir"], "alignment.phy")
    return f"{params['mafft_bin']} --auto --phylipout --reorder \"{seq_file}\" > \"{align_file}\""

@python_app
def CreateFolders(params, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os, shutil
    from pathlib import Path
    work_dir = params['dir']
    try:
        Path(os.path.join(work_dir, "bootstrap")).mkdir(exist_ok=True)
    except:
        pass
    for i in range(1, (params['bootstrap'] + 1)):
        full_path = os.path.join(work_dir, "bootstrap")
        full_path = os.path.join(full_path, f"{i}")
        if(os.path.exists(full_path)):
            try:
                shutil.rmtree(full_path, ignore_errors=True)
            except:
                pass
    for i in range(1, (params['bootstrap'] + 1)):
        full_path = os.path.join(work_dir, "bootstrap")
        full_path = os.path.join(full_path, f"{i}")
        try:
            Path(full_path).mkdir(exist_ok=True)
        except:
            pass
    return

@bash_app
def SeqBoot(params, path, seed, inputs = [], outputs = []):
    import os
    os.system(f"printf \"{os.path.join(params['dir'], 'alignment.phy')}\nR\n1\nY\n{seed}\" > {os.path.join(path, 'input_seq.txt')}")
    return f"cd {path};{params['seqboot_bin']} < {os.path.join(path, 'input_seq.txt')}"

@python_app
def CreateBootstrap(params, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os
    from Bio import Phylo, AlignIO
    from Bio.Phylo.Consensus import bootstrap
    msa = AlignIO.read(os.path.join(params['dir'], "alignment.phy"), "phylip")
    samples = bootstrap(msa, params['bootstrap'])
    i = 1
    for s in list(samples):
        folder = os.path.join(params['dir'], os.path.join("bootstrap", str(i)))
        AlignIO.write(s, os.path.join(folder, "alignment.phy"), "phylip")
        i+=1
    return 

@bash_app
def DnaDist(params, path, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os
    os.system(f"printf \"{os.path.join(path, 'alignment.phy')}\nY\" > {os.path.join(path, 'input.txt')}")
    return f"cd {path};{params['dnadist_bin']} < {os.path.join(path, 'input.txt')}"

@python_app
def Rename(params, path, old_name, new_name, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os
    os.rename(os.path.join(path, old_name), os.path.join(path, new_name))
    return

@bash_app
def CreateTree(params, path, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os
    if params["method"] == "NJ":
        os.system(f"printf \"{os.path.join(path, 'dna_dist.data')}\n123\nY\" > {os.path.join(path, 'input_tree.txt')}")
        return f"cd {path};{params['neighbor_bin']} < input_tree.txt"
    return

@python_app
def ConcatenateBSTrees(params, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os
    work_dir = params['dir']
    files = list()
    for i in range(1, params['bootstrap'] + 1):
        folder = os.path.join("bootstrap", str(i))
        full_path = os.path.join(work_dir, folder)
        files.append(os.path.join(full_path, "bstree.newick"))
    with open(os.path.join(work_dir, "bstreeslist.newick"), 'w') as outfile:
        for fname in files:
            with open(fname) as infile:
                outfile.write(infile.read())
    return

@bash_app
def ConsenseTree(params, inputs = [], outputs = [], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME):
    import os
    os.system(f"printf \"{os.path.join(params['dir'], 'bstreeslist.newick')}\nY\" > {os.path.join(params['dir'], 'input_consense.txt')}")
    return f"cd {params['dir']};{params['consense_bin']} < input_consense.txt"