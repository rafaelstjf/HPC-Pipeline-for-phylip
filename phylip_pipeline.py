import parsl
import os
import random
import json
import argparse
import logging
from parsl.channels import LocalChannel
from parsl.launchers import SrunLauncher, SingleNodeLauncher
from parsl.addresses import address_by_interface, address_by_hostname
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider, SlurmProvider
from apps import *


def GenerateSettings(params):
    name = "HPC-Phylip-Pipeline"
    interval = 30
    monitor = params['monitor']
    parsl.set_file_logger(os.path.join(params['dir'], 'log_script.output'))
    logging.info('Configuring Parsl Workflow Infrastructure')
    env_str = params['environment_variables']
    env_str += f';export PYTHONPATH=$PYTHONPATH:{params["dir"]}'
    logging.info(f'Task Environment {env_str}')
    mon_hub = parsl.monitoring.monitoring.MonitoringHub(
        workflow_name=name,
        hub_address=address_by_hostname(),
        resource_monitoring_enabled=True,
        monitoring_debug=False,
        resource_monitoring_interval=interval,
    ) if monitor else None
    config = None
    if (params['environment_kind']).lower() == "slurm":
        config = parsl.config.Config(
            executors=[
                HighThroughputExecutor(
                    address=address_by_hostname(),
                    cores_per_worker=1,
                    max_workers=params['num_cores'],
                    worker_debug=False,
                    provider=LocalProvider(
                        nodes_per_block=1,
                        channel=LocalChannel(script_dir=params['dir']),
                        parallelism=1,
                        init_blocks=params['num_blocks'],
                        max_blocks=params['num_blocks'],
                        worker_init=env_str,
                        launcher=SrunLauncher(
                            overrides=f"-c {params['num_cores']}")
                    ),
                ),
            ],
            monitoring=mon_hub,
            strategy=None,
        )
    elif (params['environment_kind']).lower() == "local":
        config = parsl.config.Config(
            executors=[
                HighThroughputExecutor(
                    worker_debug=False,
                    cores_per_worker=1,
                    max_workers=params['num_cores']
                    provider=LocalProvider(
                        channel=LocalChannel(),
                        init_blocks=1,
                        max_blocks=1,
                    ),
                )
            ],
            monitoring=mon_hub,
            strategy=None,
        )
    return config


def main(params):
    parsl.load(GenerateSettings(params))
    r_mafft = Mafft(params)
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
        r_dnadist = DnaDist(params, path, inputs=[r_ren0])
        r_ren1 = Rename(params, path, "outfile",
                        "dna_dist.data", inputs=[r_dnadist])
        r_tree = CreateTree(params, path, inputs=[r_ren1])
        r_ren2 = Rename(params, path, "outtree",
                        "bstree.newick", inputs=[r_tree])
        r_ren3 = Rename(params, path, "outfile", "tree.log", inputs=[r_ren2])
        r_trees.append(r_ren3)
    r_conc = ConcatenateBSTrees(params, inputs=r_trees)
    r_consense = ConsenseTree(params, inputs=[r_conc])
    return r_consense.result()


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
    main(json_data)
