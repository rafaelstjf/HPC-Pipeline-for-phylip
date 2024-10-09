import parsl
import logging
import os
from parsl.launchers import SrunLauncher, SingleNodeLauncher
from parsl.adresses import address_by_hostname
from parsl.executors import HighThroughputExecutor, ThreadPoolExecutor
from parsl.providers import LocalProvider, SlurmProvider


def GenerateSettings(params):
    name = "HPC-Phylip-Pipeline"
    interval = 30
    monitor = params["monitor"]
    parsl.set_file_logger(os.path.join(params['dir'], 'log_script.output'))
    logging.info('Configuring Parsl Workflow Infrastructure')
    env_str = params['environment_variables']
    env_str += f';export PYTHONPATH=$PYTHONPATH:{os.getcwd()}'
    logging.info(f'Task Environment {env_str}')
    mon_hub = parsl.monitoring.monitoring.MonitoringHub(
        workflow_name=name,
        hub_address=address_by_hostname(),
        resource_monitoring_enabled=True,
        monitoring_debug=False,
        resource_monitoring_interval=interval,
    ) if monitor else None
    config = None
    if (params['environment_kind']).lower() == "slurm-inside-job":
        config = parsl.config.Config(
            executors=[
                HighThroughputExecutor(
                    address=address_by_hostname(),
                    cores_per_worker=1,
                    max_workers_per_node=params['num_cores'],
                    worker_debug=False,
                    provider=LocalProvider(
                        nodes_per_block=1,
                        parallelism=1,
                        init_blocks=params['num_nodes'],
                        max_blocks=params['num_nodes'],
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
            executors=[ThreadPoolExecutor(
                max_threads=params["num_cores"]
            )
            ],
            monitoring=mon_hub,
            strategy=None,
        )
    elif (params["environment_kind"]).lower() == "slurm":
        config = parsl.config.Config(
            executors[HighThroughputExecutor(
                address=address_by_hostname(),
                max_workers_per_node=params['num_cores'],
                worker_debug=False,
                provider=SlurmProvider(
                    partition=params["partition"],
                    nodes_per_block=1,
                    cores_per_node=params["num_cores"],
                    parallelism=1,
                    max_blocks=params['num_nodes'],
                    walltime=params["walltime"],
                    worker_init=params["env_variables"],
                    move_files=False,
                    launcher=SrunLauncher(
                        overrides=f"-c {params['num_cores']}")

                )
            )]
        )
    return config
