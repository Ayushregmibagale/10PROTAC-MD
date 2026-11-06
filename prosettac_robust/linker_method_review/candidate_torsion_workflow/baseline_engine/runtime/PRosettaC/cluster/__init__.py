import sys

from .Cluster import Cluster
from .PBS.PBS import PBS
from .SGE.SGE import SGE
from .SLURM.SLURM import SLURM
from .LOCAL.LOCAL import LOCAL


def getCluster(name: str) -> Cluster:
        """Get a handle to the cluster given by name."""
        if name == 'SGE':
                return SGE()
        elif name == 'PBS':
                return PBS()
        elif name == 'SLURM':
                return SLURM()
        elif name == 'LOCAL':
                return LOCAL()
        else:
                sys.exit(f'Queue "{name}" not supported.')
