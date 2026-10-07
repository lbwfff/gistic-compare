"""Plotting for gistic-compare, split by analysis module.

``theme`` handles the publication style; each submodule here only draws the
figures for its module.  Public entry points (``figures_*``) are re-exported
below with the same names the CLI expects.
"""

from ._core import *  # noqa: F401,F403
from .benchmark_a import *  # noqa: F401,F403
from .benchmark_b import *  # noqa: F401,F403
from .benchmark_c import *  # noqa: F401,F403
from .benchmark_d import *  # noqa: F401,F403
from .benchmark_e import *  # noqa: F401,F403
from .benchmark_f import *  # noqa: F401,F403
from .benchmark_g import *  # noqa: F401,F403
from .confounding import *  # noqa: F401,F403
from .confounding import figures_confounding  # noqa: F401
from .entry_benchmark import (  # noqa: F401
    figures_a,
    figures_b,
    figures_c,
    figures_d,
    figures_e,
    figures_f,
)
from .entry_g import figures_g  # noqa: F401
from .gs import *  # noqa: F401,F403
from .gs import figures_gs  # noqa: F401
from .oncoprint import figures_oncoprint, plot_oncoprint  # noqa: F401
from .pairwise_p1 import *  # noqa: F401,F403
from .pairwise_p1 import figures_p1  # noqa: F401
from .pairwise_p2 import *  # noqa: F401,F403
from .pairwise_p2 import figures_p2  # noqa: F401
from .pairwise_p3 import *  # noqa: F401,F403
from .pairwise_p3 import figures_p3  # noqa: F401
from .pairwise_p4 import *  # noqa: F401,F403
from .pairwise_p4 import figures_p4, figures_p5  # noqa: F401
