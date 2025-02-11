from .base import dispatch
from .parse_input import InputParseCLI
from .taxa_filter import TaxaFilterCLI
from .render_phylogeny import PhylogenyRenderCLI
from .inference import InferenceCLI
from .infer_negbin import NegBinCLI
from .visualize_negbin import NegBinVisualizationCLI
from .interactions_to_cytoscape import InteractionToCytoscapeCLI
from .visualize_posterior import PosteriorVisualizationCLI
from .forward_simulate import ForwardSimulationCLI
from .compute_keystoneness import KeystonenessCLI


def main():
    # ========= Mapping of subcommands to cli modules.
    clis = [
        InputParseCLI(subcommand="input"),
        InferenceCLI(subcommand="infer"),
        TaxaFilterCLI(subcommand="filter"),
        PhylogenyRenderCLI(subcommand="render-phylogeny"),
        NegBinCLI(subcommand="infer-negbin"),
        NegBinVisualizationCLI(subcommand="visualize-negbin"),
        InteractionToCytoscapeCLI(subcommand="interaction-to-cytoscape"),
        PosteriorVisualizationCLI(subcommand="visualize-posterior"),
        ForwardSimulationCLI(subcommand="forward-simulate"),
        KeystonenessCLI(subcommand="evaluate-keystoneness")
    ]

    dispatch({
        cli.subcommand: cli for cli in clis
    })
