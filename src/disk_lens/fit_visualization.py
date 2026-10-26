"""Handle invalid intermediate diagnostic plots without changing fit likelihoods."""

import logging

from autoarray.exc import InversionException
import autolens as al
from autolens.imaging.model.visualizer import VisualizerImaging


logger = logging.getLogger(__name__)


class PlotSafeVisualizerImaging(VisualizerImaging):
    """Skip only intermediate plots whose best-fit inversion fails a mesh check.

    PyAutoFit evaluates likelihoods separately. A model used for an update plot
    can fail the same mesh check even after the search has accepted other valid
    samples. Keep the check in the likelihood and retain final-result failures.
    """

    @staticmethod
    def visualize(analysis, paths, instance, during_analysis):
        try:
            return VisualizerImaging.visualize(
                analysis=analysis,
                paths=paths,
                instance=instance,
                during_analysis=during_analysis,
            )
        except InversionException:
            if not during_analysis:
                raise
            logger.warning(
                "Skipping intermediate diagnostic plot: selected model failed "
                "the inversion mesh check; likelihood and search settings are unchanged."
            )


class DiagnosticSafeAnalysisImaging(al.AnalysisImaging):
    """Keep Pix2 inference unchanged if an optional update diagnostic is invalid."""

    Visualizer = PlotSafeVisualizerImaging

    def profile_log_likelihood_function(self, instance, paths=None):
        try:
            return super().profile_log_likelihood_function(
                instance=instance, paths=paths
            )
        except InversionException:
            logger.warning(
                "Skipping optional runtime profile: selected model failed "
                "the inversion mesh check; likelihood and search settings are unchanged."
            )
            return None
