# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Various utility functions
# Kristian Jensen, March 2015
#

from __future__ import print_function

import math
import cameo
import colorsys
from cameo.util import TimeMachine
from functools import partial
from cameo import config
import re
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
import bokeh
import sys
from scipy import stats
from cameo.exceptions import SolveError


class FeatureMap(object):
    def __init__(self, seq_record, feature_id=None):
        self.genome = seq_record
        self._cds_features = [feature for feature in seq_record.features if feature.type == "CDS"]
        self._feature_dict = self._index_features(self._cds_features, feature_id)

    def _index_features(self, features, feature_id):
        if feature_id is None:
            feature_id = lambda feat: feat.qualifiers["locus_tag"][0]
        feature_dict = {feature_id(feature): feature for feature in features}
        return feature_dict

    def features(self, pos=None, end=None):
        if pos is None:
            return self._cds_features
        elif end is None:
            return list(f for f in self._cds_features if f.location.start <= pos < f.location.end)
        else:
            return [f for f in self._cds_features if f.location.start < end and pos < f.location.end]

    def mutate(self, mutation): # TODO several mutations at a time
        pos = mutation.position
        mutated_component = self._component.mutate([mutation])
        mutated_features = [f for f in mutated_component.features if f.location.start <= pos <= f.location.end and f.type == "CDS"]
        return mutated_features

    def get_by_id(self, value):
        return self._feature_dict[value]


def round_fva_bound(bound, ndecimals):
    """
    Rounds a bound away from zero to ndecimals precision.
    """
    if bound > 0:
        bound = math.ceil(bound * 10**ndecimals)/10**ndecimals
    elif bound < 0:
        bound = math.floor(bound * 10**ndecimals)/10**ndecimals
    return bound


def apply_fva_constraints(model, fva_result, skip_zeros=False, skip_list=[]):
    """
    Applies the reaction bounds found with an FVA analysis to the model.
    If skip_zeros is set to True, only reactions with non-zero FVA bounds are changed. Reversible
    reactions that are only used in one direction in the FVA are only limited in that direction.
    """
    from cameo import config
    ndecimals = config.ndecimals
    exchanges = model.exchanges

    for reaction_id in fva_result.index:
        if skip_zeros:
            if fva_result["upper_bound"][reaction_id] == 0 and fva_result["lower_bound"][reaction_id] == 0:
                continue
        reaction = model.reactions.get_by_id(reaction_id)
        if reaction in exchanges or reaction in skip_list or reaction_id in skip_list or reaction_id == "ATPM":
            pass
        else:
            if reaction.upper_bound > 0:
                if fva_result["upper_bound"][reaction_id] <= 0:
                    if skip_zeros is False:
                        reaction.upper_bound = 0
                else:
                    rounded_bound = round_fva_bound(fva_result["upper_bound"][reaction_id], ndecimals)
                    if reaction.upper_bound == 1000:
                        reaction.upper_bound = rounded_bound
                    else:
                        reaction.upper_bound = max(rounded_bound, reaction.upper_bound)
            if reaction.lower_bound < 0:
                if fva_result["lower_bound"][reaction_id] >= 0:
                    if skip_zeros is False:
                        reaction.lower_bound = 0
                else:
                    rounded_bound = round_fva_bound(fva_result["lower_bound"][reaction_id], ndecimals)
                    if reaction.lower_bound == -1000:
                        reaction.lower_bound = rounded_bound
                    else:
                        reaction.lower_bound = min(rounded_bound, reaction.lower_bound)


def binding_constraints(model, ignore_zeros=True):
    """
    Solves the model and returns a list of the reactions that are bound by a constraint.
    :param model: A model
    :param ignore_zeros: If True reactions that are constrained to zero will not be reported.
    :return: A list of reactions
    """
    binding = []
    fluxes = model.solve().fluxes
    for reaction in model.reactions:
        flux = fluxes[reaction.id]
        if ignore_zeros:
            if flux == 0:
                continue
        if flux == reaction.upper_bound or flux == reaction.lower_bound:
            binding.append(reaction)
    return binding


def find_mutant_growth_rate(model, relative_bounds, original_bounds=None, skip_essential=False, solve_method=cameo.fba,
                            essential_reactions=None, **simulation_kwargs):
    """
    Calculates the fluxes for a model where bounds have been modified.
    :param model: SolverBasedModel
    :param relative_bounds: A dictionary of relative bounds that will be applied to the bounds of each reaction
    :original_bounds: A pandas dataframe (as output by cameo.flux_variability_analysis) containing initial bounds.
        If None, an FVA will be run to determine them.
    :param skip_essential: Reaction bounds will not be modified if relative bound is less than this parameter.
        Set to False or None to disable skipping.
    """
    with TimeMachine() as tm:
        # Make sure bounds are reset
        for reaction in model.reactions:
            tm(do=int, undo=partial(setattr, reaction, "upper_bound", reaction.upper_bound))
            tm(do=int, undo=partial(setattr, reaction, "lower_bound", reaction.lower_bound))

        # Find essential reactions
        if skip_essential:
            if essential_reactions is None:
                essential_reactions = model.essential_reactions()

        # Perform FVA
        if original_bounds is None:
            fva_result = cameo.flux_variability_analysis(model, fraction_of_optimum=1.0, remove_cycles=False)
            fva_result["lower_bound"] = fva_result["lower_bound"].round(6)
            fva_result["upper_bound"] = fva_result["upper_bound"].round(6)
        else:
            fva_result = original_bounds

        # Set initial bounds
        for reaction_id in fva_result.index:
                reaction = model.reactions.get_by_id(reaction_id)
                if fva_result["upper_bound"][reaction_id] >= 0:
                    reaction.upper_bound = round_fva_bound(fva_result["upper_bound"][reaction_id], config.ndecimals)
                if fva_result["lower_bound"][reaction_id] <= 0:
                    reaction.lower_bound = round_fva_bound(fva_result["lower_bound"][reaction_id], config.ndecimals)

        # Modify bounds
        for reaction_id in relative_bounds:
            reaction = model.reactions.get_by_id(convert_reaction(reaction_id))
            if skip_essential and reaction in essential_reactions and round(relative_bounds[reaction_id], config.ndecimals) <= skip_essential:
                #print "Skipped "+reaction_id
                continue
            if convert_reaction(reaction_id) not in fva_result.index:
                raise Warning("Reaction "+reaction_id+" is not in the specified original bounds")
            reaction.lower_bound *= relative_bounds[reaction_id]
            reaction.upper_bound *= relative_bounds[reaction_id]

        # Solve mutated model
        #solution = model.solve()
        solution = solve_method(model, **simulation_kwargs)
    return solution


def convert_reaction(reaction):
    """
    Converts a reaction id to cameo format.
    :param reaction:
    :return:
    """
    reaction = re.sub(r"-", "_dsh_", reaction)
    reaction = re.sub(r"\(", "_lp_", reaction)
    reaction = re.sub(r"\)", "_rp_", reaction)
    return reaction


def compare_flux_distributions(flux_dist_1, flux_dist_2, reactions, ndecimals=None):
    result = []
    if ndecimals is None:
        ndecimals = config.ndecimals
    for reaction in reactions:
        try:
            reac_id = reaction.id
        except AttributeError:
            reac_id = reaction

        if round(flux_dist_1[reac_id], ndecimals) != round(flux_dist_2[reac_id], ndecimals):
            result.append( (reac_id, flux_dist_1[reac_id], flux_dist_2[reac_id], flux_dist_1[reac_id]-flux_dist_2[reac_id]) )

    if len(result) > 0:
        return pd.DataFrame(result, columns=["reaction", "flux_1", "flux_2", "difference"])
    else:
        return pd.DataFrame(columns=["reaction", "flux_1", "flux_2", "difference"])


def plot_predicted_growth_rates(predicted_solutions, experimental_growth_rate_filename, standardize=True, show_correlation=True):

    gr_df = pd.read_csv(experimental_growth_rate_filename)

    calc_growth_rates = {strain.split()[0]: solution.fluxes["Ec_biomass_iJO1366_core_53p95M"] for strain, solution in predicted_solutions.items()}

    plot_df = gr_df.copy()
    plot_df["calc"] = plot_df["strain"].map(calc_growth_rates.get)
    plot_df = plot_df[plot_df["calc"].notnull()]

    calc_mean = plot_df["calc"].values.mean()
    calc_std = plot_df["calc"].values.std()

    mean_mean = plot_df["mean"].values.mean()
    mean_std = plot_df["mean"].values.std()

    def standardize_calc(num):
        if standardize:
            return (num - calc_mean) / calc_std
        else:
            return num

    def standardize_mean(num):
        if standardize:
            return (num - mean_mean) / mean_std
        else:
            return num

    cor_coef = calculate_growth_correlation_coefficient(predicted_solutions, experimental_growth_rate_filename)[0]

    if standardize:
        x_range = [-2, 3]
        y_range = [-2, 3]
        text_spot = [-2.5, 2.5]
    else:
        x_range = [0, 1]
        y_range = [-0.05, 1]
        text_spot = [0.1, 0.9]

    #x_range = None
    #y_range = None

    slope, intercept, pearson, p_val, stderr = stats.linregress(plot_df["mean"].values, plot_df["calc"].values)

    fig = bokeh.plotting.figure(title=None, x_range=x_range, y_range=y_range)

    fig.scatter(standardize_mean(plot_df["mean"].values), standardize_calc(plot_df["calc"].values))
    #fig.line(np.arange(0, 1, 0.1), (6.83)*np.arange(0, 1, 0.1)-3.74, color="red")
    #fig.line(*([np.arange(-3,3,0.1)]*2), color="red")
    for idx in plot_df.index:
        mean = plot_df["mean"][idx]
        std = plot_df["std"][idx]
        n = plot_df["n"][idx]
        sem = std / math.sqrt(n)
        calc = standardize_calc(plot_df["calc"][idx])
        low, high = standardize_mean(mean-sem), standardize_mean(mean+sem)
        fig.line([low, high], [calc, calc])

    if show_correlation:
        fig.text(*text_spot, text="Corr. coefficient: "+str(round(cor_coef, 3)))
    fig.xaxis.axis_label = "Exp"
    fig.yaxis.axis_label = "Calc"

    return fig


def calculate_growth_correlation_coefficient(predicted_solutions, experimental_growth_rate_filename):
    exps = []
    calcs = []

    gr_df = pd.read_csv(experimental_growth_rate_filename)
    for idx, row in gr_df.iterrows():
        strain = row[0]+" v1"
        if strain in predicted_solutions:
            exp = row[1]
            calc = predicted_solutions[strain].fluxes["Ec_biomass_iJO1366_core_53p95M"]
            exps.append(exp)
            calcs.append(calc)

    cor_coef = np.corrcoef(np.array(exps), np.array(calcs))[0, 1]
    p_val = stats.pearsonr(np.array(exps), np.array(calcs))[1]
    return cor_coef, p_val


def calculate_mean_squared_error(predicted_solutions, experimental_growth_rate_filename):
    exps = []
    calcs = []

    gr_df = pd.read_csv(experimental_growth_rate_filename)
    for idx, row in gr_df.iterrows():
        strain = row[0]+" v1"
        if strain in predicted_solutions:
            exp = row[1]
            calc = predicted_solutions[strain].fluxes["Ec_biomass_iJO1366_core_53p95M"]
            exps.append(exp)
            calcs.append(calc)

    sqsum = sum((np.array(exps)-np.array(calcs))**2)
    root_mean_sq = (sqsum / len(exps))**0.5
    return root_mean_sq


def calculate_growth_spearman_rank(predicted_solutions, experimental_growth_rate_filename):
    exps = []
    calcs = []

    gr_df = pd.read_csv(experimental_growth_rate_filename)
    for idx, row in gr_df.iterrows():
        strain = row[0]+" v1"
        if strain in predicted_solutions:
            exp = row[1]
            calc = predicted_solutions[strain].fluxes["Ec_biomass_iJO1366_core_53p95M"]
            exps.append(exp)
            calcs.append(calc)

    cor_coef = stats.spearmanr(np.array(exps), np.array(calcs))
    return cor_coef


def find_strain_solutions(model, strains, rel_bounds, solve_method, original_bounds, **kwargs):
    """
    :param model:
    :param strains:
    :param rel_bounds:
    :param solve_method:
    :param original_bounds:
    :param kwargs:
    :return:
    """
    strain_solutions = {}
    for strain in strains:
        strain_rel_bounds = rel_bounds[strain]
        try:
            solution = find_mutant_growth_rate(model, strain_rel_bounds, original_bounds=original_bounds,
                                               solve_method=solve_method, **kwargs)
        except SolveError as e:
            sys.stderr.write(strain+" could not be solved ("+str(e)+")\n")
            continue
        strain_solutions[strain] = solution

    return strain_solutions


def find_initial_bounds_from_c13(model, c13_path, activity_level, reference=None):
    """
    :param model:
    :param c13_path:
    :param activity_level:
    :return: Reference flux distribution result and initial_bounds
    """
    # Reactions to skip
    exchanges = [reac.id for reac in model.exchanges]
    skip_reacs = exchanges + ["ATPM", "Ec_biomass_iJO1366_core_53p95M"]

    c13_bounds = pd.DataFrame.from_csv(c13_path)
    if reference is None:
        with TimeMachine() as tm:
            for r_id in c13_bounds.index:
                lower = round(c13_bounds["lower_bound"][r_id], 6)
                upper = round(c13_bounds["upper_bound"][r_id], 6)
                r = model.reactions.get_by_id(r_id)
                tm(do=partial(setattr, r, "lower_bound", lower),
                   undo=partial(setattr, r, "lower_bound", r.lower_bound))
                tm(do=partial(setattr, r, "upper_bound", upper),
                   undo=partial(setattr, r, "upper_bound", r.upper_bound))

            reference = cameo.pfba(model)

    with TimeMachine() as tm:
        for r_id in c13_bounds.index:
            if r_id in skip_reacs:
                continue
            lower = round(c13_bounds["lower_bound"][r_id], 6)
            upper = round(c13_bounds["upper_bound"][r_id], 6)
            reac = model.reactions.get_by_id(r_id)
            if lower <= 0:
                tm(do=partial(setattr, reac, "lower_bound", lower/activity_level),
                   undo=partial(setattr, reac, "lower_bound", reac.lower_bound))
            if upper >= 0:
                tm(do=partial(setattr, reac, "upper_bound", upper/activity_level),
                   undo=partial(setattr, reac, "upper_bound", reac.upper_bound))

        biomass = model.reactions.Ec_biomass_iJO1366_core_53p95M
        #tm(do=partial(setattr, r, "upper_bound", reference["Ec_biomass_iJO1366_core_53p95M"]),
        #   undo=partial(setattr, r, "upper_bound", biomass.upper_bound))
        wt_bounds = cameo.flux_variability_analysis(model, fraction_of_optimum=1, remove_cycles=False).data_frame
        wt_bounds["lower_bound"] = wt_bounds["lower_bound"].round(6)
        wt_bounds["upper_bound"] = wt_bounds["upper_bound"].round(6)
        for r_id in c13_bounds.index:
            if r_id in skip_reacs:
                continue
            lower = round(c13_bounds["lower_bound"][r_id], 6)
            upper = round(c13_bounds["upper_bound"][r_id], 6)
            if lower <= 0:
                wt_bounds.loc[r_id, "lower_bound"] = lower / activity_level
            if upper >= 0:
                wt_bounds.loc[r_id, "upper_bound"] = upper / activity_level
        wt_result = reference

        return wt_result, wt_bounds