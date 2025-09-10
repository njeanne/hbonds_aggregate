#!/usr/bin/env python3

"""
Created on 01 Feb. 2024
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.2.0"

import argparse
import logging
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import numpy
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import seaborn as sns
from statannotations.Annotator import Annotator


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def get_domains(domains_file_path):
    """
    Extract the domains in order as they are set in the CSV file.
    
    :param domains_file_path: the path to the protein domains CSV file.
    :type domains_file_path: str
    :return: the ordered domains.
    :rtype: list
    """
    domains = None
    if domains_file_path:
        df = pd.read_csv(domains_file_path, sep=",")
        domains = list(df["domain"])
    return domains


def extract_colors(path, grouped):
    """
    Extract the colors by conditions for the boxplots and the dots.

    :param path: the path to the CSV file.
    :type path: str
    :param grouped: the grouped conditions.
    :type grouped: list
    :return: the colors to use.
    :rtype: dict
    """
    colors_by_condition = {"boxplots": {}, "dots": {}}
    first_item_group_not_stored = True
    df = pd.read_csv(path, sep=",", header=0)
    for _, row in df.iterrows():
        if grouped and row["condition"] in grouped:
            if first_item_group_not_stored:
                colors_by_condition["boxplots"]["/".join(grouped)] = row["boxplot color"]
                colors_by_condition["dots"]["/".join(grouped)] = row["dot color"]
                first_item_group_not_stored = False
        else:
            colors_by_condition["boxplots"][row["condition"]] = row["boxplot color"]
            colors_by_condition["dots"][row["condition"]] = row["dot color"]
    return colors_by_condition


def aggregate_contacts(conditions, md_time, dir_path, grouped):
    """
    Extract the number of contacts by region for each sample in a condition.

    :param conditions: the conditions dataframe.
    :type conditions: pandas.DataFrame
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param grouped: the grouped conditions.
    :type grouped: list
    :return: the aggregated data for each frame and the conditions (in case ony condition is removed) and the domain of
    interest.
    :rtype: pandas.DataFrame, str
    """
    pattern_sample = re.compile("outliers_(.+).csv")
    raw_dict = {}
    whole_domains = set()
    conditions_to_remove = []
    roi_set = set()
    for _, row_condition in conditions.iterrows():

        by_condition = []
        try:
            for fn in os.listdir(row_condition["path"]):
                if fn.startswith("outliers") and fn.endswith(".csv"):
                    by_condition.append(fn)
        except FileNotFoundError as exc:
            logging.error(exc, exc_info=True)
            sys.exit(1)
        if len(by_condition) == 0:
            conditions_to_remove.append(row_condition["condition"])
            logging.warning(f"Condition {row_condition['condition']}: no RMSD files, this condition is skipped.")
            continue
        logging.info(f"Aggregating {len(by_condition)} file{'s' if len(by_condition) > 1 else ''} data for condition: "
                     f"{row_condition['condition']}")

        # check if the condition belongs to the grouped conditions
        if grouped and row_condition.iloc[0] in grouped:
            condition = "/".join(grouped)
        else:
            condition = row_condition.iloc[0]
        if condition not in raw_dict:
            raw_dict[condition] = {}

        for item in sorted(by_condition):
            logging.info(f"\t\t- {item}")
            match_sample = pattern_sample.match(item)
            if match_sample:
                sample = match_sample.group(1)
            else:
                logging.error(f"No sample found with the pattern \"{pattern_sample.pattern}\" in the file {item}")
                sys.exit(1)
            raw_dict[condition][sample] = {}
            df_current = pd.read_csv(os.path.join(row_condition["path"], item), sep=",")
            for _, row in df_current.iterrows():
                whole_domains.add(row["second partner domain"])
                roi_set.add(row["ROI partner domain"])
                if row["second partner domain"] not in raw_dict[condition][sample]:
                    raw_dict[condition][sample][row["second partner domain"]] = row["number atoms contacts"]
                else:
                    raw_dict[condition][sample][row["second partner domain"]] += row["number atoms contacts"]
            # check if there is only one region of interest making contacts in each of the files used
            if len(roi_set) > 1:
                logging.error(f"For {sample}: more than one domain in the columns 'ROI partner domain' "
                              f"({', '.join(list(roi_set))}) of the outliers contact CSV file.")
                sys.exit(1)

    roi = list(roi_set)[0]

    # complete missing data in some domains
    for condition in raw_dict:
        for smp in raw_dict[condition]:
            for domain in whole_domains:
                if domain not in raw_dict[condition][smp]:
                    raw_dict[condition][smp][domain] = 0

    # reorganize the data
    reorganized_dict = {"sample": [], "conditions": [], "domains": [], "contacts": []}
    for condition in raw_dict:
        for smp in raw_dict[condition]:
            for domain in raw_dict[condition][smp]:
                reorganized_dict["sample"].append(smp)
                reorganized_dict["conditions"].append(condition)
                reorganized_dict["domains"].append(domain)
                reorganized_dict["contacts"].append(raw_dict[condition][smp][domain])

    df_out = pd.DataFrame.from_dict(reorganized_dict)
    out_path = os.path.join(dir_path, f"hydrogen-bonds_aggregated_{roi.lower().replace(' ', '-')}_{md_time}-ns.csv")
    df_out.to_csv(out_path, index=False)
    logging.info(f"Aggregated CSV file saved: {os.path.abspath(out_path)}")
    return df_out, roi


def compute_stats(src, domains, out_path):
    """
    Test the different domains contacts with the region of interest between the different conditions.
    A Mann-Whitney U test is performed with the null hypothesis is that the condition 1 group is greater than the
    condition 2 group.

    :param src: the contacts dataframe.
    :type src: pandas.DataFrame
    :param domains: the domains.
    :type domains: list
    :param out_path: the output file path.
    :type out_path: str
    """
    logging.info("Computing Mann-Whitney U test with a null hypothesis group 1 is greater than group 2:")
    data = {"contact with": [], "group 1": [], "group 2": [], "p-value": [], "statistic": [], "test":[], "H0": [],
            "comment": []}
    # get the conditions as a list, for loop performed to keep the conditions' order
    conditions = []
    for condition in src["conditions"]:
        if condition not in conditions:
            conditions.append(condition)
    # extract the data and compute the statistic test
    for domain in domains:
        domain_rows = src[src["domains"] == domain]
        for i in range(0, len(conditions) - 1):
            for j in range(i + 1, len(conditions)):
                data["contact with"].append(domain)
                data["test"].append("Mann-Whitney U")
                data["group 1"].append(conditions[i])
                data["group 2"].append(conditions[j])
                try:
                    test = mannwhitneyu(x=domain_rows["contacts"][domain_rows["conditions"] == conditions[i]],
                                        y=domain_rows["contacts"][domain_rows["conditions"] == conditions[j]],
                                        alternative="greater")
                    data["p-value"].append(test.pvalue)
                    data["statistic"].append(test.statistic)
                    data["comment"].append("")
                except ValueError as exc:
                    txt = "All the numbers are identical in the Mann-Whitney U test"
                    data["p-value"].append("N/A")
                    data["statistic"].append("N/A")
                    data["comment"].append(txt)
                    logging.warning(f"\t{txt} for the domain {domain} between {conditions[i]} and {conditions[j]}. "
                                    f"The test output is set to N/A.")
                data["H0"].append(f"{conditions[i]} is greater than {conditions[j]}")
    df = pd.DataFrame.from_dict(data)
    df.to_csv(out_path, index=False)
    logging.info(f"\tStatistics file saved: {out_path}")


def update_domains_order(labels, domains_ordered):
    """
    Update and order the domains by adding before, between and after annotations if some contacts are present outside
    the domains.

    :param labels: the labels.
    :type labels: list
    :param domains_ordered: the ordered domains on the protein.
    :type domains_ordered: list
    :return: the X axis labels ordered.
    :rtype: list
    """
    updated_labels = []
    # get the annotation before any domains
    for i in range(len(labels)):
        if labels[i].startswith("before"):
            updated_labels.append(labels[i])
            break
    labels[:] = [x for x in labels if not x.startswith("before")]
    logging.debug("Updating the labels:")
    logging.debug(f"\tinitial labels: {labels}")
    logging.debug(f"\tbefore domains:")
    logging.debug(f"\t\tnew labels:\t\t{updated_labels}")
    logging.debug(f"\t\tinitial labels:\t{labels}")
    # get the domains as in the ordered domains and add the between domains
    for dom in domains_ordered:
        domain_index_in_labels = {}
        logging.debug(f"\tdomain {dom}:")
        for i in range(len(labels)):
            if dom == labels[i]:
                domain_index_in_labels["dom"] = i
            elif labels[i].startswith(f"between {dom}"):
                domain_index_in_labels["between"] = i
        if "dom" in domain_index_in_labels:
            updated_labels.append(labels[domain_index_in_labels["dom"]])
        if "between" in domain_index_in_labels:
            updated_labels.append(labels[domain_index_in_labels["between"]])
        if "dom" in domain_index_in_labels:
            labels[:] = [x for x in labels if not x == dom]
        if "between" in domain_index_in_labels:
            labels[:] = [x for x in labels if not x.startswith(f"between {dom}")]
        logging.debug(f"\t\tnew labels:\t\t{updated_labels}")
        logging.debug(f"\t\tinitial labels:\t{labels}")
    # get the annotation after all the domains
    for i in range(len(labels)):
        if labels[i].startswith("after"):
            updated_labels.append(labels[i])
            break
    labels[:] = [x for x in labels if not x.startswith("after")]
    logging.debug(f"\tafter domains:")
    logging.debug(f"\t\tnew labels:\t\t{updated_labels}")
    logging.debug(f"\t\tinitial labels:\t{labels}")
    
    return updated_labels


def all_values_equals_correction(df, pairs_list):
    """
    If the values between two conditions are the sames, the Mann-Whitney test cannot be performed, a small variation is
    added to the first value of the first condition.

    :param df: the contacts dataframe.
    :type df: pandas.DataFrame
    :param pairs_list: the list of two tuples (domains and condition) to test with the Mann-Whitney test.
    :type pairs_list: list
    :return: the updated contacts dataframe.
    :rtype: pandas.DataFrame
    """
    for pairs in pairs_list:
        pair_1 = set(df["contacts"][(df["domains"] == pairs[0][0]) & (df["conditions"] == pairs[0][1])])
        pair_2 = set(df["contacts"][(df["domains"] == pairs[1][0]) & (df["conditions"] == pairs[1][1])])
        if len(pair_1) == 1 and len(pair_1) == len(pair_2):
            # add 0.00000001 to the first value of the first condition
            # to be able to perform the Mann-Withney test if all the values are the same, see:
            # https://stackoverflow.com/questions/54212583/change-1st-row-of-a-dataframe-based-on-a-condition-in-pandas
            # for explanations
            mask = (df["domains"] == pairs[0][0]) & (df["conditions"] == pairs[0][1])
            idx = mask.idxmax() if mask.any() else numpy.repeat(False, len(df))
            previous_value = df.loc[idx, "contacts"]
            df.loc[idx, "contacts"] = previous_value + 0.00000001
            logging.warning(f"\tdomain \"{pairs[0][0]}\" conditions \"{pairs[0][1]}\" and \"{pairs[1][1]}\" have the same values "
                            f"{previous_value}. The first value of the condition \"{pairs[0][1]}\" is set to "
                            f"{df.loc[idx, 'contacts']} to perform the Mann-Withney test.")
    return df


def boxplot_aggregated(src, roi, colors_plot, md_time, dir_path, fmt, domains, subtitle_arg):
    """
    Create the boxplots by conditions.

    :param src: the data source.
    :type src: pandas.DataFrame
    :param roi: region of interest, the region in contact with the other domains.
    :type roi: str
    :param colors_plot: the colors to use.
    :type colors_plot: dict
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param domains: the updated and ordered list of domains.
    :type domains: list
    :param subtitle_arg: the subtitle of the plot.
    :type subtitle_arg: str
    """
    logging.info("Plotting the aggregated hydrogen bonds by condition:")
    plt.figure(figsize=(15, 15))
    # create the statistical pairs annotations
    boxplot_pairs = []
    conditions = list(set(src["conditions"]))
    for domain in domains:
        for i in range(0, len(conditions) - 1):
            for j in range(i + 1, len(conditions)):
                boxplot_pairs.append(((domain, conditions[i]), (domain, conditions[j])))
    # for a domain, if all the values for both tested conditions are equals, add a little variation to the first value
    # of the first condition to be able to perform a Mann-Whitney test
    src = all_values_equals_correction(src, boxplot_pairs)
    # creating the plotting parameters
    plotting_parameters = {
        "data": src,
        "x": "domains",
        "y": "contacts",
        "hue": "conditions",
        "order": domains
    }

    # create the plot
    with sns.plotting_context():
        ax = sns.boxplot(**plotting_parameters, palette=colors_plot["boxplots"])
        sns.stripplot(**plotting_parameters, size=8, marker="o", linewidth=2, dodge=True, palette=colors_plot["dots"])
        annotator = Annotator(ax, boxplot_pairs, **plotting_parameters)
        annotator.configure(test="Mann-Whitney", text_format="star", hide_non_significant=True)
        annotator.apply_test(alternative="greater")
        annotator.annotate()

        # add separators between conditions
        [ax.axvline(x + 0.5, alpha=0.2) for x in ax.get_xticks()]

        # modify the ticks labels for the X axis by adding new lines every 3 words
        modified_x_labels = [re.sub(r'(\w+ \w+ \w+)( )',
                                    r'\1\n', x_label.get_text()) for x_label in ax.get_xticklabels()]
        # set the number of ticks for the X axis to avoid a matplotlib warning
        ax.set_xticks([num_tick for num_tick in range(len(modified_x_labels))])
        ax.set_xticklabels(modified_x_labels, rotation=45, horizontalalignment="right")

        # remove extra legend handles and add the count of samples by condition
        handles, labels = ax.get_legend_handles_labels()
        custom_labels = []
        number_of_labels = len(colors_plot["boxplots"])
        for label in labels[:number_of_labels]:
            sample_set = set(src[src["conditions"] == label]["sample"])
            custom_labels.append(f"{label} ({len(sample_set)})")
        ax.legend(handles[:3], custom_labels, title="Condition")

        plt.suptitle(f"Hydrogen bonds by domain with the {roi} at {md_time} ns of molecular dynamics",
                     fontsize="large", fontweight="bold")
        subtitle = "Mann-Withney H0: first condition greater than the second."
        if subtitle_arg:
            subtitle = f"{subtitle_arg}, {subtitle}"
        plt.title(subtitle)
        plt.xlabel("Domains", fontweight="bold")
        plt.ylabel(f"Hydrogen bonds", fontweight="bold")
        plot = ax.get_figure()
        out_path_plot = os.path.join(dir_path, f"hydrogen-bonds_aggregated_{roi.lower().replace(' ', '-')}_"
                                               f"{md_time}-ns.{fmt}")
        plot.savefig(out_path_plot)
    logging.info(f"\tAggregated hydrogen bonds by condition: {os.path.abspath(out_path_plot)}")


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    Aggregate the outliers hydrogen bonds by regions in one plot to compare between various conditions.

    The input is a comma separated file without header which first column is the condition, the second column the path 
    of the directory containing the contacts analysis files and the third column the color in hexadecimal format. i.e:

    insertions,tests/inputs/insertions,#fc030b
    WT,tests/inputs/WT,#0303fc

    The output is a plot with the aggregated contacts data by region and by condition.
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-t", "--md-time", required=True, type=int,
                        help="the molecular dynamics duration in nanoseconds.")
    parser.add_argument("-d", "--domains", required=True, type=str,
                        help="a sample CSV domains annotation file, to set the order of the protein domains on the X "
                             "axis. If this option is not used, the domains will be displayed randomly.")
    parser.add_argument("-g", "--group", required=False, nargs="+", type=str,
                        help="a list of conditions, separated by spaces, to group as they appear in the first column "
                             "of the input file. The color used will be the color of the first condition.")
    parser.add_argument("-s", "--subtitle", required=False, type=str,
                        help="Free text used as a subtitle for the boxplots.")
    parser.add_argument("-x", "--format", required=False, default="svg",
                        choices=["eps", "jpg", "jpeg", "pdf", "pgf", "png", "ps", "raw", "svg", "svgz", "tif", "tiff"],
                        help="the output plots format: 'eps': 'Encapsulated Postscript', "
                             "'jpg': 'Joint Photographic Experts Group', 'jpeg': 'Joint Photographic Experts Group', "
                             "'pdf': 'Portable Document Format', 'pgf': 'PGF code for LaTeX', "
                             "'png': 'Portable Network Graphics', 'ps': 'Postscript', 'raw': 'Raw RGBA bitmap', "
                             "'rgba': 'Raw RGBA bitmap', 'svg': 'Scalable Vector Graphics', "
                             "'svgz': 'Scalable Vector Graphics', 'tif': 'Tagged Image File Format', "
                             "'tiff': 'Tagged Image File Format'. Default is 'svg'.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("input", type=str,
                        help="the path to the CSV (comma separated without header) file which first column is the "
                             "condition, the second column the path of the directory containing the plots_contacts "
                             "script CSV output files and the third column the color.")
    args = parser.parse_args()

    # create output directory if necessary
    os.makedirs(args.out, exist_ok=True)
    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(args.out, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")
    logging.info(f"MD simulation time: {args.md_time} ns")

    ordered_domains = get_domains(args.domains)
    data_conditions = pd.read_csv(args.input, sep=",", header=0)
    colors = extract_colors(args.input, args.group)
    df_contacts, domain_of_interest = aggregate_contacts(data_conditions, args.md_time, args.out, args.group)
    updated_ordered_domains = update_domains_order(list(set(df_contacts["domains"])), ordered_domains)
    compute_stats(df_contacts, updated_ordered_domains,
                  os.path.join(args.out, f"hydrogen-bonds_aggregated_statistics_{domain_of_interest.replace(' ', '-')}_"
                                         f"{args.md_time}-ns.csv"))
    boxplot_aggregated(df_contacts, domain_of_interest, colors, args.md_time, args.out, args.format,
                       updated_ordered_domains, args.subtitle)
