#!/usr/bin/env python3

"""
Created on 01 Feb. 2024
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
import logging
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
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



def get_conditions(path):
    """
    Extract the conditions, the paths and the colors.

    :param path: the path to the CSV file.
    :type path: str
    :return: the conditions.
    :rtype: pd.DataFrame
    """
    df = pd.read_csv(path, sep=",", header=None)
    df.columns = ["condition", "path", "color"]
    return df


def aggregate_contacts(conditions, md_time, dir_path):
    """
    Extract the number of contacts by region for each sample in a condition.

    :param conditions: the conditions dataframe.
    :type conditions: pandas.DataFrame
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :return: the aggregated data for each frame and the conditions (in case ony condition is removed) and the domain of
    interest.
    :rtype: pandas.DataFrame, str
    """
    pattern_sample = re.compile("outliers_(.+)_ORF1.csv")
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
        raw_dict[row_condition.iloc[0]] = {}
        for item in sorted(by_condition):
            logging.info(f"\t\t- {item}")
            match_sample = pattern_sample.match(item)
            if match_sample:
                sample = match_sample.group(1)
            else:
                logging.error(f"No sample found with the pattern \"{pattern_sample.pattern}\" in the file {item}")
                sys.exit(1)
            raw_dict[row_condition.iloc[0]][sample] = {}
            df_current = pd.read_csv(os.path.join(row_condition["path"], item), sep=",")
            for _, row in df_current.iterrows():
                whole_domains.add(row["second partner domain"])
                roi_set.add(row["ROI partner domain"])
                if row["second partner domain"] not in raw_dict[row_condition.iloc[0]][sample]:
                    raw_dict[row_condition.iloc[0]][sample][row["second partner domain"]] = row["number atoms contacts"]
                else:
                    raw_dict[row_condition.iloc[0]][sample][row["second partner domain"]] += row["number atoms contacts"]

    # check if there is only one region of interest making contacts in all the files used
    if len(roi_set) == 1:
        roi = list(roi_set)[0]
    else:
        logging.error(f"More than one domain in the columns 'ROI partner domain' ({', '.join(list(roi_set))}) of the "
                      f"files in the directories provided by the input CSV file.")
        sys.exit(1)
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
    out_path = os.path.join(dir_path, f"contacts_aggregated_{roi}_{md_time}-ns.csv")
    df_out.to_csv(out_path, index=False)
    logging.info(f"Aggregated CSV file saved: {os.path.abspath(out_path)}")
    return df_out, roi


def order_x_axis(labels, domains_ordered):
    """
    Order the x axis values depending on the domains order.

    :param labels: the labels.
    :type labels: list
    :param domains_ordered: the ordered domains on the protein.
    :type domains_ordered: list
    :return: the X axis labels ordered.
    :rtype: list
    """
    ordered_labels = []
    # get the annotation before any domains
    for i in range(len(labels)):
        if labels[i].startswith("before"):
            ordered_labels.append(labels[i])
            break
    labels[:] = [x for x in labels if not x.startswith("before")]
    logging.debug("Reordering the X axis:")
    logging.debug(f"\tinitial X labels: {labels}")
    logging.debug(f"\tbefore domains:")
    logging.debug(f"\t\tnew X labels:\t\t{ordered_labels}")
    logging.debug(f"\t\tinitial X labels:\t{labels}")
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
            ordered_labels.append(labels[domain_index_in_labels["dom"]])
        if "between" in domain_index_in_labels:
            ordered_labels.append(labels[domain_index_in_labels["between"]])
        if "dom" in domain_index_in_labels:
            labels[:] = [x for x in labels if not x == dom]
        if "between" in domain_index_in_labels:
            labels[:] = [x for x in labels if not x.startswith(f"between {dom}")]
        logging.debug(f"\t\tnew X labels:\t\t{ordered_labels}")
        logging.debug(f"\t\tinitial X labels:\t{labels}")
    # get the annotation after all the domains
    for i in range(len(labels)):
        if labels[i].startswith("after"):
            ordered_labels.append(labels[i])
            break
    labels[:] = [x for x in labels if not x.startswith("after")]
    logging.debug(f"\tafter domains:")
    logging.debug(f"\t\tnew X labels:\t\t{ordered_labels}")
    logging.debug(f"\t\tinitial X labels:\t{labels}")
    
    return ordered_labels


def boxplot_aggregated(src, doi, md_time, dir_path, fmt, domains, subtitle):
    """
    Create the boxplots by conditions.

    :param src: the data source.
    :type src: pandas.DataFrame
    :param doi: domain of interest, the domain in contacts with the other domains.
    :type doi: str
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param domains: the ordered list of domains.
    :type domains: list
    :param subtitle: the subtitle of the plot.
    :type subtitle: str
    """
    plt.figure(figsize=(15, 15))
    # reorder the boxplots X axis if the domains argument is not null
    if domains:
        x_order = order_x_axis(list(set(src["domains"])), domains)
        logging.info(f"The boxplots were ordered as in the domains file: {', '.join(domains)}.")
    else:
        x_order = list(set(src["domains"]))
        logging.warning(f"No specific order provided for the boxplots.")

    # create the statistical pairs annotations
    boxplot_pairs = []
    conditions = list(set(src["conditions"]))
    for domain in x_order:
        for i in range(0, len(conditions) - 1):
            for j in range(i + 1, len(conditions)):
                boxplot_pairs.append(((domain, conditions[i]), (domain, conditions[j])))

    # creating the plotting parameters
    plotting_parameters = {
       "data": src,
        "x": "domains",
        "y": "contacts",
        "hue": "conditions",
        "order": x_order
    }

    # create the plot
    with sns.plotting_context():
        ax = sns.boxplot(**plotting_parameters, palette={"insertions": "red", "duplications": "orange", "WT": "cyan"})
        sns.stripplot(**plotting_parameters, size=8, marker="o", linewidth=2, dodge=True,
                      palette={"insertions": "darkred", "duplications": "chocolate", "WT": "blue"})
        annotator = Annotator(ax, boxplot_pairs, **plotting_parameters)
        annotator.configure(test="t-test_welch", text_format="star", hide_non_significant=True)
        annotator.apply_and_annotate()

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
        for label in labels[:3]:
            sample_set = set(src[src["conditions"] == label]["sample"])
            custom_labels.append(f"{label} ({len(sample_set)})")
        ax.legend(handles[:3], custom_labels, title="Condition")

        plt.suptitle(f"Contacts by domain with the {doi} at {md_time} ns of molecular dynamics", fontsize="large",
                     fontweight="bold")
        if subtitle:
            plt.title(subtitle)
        plt.xlabel("Domains", fontweight="bold")
        plt.ylabel(f"Number of contacts", fontweight="bold")
        plot = ax.get_figure()
        out_path_plot = os.path.join(dir_path, f"contacts_aggregated_{doi.lower().replace(' ', '-')}_"
                                               f"{md_time}-ns.{fmt}")
        plot.savefig(out_path_plot)
        logging.info(f"Aggregated contacts by condition: {os.path.abspath(out_path_plot)}")


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    Aggregate the outliers contacts by regions in one plot to compare between various conditions.

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
    parser.add_argument("-d", "--domains", required=False, type=str,
                        help="a sample CSV domains annotation file, to set the order of the protein domains on the X "
                             "axis. If this option is not used, the domains will be displayed randomly.")
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
    print(ordered_domains)
    data_conditions = get_conditions(args.input)
    df_contacts, domain_of_interest = aggregate_contacts(data_conditions, args.md_time, args.out)
    boxplot_aggregated(df_contacts, domain_of_interest, args.md_time, args.out, args.format, ordered_domains,
                       args.subtitle)
