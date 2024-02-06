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
    :return: the aggregated data for each frame and the conditions (in case ons condition is removed).
    :rtype: pandas.DataFrame, pandas.DataFrame
    """
    pattern_sample = re.compile("outliers_(.+)_ORF1.csv")
    raw_dict = {}
    whole_domains = set()
    conditions_to_remove = []
    for _, row_condition in conditions.iterrows():
        by_condition = [fn for fn in os.listdir(row_condition["path"]) if fn.startswith("outliers") and fn.endswith(".csv")]
        if len(by_condition) == 0:
            conditions_to_remove.append(row_condition["condition"])
            logging.warning(f"Condition {row_condition['condition']}: no RMSD files, this condition is skipped.")
            continue
        logging.info(f"Aggregating {len(by_condition)} file{'s' if len(by_condition) > 1 else ''} data for condition: "
                     f"{row_condition['condition']}")
        raw_dict[row_condition[0]] = {}
        for item in sorted(by_condition):
            logging.info(f"\t\t- {item}")
            match_sample = pattern_sample.match(item)
            if match_sample:
                sample = match_sample.group(1)
            else:
                logging.error(f"No sample found with the pattern \"{pattern_sample.pattern}\" in the file {item}")
                sys.exit(1)
            raw_dict[row_condition[0]][sample] = {}
            df_current = pd.read_csv(os.path.join(row_condition["path"], item), sep=",")
            for _, row in df_current.iterrows():
                whole_domains.add(row["second partner domain"])
                if row["second partner domain"] not in raw_dict[row_condition[0]][sample]:
                    raw_dict[row_condition[0]][sample][row["second partner domain"]] = row["number atoms contacts"]
                else:
                    raw_dict[row_condition[0]][sample][row["second partner domain"]] += row["number atoms contacts"]

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
    out_path = os.path.join(dir_path, f"contacts_aggregated_{md_time}-ns.csv")
    df_out.to_csv(out_path, index=False)
    logging.info(f"Aggregated CSV file saved: {os.path.abspath(out_path)}")
    return df_out

def boxplot_aggregated(src, md_time, dir_path, fmt, subtitle):
    """
    Create the boxplots by conditions.

    :param src: the data source.
    :type src: pandas.DataFrame
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param subtitle: the subtitle of the plot.
    :type subtitle: str
    """
    plt.figure(figsize=(15, 15))
    ax = sns.boxplot(data= src, x="domains", y="contacts", hue="conditions",
                             palette={"insertions": "red", "duplications": "orange", "WT": "cyan"})
    sns.stripplot(data= src, x="domains", y="contacts", size=8, hue="conditions", marker="o",
                                 linewidth=2, dodge=True, edgecolor="gray",
                                 palette={"insertions": "darkred", "duplications": "chocolate", "WT": "blue"})
    x_labels = ax.get_xticklabels()
    new_x_labels = [re.sub(r'(\w+ \w+ \w+)( )',r'\1\n', x.get_text()) for x in x_labels]
    ax.set_xticklabels(new_x_labels)
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment="right")

    # remove extra legend handles and add the count of samples by condition
    handles, labels = ax.get_legend_handles_labels()
    custom_labels = []
    for label in labels[:3]:
        sample_set = set(src[src["conditions"] == label]["sample"])
        custom_labels.append(f"{label} ({len(sample_set)})")
    ax.legend(handles[:3], custom_labels, title="Condition")

    plt.suptitle(f"Contacts by domain at {md_time} ns", fontsize="large", fontweight="bold")
    if subtitle:
        plt.title(subtitle)
    plt.xlabel("Domains", fontweight="bold")
    plt.ylabel(f"Number of contacts", fontweight="bold")
    out_path_plot = os.path.join(dir_path, f"contacts_aggregated_{md_time}-ns.{fmt}")
    plot = ax.get_figure()
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

    data_conditions = get_conditions(args.input)
    df = aggregate_contacts(data_conditions, args.md_time, args.out)
    boxplot_aggregated(df, args.md_time, args.out, args.format, args.subtitle)
