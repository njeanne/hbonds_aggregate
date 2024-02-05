# Aggregate the contacts

From the CSV files describing the amino acids contacts between a Region of Interest and the other regions of a protein 
during a Molecular Dynamics simulation, the script will produce a plot of the boxplots representing those contacts by 
domains and conditions.

The input CSV data are produced by the [plot_contacts](https://github.com/njeanne/plot_contacts/tree/main) script.

## Conda environment

A [conda](https://docs.conda.io/projects/conda/en/latest/index.html) YAML environment file is provided: 
`conda_env/python3_env.yml`. The file contains all the dependencies to run the script.
The conda environment is generated using the command:
```shell script
# create the environment
conda env create -f conda_env/python3_env.yml

# activate the environment
conda activate python3
```

## Usage

The script can be tested with the test data provided in the `data` directory, which contains a CSV file describing the 
different conditions and the location of the directory containing the CSV output files from the `plot_contacts.py` 
script. 
The commands to use the script are:

```shell script
conda activate python3

./contacts_aggregate.py --md-time 1002 --subtitle "annotations Koonin" --out results data/conditions.csv

conda deactivate
```

## Outputs

The script outputs are:
- boxplots of the contacts by conditions and domains:

![contacts heatmap](doc/_static/boxplots.svg)

- a CSV file of the contacts by condition and domain.
