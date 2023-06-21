# Gene-based rare variant burden test

This repository hosts the code and related data used for gene-based rare variant burden test。

## Repository Structure

This repository is divided into six main folders, each serving a specific purpose in the analysis pipeline:

1. `Data`: This folder contains the raw and processed data used in the study.

2. `Genotype_filtering`: This folder contains scripts for filtering genotype data. It ensures the inclusion of high-quality data by filtering out low-quality or inconsistent genotypes.

3. `Individual_filtering`: The scripts here are used for filtering individuals in the study. This helps remove outliers and those individuals who do not meet the specific inclusion criteria for the study.

4. `Qualifying_variant`: The code in this folder is designed to identify qualifying variants for the study. These are the genetic variants that are believed to contribute to congenital vertebral malformation.

5. `Burden_test`: This folder contains scripts for conducting the gene-based burden test. The burden test is used to detect the accumulation of rare variants in a given gene across individuals.

6. `Visualization`: The scripts in this folder are used for visualizing the results of the study. This includes generating graphs, charts, and other visual representations to better understand and interpret the data.

In addition, a shell file `pipeline.sh` is included at the root of the repository. This file serves as an executable script to run the entire analysis pipeline in sequence.

## Setup

To set up and run the scripts, you need to:

1. Clone the repository.
2. Install the required dependencies mentioned in the code.
3. Run the `pipeline.sh` script to perform all the analysis in a streamlined manner.

## Contributing

We welcome contributions to this repository. If you find an error or a way to improve the scripts, please create an issue or a pull request.

## Contact

For any questions or concerns related to the code, please feel free to reach out to:

- Sen Zhao: zhaosen830@163.com
- Hengqiang Zhao: zhaohq921@gmail.com

Please mention the title of the study in your subject line for a prompt response.
