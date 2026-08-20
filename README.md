# sturgeon-v2
This repository presents a pip installable software package to easily run Sturgeon V2 in a clinical setting. Sturgeon V2 was trained using a relabeling strategy to make use of unlabeled data. Sturgeon V2 makes use of an ensemble of submodels. The average of these submodels is calibrated based on the sequencing depth.

## Installation
* Install using `pip install git+ssh://git@github.com/UMCUGenetics/sturgeon-v2.git`
* Download the cns-v2.4b [zip file](https://drive.google.com/file/d/18cHuUTN_7TyOZzcbaIGeVrugoYCHNkOE/view?usp=share_link)
* Download the cns-v2.4a (deprecated) [zip file](https://drive.google.com/file/d/1opEUbGWtUCONWhlfrv9NfSuGM8pXj-D8/view?usp=share_link)

## Usage
<pre>
usage: sturgeon-v2 [-h] {single,multi} ...

Sturgeon v2: Intraoperative methylation classification

positional arguments:
  {single,multi}
    single        Predict on a single .bed or .csv file
    multi         Predict on a directory of .bed or .csv files

Usage: sturgeon-v2 single/multi [-h] [-i INPUT_FILE] [-m MODEL] [-c NUM_CPGS] [-o OUTPUT_FILE] [-s NUM_SAMPLINGS]
                   [-f FILE_TYPE]

Options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        path to input file or directory
  -m MODEL, --model MODEL
                        specify path to model zip you want to use.
  -c NUM_CPGS, --num_CpGs NUM_CPGS
                        specify the number of samples CpG sites (default is to use all available sites).
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        input
  -s NUM_SAMPLINGS, --num_samplings NUM_SAMPLINGS
                        Specify the number of random samples of size num_CpGs. Default is 1 random
                        sampling.
  -f FILE_TYPE, --file_type FILE_TYPE
                        input file type 'bed' or 'csv'
</pre>
