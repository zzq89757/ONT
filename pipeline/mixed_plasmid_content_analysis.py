import argparse
import re
import subprocess
import time
from os import system
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
from docx.shared import Mm
from docxtpl import DocxTemplate, InlineImage
from pysam import index, depth, flagstat, FastqFile, VariantFile

from sequencing_analysis import quality_check, map2reference



def main() -> None:
    quality_check("/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/run/mixed/AX09-1.fastq","/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/run/mixed/")
    map2reference("/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/run/mixed/")


if __name__ == "__main__":
    main()