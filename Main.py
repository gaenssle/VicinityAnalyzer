#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

import os
import re
from multiprocessing import Pool
import pandas as pd
import argparse

# Own modules
import Import_Export as IE
import Download_GenomeJP as Genome
import Download_KEGG as KEGG

