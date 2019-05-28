from scipy.stats import zscore
import pandas as pd
import os
import argparse
import itertools
import getpass
import uuid
import datetime
import logging
import uuid
import numpy as np
from Bio.Seq import Seq
from joblib import Parallel, delayed
from StringIO import StringIO
import pandas as pd
import sys
from PIL import Image
import numpy as np
import os.path
import re

def write_fasta(file_name,myDict):
	out = open(file_name,"wt")
	for k in myDict:
		out.write(">"+k+"\n")
		out.write(myDict[k]+"\n")
	out.close()
	
	
def setup_custom_logger(jid):
	logFormatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')
	rootLogger = logging.getLogger("root")
	fileHandler = logging.FileHandler(jid+".log")
	fileHandler.setFormatter(logFormatter)
	rootLogger.addHandler(fileHandler)
	consoleHandler = logging.StreamHandler()
	consoleHandler.setFormatter(logFormatter)
	rootLogger.addHandler(consoleHandler)
	rootLogger.setLevel(logging.INFO)
	return rootLogger

