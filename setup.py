#!/usr/bin/env python 
# Copyright (c) 2015 
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License.
import shutil
import sys
import os

try:
	if sys.argv[1] == "install":
		command1 = "cp -rf Emotif_alpha.py /usr/local/bin/Emotif_alpha"
		command2 = "cp -rf Emotif /usr/local/lib/python2.7/dist-packages/"
		os.system(command1)
		os.system(command2)
except:
	print '''
	Usage: sudo python setup.py install
	
	'''




















