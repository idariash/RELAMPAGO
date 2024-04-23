#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 18:50:37 2020

@author: idariash
"""

import tempfile
import shutil

dirpath = tempfile.mkdtemp()
# ... do stuff with dirpath

shutil.rmtree(dirpath)

