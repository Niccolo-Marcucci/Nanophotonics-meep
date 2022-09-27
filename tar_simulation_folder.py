#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys

hashtag = sys.argv[1]

files = os.listdir("./data")
found = False
for file in files :
    if file.find( hashtag ) >= 0 and os.path.isdir('data/' + file):
        folder = 'data/' + file
        found = True

if not found :
    raise FileExistsError(f"Missing simulation folder for hash {hashtag}")

os.system(f"tar -czvf {folder}.tar.gz -C {folder} .")


