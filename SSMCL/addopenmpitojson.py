import os
from os.path import exists, join 
import json
import sys
if exists('build/compile_commands.json'):
    outjs = []
    with open('build/compile_commands.json', 'r') as f:
        js = json.load(f)
        for i in range(len(js)):
            sl = js[i]['command'].split(" ")
            sl.insert(1, "-I/global/common/software/nersc/pe/gnu/openmpi/5.0.0rc12/")
            js[i]['command'] = " ".join(sl)
    with open('build/compile_commands.json', 'w') as f:
        json.dump(js,f)
