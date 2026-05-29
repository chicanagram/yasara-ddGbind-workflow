# **********************************************************
# *                                                        *
# *                    Y  A  S  A  R  A                    *
# *                                                        *
# * Yet Another Scientific Artificial Reality Application  *
# *                                                        *
# **********************************************************
# *  yasara.py - The loader for the YASARA Python module   *
# * You can find the complete documentation of this module *
# *    in YASARA's HTML manual yasara/doc/index.html at    *
# *       'Scripts - Use YASARA as a Python module'.       *
# *      To get a list of all defined functions, type      *
# *                 help('yasaramodule')                   *
# *          License for this Python module: BSD           *
# **********************************************************

# IMPORTANT: Copy this module to a place where Python can find
# it, e.g. the directory where you keep your own Python modules.
# If you move YASARA somewhere else, please adapt the path below:

import platform
opsys = platform.system()
version_suffix = '2026'

if opsys=='Linux':
    yasaradir = f'/home/cchiajt-bii/yasara{version_suffix}/'
elif opsys=='Windows':
    yasaradir = f"C:/Users/cchiajt-bii/Documents/projects/yasara{version_suffix}/"
elif opsys=='Darwin':
    yasaradir = '/Applications/YASARA.app/Contents/yasara/'
print('YASARA directory:', yasaradir)

import sys,os
sys.path.append(os.path.join(yasaradir,'pym'))
sys.path.append(os.path.join(yasaradir,'plg'))
from yasaramodule import *
