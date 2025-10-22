import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'generate_flarelist_python'))

from generate_flarelist_python.flarelist_generate import get_flares

flares = get_flares('2023-01-01', '2023-12-31', 'raw_fits')
