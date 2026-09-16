"""Compile and run the original PSM examples in separate build directories."""

from pathlib import Path
import shutil
import subprocess
import sys
from split_adgj import split

root = Path(__file__).resolve().parent
kinds = sys.argv[1:] or ['EE', 'EO', 'OE', 'OO']
if any(kind not in ('EE', 'EO', 'OE', 'OO') for kind in kinds):
    sys.exit('Usage: python PSM/run_examples.py [EE EO OE OO]')
split(root / 'adgj', root / 'split')
for kind in kinds:
    source = root / 'split' / kind
    work = root / 'build' / kind
    work.mkdir(parents=True, exist_ok=True)
    executable = work / f'PSM_{kind}.exe'
    subprocess.run(['gfortran', '-O2', '-std=legacy',
                    '-fallow-argument-mismatch', '-I', str(source),
                    str(source / f'PSM_{kind}.f'), '-o', str(executable)],
                   check=True)
    shutil.copyfile(source / f'{kind}_DATA', work / f'{kind}_DATA')
    with (work / 'run.log').open('w') as log:
        subprocess.run([str(executable)], cwd=work, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    print(f'{kind}: process completed; outputs in {work}', flush=True)
