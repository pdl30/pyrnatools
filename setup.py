import os
from setuptools import setup, find_packages

setup(name='pyrnatools',
      version='0.0.1',
      description='pyrnatools is a Python module to analyze ChIP-seq NGS data',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=find_packages(),
      package_data={"pyrnatools":['data/*']},
      scripts=['scripts/pyrna_align.py', 'scripts/pyrna_assembly.py', 'scripts/pyrna_diff.py', 'scripts/pyrna_count.py','scripts/pyrna_pair_tool.py','scripts/pyrna_ucsc.py', 'scripts/pyrna_download.py'],
      install_requires=['pysam', 'pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Environment :: Console',
      ],
      long_description="""

pyrnatools is a Python module to analyze RNA-seq NGS data

 Contact
=============

If you have any questions or comments about pyrnatools, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
