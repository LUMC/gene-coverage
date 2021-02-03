# Copyright (C) 2020 Leiden University Medical Center
# This file is part of gene-coverage
#
# gene-coverage is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# gene-coverage is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with gene-coverage.  If not, see <https://www.gnu.org/licenses/

from pathlib import Path

from setuptools import setup, find_packages

setup(name="gene-coverage",
      version="1.1.0-dev",
      description="Chunk and scatter the regions in a bed or sequence dict "
                  "file",
      long_description=Path("README.rst").read_text(),
      long_description_content_type='text/x-rst',
      classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: "
        "GNU Affero General Public License v3 or later (AGPLv3+)",
      ],
      python_requires=">=3.6",
      # This pysam version supports vcf and not much changes to the interface
      # after this release.
      # pysam has much less dependencies than cyvcf2
      install_requires=["pybedtools", "v-tools", "numpy"],
      keywords="bioinformatics gene coverage transcript exon",
      url="https://github.com/lumc/gene-coverage",
      author="Leiden University Medical Center",
      author_email="sasc@lumc.nl",
      zip_safe=False,
      license="AGPL-3.0-or-later",
      packages=find_packages('src'),
      package_dir={'': 'src'},
      entry_points={
          "console_scripts":
              ["gene-coverage=gene-coverage:main"]
      })
