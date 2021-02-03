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
from typing import NamedTuple, Generator

from pybedtools.bedtool import BedTool
from pybedtools.bedtool import Interval
from vtools.gcoverage import Region as OneBasedRegion
from vtools.gcoverage import RefRecord


def refflat_to_bedtools(filename: str) -> Generator[BedTool, None, None]:
    with open(filename, "rt") as file_h:
        for line in file_h:
            ref_record = RefRecord.from_line(line)
            bed_regions = ((
                f"{exon.chr}\t{exon.start}\t{exon.end}\t"
                f"{ref_record.transcript}\t.\t"  # name, score
                f"{'+' if ref_record.forward else '-'}\n")  # strand
                for exon in ref_record.exons)
            yield BedTool("".join(bed_regions), from_string=True)
