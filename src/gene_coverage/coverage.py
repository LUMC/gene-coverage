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
from typing import Generator

from pybedtools.bedtool import BedTool
from vtools.gcoverage import RefRecord


def file_to_refflat_records(filename: str) -> Generator[RefRecord, None, None]:
    with open(filename, "rt") as file_h:
        for line in file_h:
            yield RefRecord.from_line(line)


def refflat_record_to_bedtool(refflat_record: RefRecord) -> BedTool:
    bed_regions = ((
        f"{exon.chr}\t{exon.start}\t{exon.end}\t"
        f"{refflat_record.transcript}\t.\t"  # name, score
        f"{'+' if refflat_record.forward else '-'}\n")  # strand
        for exon in refflat_record.exons)
    return BedTool("".join(bed_regions), from_string=True)
