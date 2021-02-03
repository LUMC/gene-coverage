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
import itertools
import os
from typing import Generator, Union, Iterable

import numpy as np
from pybedtools.bedtool import BedTool
from vtools.gcoverage import RefRecord


def file_to_refflat_records(filename: Union[str, os.PathLike]
                            ) -> Generator[RefRecord, None, None]:
    """Generate a Refflat record for each line in the refflat file"""
    with open(filename, "rt") as file_h:
        for line in file_h:
            yield RefRecord.from_line(line)


def refflat_record_to_exons_bedtool(refflat_record: RefRecord) -> BedTool:
    """
    Convert a refflat record to a BedTool object with all the exons
    """
    bed_regions = ((
        f"{exon.chr}\t{exon.start}\t{exon.end}\t"
        f"{refflat_record.transcript}\t.\t"  # name, score
        f"{'+' if refflat_record.forward else '-'}\n")  # strand
        for exon in refflat_record.exons)
    return BedTool("".join(bed_regions), from_string=True)


def feature_coverage_multisample_to_array(
        feature: BedTool, samples: Iterable[BedTool]) -> np.ndarray:
    """

    :param feature: The BedTool with the features of interest
    :param samples: Sample files (as BedTools) that the coverage should be
                    calculated over.
    :return: A numpy array with coverage at all the positions for each sample.
             For example if the feature has 800 bp in total and there are 3
             samples, the numpy array will contain 3 x 800 = 2400 positions.
    """
    # This way the true median, mode and percentage histogram can be calculated
    # correctly. And the numpy array only has to be constructed once. (Not
    # merging multiple independently created numpy arrays.)

    # Get the counts from the bedtools coverage result and mash all the samples
    # together in one array.
    counts = (pos.count for pos in
              (feature.coverage(sample, d=True, stream=True)
               for sample in samples))
    # np.uint16 has a max of 65535. That should be enough. Weird stuff may
    # happen if coverages is higher than that though.
    # TODO: If necessary benchmark the performance difference between uint16
    # TODO: and uint32. Choose uint32 when there is a negligible difference.
    return np.fromiter(counts, dtype=np.uint16)


def coverage_stats_from_base_coverages(base_coverages: np.ndarray):
    """
    Calculate the median, fractions of at least 10, 20, 30 and 50 coverage
    from a numpy array with coverage values.
    :param base_coverages: A numpy array with coverage for each base of the
                           feature.
    :return: Tuple with median, fraction >= 10, >= 20, >= 30, >=50.
    """
    total = base_coverages.size
    median = np.median(base_coverages)
    max_value = np.iinfo(base_coverages.dtype).max
    counts, _ = np.histogram(base_coverages, [0, 10, 20, 30, 50, max_value])
    # >= 50 includes only the last bin. >= 30 includes 50-maxvalue and 30-50.
    # by using accumulate  and reversed we calculate the sums of these bins at
    # once.
    reverse_cumulative_counts = list(itertools.accumulate(reversed(counts)))
    fraction_gt_50 = reverse_cumulative_counts[0] / total
    fraction_gt_30 = reverse_cumulative_counts[1] / total
    fraction_gt_20 = reverse_cumulative_counts[2] / total
    fraction_gt_10 = reverse_cumulative_counts[3] / total
    return (median, fraction_gt_10, fraction_gt_20, fraction_gt_30,
            fraction_gt_50)
