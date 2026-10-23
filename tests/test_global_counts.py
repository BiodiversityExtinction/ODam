#!/usr/bin/env python3

import argparse
import contextlib
import csv
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pysam


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import ODam  # noqa: E402


class GlobalCountTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory(prefix="odam_test_")
        self.work = Path(self.tempdir.name)
        self.reference = self.work / "reference.fa"
        self.bam = self.work / "reads.bam"
        self._write_fixture()

    def tearDown(self):
        self.tempdir.cleanup()

    def _write_fixture(self):
        reference_sequence = "GC" * 1500
        self.reference.write_text(f">chr1\n{reference_sequence}\n", encoding="ascii")
        pysam.faidx(str(self.reference))

        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": len(reference_sequence)}],
        }
        records = []

        def add_read(name, length, start, flag=0, mapq=60, low_quality_test=False):
            segment = pysam.AlignedSegment()
            segment.query_name = name
            segment.flag = flag
            segment.reference_id = 0
            segment.reference_start = start
            segment.mapping_quality = mapq
            segment.cigartuples = [(0, length)]
            sequence = list(reference_sequence[start : start + length])
            # One of each tracked transversion per eligible read.
            sequence[0] = "T"  # G>T
            sequence[1] = "A"  # C>A
            sequence[2] = "C"  # G>C
            sequence[3] = "G"  # C>G
            if low_quality_test:
                sequence[4] = "T"  # Excluded G>T because base quality is low.
            segment.query_sequence = "".join(sequence)
            qualities = [40] * length
            if low_quality_test:
                qualities[4] = 10
            segment.query_qualities = qualities
            records.append(segment)

        # Two forward and two reverse alignments with different read lengths.
        add_read("read50", 50, 0, low_quality_test=True)
        add_read("read100_reverse", 100, 250, flag=16)
        add_read("read150", 150, 500)
        add_read("read200_reverse", 200, 750, flag=16)

        # These would strongly alter the result if existing read filters changed.
        add_read("low_mapq", 50, 1000, mapq=10)
        add_read("duplicate", 50, 1100, flag=1024)
        add_read("secondary", 50, 1200, flag=256)
        add_read("supplementary", 50, 1300, flag=2048)
        add_read("qc_fail", 50, 1400, flag=512)

        records.sort(key=lambda record: record.reference_start)
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as output:
            for record in records:
                output.write(record)
        pysam.index(str(self.bam))

    def _compute(self, max_pos=10, max_reads=0, random_seed=1, normalize_ends=True):
        return ODam.compute_counts(
            bam_path=str(self.bam),
            ref_fasta_path=str(self.reference),
            max_pos=max_pos,
            min_mapq=30,
            min_baseq=30,
            max_reads=max_reads,
            threads=1,
            normalize_ends=normalize_ends,
            random_seed=random_seed,
        )

    @staticmethod
    def _args(**overrides):
        values = {
            "reuse_existing_pos": False,
            "max_pos": 10,
            "min_mapq": 30,
            "min_baseq": 30,
            "max_reads": 0,
            "threads": 1,
            "normalize_ends": True,
            "random_seed": 1,
            "region": None,
            "pseudocount": 0.5,
            "window": 3,
            "plot_max_pos": 5,
            "plot_log_y": False,
            "plot_y_max": None,
            "sample_list": None,
        }
        values.update(overrides)
        return argparse.Namespace(**values)

    def test_each_eligible_aligned_base_is_counted_once_globally(self):
        profiles, global_counts, reads_used, diagnostics = self._compute(max_pos=100)

        self.assertEqual(reads_used, 4)
        self.assertEqual(diagnostics["missing_reference_contigs"], 0)
        self.assertEqual(
            global_counts,
            {"G": 249, "C": 250, "G>T": 4, "G>C": 4, "C>A": 4, "C>G": 4},
        )
        self.assertAlmostEqual(ODam.mismatch_proportion(global_counts["G>T"], global_counts["G"]), 4 / 249)
        self.assertAlmostEqual(ODam.mismatch_proportion(global_counts["C>A"], global_counts["C"]), 4 / 250)

        # With max_pos=100, short-read positions occur in both end profiles.
        # Their aggregate therefore must not be used as the global count.
        positional_total = ODam.aggregate_counts(
            list(profiles["3p"].values()) + list(profiles["5p"].values())
        )
        self.assertNotEqual(positional_total["G"], global_counts["G"])

    def test_global_counting_respects_cigar_alignment_operations(self):
        reference_sequence = "GC" * 100
        reference = self.work / "cigar_reference.fa"
        bam = self.work / "cigar.bam"
        reference.write_text(f">chr1\n{reference_sequence}\n", encoding="ascii")
        pysam.faidx(str(reference))

        segment = pysam.AlignedSegment()
        segment.query_name = "cigar_read"
        segment.flag = 0
        segment.reference_id = 0
        segment.reference_start = 0
        segment.mapping_quality = 60
        # 2S 4M 1I 4M 1D 4M 2S: soft clips, insertion, and deletion must
        # not contribute to whole-read reference-base opportunities.
        segment.cigartuples = [(4, 2), (0, 4), (1, 1), (0, 4), (2, 1), (0, 4), (4, 2)]
        sequence = ["T"] * 17
        for qpos, rpos in segment.get_aligned_pairs(matches_only=True):
            sequence[qpos] = reference_sequence[rpos]
        aligned_pairs = segment.get_aligned_pairs(matches_only=True)
        for (qpos, _), mismatch in zip(aligned_pairs[:4], ("T", "A", "C", "G")):
            sequence[qpos] = mismatch
        segment.query_sequence = "".join(sequence)
        segment.query_qualities = [40] * len(sequence)

        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": len(reference_sequence)}],
        }
        with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
            output.write(segment)
        pysam.index(str(bam))

        _, global_counts, reads_used, _ = ODam.compute_counts(
            bam_path=str(bam),
            ref_fasta_path=str(reference),
            max_pos=20,
            min_mapq=30,
            min_baseq=30,
            max_reads=0,
            threads=1,
            normalize_ends=True,
            random_seed=1,
        )
        self.assertEqual(reads_used, 1)
        self.assertEqual(
            global_counts,
            {"G": 6, "C": 6, "G>T": 1, "G>C": 1, "C>A": 1, "C>G": 1},
        )

    def test_global_counts_ignore_max_pos_and_existing_profiles_are_stable(self):
        profiles5, global5, _, _ = self._compute(max_pos=5)
        profiles100, global100, _, _ = self._compute(max_pos=100)

        self.assertEqual(global5, global100)
        for end in ("5p", "3p"):
            self.assertEqual(dict(profiles5[end]), {p: profiles100[end][p] for p in range(1, 6)})

    def test_reverse_reads_retain_strand_normalized_end_assignment(self):
        profiles, _, _, _ = self._compute(max_pos=5, normalize_ends=True)
        raw_profiles, _, _, _ = self._compute(max_pos=5, normalize_ends=False)

        self.assertEqual(profiles["5p"][1]["G>T"], 2)
        self.assertEqual(profiles["3p"][1]["G>T"], 2)
        self.assertEqual(raw_profiles["5p"][1]["G>T"], 4)
        self.assertEqual(raw_profiles["3p"][1]["G>T"], 0)

    def test_random_sampling_is_reproducible(self):
        first = self._compute(max_pos=10, max_reads=2, random_seed=27)
        second = self._compute(max_pos=10, max_reads=2, random_seed=27)

        self.assertEqual(first[1:], second[1:])
        for end in ("5p", "3p"):
            self.assertEqual(dict(first[0][end]), dict(second[0][end]))

    def test_plot_max_pos_does_not_change_global_summary(self):
        with contextlib.redirect_stdout(io.StringIO()):
            row2 = ODam.analyze_sample(
                "sample", str(self.bam), str(self.reference), self._args(plot_max_pos=2)
            )
            row5 = ODam.analyze_sample(
                "sample", str(self.bam), str(self.reference), self._args(plot_max_pos=5)
            )
        global_fields = [
            "g_bases_observed",
            "c_bases_observed",
            "gt_count",
            "gc_count",
            "ca_count",
            "cg_count",
            "overall_gt_over_gc",
            "overall_ca_over_cg",
            "overall_combined",
            "overall_p_gt_over_g",
            "overall_p_gc_over_g",
            "overall_p_ca_over_c",
            "overall_p_cg_over_c",
        ]
        self.assertEqual({k: row2[k] for k in global_fields}, {k: row5[k] for k in global_fields})

    def test_reuse_position_file_rescans_bam_for_global_counts(self):
        pos_path = self.work / "cached_pos.tsv"
        with contextlib.redirect_stdout(io.StringIO()):
            original = ODam.analyze_sample(
                "sample", str(self.bam), str(self.reference), self._args(), str(pos_path)
            )

        rows = ODam.read_pos_tsv(str(pos_path))
        for row in rows:
            for key in ("G", "C", "G>T", "G>C", "C>A", "C>G"):
                row[key] = 999
        ODam.write_pos_tsv(str(pos_path), rows)

        with contextlib.redirect_stdout(io.StringIO()):
            reused = ODam.analyze_sample(
                "sample",
                str(self.bam),
                str(self.reference),
                self._args(reuse_existing_pos=True),
                str(pos_path),
            )
        self.assertEqual(reused["status"], "ok_reused_pos")
        for key in (
            "g_bases_observed",
            "c_bases_observed",
            "gt_count",
            "gc_count",
            "ca_count",
            "cg_count",
            "overall_combined",
        ):
            self.assertEqual(reused[key], original[key])

        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises((FileNotFoundError, OSError)):
                ODam.analyze_sample(
                    "sample",
                    str(self.work / "missing.bam"),
                    str(self.reference),
                    self._args(reuse_existing_pos=True),
                    str(pos_path),
                )

    def test_batch_and_single_summary_outputs_are_consistent(self):
        single_summary = self.work / "single.tsv"
        batch_summary = self.work / "batch.tsv"
        sample_list = self.work / "samples.tsv"
        sample_list.write_text(
            f"sample\tbam\treference\nsynthetic\t{self.bam}\t{self.reference}\n",
            encoding="utf-8",
        )
        common = [
            "--max-reads",
            "0",
            "--max-pos",
            "10",
            "--plot-max-pos",
            "5",
            "--min-mapq",
            "30",
            "--min-baseq",
            "30",
        ]
        subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "ODam.py"),
                "--bam",
                str(self.bam),
                "--reference",
                str(self.reference),
                "--sample-label",
                "synthetic",
                "--summary-tsv-out",
                str(single_summary),
                *common,
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "ODam.py"),
                "--sample-list",
                str(sample_list),
                "--batch-summary-out",
                str(batch_summary),
                *common,
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )

        with single_summary.open(newline="", encoding="utf-8") as handle:
            single_row = next(csv.DictReader(handle, delimiter="\t"))
        with batch_summary.open(newline="", encoding="utf-8") as handle:
            batch_row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(single_row, batch_row)


if __name__ == "__main__":
    unittest.main()
