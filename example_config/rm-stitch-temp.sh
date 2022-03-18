#!/bin/sh
# This script removes all temporary files created by the non-rigid stitching program.
# WARNING: Run only on directories that contain stitching projects, otherwise some
# unrelated files may be deleted (see rm commands below)!

rm ./*_refpoints.txt
rm ./*-err.txt
rm ./*-out.txt
rm ./*_transformation.txt
rm ./*_gof_*.raw
rm ./*_defpoints_*.raw
rm ./*_world_to_local_shifts_*.raw
rm ./*_index.txt
rm ./*_done.tif
rm -r ./slurm-io-files