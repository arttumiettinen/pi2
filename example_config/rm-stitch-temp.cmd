@echo off
rem This script removes all temporary files created by the non-rigid stitching program.
rem WARNING: Run only on directories that contain stitching projects, otherwise some
rem unrelated files may be deleted (see del commands below)!

del *_refpoints.txt
del *_refpoints_settings.txt
del *_transformation.txt
del *_gof_*.raw
del *_defpoints_*.raw
del *_world_to_local_shifts_*.raw
del *_index.txt
del *_done.tif
del *_mosaic_settings.txt
del *_transformations_settings.txt
del *_position_settings.txt
rem Slurm I/O folder is not deleted here as Slurm is never available in Windows anyway.