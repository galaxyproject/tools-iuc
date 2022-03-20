. <(planemo conda_env pyGenomeTracks.xml)
# After
# planemo t  --no_cleanup --galaxy_python_version 3.8 --galaxy_root $galaxyF --job_output_files outputs
python test-data/compare_images.py
conda_env_deactivate
