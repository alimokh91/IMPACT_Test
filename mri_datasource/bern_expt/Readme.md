## How to convert datasets into the hpc-predic-io format

### Bern experimental dataset

Download Bern experimental data and convert it to hpc-predict-io format
```
source python/venv/bin/activate
cd mri_datasource
../../../data/fetch_scripts/fetch_bern_experiments_assimilation_in.sh
python generate_bern_exp_meta_info.py --data bern_data_experiments_source --output bern_exp_metadata.json
PYTHONPATH=../../python/ python  convert_bern_expt_to_hpc_predict.py --input bern_data_experiments_source/bern_exp_metadata.json --output bern_data_experiments_hpc_predict/bern_experimental_dataset_flow_mri.h5 --downsample 1 --single_rep -1

Downsample factor and option to save only an instantaneous velocity field can be changed via the arguments
```

To generate an IMPACT config.txt file, then use the following command (feel free to adjust the numerical parameters)

```
cd bern_data_experiments_hpc_predict/
PYTHONPATH=../../../ python ../../../python/mr_io_impact_config.py --input-mri bern_experimental_dataset_flow_mri.h5 --output-mri  bern_experimental_dataset_assimilation_results.h5 --sr 2 2 2 --padding 0.5 0.5 0.5 --tr 2 --config ../../../python/config.txt.j2 --output config.txt --np 4
```

As of the first version, you can then copy both config.txt and bern_experimental_dataset_flow_mri.h5 to the prog directory of IMPACT and start impact_debug.exe from there or start impact_debug.exe directly in the bern_data_experiments_hpc_predict directory (make sure it is not spinning in the gdb-debug-loop) using

```
mpiexec -np 4 ../../../IMPACT/prog/impact_debug.exe
```
