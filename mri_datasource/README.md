## How to convert datasets into the hpc-predic-io format

### Bern experimental dataset

Download Bern experimental data and convert it to hpc-predict-io format
```
source python/venv/bin/activate
cd mri_datasource
./fetch_bern_experimental_data.sh
python generate_bern_expt_metadata.py --input bern_data_experiments_source/ --output metadata.json
PYTHONPATH=../python/ python  convert_bern_expt_to_hpc_predict.py --input metadata.json --output bern_data_experiments_hpc_predict/
```

To generate an IMPACT config.txt file, then use the following command (feel free to adjust the numerical parameters)

```
cd bern_data_experiments_hpc_predict/
PYTHONPATH=../../ python ../../python/mr_io_impact_config.py --mri bern_experimental_dataset_flow_mri.h5 --sr 2 2 2 --padding 0.5 0.5 0.5 --tr 2 --config ../../python/config.txt.j2 --output config.txt --np 4
```

As of the first version, you can then copy both config.txt and bern_experimental_dataset_flow_mri.h5 to the prog directory of IMPACT and start impact_debug.exe from there.
