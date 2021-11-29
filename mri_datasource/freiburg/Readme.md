# DICOM interface for HPC-PREDICT

This documents the code to integrate clinical flow MRI data into the HPC-PREDICT pipeline by providing a traceable DICOM-to-HDF5 conversion and vice-versa with dataset-wide analyses to guide data cleaning steps.  

## Setup

To run the Jupyter notebooks implemented here, install

```
pip install -r requirements.txt
```

## DICOM Reader 

Dataset-wide analysis is available from `freiburg_full_dataset_analysis.ipynb`.

Single-study conversion is implemented in Jupyter notebook `freiburg_dicom_reader.ipynb`. To run, use

```
../../docker/dicom_to_hpc_predict_io_conversion.sh
```

or in a docker container

```
../../docker/run_dicom_to_hpc_predict_io_conversion.sh
```

and to generate the DVC stage use

```
../../docker/dvc/dvc_run_dicom_to_hpc_predict_io_conversion.sh
```

The idea is to maintain the notebook with only the code but no data here and run it on data in the dvc-managed repository with papermill (outputting the resulting notebook there).

The HDF5 results can be investigated with the viewer in 

```
../../docker/run_visualization_notebook.sh
```

To adapt the DICOM interface to a new setting, proceed as follows 

1. Run the full dataset analysis notebook to get an intuition on the headers
2. Develop single-study header visualizations analogous to the one for flow MRI
3. Once a single study is sufficiently understood, implement the conversion
4. Run the resulting notebook on each individual study using papermill and output it to the data repo alongside the conversion result

## DICOM Writer

We use [highdicom](https://github.com/MGHComputationalPathology/highdicom) to write DICOM segmentation masks. A proof-of-concept for how to output a 4D flow MRI with segmentation and abnormality detection masks is available in `segmentation_writer.py`. A Jupyter notebook implementation analogous to the reader has been started in `freiburg_dicom_writer.ipynb`, but remains to be finished with the highdicom output functionality.

## Testing

Currently, there are no tests available. Once the Writer is completed, the functionality can be demonstrated with a round-trip test `DICOM -> HDF5 -> DICOM`.

## Further documentation

The latest presentation of the DICOM interface is available [here](https://github.com/HPC-PREDICT/hpc-predict/blob/master/doc/presentations/HPC-PREDICT%20meeting%2030.04.21%20-%20Slides%20CSCS.pdf). Development issues are discussed [here](https://github.com/HPC-PREDICT/hpc-predict/issues?q=dicom).
