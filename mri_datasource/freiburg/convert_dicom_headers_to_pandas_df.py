import os
import pydicom
pydicom.config.datetime_conversion = True

import itertools
import pandas as pd
from dicom_visitors import visit_dicom_files, visit_dicom_tar
from collections import OrderedDict
from types import SimpleNamespace
import json

# From https://github.com/pydicom/pydicom/issues/319#issuecomment-283003834
def dicom_dataset_to_dict(dicom_header):
    dicom_headers = OrderedDict()

    # https://pydicom.github.io/pydicom/stable/old/base_element.html#dataelement
    # DataElement is a simple object which stores the following things:
    # tag – the element’s tag (as a BaseTag object)
    # VR – the element’s Value Representation – a two letter str that describes to the format of the stored value
    # VM – the element’s Value Multiplicity as an int. This is automatically determined from the contents of the value.
    # value – the element’s actual value. A regular value like a number or string (or list of them if the VM > 1), or a Sequence.
    for dicom_elem in itertools.chain(dicom_header.file_meta, dicom_header) if hasattr(dicom_header, 'file_meta') else dicom_header:
        if dicom_elem.tag == (0x7fe0, 0x0010):
            # discard pixel data
            continue
        if dicom_elem.VR == 'SQ':
            dicom_headers[dicom_elem.tag] = data_element_to_table_entry(dicom_elem)
        else:
            # Could handle other value representations explicitly in the following
            # http://dicom.nema.org/medical/dicom/current/output/chtml/part05/sect_6.2.html#table_6.2-1
            # Cf. conversion of bytes to python objects during pydicom parsing in
            # https://github.com/pydicom/pydicom/blob/22cfa15dae972d82fb38fcb2925b3daf92eaf1d4/pydicom/values.py#L695
            if type(dicom_elem.value) == pydicom.dataset.Dataset: # when is this path actually taken?
                raise RuntimeError("Processing DICOM element of type Dataset: {}".format(dicom_elem.value))
                dicom_headers[dicom_elem.tag] = dicom_dataset_to_dict(dicom_elem.value)
            else:
                #print(" "*4 + "{} ({})".format(dicom_elem, type(dicom_elem.value)))
                dicom_headers[dicom_elem.tag] = dicom_elem
    return dicom_headers


def data_element_to_table_entry(data_elem):
    # Could also use vars(data_elem) to fill in the SimpleNamespace
    data_elem_value = data_elem.value if data_elem.VR != 'SQ' else \
                      [dicom_dataset_to_dict(item) for item in data_elem]
    return SimpleNamespace(tag=data_elem.tag,
                           VR=data_elem.VR,
                           VM=data_elem.VM,
                           name=data_elem.name,
                           keyword=data_elem.keyword,
                           private_creator=data_elem.private_creator,
                           value=data_elem_value)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate pandas dataframe from DICOM headers for exploratory data analysis.')
    parser.add_argument('--hpc-predict-data-root', type=str, default="/home/lukasd/src/hpc-predict/data/v0/",
                        help='(possibly encrypted) HPC-PREDICT data root directory, e.g. data/v1 or data/v0')
    parser.add_argument('--mri-data-root', type=str, default="input_data/original/mri/MRT Daten Bern/",
                        help='DICOM root directory (relative to HPC-PREDICT data root directory)')
    parser.add_argument('--mri-samples', type=int, nargs='+', #default=[3],
                        help='Sample directories to process')
    parser.add_argument('--tar', dest='tar', action='store_true')
    parser.set_defaults(tar=False)
    parser.add_argument('--output-root', type=str, default="input_data/preprocessed/mri/MRT Daten Bern DICOM Header/",
                        help='Output root directory for serialized DataFrames (relative to HPC-PREDICT data root directory)')
    return parser.parse_args()

### FIXME: replace all occurrences of SimpleNamespace

if __name__ == '__main__':
    args = parse_args()

    assert os.path.exists(args.hpc_predict_data_root)
    if os.path.exists(os.path.join(args.hpc_predict_data_root, 'encrypt')):
        abs_hpc_predict_data_root = os.path.join(args.hpc_predict_data_root, 'decrypt')
    else:
        abs_hpc_predict_data_root = args.hpc_predict_data_root
        
    assert args.mri_data_root.startswith('input_data')
    abs_mri_data_root = os.path.join(abs_hpc_predict_data_root, args.mri_data_root)
    assert os.path.exists(abs_mri_data_root)

    # TODO: More elegant would be to 1. analyze all DICOM files jointly to partition them and run the conversion
    #  on partitioned sets, but not yet needed for Freiburg dataset
    for sample in args.mri_samples:
        assert os.path.exists(os.path.join(abs_mri_data_root, str(sample) + ('.tar' if args.tar else '/') ))

    abs_output_root = os.path.join(abs_hpc_predict_data_root, args.output_root)
    if not os.path.exists(abs_output_root):
        print(f"Creating output root directory {abs_output_root}.")
        os.makedirs(abs_output_root, exist_ok=True)

    for sample in args.mri_samples:
        header_sample_table = []

        abs_mri_sample_path = os.path.join(abs_mri_data_root, str(sample) + ('.tar' if args.tar else '/'))
        abs_dataframe_sample_output_path = os.path.join(abs_output_root, str(sample), 'dicom_header_df.pkl')

        if os.path.exists(abs_dataframe_sample_output_path):
            print(f"Skipping   {os.path.relpath(abs_mri_sample_path, start=abs_hpc_predict_data_root)} (output {os.path.relpath(abs_dataframe_sample_output_path, start=abs_hpc_predict_data_root)} already exists).")
        else:
            print(f"Converting {os.path.relpath(abs_mri_sample_path, start=abs_hpc_predict_data_root)} -> "
                  f"{os.path.relpath(abs_dataframe_sample_output_path, start=abs_hpc_predict_data_root)}")
            os.makedirs(os.path.dirname(abs_dataframe_sample_output_path), exist_ok=False)

            def header_table_from_dicom_dataset(path, dataset):
                dicom_headers = dicom_dataset_to_dict(dataset)
                dicom_headers["FileCollectionID"] = SimpleNamespace(tag="FileCollectionID",
                                                            VR='CS',
                                                            VM=1,
                                                            name="FileCollectionID",
                                                            keyword="FileCollectionID",
                                                            private_creator=None,
                                                            value=sample)
                file_path = os.path.relpath(abs_mri_sample_path, start=abs_hpc_predict_data_root) + ':' + path if args.tar else os.path.relpath(path, start=abs_hpc_predict_data_root)
                dicom_headers["FilePath"] = SimpleNamespace(tag="FilePath",
                                                            VR='CS',
                                                            VM=1,
                                                            name="FilePath",
                                                            keyword="FilePath",
                                                            private_creator=None,
                                                            value=file_path)
                for i in range(1,5):
                    split_file_path = file_path.rsplit(':', maxsplit=1)
                    split_file_path[-1] = '/'.join(split_file_path[-1].rsplit('/')[-5:-5+i])
                    dicom_headers[f'FilePath_{i}'] = SimpleNamespace(tag=f'FilePath_{i}',
                                                                     VR='CS',
                                                                     VM=1,
                                                                     name=f'FilePath_{i}',
                                                                     keyword=f'FilePath_{i}',
                                                                     private_creator=None,
                                                                     value=':'.join(split_file_path))

                header_sample_table.append(dicom_headers)

            if not args.tar:
                visit_dicom_files(abs_mri_sample_path, header_table_from_dicom_dataset)
            else:
                visit_dicom_tar(abs_mri_sample_path, header_table_from_dicom_dataset)                
            header_sample_df = pd.DataFrame(header_sample_table)
            header_sample_df.to_pickle(abs_dataframe_sample_output_path, compression=None)

            header_sample_table.clear()
