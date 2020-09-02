import os
import pydicom
pydicom.config.datetime_conversion = True

import pandas as pd
from dicom_visitors import visit_dicom_files
from collections import OrderedDict
from types import SimpleNamespace

# From https://github.com/pydicom/pydicom/issues/319#issuecomment-283003834
def dicom_dataset_to_dict(dicom_header):
    dicom_headers = OrderedDict()

    # https://pydicom.github.io/pydicom/stable/old/base_element.html#dataelement
    # DataElement is a simple object which stores the following things:
    # tag – the element’s tag (as a BaseTag object)
    # VR – the element’s Value Representation – a two letter str that describes to the format of the stored value
    # VM – the element’s Value Multiplicity as an int. This is automatically determined from the contents of the value.
    # value – the element’s actual value. A regular value like a number or string (or list of them if the VM > 1), or a Sequence.
    for dicom_elem in dicom_header:
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
    parser.add_argument('--mri-data-root', type=str, default="/home/lukasd/src/hpc-predict/data/v0/input_data/original/mri/MRT Daten Bern/",
                        help='DICOM root directory')
    parser.add_argument('--mri-samples', type=int, nargs='+', default=[3, 5],
                        help='Sample directories to process')
    parser.add_argument('--output-root', type=str, default="/home/lukasd/src/hpc-predict/data/v0/input_data/preprocessed/mri/MRT Daten Bern DICOM Header/",
                        help='Output root directory for serialized DataFrames')
    return parser.parse_args()

### FIXME: replace all occurrences of SimpleNamespace

if __name__ == '__main__':
    args = parse_args()

    for sample in args.mri_samples:
        assert os.path.exists(os.path.join(args.mri_data_root, str(sample) + '/'))

    assert os.path.exists(args.output_root)

    for sample in args.mri_samples:
        header_sample_table = []

        mri_data_path = os.path.join(args.mri_data_root, str(sample) + '/')
        dataframe_output_path = os.path.join(args.output_root, str(sample) + '.pkl')

        if os.path.exists(dataframe_output_path):
            print("File already exists at {} - skipping sample {}.".format(
                dataframe_output_path, sample))
            continue
        print("Converting sample {}...".format(
            sample))

        def header_table_from_dicom_dataset(path, dataset):
            dicom_headers = dicom_dataset_to_dict(dataset)
            dicom_headers["FileCollectionID"] = SimpleNamespace(tag="FileCollectionID",
                                                        VR='CS',
                                                        VM=1,
                                                        name="FileCollectionID",
                                                        keyword="FileCollectionID",
                                                        private_creator=None,
                                                        value=sample)
            dicom_headers["FilePath"] = SimpleNamespace(tag="FilePath",
                                                        VR='CS',
                                                        VM=1,
                                                        name="FilePath",
                                                        keyword="FilePath",
                                                        private_creator=None,
                                                        value=path)

            header_sample_table.append(dicom_headers)

        visit_dicom_files(mri_data_path, header_table_from_dicom_dataset)
        header_sample_df = pd.DataFrame(header_sample_table)
        header_sample_df.to_pickle(dataframe_output_path, compression=None)

        header_sample_table.clear()
