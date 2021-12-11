import os
import pydicom
pydicom.config.datetime_conversion = True

import pandas as pd
from collections import namedtuple, OrderedDict

pydicom_tuples = dict()

treated_as_builtin = set()
def get_builtin_from_pydicom(val):
    from pydicom.valuerep import DA, DT, TM, DSfloat, DSdecimal, IS, PersonName
    from pydicom.uid import UID
    from pydicom.multival import MultiValue
    from datetime import (date, datetime, time, timedelta, timezone)
    from decimal import Decimal

    if isinstance(val, DA):
        return date(val.year, val.month, val.day)
    elif isinstance(val, DT):
        return datetime(val.year, val.month, val.day,
                        val.hour, val.minute, val.second,
                        val.microsecond, val.tzinfo)
    elif isinstance(val, TM):
        return time(val.hour, val.minute, val.second,
                    val.microsecond)
    elif isinstance(val, DSfloat):
        return float(val)
    elif isinstance(val, DSdecimal):
        return Decimal(val)
    elif isinstance(val, IS):
        return int(val)
    elif isinstance(val, MultiValue):
        return tuple([get_builtin_from_pydicom(el) for el in val])
    elif isinstance(val, UID):
        return str(val)
    elif isinstance(val, (PersonName, str, int, float)):
        return val
    elif isinstance(val, list):
        return tuple(val)
    else:
        treated_as_builtin.add(type(val))
        return val


def get_vars_pydicom(value):
    from pydicom.valuerep import DSfloat, DA, DT, TM
    from pydicom.uid import UID
    from datetime import date, time

    if isinstance(value, DSfloat):
        return value.__getstate__()
    elif isinstance(value, UID):
        return {'uid': str(value)}
    elif isinstance(value, DA) or isinstance(value, DT) or isinstance(value, TM):
        return {'original_string': value.original_string}
    else:
        return vars(value)


def hashify_pydicom(value):
    if hasattr(type(value), '__module__') and type(value).__module__.startswith('pydicom'):
        # if not type(value) in pydicom_tuples:
        #     pydicom_tuples[type(value)] = namedtuple(type(value).__name__,
        #                                              [k.lstrip("_") for k in get_vars_pydicom(value).keys()])
        # return pydicom_tuples[type(value)](
        #     **{k.lstrip("_"): hashify_pydicom(v) if not type(v).__module__.endswith('valuerep') else v for
        #        (k, v) in get_vars_pydicom(value).items()})
        return get_builtin_from_pydicom(value)
    elif isinstance(value, list):
        return tuple([hashify_pydicom(v) for v in value])
    elif isinstance(value, OrderedDict):
        return tuple([hashify_pydicom(v.value) for v in value.values()])
    elif isinstance(value, dict):
        return tuple(value.items())
    else:
        return value


def agg_unique(group):
    labels = []
    row = []
    for (name, col) in group[group.columns[:]].iteritems():
        labels.append(name)
        row.append(col.dropna().unique())

    return pd.DataFrame([row], columns=labels)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate pandas dataframe from DICOM headers for exploratory data analysis.')
    parser.add_argument('--headers-root', type=str, default="/home/lukasd/src/hpc-predict/data/v0/input_data/preprocessed/mri/MRT Daten Bern DICOM Header/",
                        help='DICOM headers directory')
    parser.add_argument('--mri-samples', type=str, nargs='+', default=["3.pkl"],
                        help='Sample directories to process')
    parser.add_argument('--output-root', type=str, default="/home/lukasd/src/hpc-predict/data/v0/input_data/preprocessed/mri/MRT Daten Bern DICOM Header Unique/",
                        help='Output root directory for serialized DataFrames with unique tags')
    return parser.parse_args()

### FIXME: replace all occurrences of SimpleNamespace by corresponding namedtuples to make the objects hashable 


if __name__ == '__main__':
    args = parse_args()

    for sample in args.mri_samples:
        assert os.path.exists(os.path.join(args.headers_root, sample))

    assert os.path.exists(args.output_root)

    for sample in args.mri_samples:

        input_path = os.path.join(args.headers_root, sample)
        output_path = os.path.join(args.output_root, sample)

        if os.path.exists(output_path):
            print("File already exists at {} - skipping sample {}.".format(
                output_path, sample))
            continue

        df = pd.read_pickle(input_path)

        tag_names = df.apply(lambda c: c.dropna().apply(lambda x: x.name).unique(), axis=0)
        assert tag_names.shape[0] == 1
        tag_names = tag_names.loc[0]
        df_renamed = df.rename(columns=tag_names.to_dict())


        def remove_ordered_dict(e):
            if isinstance(e, OrderedDict):
                return tuple([e_d.value for e_d in e.values()])
            else:
                return e

        def extract_value(x):
            if hasattr(x,'value'):
                if isinstance(x.value,list):
                    return tuple([ remove_ordered_dict(e) for e in x.value ])
                else:
                    return x.value
            else:
                return x



        # df_renamed_hashified = df_renamed.drop(['FileCollectionID', 'FilePath'], axis=1).applymap(
        #     lambda x: extract_value(x))
        df_renamed_hashified = df_renamed.applymap(
            lambda x: hashify_pydicom(x.value) if pd.notnull(x) else x)
        print("Treated as built-ins: {}".format(treated_as_builtin))
        # df_grouped_unique = df_renamed_hashified.groupby(
        #     ['FileCollectionID', 'Instance Creation Date',
        #      'Patient ID', 'Study Instance UID', 'Series Instance UID',
        #      'Image Type', 'Sequence Name']).apply(lambda x: agg_unique(x))
        group_by_keys = ['FileCollectionID', 'Instance Creation Date', 'Patient ID', 'Study Instance UID',
                         'Sequence Name', 'Series Instance UID', 'Image Type']
        df_grouped_unique = df_renamed_hashified.groupby(group_by_keys).apply(lambda x: agg_unique(x))
        df_grouped_unique.drop(labels=['FilePath'], axis=1).to_pickle(output_path, compression=None)

        # df_grouped_unique_singleton = \
        #     df_grouped_unique.drop(df_grouped_unique.columns[df_grouped_unique.applymap(lambda x: len(x) > 1).any()],
        #                        axis=1).applymap(lambda x: x[0] if len(x) == 1 else x)
        # df_grouped_unique_singleton.to_pickle(output_path, compression=None)

        print("Done writing " + os.path.join(args.output_root, sample), flush=True)




