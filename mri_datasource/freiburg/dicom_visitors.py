import os
import tarfile
import pydicom

# TODO: Pandas reader that converts DICOM file set into DataFrame table (one column per header)

# Example usage can be found in the explore_freiburg_data Jupyter notebook
# Callbacks are meant to modify some external state based on the images they're invoked on according to the visitor pattern

# DICOM file visitor
# image_callback signature:
#   (image_filename: str,
#   dcm_image: pydicom.dataset.FileDataset (DICOM image read with dicom.dcm_read))
#    -> None

def visit_dicom_files(directory, image_callback, filter_images=lambda f: True):
    local_image_filenames = []
    for f in os.scandir(directory):
        if f.is_dir():
            visit_dicom_files(f.path, image_callback, filter_images)
        else:
            local_image_filenames.append(f.path)
    for f_path in local_image_filenames:
        if os.path.basename(f_path) != "DICOMDIR" and filter_images(f_path):
            image_callback(f_path, pydicom.dcmread(f_path))


def visit_dicom_tar(tar_file, image_callback, filter_images=lambda f: True):
    
    tar = tarfile.open(tar_file)    

    d = {f.name: list() for f in tar.getmembers() if f.isdir()}
    for f in tar.getmembers():
        if f.isfile():
            d[os.path.dirname(f.name)].append(f.name)
            
    for d, filenames in d.items(): # possibly only iterate over a selection of directories
        for filename in filenames:
            file = tar.extractfile(filename)
            if os.path.basename(filename) != "DICOMDIR" and filter_images(filename):
                image_callback(filename, pydicom.dcmread(file))
                                    #pydicom.filereader.dcmread(file)   

def visit_dicom_file_collection(mri_data_root, image_callback, max_dirs=2):
    mri_data_collection = sorted([s.path for s in os.scandir(mri_data_root) if s.is_dir()], key=lambda x: int(os.path.basename(x)))

    for mri_data_sample in mri_data_collection if max_dirs is None else mri_data_collection[:max_dirs]:
        mri_data_path = os.path.join(mri_data_root, mri_data_sample)

        visit_dicom_files(mri_data_path, image_callback)


# DICOMDIR visitor
# image_callback signature:
#   (image_filename: str,
#   dcm_image: pydicom.dataset.FileDataset (DICOM image read with dicom.dcm_read),
#   dcmdir_entry: dict with entries for 'patient', 'study', 'series' with values corresponding to the records of the DICOM image)
#    -> None

def visit_all_series(base_dir, patient_record, study, image_callback):
    print(" " * 4 + "Study {}: {}".format(study.StudyID,
                                          study.StudyDate))
    all_series = study.children
    # go through each serie
    for series in all_series:
        image_count = len(series.children)
        plural = ('', 's')[image_count > 1]

        # Write basic series info and image count

        # Put N/A in if no Series Description
        if 'SeriesDescription' not in series:
            series.SeriesDescription = "N/A"
        print(" " * 8 + "Series {}: {}: {} ({} image{})".format(
            series.SeriesNumber, series.Modality, series.SeriesDescription,
            image_count, plural))

        # Open and read something from each image, for demonstration
        # purposes. For simple quick overview of DICOMDIR, leave the
        # following out
        print(" " * 12 + "Reading images...")
        image_records = series.children
        for image in [(image_filename_rec[0], pydicom.dcmread(image_filename_rec[0]), image_filename_rec[1])
                      for image_filename_rec in [(os.path.join(base_dir, *image_rec.ReferencedFileID), image_rec)
                                                 for image_rec in image_records]]:
            image_callback(*image[:2], dcmdir_entry={'patient': patient_record, 'study': study, 'series': series,
                                                     'image': image[2]})


def visit_dicom_dir(base_dir, dicom_dir, image_callback):
    # go through the patient record and print information
    for patient_record in dicom_dir.patient_records:
        if (hasattr(patient_record, 'PatientID') and
                hasattr(patient_record, 'PatientName')):
            print("Patient: {}: {}".format(patient_record.PatientID,
                                           patient_record.PatientName))
        studies = patient_record.children
        # go through each series
        for study in studies:
            visit_all_series(base_dir, patient_record, study, image_callback)


def visit_dicom_dir_collection(mri_data_root, image_callback, max_dirs=2):
    mri_data_collection = sorted([s.path for s in os.scandir(mri_data_root) if s.is_dir()],
                                 key=lambda x: int(os.path.basename(x)))

    for mri_data_sample in mri_data_collection if max_dirs is None else mri_data_collection[:max_dirs]:
        # fetch the path to the test data
        mri_data_path = os.path.join(mri_data_root, mri_data_sample)
        filepath = os.path.join(mri_data_path, 'DICOMDIR')
        print('Visiting {}...'.format(mri_data_sample))
        # load the data
        dicom_dir = pydicom.filereader.read_dicomdir(filepath)
        base_dir = os.path.dirname(filepath)

        visit_dicom_dir(base_dir, dicom_dir, image_callback)
