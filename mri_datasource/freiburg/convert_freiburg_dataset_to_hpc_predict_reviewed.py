import os
import numpy
import pydicom

from dicom_visitors import visit_dicom_files

# Code review for Nur's pull request (#27)
# TODO: - DICOMDIR seem to incompletely reference DICOM datasets, therefore singleing out flow MRIs manually
#       seems necessary
#       - Validate groupings of MRIs: Headers being used to single out an instance of a flow MRI (currently only
#       MRAcquisitionType == '3D', ImageType == 'M' or 'P', Modality == 'MR', but not series/study/patient etc. The
#       current implementation seems to assume every numeric folder corresponds to a single flow MRI, but this seems
#       wrong as the data sometimes comes from different patients (and the CardiacNumberOfImages also suggest a
#       different number of slices along the time axis).
#       - geometric validation (MRI coordinate system): orthogonality of Row/Col vectors (ImageOrientationPatient)
#       with inter-slice origin offset (ImagePositionPatient), consistency of these across all images of a single
#       flow MRI instance, i.e. that the spatial coordinate grid of the DICOM images matches the assumptions made in
#       HPC-PREDICT-IO (and is consistent across all images in a single flow MRI)
#       - understand and decide what to do about staggered time grids for x-/y-/z-velocity components (also shifted
#       w.r.t. magnitude values)
#       - many of these tasks that require to understand the data are better suited to a Jupyter notebook rather than a
#       plain python script (for reproducibility of the findings)

def create_image(file_list, seq_name, array_type, row, column, ImageType):
    l_file_list = []
    l_slice_locations = []
    l_acquisition_times = []
    l_rows = []
    l_columns = []

    # Review: - This completes filtering by most frequent SequenceName and does filtering by ImageType[2] (where probably
    #         M or P is expected - is this everywhere the case?).
    #         - All these filtering steps should be done in one single step and rereading the files should be
    #         avoided. Also, consider expressing these filtering operations etc. with a functional programming model,
    #         e.g. with pandas might be an elegant solution.
    for i in range(len(file_list)):
        dcm_image = pydicom.dcmread(file_list[i])

        if dcm_image.ImageType[2] == ImageType and dcm_image.SequenceName == seq_name:
            l_file_list.append(file_list[i])
            l_slice_locations.append(dcm_image.get("SliceLocation", None))
            l_acquisition_times.append(dcm_image.get("AcquisitionTime", None))
            l_rows.append(dcm_image.get("Rows", None))
            l_columns.append(dcm_image.get("Columns", None))

    # Review: - the list of row sizes/column sizes of the images (# images of these) is shrinked down to
    #         the row size/column size to be used for the geometry parameters. The length of the expected 1d geometry
    #         coordinate arrays correspond to the number of rows, but to compute them we need to get the PixelSpacing
    #         as well as the coordinates from ImagePositionPatient between adjacent images (to be validated with
    #         ImageOrientationPatient, SliceThickness as used below is not needed).
    l_rows = numpy.array(l_rows)
    unique_rows = numpy.resize(l_rows, (l_rows[0],))
    l_columns = numpy.array(l_columns)
    unique_columns = numpy.resize(l_columns, (l_columns[0],))

    # Review: - Are the slice locations the same at all time steps? For unique acquisition for magnitude image type times
    #         I get 46 on the '3' folder and 138 for the phase image type, whereas e.g.
    #         '3/10000000/10000001/10000007/10002597' (magnitude) has dcm_img.CardiacNumberOfImages == "23" and so does
    #         '3/10000000/10000001/10000002/10002544' (phase). Is it possible that two flow MRIs are mixed up here?
    unique_slice_loc = numpy.sort(numpy.unique(l_slice_locations))
    l_acquisition_times = numpy.array(l_acquisition_times)
    unique_acq_times = numpy.sort(numpy.unique(l_acquisition_times))

    print("List of ", ImageType, " Dicom Files", len(l_file_list))
    PixelDims = (row, column, unique_slice_loc.shape[0], unique_acq_times.shape[0])
    print(ImageType, " Pixel Dimension", PixelDims)
    img = numpy.zeros(PixelDims, dtype=array_type)

    # Review: - This relies on different AcquisitionTimes for phases in different directions (is this consistent across
    #         the entire dataset?), which is then corrected again in build_4D image. It would be better to not mix up
    #         different components from the start and directly build the final 4d-cube of image data with np.stack.
    #         - Since the component time grids are staggered, the HPC-PREDICT-IO format needs to be modified.
    #         - Also, there is no validation that the spatial coordinates of the different component measurements are the
    #         same - i.e. that ImageOrientationPatient and ImagePositionPatient (and possibly also PixelSpacing) are the
    #         same? Again, testing on the two studies in folder 3 showed one and two, respectively, ImageOrientationPatient
    #         values per study, while the 30 values of ImagePositionPatient per study were consistent up to 1 strong
    #         outlier
    #         - Rereading the files should be avoided
    for i in range(len(l_file_list)):
        dcm_image = pydicom.dcmread(l_file_list[i])
        slice_idx = numpy.where(dcm_image.SliceLocation == unique_slice_loc)[0][0]
        time_idx = numpy.where(dcm_image.AcquisitionTime == unique_acq_times)[0][0]
        img[:, :, slice_idx, time_idx] = dcm_image.pixel_array
        # PixelSpacing = (float(dcm_image.PixelSpacing[0]), float(dcm_image.PixelSpacing[1]), float(dcm_image.SliceThickness))
        # print('Pixel spacing in mm: ', PixelSpacing)
        # print('Pixel volume in mm^3:', PixelSpacing[0] * PixelSpacing[1] * PixelSpacing[2])

    return img, unique_slice_loc, unique_acq_times, unique_rows, unique_columns


def check_sequence(sequence_name_list):
    unique_sequence_name = numpy.unique(sequence_name_list)
    num_of_sequences = len(sequence_name_list)
    unique_sequences = len(unique_sequence_name)
    array_seq = numpy.zeros(unique_sequences)
    position_to_pick = int(0)
    if unique_sequences < 1:
        raise ValueError("Error: no DICOM images with sequence tag setted!")
    for j in range(num_of_sequences):
        for i in range(unique_sequences):
            if sequence_name_list[j] == unique_sequence_name[i]:
                array_seq[i] = array_seq[i] + 1
    sequence_name = unique_sequence_name[numpy.argmax(array_seq)]
    for i in range(num_of_sequences):
        if sequence_name == sequence_name_list[i]:
            position_to_pick = i
            break
        else:
            position_to_pick = None
    return position_to_pick, sequence_name


def build_4Dimage(image_mag, image_pha):
    image_phx = image_pha[:, :, :, 0::3]
    image_phy = image_pha[:, :, :, 1::3]
    image_phz = image_pha[:, :, :, 2::3]
    image = numpy.stack((image_mag, image_phx, image_phy, image_phz), axis=-1)
    return image


def read_dicom(root_path):
    # Review: - Functional programming would allow to reduce the number of data structures (lists, etc.), avoid repeated
    #         reading of files and to express the data transformations in a much cleaner way. Maybe pandas might be
    #         worth giving a try? (it would probably also be useful for data set exploration)
    listFilesDicom = []
    sequence_name_list = []
    pixel_array_dtype = []
    rows = []
    columns = []

    def read_dicom_files_visitor(image_filename, dcm_image):
        # Select the DICOM data that has only 3D MRI with magnitude or phase
        if dcm_image.MRAcquisitionType == '3D' and \
                ('M' in dcm_image.ImageType or 'P' in dcm_image.ImageType) and \
                dcm_image.Modality == 'MR':
            listFilesDicom.append(image_filename)
            sequence_name_list.append(dcm_image.SequenceName)  # why sequence name and not SeriesInstanceUID e.g.?
            rows.append(dcm_image.Rows)
            columns.append(dcm_image.Columns)
            pixel_array_dtype.append(dcm_image.pixel_array.dtype)
        else:
            print("Skipping {} with attributes (MRAcquisitionType) {} (ImageType) {} (Modality) {}".format(
                image_filename, dcm_image.MRAcquisitionType, dcm_image.ImageType, dcm_image.Modality))

    visit_dicom_files(root_path, read_dicom_files_visitor)

    # Review: - Below we're creating one single MRI per folder, but the Jupyter notebook shows that these images
    #         are from different patients. How can we be sure the computed flow MRI is correct then?
    #         - check_sequence finds the most frequent SequenceName in the list of DICOM files and
    #         returns the index of its first occurrence in the above list (position_to_pick)
    #         - The filtering is then completed in create_image (cf. comment there)
    position_to_pick, seq_name = check_sequence(sequence_name_list)

    row_tp = int(rows[position_to_pick])
    col_tp = int(columns[position_to_pick])
    px_array_tp = pixel_array_dtype[position_to_pick]

    # Review: use imgMagnitude, _, _, _, _ to discard return values after the first one
    imgMagnitude = create_image(listFilesDicom, seq_name, px_array_tp, row_tp, col_tp, 'M')

    imgPhase, unique_slice_loc, unique_acq_times, unique_rows, unique_columns = create_image(listFilesDicom, seq_name,
                                                                                             px_array_tp, row_tp,
                                                                                             col_tp, 'P')

    # Review: - unique_rows and unique_columns contain the number of pixels along their direction
    #         in all entries - they should be replaced by something like np.linspace with number
    #         of rows/columns multiplied by the PixelSpacing of the images. unique_slice_loc should be
    #         validated - in a test I've found that axis to be the reverted z-axis compared to that defined by
    #         ImageOrientationPatient (by completing the right-handed system). This means that the data is reflected
    #         along the z-axis, which would anatomically be wrong.
    #         - geometry should just be a list of np.arrays with the coordinates of the voxel centers along the
    #         separate axes.
    #         - generally, I'd prefer having the HDF5-writing code in a separate function from that that computes the
    #         data structures from DICOM files

    # geometry = numpy.stack((unique_rows, unique_columns, unique_slice_loc), axis=-1)
    #    hpc_predict_mri = mr_io.FlowMRI(geometry,  # ...,
    #                                    time=unique_acq_times,  # ...,
    #                                    time_heart_cycle_period=60/70,
    #                                    intensity=imgMagnitude[0],  # ...,
    #                                    velocity_mean=imgPhase,  # ...,
    #                                    velocity_cov=numpy.zeros(shape=imgPhase.shape+(3,)))  # ...)
    #    hpc_predict_mri.write_hdf5(path + '/freiburg_dataset_mri.h5')

    image = build_4Dimage(imgMagnitude[0], imgPhase)

    return image


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate hpc-predict-io HDF5-message from Freiburg dataset sample (DICOM).')
    parser.add_argument('--mri-data-root', type=str, default="/home/lukasd/src/hpc-predict/data/v0/input_data/original/mri/MRT Daten Bern/",
                        help='DICOM directory containing the DICOMDIR')
    parser.add_argument('--mri-samples', type=int, nargs='+', default=[3,5],
                        help='Sample directories to process')
    # parser.add_argument('--output', type=str,
    #                     help='Output directory for HPC-PREDICT-IO files')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    for sample in args.mri_samples:
        mri_data_root = "/home/lukasd/src/hpc-predict/data/v0/input_data/original/mri/MRT Daten Bern/"
        mri_data_path = os.path.join(mri_data_root, str(sample) + '/')
        read_dicom(mri_data_path)
        # ... (read_dicom should probably just read the data and return the HPC-PREDICT-IO data structures)

