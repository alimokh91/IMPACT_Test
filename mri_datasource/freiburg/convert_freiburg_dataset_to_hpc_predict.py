import argparse
import logging
import pydicom
import numpy
import fleep
import zipfile
import time
import os
import h5py
import matplotlib.pyplot as plt
import imageio
import mr_io # Requires adding ../python to PYTHONPATH
from pydicom.filereader import read_dicomdir

# Parse data input and output directories
def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Generate hpc-predict-io HDF5-message from Freiburg dataset sample (DICOM).')
    parser.add_argument('--input', type=str,
                        help='Directory containing the dataset')
    parser.add_argument('--input_name', type=str, default=None,
                        help='Name of the dataset')
    parser.add_argument('--password', type=str, default=None,
                        help='Password dataset')
    parser.add_argument('--output', type=str,
                        help='Output directory for HDF5 files')
    parser.add_argument('--log', type=str, default="warn", help="Logging level")
    return parser.parse_args()

args = parse_args()
logging.basicConfig(level=args.log.upper())

def read_dicom(path, filename=None, remove_temp_files=False, password=None):
    ticks = time.clock()
    print("Initial Time: ", ticks)
    listFilesDicom = []
    sequence_name_list = []
    pixel_array_dtype = []
    rows = []
    columns = []

    if filename is not None:
        full_path = path + filename
        temp_path = './temp_dicom_' + filename
        if is_zip(full_path):
            extract_zip(full_path, filename, 'r', password)
    else:
        temp_path = path

    if os.path.isdir(temp_path):
        # browser folder (or zipped folder) and subdirectories
        for dirName, subdirList, fileList in os.walk(temp_path):  # browser zip file and subdirectories
            for path in fileList:
                l_file = os.path.join(dirName, path)
                #check if it is DICOM file
                if is_dicom(l_file):

                    RefDs = pydicom.read_file(l_file, force=True)
                    #select the DICOM data that has only 3D MRI with magnitude or phase
                    if RefDs.get("MRAcquisitionType", None) == '3D' and (
                            'M' in RefDs.get("ImageType", None) or 'P' in RefDs.get("ImageType", None)) and RefDs.get(
                            "Modality", None) == 'MR':
                        listFilesDicom.append(l_file)
                        sequence_name_list.append(RefDs.SequenceName)
                        rows.append(RefDs.get("Rows", None))
                        columns.append(RefDs.get("Columns", None))
                        pixel_array_dtype.append(RefDs.pixel_array.dtype)

        #evaluate which MRI dataset we need
        position_to_pick, seq_name = check_sequence(sequence_name_list)
        row_tp = int(rows[position_to_pick])
        col_tp = int(columns[position_to_pick])
        px_array_tp = pixel_array_dtype[position_to_pick]
        del sequence_name_list, rows, columns, pixel_array_dtype

        imgMagnitude = create_image(listFilesDicom, seq_name, px_array_tp, row_tp, col_tp, 'M')
        #print ('m ' ,imgMagnitude[0].size)
        
        imgPhase, unique_slice_loc, unique_acq_times, unique_rows, unique_columns = create_image(listFilesDicom, seq_name, px_array_tp, row_tp, col_tp, 'P')
        #print('p ', imgPhase.size)
        del listFilesDicom

        #if filename is not None:
        #    h5f = h5py.File(filename + '.h5', 'w')
        #else:
        #    h5f = h5py.File('dicom_image.h5', 'w')
        #grp = h5f.create_group('DICOM_data')
        #print('dim ', imgMagnitude[0].ndim)
        #print('len ', len(imgMagnitude))
        #print('size ', imgMagnitude[0].size)
        #print('len2 ', imgMagnitude[0].len())
        #print('shape ', imgMagnitude[0].shape)
        
        #grp.create_dataset('Magnitude', imgMagnitude[0].shape, data=imgMagnitude[0], dtype=imgMagnitude[0].dtype)
        #dataset = grp.create_dataset('Phase', imgPhase.shape, data=imgPhase, dtype=imgPhase.dtype)
        
        #print("Dataset dataspace is", dataset.shape)
        #print("Dataset Numpy datatype is", dataset.dtype)
        #print("Dataset name is", dataset.name)
        #print("Dataset is a member of the group", dataset.parent)
        #print("Dataset was created in the file", dataset.file)
        #
        #h5f.close()
        #mat_data
        #i*spatial_voxel_width[j] for i in range(mat_data['results_v'].shape[j])
        #geometry = [ np.array([ i*spatial_voxel_width[j] for i in range(mat_data['results_v'].shape[j]) ]) for j in range(3) ]
        #geometry = [ numpy.array([ 1,2,3])] # for j in range(3) ]  #TODO
        geometry = numpy.stack((unique_rows, unique_columns, unique_slice_loc), axis=-1)
 
        hpc_predict_mri = mr_io.FlowMRI(geometry,  # ...,
                                        time=unique_acq_times,  # ...,
                                        time_heart_cycle_period=60/70,
                                        intensity=imgMagnitude[0],  # ...,
                                        velocity_mean=imgPhase,  # ...,
                                        velocity_cov=numpy.zeros(shape=imgPhase.shape+(3,)))  # ...)
        hpc_predict_mri.write_hdf5(path + '/freiburg_dataset_mri.h5')


        image = build_4Dimage(imgMagnitude[0], imgPhase)
        del imgMagnitude, imgPhase

        if remove_temp_files:
            shutil.rmtree(os.getcwd() + '/temp_dicom_' + filename)

    else:
        if is_dicom(full_path):
            RefDs = pydicom.read_file(full_path, force=True)  # assuming that one zip file doesn't have only one file
            # TODO extract image from a single file
            #     should be like:
            image = RefDs.pixel_array
        else:
            raise ValueError("Error: not a DICOM file!")

    ticks2 = time.clock()
    print("Final Time: ", ticks2)
    output = ticks2 - ticks
    print("Total Time: ", output)
    return image

#... # read dicom file from input directory
#... # and create the arrays required below to create an FlowMRI object
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
                array_seq[i] = array_seq[i] +1
    sequence_name = unique_sequence_name[numpy.argmax(array_seq)]
    for i in range(num_of_sequences):
        if sequence_name == sequence_name_list[i]:
            position_to_pick = i
            break
        else:
            position_to_pick = None
    return position_to_pick,sequence_name

def create_image(file_list, seq_name, array_type, row, column, Type):
    l_file_list = []
    l_slice_locations = []
    l_acquisition_times = []
    l_rows = []
    l_columns = []

    for i in range(len(file_list)):
        RefDs = pydicom.read_file(file_list[i], force=True)

        if RefDs.ImageType[2] == Type and RefDs.SequenceName == seq_name:
            l_file_list.append(file_list[i])
            l_slice_locations.append(RefDs.get("SliceLocation", None))
            l_acquisition_times.append(RefDs.get("AcquisitionTime", None))
            l_rows.append(RefDs.get("Rows", None))
            l_columns.append(RefDs.get("Columns", None))

    l_slice_locations = numpy.array(l_slice_locations)
    l_acquisition_times = numpy.array(l_acquisition_times)
    l_rows = numpy.array( l_rows )
    l_columns = numpy.array(l_columns)
    
    print("l_rows[0] ", l_rows[0])
    unique_rows = l_rows.resize(l_rows[0])
    unique_columns = l_columns.resize(l_columns[0])

    #print("l_acquisition_times", l_acquisition_times, "lunghezza", len(l_acquisition_times) )
    #print("l_slice_locations", l_slice_locations.shape, "lunghezza", len(l_slice_locations) )
    ### sort and unique arrays ###
    print("unique_rows", unique_rows.shape)

    unique_slice_loc = numpy.sort(numpy.unique(l_slice_locations))
    unique_acq_times = numpy.sort(numpy.unique(l_acquisition_times))
    #print("unique_acquisition_times", unique_acq_times, "lunghezza", len(unique_acq_times) )
    #print("l_slice_locations", unique_slice_loc.shape, "lunghezza", len(unique_slice_loc) )
    del l_slice_locations, l_acquisition_times

    print("List of ", Type, " Dicom Files", len(l_file_list))
    PixelDims = (row, column, unique_slice_loc.shape[0] , unique_acq_times.shape[0])
    print(Type, " Pixel Dimension", PixelDims)
    img = numpy.zeros(PixelDims, dtype = array_type)

    for i in range(len(l_file_list)):
        RefDs = pydicom.read_file(l_file_list[i], force=True)
        slice_idx = numpy.where(RefDs.get("SliceLocation", None) == unique_slice_loc)[0][0]
        time_idx = numpy.where(RefDs.get("AcquisitionTime", None) == unique_acq_times)[0][0]
        img[:, :,slice_idx, time_idx]= RefDs.pixel_array
        #PixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))
        #print('Pixel spacing in mm: ', PixelSpacing)
        #print('Pixel volume in mm^3:', PixelSpacing[0] * PixelSpacing[1] * PixelSpacing[2])
    
    return img, unique_slice_loc, unique_acq_times, unique_rows, unique_columns

def build_4Dimage(image_mag,image_pha):
    image_phx = image_pha[:,:,:,0::3]
    image_phy = image_pha[:,:,:,1::3]
    image_phz = image_pha[:,:,:,2::3]
    image = numpy.stack((image_mag, image_phx, image_phy, image_phz), axis=-1)
    return image

# ================================
# This function takes a 4d image array [imgx, imgy, imgz, time] and saves a gif of all xy slices over all times.
# ================================
def save_gif(im, save_path):
    num_times = im.shape[3]
    num_slices = im.shape[2]
    print('num_times',num_times)
    print('num_slices',num_slices)
    # all slices of each time-index as png
    for t in range(num_times):
        plt.figure(figsize=[30,30])
        for j in range(num_slices):
            plt.subplot(8, 8, j+1)
            plt.imshow(im[:,:,j,t], cmap='gray')
            plt.xticks([], []); plt.yticks([], [])
            plt.title('slice: ' + str(j) + ', time:' + str(t//2))
        plt.savefig(save_path + '_time' + str(t) + '.png')
        plt.close()
    # all time indices in gif
    gif = []
    for t in range(num_times):
        gif.append(imageio.imread(save_path + '_time' + str(t) + '.png'))
    imageio.mimsave(save_path + 'image.gif', gif, format='GIF', duration=0.5)

def write_dicom(path='./', filename='DefautlDICOMfile', write_compress=False,password=None):
    ticks = time.clock()
    print("Initial Time: ", ticks)

    fullpath = path + filename + '.dcm'
    #for i in range(len(RefDs)):
    #    RefDs[i] = pydicom.write_file(filename), RefDs , force=True)
        #RefDs[i] = pydicom.save_as(fullpath,False)#, force=True)
#        RefDs = pydicom.write_file(fullpath)# , force=True)

    if len(RefDs) is 1: #TODO to check
        #####################################################
        ## single DICOM case
        fullpath = path + filename + '.dcm'
        pydicom.write_file(fullpath, RefDs)

        if write_compress:
            write_zip(fullpath,filename,password)
            os.remove(fullpath)
        #####################################################
    else:
        #####################################################
        ## folder DICOM case
        fullpath = path + filename
        for i in range(len(RefDs)):
            iterativeName= fullpath + "_" + str(i)
            pydicom.write_file(iterativeName, RefDs[i])
        #####################################################

    ticks2 = time.clock()
    print("Final Time: ", ticks2)
    output = ticks2 - ticks
    print("Total Time: ", output)

def is_zip(path):
    with open(path, "rb") as file:
        info = fleep.get(file.read(128))
    if (info.mime_matches("application/zip")):
        print("It is a zip file")
        return True
    else:
        return False

def is_dicom(path):
    with open(path, "rb") as file:
        info = fleep.get(file.read(128))
        if (info.mime_matches("application/dicom")):
            print("It is a DICOM file")
            return True
        else:
            n = file.read(132)
            if n.find(bytes('DICM', 'utf-8')) is -1 :
                raise ValueError("The file {} is NOT a DICOM file!!".format(path))
                return False
            else:
               # print("It is a DICOM file")
                return True

def extract_zip(path, filename='filename', fileMode='r', password=None):
    with zipfile.ZipFile(path, fileMode) as zip:
        # printing all the contents of the zip file
        #zip.printdir()
        print('Extracting all the files now...')
        for zinfo in zip.infolist():
            is_encrypted = zinfo.flag_bits & 0x1
        if is_encrypted:
            print('Extracting the encrypted files now...')
            #zip.extractall(pwd=bytes(password,"utf-8"))
            #os.system("7z x "+ path + " -p" + password)
            #TODO modify in order to run in a encrypted environment, for example enforcing local filesystem
            os.system('unzip -qd ./temp_dicom_' + filename + ' -P ' + password + ' ' + path)
        else:
            # extracting all the files
            #zip.extractall()
            os.system('unzip -qd ./temp_dicom_' + filename + ' ' + path)
        print('Done!')

def get_all_file_paths(directory):
    # initializing empty file paths list
    file_paths = []
    # crawling through directory and subdirectories
    for root, directories, files in os.walk(directory):
        for filename in files:
            # join the two strings in order to form the full filepath
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)
    # returning all file paths
    return file_paths

def write_zip(pathToZip, outputName='defaultiDicomName', password=None):
    ticks = time.clock()
    print("Initial Time: ", ticks)
    outputName = outputName + '.zip'
    if password is None:
        if os.path.isdir(pathToZip):
            allFilePaths = get_all_file_paths(pathToZip)
            # printing the list of all files to be zipped
            #print('Following files will be zipped:')
            #for file_name in allFilePaths:
            #     print(file_name)
            #writing files to a zipfile
            with zipfile.ZipFile(outputName,'w') as zip:
                # writing each file one by one
                for file in allFilePaths:
                    zip.write(file)
        if os.path.isfile(pathToZip):
            with zipfile.ZipFile(outputName,'w') as zip:
                zip.write(pathToZip)
    else:
        os.system('zip -rq ' + outputName + ' -P ' + password + ' ' + pathToZip) #-Zb and -9

    print('All files zipped successfully!')
    ticks2 = time.clock()
    print("Final Time: ", ticks2)
    output = ticks2 - ticks
    print("Total Time: ", output)


#dicom_dir = read_dicomdir(args.input + "/DICOMDIR")
read_dicom(args.input)
#read_dicom("./3")
# hpc_predict_mri = mr_io.FlowMRI(geometry=None, #...,
#                                 time=None, #...,
#                                 intensity=None, #...,
#                                 velocity_mean=None, #...,
#                                 velocity_cov=None) #...)
# hpc_predict_mri.write_hdf5(args.output + '/freiburg_dataset_mri.h5')
