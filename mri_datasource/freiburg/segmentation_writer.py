from pathlib import Path

import numpy as np
import mr_io
import os
import datetime

from pydicom.sr.codedict import codes
from pydicom.filereader import dcmread
from pydicom.uid import generate_uid
from pydicom import config
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset

from highdicom.content import AlgorithmIdentificationSequence
from highdicom.seg.content import SegmentDescription
from highdicom.seg.enum import (
    SegmentAlgorithmTypeValues,
    SegmentationTypeValues
)
from highdicom.seg.sop import Segmentation

#enforce standard DICOM
config.enforce_valid_values = True
config.future_behavior(False)

folder = '2021-04-20_00-33-21_copper_patient1_R10'
name = 'kspc_R10_vn.mat_segmented.h5'

name_file =  folder + '/output/' + name

path = '/home/nfadel/Desktop/sw/hpc-predict/data/v1/decrypt/segmenter/cnn_segmenter/hpc_predict/v1/inference/' + name_file

#path = '/home/nfadel/Desktop/v4_r20.h5'
mri = mr_io.SegmentedFlowMRI.read_hdf5(path)

# Path to directory containing single-frame legacy CT Image instances
# stored as PS3.10 files
series_dir = Path('/home/nfadel/Desktop/sw/hpc-predict/hpc-predict-io/mri_datasource/freiburg/MRT_Daten_Bern/127/10005583/10005584/10005609')
image_files = series_dir.glob('100096CC')#, '100096CC', '100096E1')#'1000','10005586','100055D8','100055C2')

# Read MRI Image data sets from PS3.10 files
image_datasets = [dcmread(str(f)) for f in image_files] #file dataset


# Create a boolean segmentation mask
#mask = np.zeros(
#    shape=(
#        len(image_datasets),
#        image_datasets[0].Rows,
#        image_datasets[0].Columns
#    ),
#    dtype=np.bool
#)
#mask[1:-1, 10:-10,100:-100] = 1

# Describe the algorithm that created the segmentation
algorithm_identification = AlgorithmIdentificationSequence(
    name='test',
    version='v1.0',
    family=codes.cid7162.ArtificialIntelligence
)

# Describe the segment
description_segment_1 = SegmentDescription(
    segment_number=1,
    segment_label='first segment',
    segmented_property_category=codes.cid7150.Tissue,
    segmented_property_type=codes.cid7166.ConnectiveTissue,
    algorithm_type=SegmentAlgorithmTypeValues.AUTOMATIC,
    algorithm_identification=algorithm_identification,
    tracking_uid=generate_uid(),
    tracking_id='test segmentation of computed tomography image'
)

positions= np.zeros(mri.segmentation_prob.shape[0])


## METADATA
file_meta = FileMetaDataset()
####################
ds = FileDataset('pota.dcm', {}, file_meta=file_meta, preamble=b"\0" * 128)
####################
## DICOM TAGS
ds = image_datasets
ds[0].PatientName = "Test^Firstname"
ds[0].PatientID = "123456"
## Set the transfer syntax
ds[0].is_little_endian = True
ds[0].is_implicit_VR = True
# Set creation date/time
#dt = datetime.datetime.now()
#ds.ContentDate = dt.strftime('%Y%m%d')
#timeStr = dt.strftime('%H%M%S.%f')  # long format with micro seconds
#ds.ContentTime = timeStr

#ds.save_as('output.dcm')


name_output_file = "dicom"
#print('all ', mri.intensity[:,:,:,:].size) 
#print('x ', mri.intensity[:,1,1,1].size)
#print('y ', mri.intensity[1,:,1,1].size)
#print('time ', mri.intensity[1,1,1,:].size)
#print('z ', mri.intensity[1,1,:,1].size)
os.makedirs(name_output_file + "/intensity", exist_ok = True)
series_instance_uid=generate_uid()
for mri_time in range(mri.intensity[1,1,1,:].size):
    for z_axis in range(mri.intensity[1,1,:,1].size): 
        # Create the Segmentation instance ## TODO remove mask
        #image_datasets[z_axis].image_type = 'M'
        ds[0].image_type = 'M',
        ds[0].pixel_data=mri.intensity[:,:,z_axis,mri_time],
        ds[0].sop_instance_uid=generate_uid(),
        ds[0].fix_meta_info(True)
        ds[0].save_as(f"dicom/intensity/TS_{mri_time}_INT_{z_axis}.dcm",write_like_original=False)

#print('all ', mri.velocity_mean[:,:,:,:].size)
#print('time ', mri.velocity_mean[1,1,1,:].size)
#print('x ', mri.velocity_mean[:,1,1,1].size)
#print('y ', mri.velocity_mean[1,:,1,1].size)
#print('z ', mri.velocity_mean[1,1,:,1].size)
os.makedirs(name_output_file + "/velocity_mean", exist_ok = True)
series_instance_uid=generate_uid()
for mri_time in range(int(mri.velocity_mean[1,1,1,:].size/3)):
    for z_axis in range(int(mri.velocity_mean[1,1,:,1].size/3)): 
#        ds.append(ds[0])
        ds[0].image_type = 'P'
        ds[0].pixel_data=mri.velocity_mean[:,:,z_axis,mri_time],
        ds[0].sop_instance_uid=generate_uid(),
        ds[0].fix_meta_info(True)
        ds[0].save_as(f"dicom/velocity_mean/TS_{mri_time}_VELMEAN_{z_axis}.dcm",write_like_original=False)


image_datasets[0].StudyInstanceUID = generate_uid()
os.makedirs(name_output_file + "/segmentation", exist_ok = True) 
series_instance_uid=generate_uid()
for mri_time in range(mri.segmentation_prob[1,1,1,:].size):
    for z_axis in range(mri.segmentation_prob[1,1,:,1].size): 
        # Create the Segmentation instance ##
        seg_dataset = Segmentation(
                source_images=ds,
                pixel_array=mri.segmentation_prob[:,:,z_axis,mri_time],
                segmentation_type=SegmentationTypeValues.FRACTIONAL,
                segment_descriptions=[description_segment_1],
                series_instance_uid=series_instance_uid,
                series_number=2,
                sop_instance_uid=generate_uid(),
                instance_number=z_axis,
                manufacturer='Manufacturer',
                manufacturer_model_name='Model',
                software_versions='v1',
                device_serial_number='Device XYZ',
                #plane_positions=positions#works only if the first dimension of pixel_array = image_datasets
                )
        seg_dataset.fix_meta_info(True)
        seg_dataset.save_as(f"dicom/segmentation/TS_{mri_time}_SEG_{z_axis}.dcm",write_like_original=False)
