from pathlib import Path

import numpy as np
import mr_io

from pydicom.sr.codedict import codes
from pydicom.filereader import dcmread
from pydicom.uid import generate_uid
from pydicom import config

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

mri = mr_io.SegmentedFlowMRI.read_hdf5('/home/nfadel/Desktop/sw/hpc-predict/hpc-predict-io/mri_datasource/recon_volN1_vn.mat.h5')#'/home/nfadel/Desktop/sw/hpc-predict/data/v1/decrypt/segmenter/cnn_segmenter/hpc_predict/v1/inference/2021-01-29_14-47-59_daint102/output/recon_volN1_vn.mat_segmented.h5')

# Path to directory containing single-frame legacy CT Image instances
# stored as PS3.10 files
series_dir = Path('/home/nfadel/Desktop/sw/hpc-predict/hpc-predict-io/mri_datasource/freiburg/MRT_Daten_Bern/127/10005583/10005584/10005585/')
image_files = series_dir.glob('10005605')#,'10005586','100055D8','100055C2')

# Read MRI Image data sets from PS3.10 files
image_datasets = [dcmread(str(f)) for f in image_files]

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
#print('&&&&&&&&&&&&&&&&&&&&')
#print(positions.shape)
positions=mri.segmentation_prob[:,1,:,1]
#print('################')
#print(positions.ndim)
#print(positions.shape)
#print('@@@@@@@@@@@@@@@@@@@@')
#print(mri.segmentation_prob[:,:,:,1].shape[0])#92
#print(len(positions))#96
#print('@@@@@@@@@@@@@@@@@@@@@')

#mri.segmentation_prob[:,:,:,mri_time].shape[0] != len(positions)

series_instance_uid=generate_uid()
#TODO check and in case create folder
for mri_time in range(mri.segmentation_prob[1,1,1].size): 
    # Create the Segmentation instance ##
    seg_dataset = Segmentation(
            source_images=image_datasets,
            pixel_array=mri.segmentation_prob[:,:,:,mri_time],
            segmentation_type=SegmentationTypeValues.FRACTIONAL, #BINARY,
            segment_descriptions=[description_segment_1],
            series_instance_uid=series_instance_uid,
            series_number=2,
            sop_instance_uid=generate_uid(),
            instance_number=mri_time,
            manufacturer='Manufacturer',
            manufacturer_model_name='Model',
            software_versions='v1',
            device_serial_number='Device XYZ',
            #plane_positions=positions#works only if the first dimension of pixel_array = image_datasets
            )
    seg_dataset.fix_meta_info(True)
    seg_dataset.save_as(f"dicom/SEG{mri_time}",write_like_original=False)
#seg_dataset.save_as("seg.dcm")
#TODO create dicomdir with dcmgpdir
