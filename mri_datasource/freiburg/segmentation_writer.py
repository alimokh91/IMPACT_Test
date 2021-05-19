from pathlib import Path

import numpy as np
import mr_io
import os
import datetime as d

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

enable_intensity = True
enable_velocity_mean = True
enable_segmentation = True
enable_abnormality = True
enable_debug = False

def debug_print_tags(ds):
    print('AccessionNumber ',ds[0].AccessionNumber)
    print('AcquisitionDate ',ds[0].AcquisitionDate)
    print('AcquisitionMatrix ',ds[0].AcquisitionMatrix)
    print('AcquisitionNumber ',ds[0].AcquisitionNumber)
    print('AcquisitionTime ',ds[0].AcquisitionTime)
    #print('AdditionalPatientHistory ',ds[0].AdditionalPatientHistory)
    #print('AdmittingDiagnosesDescription ',ds[0].AdmittingDiagnosesDescription)
    print('AngioFlag ',ds[0].AngioFlag)
    print('BitsAllocated ',ds[0].BitsAllocated)
    print('BitsStored ',ds[0].BitsStored)
    print('CardiacNumberOfImages ',ds[0].CardiacNumberOfImages)
    print('Columns ',ds[0].Columns)
    print('ContentDate ',ds[0].ContentDate)
    print('ContentTime ',ds[0].ContentTime)
    #print('DerivationCodeSequence ',ds[0].DerivationCodeSequence)
    print('DerivationDescription ',ds[0].DerivationDescription)
    print('DeviceSerialNumber ',ds[0].DeviceSerialNumber)
    print('EchoNumbers ',ds[0].EchoNumbers)
    print('EchoTime ',ds[0].EchoTime)
    print('EchoTrainLength ',ds[0].EchoTrainLength)
    print('FlipAngle ',ds[0].FlipAngle)
    print('FrameOfReferenceUID ',ds[0].FrameOfReferenceUID)
    print('HighBit ',ds[0].HighBit)
    print('ImageComments ',ds[0].ImageComments)
    print('ImageOrientationPatient ',ds[0].ImageOrientationPatient)
    print('ImagePositionPatient ',ds[0].ImagePositionPatient)
    print('ImageType ',ds[0].ImageType)
    print('ImagedNucleus ',ds[0].ImagedNucleus)
    print('ImagingFrequency ',ds[0].ImagingFrequency)
    print('InPlanePhaseEncodingDirection ',ds[0].InPlanePhaseEncodingDirection)
    print('InstanceCreationDate ',ds[0].InstanceCreationDate)
    print('InstanceCreationTime ',ds[0].InstanceCreationTime)
    print('InstanceNumber ',ds[0].InstanceNumber)
    print('InstitutionAddress ',ds[0].InstitutionAddress)
    print('InstitutionName ',ds[0].InstitutionName)
    print('InstitutionalDepartmentName ',ds[0].InstitutionalDepartmentName)
    #print('IssuerOfAccessionNumberSequence ',ds[0].IssuerOfAccessionNumberSequence)
    print('IssuerOfPatientID ',ds[0].IssuerOfPatientID)
    print('LargestImagePixelValue ',ds[0].LargestImagePixelValue)
    print('MRAcquisitionType ',ds[0].MRAcquisitionType)
    print('MagneticFieldStrength ',ds[0].MagneticFieldStrength)
    print('Manufacturer ',ds[0].Manufacturer)
    print('ManufacturerModelName ',ds[0].ManufacturerModelName)
    print('Modality ',ds[0].Modality)
    print('NominalInterval ',ds[0].NominalInterval)
    print('NumberOfAverages ',ds[0].NumberOfAverages)
    print('NumberOfPhaseEncodingSteps ',ds[0].NumberOfPhaseEncodingSteps)
    print('OperatorsName ',ds[0].OperatorsName)
    #print('OriginalAttributesSequence ',ds[0].OriginalAttributesSequence)
    print('#### ')
    print('OtherPatientIDs ',ds[0].OtherPatientIDs)
    #print('OtherPatientNames ',ds[0].OtherPatientNames)
    print('PatientAge ',ds[0].PatientAge)
    print('PatientBirthDate ',ds[0].PatientBirthDate)
    #print('PatientComments ',ds[0].PatientComments)
    print('PatientID ',ds[0].PatientID)
    print('PatientMotherBirthName ',ds[0].PatientMotherBirthName)
    print('PatientName ',ds[0].PatientName)
    print('PatientPosition ',ds[0].PatientPosition)
    print('PatientSex ',ds[0].PatientSex)
    print('PatientSize ',ds[0].PatientSize)
    print('PatientWeight ',ds[0].PatientWeight)
    print('PercentPhaseFieldOfView ',ds[0].PercentPhaseFieldOfView)
    print('PercentSampling ',ds[0].PercentSampling)
    print('PerformedProcedureStepDescription ',ds[0].PerformedProcedureStepDescription)
    print('PerformedProcedureStepID ',ds[0].PerformedProcedureStepID)
    print('PerformedProcedureStepStartDate ',ds[0].PerformedProcedureStepStartDate)
    print('PerformedProcedureStepStartTime ',ds[0].PerformedProcedureStepStartTime)
    print('PerformingPhysicianName ',ds[0].PerformingPhysicianName)
    print('PhotometricInterpretation ',ds[0].PhotometricInterpretation)
    print('PixelBandwidth ',ds[0].PixelBandwidth)
    print('PixelRepresentationpp ',ds[0].PixelRepresentation)
    print('PixelSpacing ',ds[0].PixelSpacing)
    #print('PositionReferenceIndicator ',ds[0].PositionReferenceIndicator)
    print('PregnancyStatus ',ds[0].PregnancyStatus)
    #print('ProcedureCodeSequence ',ds[0].ProcedureCodeSequence)
    print('ProtocolName ',ds[0].ProtocolName)
    #print('ReferencedImageSequence ',ds[0].ReferencedImageSequence)
    #print('ReferencedPerformedProcedureStepSequence ',ds[0].ReferencedPerformedProcedureStepSequence)
    print('ReferringPhysicianName ',ds[0].ReferringPhysicianName)
    print('RepetitionTime ',ds[0].RepetitionTime)
    #print('RequestAttributesSequence ',ds[0].RequestAttributesSequence)
    #print('RequestedProcedureCodeSequence ',ds[0].RequestedProcedureCodeSequence)
    print('RequestedProcedureDescription ',ds[0].RequestedProcedureDescription)
    print('RescaleIntercept ',ds[0].RescaleIntercept)
    print('RescaleSlope ',ds[0].RescaleSlope)
    print('RescaleType ',ds[0].RescaleType)
    print('ResponsibleOrganization ',ds[0].ResponsibleOrganization)
    print('Rows ',ds[0].Rows)
    print('SAR ',ds[0].SAR)
    print('SOPClassUID ',ds[0].SOPClassUID)
    print('SOPInstanceUID ',ds[0].SOPInstanceUID)
    print('SamplesPerPixel ',ds[0].SamplesPerPixel)
    print('ScanOptions ',ds[0].ScanOptions)
    print('ScanningSequence ',ds[0].ScanningSequence)
    print('SequenceName ',ds[0].SequenceName)
    print('SequenceVariant ',ds[0].SequenceVariant)
    print('SeriesDate ',ds[0].SeriesDate)
    print('SeriesDescription ',ds[0].SeriesDescription)
    print('SeriesInstanceUID ',ds[0].SeriesInstanceUID)
    print('SeriesNumber ',ds[0].SeriesNumber)
    print('SeriesTime ',ds[0].SeriesTime)
    print('SliceLocation ',ds[0].SliceLocation)
    print('SliceThickness ',ds[0].SliceThickness)
    print('SmallestImagePixelValue ',ds[0].SmallestImagePixelValue)
    print('SoftwareVersions ',ds[0].SoftwareVersions)
    print('SpecificCharacterSet ',ds[0].SpecificCharacterSet)
    print('StationName ',ds[0].StationName)
    print('StudyComments ',ds[0].StudyComments)
    print('StudyDate ',ds[0].StudyDate)
    print('StudyDescription ',ds[0].StudyDescription)
    print('StudyID ',ds[0].StudyID)
    print('StudyInstanceUID ',ds[0].StudyInstanceUID)
    print('StudyTime ',ds[0].StudyTime)
    print('TransmitCoilName ',ds[0].TransmitCoilName)
    print('TriggerTime ',ds[0].TriggerTime)
    print('VariableFlipAngleFlag ',ds[0].VariableFlipAngleFlag)
    print('WindowCenter ',ds[0].WindowCenter)
    print('WindowCenterWidthExplanation ',ds[0].WindowCenterWidthExplanation)
    print('WindowWidth ',ds[0].WindowWidth)
    print('dBdt ',ds[0].dBdt)


folder = '2021-04-20_00-33-21_copper_patient1_R10'
name = 'kspc_R10_vn.mat_segmented.h5'

name_file =  folder + '/output/' + name

path = '/home/nfadel/Desktop/sw/hpc-predict/data/v1/decrypt/segmenter/cnn_segmenter/hpc_predict/v1/inference/' + name_file

path = '/home/nfadel/Desktop/v4_r20.h5'
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

## METADATA
file_meta = FileMetaDataset()
####################
ds = [ FileDataset(None, {}, file_meta=file_meta, preamble=b"\0" * 128) ]
####################
## DICOM TAGS
#print('type image_datsets ', type(image_datasets))
#print('type ds[0] ', type(ds[0]))
#print('type image_datsets[0] ', type(image_datasets[0]))
#ds = image_datasets
## Set the transfer syntax
ds[0].is_little_endian = True
ds[0].is_implicit_VR = True
# Set creation date/time
#dt = datetime.datetime.now()
#ds.ContentDate = dt.strftime('%Y%m%d')
#timeStr = dt.strftime('%H%M%S.%f')  # long format with micro seconds
#ds.ContentTime = timeStr

#ds.save_as('output.dcm')



ds[0].AccessionNumber = '0027084519'
ds[0].AcquisitionDate = d.date(2014,5,20)
ds[0].AcquisitionMatrix = [0, 160, 100, 0]
ds[0].AcquisitionNumber = 0
ds[0].AcquisitionTime = d.time(13,30,45,31)
#ds[0].AdditionalPatientHistory
#ds[0].AdmittingDiagnosesDescription
ds[0].AngioFlag = 'N'
ds[0].BitsAllocated = 16
ds[0].BitsStored = 12
ds[0].CardiacNumberOfImages = 44
ds[0].Columns = 120
ds[0].ContentDate = d.date(2014,5,20)
ds[0].ContentTime = d.time(14,11,17,984)
#TODO
#print(type(image_datasets[0]))
#print(type(ds[0]))
image_datasets[0].DerivationCodeSequence.CodeValue = '121327'
#add_new(tag, VR, value)
#ds[0].add_new(DerivationCodeSequence[0].CodeValue, '121327')
#ds[0].DerivationCodeSequence[0].CodeValue = '121327'
#ds[0].DerivationCodeSequence.CodingSchemeDesignator = 'DCM'
#ds[0].DerivationCodeSequence.CodeMeaning = 'Full fidelity image, uncompressed or lossless compressed'

ds[0].DerivationDescription = 'Lossless JPEG compression, selection value 1, point transform 0, compression ratio 1.7085'
ds[0].DeviceSerialNumber = '35042'
ds[0].EchoNumbers = 1
ds[0].EchoTime = '2.573' #TODO is it OK?
ds[0].EchoTrainLength = 1
ds[0].FlipAngle = 7
ds[0].FrameOfReferenceUID = '1.3.12.2.1107.5.2.32.35042.2.20140520124746281.0.0.0'
ds[0].HighBit = 11
ds[0].ImageComments = 'phase difference images\nAPST_HE_C_160'
ds[0].ImageOrientationPatient = [0.2709607843302, 0.96259038710926, -5.22059898e-008, -0.0232690420399, 0.00654997782592, -0.9997077820408]
ds[0].ImagePositionPatient = [-57.148255351844, -136.43397068637, 168.38495554621]
ds[0].ImageType = ['ORIGINAL', 'PRIMARY', 'P', 'DIS2D']
ds[0].ImagedNucleus = '1H'
ds[0].ImagingFrequency = 123.24119
ds[0].InPlanePhaseEncodingDirection = 'ROW'
ds[0].InstanceCreationDate = d.date(2014,5,20)
ds[0].InstanceCreationTime = d.time(14,11,17,984)
ds[0].InstanceNumber = 2275
ds[0].InstitutionAddress = 'Albertiner strasse, 515, Rafalburg /c528a3/,BW,CH,79106'
ds[0].InstitutionName = 'Rad. Univ. Klinik Rafaelburg'
ds[0].InstitutionalDepartmentName = 'Department'
#TODO
#ds[0].IssuerOfAccessionNumberSequence.LocalNamespaceEntityID = '6229UKF'
ds[0].IssuerOfPatientID = '6229UKF'
ds[0].LargestImagePixelValue = 4094
ds[0].MRAcquisitionType = '3D'
ds[0].MagneticFieldStrength = 3
ds[0].Manufacturer = 'SIEMENS'
ds[0].ManufacturerModelName = 'TrioTim'
ds[0].Modality = 'MR'
ds[0].NominalInterval = 1111
ds[0].NumberOfAverages = 1
ds[0].NumberOfPhaseEncodingSteps = 100
ds[0].OperatorsName = 'MS'
#TODO
#ds[0].OriginalAttributesSequence.ModifiedAttributesSequence = '124746.062000'
#ds[0].OriginalAttributesSequence.StudyTime = d.time(12,47,46,062)
#ds[0].OriginalAttributesSequence.AttributeModificationDateTime = d.datetime(2019,07,13,16,19,07)
#ds[0].OriginalAttributesSequence.ModifyingSystem = 'MERLIN'
#ds[0].OriginalAttributesSequence.SourceOfPreviousValues = 'Rad. Univ. Klinik Freiburg'
#ds[0].OriginalAttributesSequence.ReasonForTheAttributeModificati = 'COERCE'
ds[0].OtherPatientIDs = ['17067214', 'f62f9fd1-4ab9-4de8-afdc-3a309ca6d71f']
#ds[0].OtherPatientNames
ds[0].PatientAge = '070Y'
ds[0].PatientBirthDate = d.date(1944,4,3)
#ds[0].PatientComments
ds[0].PatientID = '666666'
ds[0].PatientMotherBirthName = 'Prosdocimo'
ds[0].PatientName = 'Goofy'
ds[0].PatientPosition = 'HFS'
ds[0].PatientSex = 'F'
ds[0].PatientSize = 244
ds[0].PatientWeight = 45
ds[0].PercentPhaseFieldOfView = 75
ds[0].PercentSampling = 83.3333
ds[0].PerformedProcedureStepDescription = 'MR Thorax'
ds[0].PerformedProcedureStepID = 'R05617.27084519'
ds[0].PerformedProcedureStepStartDate = d.date(2014,5,20)
ds[0].PerformedProcedureStepStartTime = d.time(12,48,31,295)
ds[0].PerformingPhysicianName = 'Frankestein'
ds[0].PhotometricInterpretation = 'MONOCHROME2'
ds[0].PixelBandwidth = 453
ds[0].PixelRepresentation = 0
ds[0].PixelSpacing = [2.125, 2.125]
#ds[0].PositionReferenceIndicator
ds[0].PregnancyStatus = 4
#TODO
#ds[0].ProcedureCodeSequence.CodeValue = '6010'
#ds[0].ProcedureCodeSequence.CodingSchemeDesignator ='99GAP'
#ds[0].ProcedureCodeSequence.CodeMeaning = 'MR Thorax'
ds[0].ProtocolName = '4DFlow_Pat5_Tres20'
#TODO
#ds[0].ReferencedImageSequence.ReferencedSOPClassUID = 'MR Image Storage'
#ds[0].ReferencedImageSequence.ReferencedSOPInstanceUID = '1.3.12.2.1107.5.2.32.35042.201405201323376824961820'
#(0008, 1150) Referenced SOP Class UID            UI: MR Image Storage
#(0008, 1155) Referenced SOP Instance UID         UI: 1.3.12.2.1107.5.2.32.35042.201405201323376824961820(0008, 1150) Referenced SOP Class #UID            UI: MR Image Storage
#(0008, 1155) Referenced SOP Instance UID         UI: 1.3.12.2.1107.5.2.32.35042.2014052012483887491258809(0008, 1150) Referenced SOP Class #UID            UI: MR Image Storage
#(0008, 1155) Referenced SOP Instance UID         UI: 1.3.12.2.1107.5.2.32.35042.2014052012513677902359453]
#TODO
#ds[0].ReferencedPerformedProcedureStepSequence.ReferencedSOPClassUID = 'Modality Performed Procedure Step SOP Class'
#ds[0].ReferencedPerformedProcedureStepSequence.ReferencedSOPClassUID = '1.2.40.0.13.1.1.1.193.196.214.172.20140515182142778.7737'
ds[0].ReferringPhysicianName = 'Amb. Neurologie'
ds[0].RepetitionTime = 20
#TODO
#ds[0].RequestAttributesSequence.AccessionNumber = '0027084519'
#ds[0].RequestAttributesSequence.StudyInstanceUID = '1.2.276.0.38.1.1.1.3150.20140520105617.27084519'
#ds[0].RequestAttributesSequence.RequestingService = '200260Amb. Neurologie'
#ds[0].RequestAttributesSequence.RequestedProcedureDescription = 'MR-Angio Thorax'
#print('pota', ds[0].RequestAttributesSequence)
#ds[0].RequestAttributesSequence.RequestedProcedureCodeSequence.CodeValue = '6011-57ad451f'
#ds[0].RequestAttributesSequence.RequestedProcedureCodeSequence.CodingSchemeDesignator = 'GAPIT'
#ds[0].RequestAttributesSequence.RequestedProcedureCodeSequence.CodeMeaning = 'MR-Angio Thorax'
#ds[0].RequestAttributesSequence.ScheduledProcedureStepDescriptio = 'MR-Angio Thorax'
#ds[0].RequestAttributesSequence.ScheduledProtocolCodeSequence.CodeValue = '6011-57ad451f'
#ds[0].RequestAttributesSequence.ScheduledProtocolCodeSequence.CodingSchemeDesignator = 'GAPIT'
#ds[0].RequestAttributesSequence.ScheduledProtocolCodeSequence.CodeMeaning = 'MR-Angio Thorax'
#ds[0].RequestAttributesSequence.ScheduledProcedureStepID = '0027084519'
#ds[0].RequestAttributesSequence.RequestedProcedureID = '0027084519'

#[(0008, 0050) Accession Number                    SH: '0027084519'
#(0020, 000d) Study Instance UID                  UI: 1.2.276.0.38.1.1.1.3150.20140520105617.27084519
#(0032, 1032) Requesting Physician                PN: ''
#(0032, 1033) Requesting Service                  LO: '200260Amb. Neurologie'
#(0032, 1060) Requested Procedure Description     LO: 'MR-Angio Thorax'
#(0032, 1064)  Requested Procedure Code Sequence  1 item(s) ---- 
#   (0008, 0100) Code Value                          SH: '6011-57ad451f'
#   (0008, 0102) Coding Scheme Designator            SH: 'GAPIT'
#   (0008, 0104) Code Meaning                        LO: 'MR-Angio Thorax'
#   ---------
#(0040, 0007) Scheduled Procedure Step Descriptio LO: 'MR-Angio Thorax'
#(0040, 0008)  Scheduled Protocol Code Sequence  1 item(s) ---- 
#   (0008, 0100) Code Value                          SH: '6011-57ad451f'
#   (0008, 0102) Coding Scheme Designator            SH: 'GAPIT'
#   (0008, 0104) Code Meaning                        LO: 'MR-Angio Thorax'
#   ---------
#(0040, 0009) Scheduled Procedure Step ID         SH: '0027084519'
#(0040, 1001) Requested Procedure ID              SH: '0027084519'
#(0040, 1002) Reason for the Requested Procedure  LO: ''
#(0040, 1400) Requested Procedure Comments        LT: '']

#TODO
#ds[0].RequestedProcedureCodeSequence.CodeValue = '6010'
#ds[0].RequestedProcedureCodeSequence.CodingSchemeDesignator = '99GAP'
#ds[0].RequestedProcedureCodeSequence.CodeMeaning = 'MR Thorax'
ds[0].RequestedProcedureDescription = 'MR Thorax'
ds[0].RescaleIntercept = -4096
ds[0].RescaleSlope = 2
ds[0].RescaleType = 'US'
ds[0].ResponsibleOrganization = 'X'
ds[0].Rows = 160
ds[0].SAR = 0.36335781457169
ds[0].SOPClassUID = '1.2.840.10008.5.1.4.1.1.4'
ds[0].SOPInstanceUID = '1.3.12.2.1107.5.2.32.35042.2014052014111617010765431'
ds[0].SamplesPerPixel = 1
ds[0].ScanOptions = 'CT'
ds[0].ScanningSequence = 'GR'
ds[0].SequenceName = 'fl3d1'
ds[0].SequenceVariant = ['SP', 'OSP']
ds[0].SeriesDate = d.date(2014,5,20)
ds[0].SeriesDescription = '4DFlow_Pat5_Tres20'
ds[0].SeriesInstanceUID = '1.3.12.2.1107.5.2.32.35042.2014052013303780789062115.0.0.0'
ds[0].SeriesNumber  = 10
ds[0].SeriesTime = d.time(14,11,12,843)
ds[0].SliceLocation = -22.107260462064
ds[0].SliceThickness = 2.5
ds[0].SmallestImagePixelValue = 4
ds[0].SoftwareVersions = 'syngo MR B17'
ds[0].SpecificCharacterSet = 'ISO_IR 100'
ds[0].StationName = 'MRC20506'
ds[0].StudyComments = 'APST_HE_C_160'
ds[0].StudyDate = d.date(2014,5,20)
ds[0].StudyDescription = 'cv_MRI^vascular'
ds[0].StudyID = '0027084519'
ds[0].StudyInstanceUID = '1.2.276.0.38.1.1.1.3150.20140520105617.27084519'
ds[0].StudyTime = d.time(12,47,46,62)
ds[0].TransmitCoilName = 'Body'
ds[0].TriggerTime = 550
ds[0].VariableFlipAngleFlag = 'N'
ds[0].WindowCenter = -162
ds[0].WindowCenterWidthExplanation = 'Algo1'
ds[0].WindowWidth = 4516
ds[0].dBdt = 0

########################################

if enable_debug:
    debug_print_tags(ds)

name_output_file = "dicom"

if hasattr(mri, 'intensity') and enable_intensity:
    if enable_debug:
        print('all ', mri.intensity[:,:,:,:].size) 
        print('x ', mri.intensity[:,1,1,1].size)
        print('y ', mri.intensity[1,:,1,1].size)
        print('time ', mri.intensity[1,1,1,:].size)
        print('z ', mri.intensity[1,1,:,1].size)
    os.makedirs(name_output_file + "/intensity", exist_ok = True)
    series_instance_uid=generate_uid()
    for mri_time in range(mri.intensity[1,1,1,:].size):
        ds[0].AcquisitionTime = d.time(ds[0].AcquisitionTime.hour,ds[0].AcquisitionTime.minute,ds[0].AcquisitionTime.second,int(mri.time[mri_time]*100))
        ds[0].ContentTime = d.time(ds[0].ContentTime.hour,ds[0].ContentTime.minute,ds[0].ContentTime.second,int(mri.time[mri_time]*100))
        ds[0].EchoTime = str(mri.time[mri_time])
        ds[0].InstanceCreationTime = d.time(ds[0].InstanceCreationTime.hour,ds[0].InstanceCreationTime.minute,ds[0].InstanceCreationTime.second,int(mri.time[mri_time]*100))
#        ds[0].OriginalAttributesSequence.StudyTime = mri.time[mri_time]
#        ds[0].OriginalAttributesSequence.AttributeModificationDateTime = mri.time[mri_time]
        ds[0].PerformedProcedureStepStartTime = d.time(ds[0].PerformedProcedureStepStartTime.hour,ds[0].PerformedProcedureStepStartTime.minute,ds[0].PerformedProcedureStepStartTime.second,int(mri.time[mri_time]*100))
        ds[0].RepetitionTime = str(mri.time[mri_time])
        ds[0].SeriesTime = d.time(ds[0].SeriesTime.hour,ds[0].SeriesTime.minute,ds[0].SeriesTime.second,int(mri.time[mri_time]*100))
        ds[0].StudyTime = d.time(ds[0].StudyTime.hour,ds[0].StudyTime.minute,ds[0].StudyTime.second,int(mri.time[mri_time]*100))
        ds[0].TriggerTime = str(mri.time[mri_time])
        for z_axis in range(mri.intensity[1,1,:,1].size): 
            ds[0].image_type = 'M',
            ds[0].pixel_data=mri.intensity[:,:,z_axis,mri_time],
            ds[0].sop_instance_uid=generate_uid(),
            ds[0].fix_meta_info(True)
            ds[0].save_as(f"dicom/intensity/TS_{mri_time}_INT_{z_axis}.dcm",write_like_original=False)

if hasattr(mri, 'velocity_mean') and enable_velocity_mean:
    if enable_debug:
        print('all ', mri.velocity_mean[:,:,:,:].size)
        print('time ', mri.velocity_mean[1,1,1,:].size)
        print('x ', mri.velocity_mean[:,1,1,1].size)
        print('y ', mri.velocity_mean[1,:,1,1].size)
        print('z ', mri.velocity_mean[1,1,:,1].size)
    os.makedirs(name_output_file + "/velocity_mean", exist_ok = True)
    series_instance_uid=generate_uid()
    for mri_time in range(int(mri.velocity_mean[1,1,1,:].size/3)):
        ds[0].AcquisitionTime = d.time(ds[0].AcquisitionTime.hour,ds[0].AcquisitionTime.minute,ds[0].AcquisitionTime.second,int(mri.time[mri_time]*100))
        ds[0].ContentTime = d.time(ds[0].ContentTime.hour,ds[0].ContentTime.minute,ds[0].ContentTime.second,int(mri.time[mri_time]*100))
        ds[0].EchoTime = str(mri.time[mri_time])
        ds[0].InstanceCreationTime = d.time(ds[0].InstanceCreationTime.hour,ds[0].InstanceCreationTime.minute,ds[0].InstanceCreationTime.second,int(mri.time[mri_time]*100))
#        ds[0].OriginalAttributesSequence.StudyTime = mri.time[mri_time]
#        ds[0].OriginalAttributesSequence.AttributeModificationDateTime = mri.time[mri_time]
        ds[0].PerformedProcedureStepStartTime = d.time(ds[0].PerformedProcedureStepStartTime.hour,ds[0].PerformedProcedureStepStartTime.minute,ds[0].PerformedProcedureStepStartTime.second,int(mri.time[mri_time]*100))
        ds[0].RepetitionTime = str(mri.time[mri_time])
        ds[0].SeriesTime = d.time(ds[0].SeriesTime.hour,ds[0].SeriesTime.minute,ds[0].SeriesTime.second,int(mri.time[mri_time]*100))
        ds[0].StudyTime = d.time(ds[0].StudyTime.hour,ds[0].StudyTime.minute,ds[0].StudyTime.second,int(mri.time[mri_time]*100))
        ds[0].TriggerTime = str(mri.time[mri_time])
        for z_axis in range(int(mri.velocity_mean[1,1,:,1].size/3)): 
            ds[0].image_type = 'P'
            ds[0].pixel_data=mri.velocity_mean[:,:,z_axis,mri_time],
            ds[0].sop_instance_uid=generate_uid(),
            ds[0].fix_meta_info(True)
            ds[0].save_as(f"dicom/velocity_mean/TS_{mri_time}_VELMEAN_{z_axis}.dcm",write_like_original=False)

if hasattr(mri, 'segmentation_prob') and enable_segmentation:
    if enable_debug:
        print('time ',mri.segmentation_prob[1,1,1,:].size)
    image_datasets[0].StudyInstanceUID = generate_uid()
    os.makedirs(name_output_file + "/segmentation", exist_ok = True) 
    series_instance_uid=generate_uid()
    for mri_time in range(mri.segmentation_prob[1,1,1,:].size):
        ds[0].AcquisitionTime = d.time(ds[0].AcquisitionTime.hour,ds[0].AcquisitionTime.minute,ds[0].AcquisitionTime.second,int(mri.time[mri_time]*100))
        ds[0].ContentTime = d.time(ds[0].ContentTime.hour,ds[0].ContentTime.minute,ds[0].ContentTime.second,int(mri.time[mri_time]*100))
        ds[0].EchoTime = str(mri.time[mri_time])
        ds[0].InstanceCreationTime = d.time(ds[0].InstanceCreationTime.hour,ds[0].InstanceCreationTime.minute,ds[0].InstanceCreationTime.second,int(mri.time[mri_time]*100))
#        ds[0].OriginalAttributesSequence.StudyTime = mri.time[mri_time]
#        ds[0].OriginalAttributesSequence.AttributeModificationDateTime = mri.time[mri_time]
        ds[0].PerformedProcedureStepStartTime = d.time(ds[0].PerformedProcedureStepStartTime.hour,ds[0].PerformedProcedureStepStartTime.minute,ds[0].PerformedProcedureStepStartTime.second,int(mri.time[mri_time]*100))
        ds[0].RepetitionTime = str(mri.time[mri_time])
        ds[0].SeriesTime = d.time(ds[0].SeriesTime.hour,ds[0].SeriesTime.minute,ds[0].SeriesTime.second,int(mri.time[mri_time]*100))
        ds[0].StudyTime = d.time(ds[0].StudyTime.hour,ds[0].StudyTime.minute,ds[0].StudyTime.second,int(mri.time[mri_time]*100))
        ds[0].TriggerTime = str(mri.time[mri_time])
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

if hasattr(mri, 'anomaly_prob') and enable_abnormality:
    image_datasets[0].StudyInstanceUID = generate_uid()
    os.makedirs(name_output_file + "/anomaly_prob", exist_ok = True) 
    series_instance_uid=generate_uid()
    for mri_time in range(mri.anomaly_prob[1,1,1,:].size):
        ds[0].AcquisitionTime = d.time(ds[0].AcquisitionTime.hour,ds[0].AcquisitionTime.minute,ds[0].AcquisitionTime.second,int(mri.time[mri_time]*100))
        ds[0].ContentTime = d.time(ds[0].ContentTime.hour,ds[0].ContentTime.minute,ds[0].ContentTime.second,int(mri.time[mri_time]*100))
        ds[0].EchoTime = str(mri.time[mri_time])
        ds[0].InstanceCreationTime = d.time(ds[0].InstanceCreationTime.hour,ds[0].InstanceCreationTime.minute,ds[0].InstanceCreationTime.second,int(mri.time[mri_time]*100))
#        ds[0].OriginalAttributesSequence.StudyTime = mri.time[mri_time]
#        ds[0].OriginalAttributesSequence.AttributeModificationDateTime = mri.time[mri_time]
        ds[0].PerformedProcedureStepStartTime = d.time(ds[0].PerformedProcedureStepStartTime.hour,ds[0].PerformedProcedureStepStartTime.minute,ds[0].PerformedProcedureStepStartTime.second,int(mri.time[mri_time]*100))
        ds[0].RepetitionTime = str(mri.time[mri_time])
        ds[0].SeriesTime = d.time(ds[0].SeriesTime.hour,ds[0].SeriesTime.minute,ds[0].SeriesTime.second,int(mri.time[mri_time]*100))
        ds[0].StudyTime = d.time(ds[0].StudyTime.hour,ds[0].StudyTime.minute,ds[0].StudyTime.second,int(mri.time[mri_time]*100))
        ds[0].TriggerTime = str(mri.time[mri_time])
        for z_axis in range(mri.anomaly_prob[1,1,:,1].size): 
            # Create the Segmentation instance
            seg_dataset = Segmentation(
                    source_images=ds,
                    pixel_array=mri.anomaly_prob[:,:,z_axis,mri_time],
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
                    )
            seg_dataset.fix_meta_info(True)
            seg_dataset.save_as(f"dicom/anomaly_prob/TS_{mri_time}_AN_{z_axis}.dcm",write_like_original=False)


