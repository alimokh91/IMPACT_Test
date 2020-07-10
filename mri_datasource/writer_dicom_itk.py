#!/usr/bin/env python
import os

import itk
import sys

if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] +
          " [DicomDirectory [outputFileName [seriesName]]]")
    print("If DicomDirectory is not specified, current directory is used\n")

# current directory by default
dirName = '.'
if len(sys.argv) > 1:
    dirName = sys.argv[1]


print('dirName', dirName)
PixelType = itk.ctype('signed short')
Dimension = 3

ImageType = itk.Image[PixelType, Dimension]

namesGenerator = itk.GDCMSeriesFileNames.New()
namesGenerator.SetUseSeriesDetails(True)
namesGenerator.AddSeriesRestriction("0080|0020")
namesGenerator.SetGlobalWarningDisplay(True)
namesGenerator.SetDirectory(dirName)

seriesUID = namesGenerator.GetSeriesUIDs()

if len(seriesUID) < 1:
    print('No DICOMs in: ' + dirName)
    sys.exit(1)

print('The directory: ' + dirName)
print('Contains the following DICOM Series: ')
for uid in seriesUID:
    print(uid)

seriesFound = False
for uid in seriesUID:
    seriesIdentifier = uid
    if len(sys.argv) > 3:
        seriesIdentifier = sys.argv[3]
        seriesFound = True
    print('Reading: ' + seriesIdentifier)
    fileNames = namesGenerator.GetFileNames(seriesIdentifier)


    writer = itk.ImageFileWriter[ImageType].New()
    outFileName = os.path.join(dirName, seriesIdentifier + '.nrrd')
    if len(sys.argv) > 2:
        outFileName = sys.argv[2]
    writer.SetFileName(outFileName)
    writer.UseCompressionOn()
    writer.SetInput(reader.GetOutput())
    print('Writing: ' + outFileName)
    writer.Update()

    if seriesFound:
        break
