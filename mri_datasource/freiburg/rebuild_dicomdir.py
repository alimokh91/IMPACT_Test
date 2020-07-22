import gdcm

self._ts = gdcm.TransferSyntax(gdcm.TransferSyntax.ExplicitVRLittleEndian)
self._image.SetTransferSyntax(self._ts )
dcmdir = gdcm.DICOMDIRGenerator();
d = gdcm.Directory()

nfiles = d.Load(os.path.abspath('./DICOMDIR'), True);
self._logger.info(nfiles)

FilenamesT = gdcm.FilenamesType()
FilenamesT = d.GetFilenames()
dcmdir.SetRootDirectory(d.GetToplevel())
dcmdir.SetFilenames(FilenamesT)
dcmdir.SetDescriptor(basename)
self._logger.info('starting dmcdir')

if not dcmdir.Generate():
    self._logger.error('Error generating dcmdir')
else:
    self._logger.info('Success gerenating dmcdir')
    gdcm.FileMetaInformation.SetSourceApplicationEntityTitle("GenerateDICOMDIR" )
    writer = gdcm.Writer()
    writer.SetFile(dcmdir.GetFile() )
    writer.SetFileName( basename )
    
    if not writer.Write():
        self._logger.error('Error saving dcmdir')
    else:
        self._logger.info('Success saving dmcdir')
