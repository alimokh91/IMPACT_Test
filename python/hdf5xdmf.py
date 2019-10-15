import h5py
import os
import sys

if len(sys.argv) < 2:
    print("Usage: {} hdf5-file".format(sys.argv[0]))
    exit(1)

group = "flow-mri"
necesary_fiels = ('x_coordinates', 'y_coordinates', 'z_coordinates', 't_coordinates', 'intensity', 'velocity_mean')
hdf5FileName = sys.argv[1]


hdf5file = h5py.File(hdf5FileName, 'r')
if not group in hdf5file:
    print("Missing group {} in file {}".format(group, hdf5FileName))
    exit(1)
hdf5file = hdf5file[group]

for field in necesary_fiels:
    if not field in hdf5file:
        print("Missing field {} in file {}".format(field, hdf5FileName))
        exit(1)

for field in ('intensity', 'velocity_mean'):
    for idx,coordinates in zip((0,1,2,3), ('z_coordinates', 'y_coordinates', 'x_coordinates', 't_coordinates')):
        if hdf5file[field].shape[idx] != hdf5file[coordinates].shape[0]:
            print("The field {} does not have the correct dimensions for {}".format(field, coordinates))
            exit(1)

xdmfFile = os.path.splitext(hdf5FileName)[0]+'.xdmf'
xdmfFile = open(xdmfFile, 'w')

zDims = hdf5file['z_coordinates'].shape[0]
yDims = hdf5file['y_coordinates'].shape[0]
xDims = hdf5file['x_coordinates'].shape[0]
tDims = hdf5file['t_coordinates'].shape[0]

xdmfFile.write("<Xdmf>\n")
xdmfFile.write("  <Domain>\n")
timeSlices = hdf5file['t_coordinates']
xdmfFile.write('    <Grid Name="GridTime" GridType="Collection" CollectionType="Temporal">\n')
for time,idx in zip(timeSlices, range(len(timeSlices))):
    xdmfFile.write('      <Grid Name="Grid" GridType="Uniform">\n')
    xdmfFile.write('        <Time Value="{}" />\n'.format(time))
    xdmfFile.write('        <Topology TopologyType="3DRectMesh" Dimensions="{} {} {}" />\n'.format(zDims, yDims, xDims))
    xdmfFile.write('        <Geometry Type="VXVYVZ">\n')
#    xdmfFile.write('        <Geometry Type="ORIGIN_DXDYDZ">\n')
    xdmfFile.write('          <DataItem Name="X" Format="HDF" NumberType="Float" Precision="8" Dimensions="{}">{}:/{}/x_coordinates</DataItem>\n'.format(xDims, os.path.basename(hdf5FileName), group))
    xdmfFile.write('          <DataItem Name="Y" Format="HDF" NumberType="Float" Precision="8" Dimensions="{}">{}:/{}/y_coordinates</DataItem>\n'.format(yDims, os.path.basename(hdf5FileName), group))
    xdmfFile.write('          <DataItem Name="Z" Format="HDF" NumberType="Float" Precision="8" Dimensions="{}">{}:/{}/z_coordinates</DataItem>\n'.format(zDims, os.path.basename(hdf5FileName), group))
    xdmfFile.write('          <DataItem Format="XML" Dimensions="3">0.0 0.0 0.0</DataItem>\n')
    xdmfFile.write('          <DataItem Format="XML" Dimensions="3">1.0 1.0 1.0</DataItem>\n')
    xdmfFile.write('        </Geometry>\n')
    xdmfFile.write('        <Attribute AttributeType="Scalar" Name="Intensity" Center="Node">\n')
    xdmfFile.write('          <DataItem ItemType="HyperSlab" Dimensions="{} {} {}" Format="XML">\n'.format(zDims, yDims, xDims))
    xdmfFile.write('            <DataItem Dimensions="3 4" Format="XML">\n')
    xdmfFile.write('              0 0 0 {} <!-- coordinate start -->\n'.format(idx))
    xdmfFile.write('              1 1 1 1 <!-- stride -->\n')
    xdmfFile.write('              {} {} {} 1 <!-- count -->\n'.format(zDims, yDims, xDims))
    xdmfFile.write('            </DataItem>\n')
    xdmfFile.write('            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{} {} {} {}">\n'.format(zDims, yDims, xDims, tDims))
    xdmfFile.write('              {}:/{}/intensity\n'.format(os.path.basename(hdf5FileName), group))
    xdmfFile.write('            </DataItem>\n')
    xdmfFile.write('          </DataItem>\n')
    xdmfFile.write('        </Attribute>\n')
    xdmfFile.write('        <Attribute AttributeType="Vector" Name="Velocity" Center="Node">\n')
    xdmfFile.write('          <DataItem ItemType="HyperSlab" Dimensions="{} {} {} 3" Format="XML">\n'.format(zDims, yDims, xDims))
    xdmfFile.write('            <DataItem Dimensions="3 5" Format="XML">\n')
    xdmfFile.write('              0 0 0 {} 0<!-- coordinate start -->\n'.format(idx))
    xdmfFile.write('              1 1 1 1 1 <!-- stride -->\n')
    xdmfFile.write('              {} {} {} 1 3 <!-- count -->\n'.format(zDims, yDims, xDims))
    xdmfFile.write('            </DataItem>\n')
    xdmfFile.write('            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{} {} {} {} 3">\n'.format(zDims, yDims, xDims, tDims))
    xdmfFile.write('              {}:/{}/velocity_mean\n'.format(os.path.basename(hdf5FileName), group))
    xdmfFile.write('            </DataItem>\n')
    xdmfFile.write('          </DataItem>\n')
    xdmfFile.write('        </Attribute>\n')
    xdmfFile.write('      </Grid>\n')
xdmfFile.write('    </Grid>\n')
xdmfFile.write('  </Domain>\n')
xdmfFile.write('</Xdmf>\n')


