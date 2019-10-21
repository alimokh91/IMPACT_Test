import h5py
import jinja2
import os
import sys


xdmfTemplate = """
<Xdmf>
  <Domain>
    <Grid Name="GridTime" GridType="Collection" CollectionType="Temporal">
      {%- for curIdx in range(timeSlices.shape[0]) %}
      <Grid Name="Grid" GridType="Uniform">
        <Time Value="{{ timeSlices[curIdx] }}" />
        <Topology TopologyType="3DRectMesh" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }}" />
        <Geometry Type="VXVYVZ">
          <DataItem Name="X" Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ xDims }}">{{ hdf5FileName }}:/{{ group }}/x_coordinates</DataItem>
          <DataItem Name="Y" Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ yDims }}">{{ hdf5FileName }}:/{{ group }}/y_coordinates</DataItem>
          <DataItem Name="Z" Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ zDims }}">{{ hdf5FileName }}:/{{ group }}/z_coordinates</DataItem>
        </Geometry>
        {%- for scalar_field in scalar_fields %}
        <Attribute AttributeType="Scalar" Name="{{ scalar_field }}" Center="Node">
          <DataItem ItemType="HyperSlab" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }}" Format="XML">
            <DataItem Dimensions="3 4" Format="XML">
              0 0 0 {{ curIdx }} <!-- coordinate start -->
              1 1 1 1 <!-- stride -->
              {{ zDims }} {{ yDims }} {{ xDims }} 1 <!-- count -->
            </DataItem>
            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }} {{ tDims }}">
              {{ hdf5FileName }}:/{{ group }}/{{ scalar_field }}
            </DataItem>
          </DataItem>
        </Attribute>
        {%- endfor %}
        {%- for vector_field in vector_fields %}
        <Attribute AttributeType="Vector" Name="{{ vector_field }}" Center="Node">
          <DataItem ItemType="HyperSlab" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }} 3" Format="XML">
            <DataItem Dimensions="3 5" Format="XML">
              0 0 0 {{ curIdx }} 0<!-- coordinate start -->
              1 1 1 1 1 <!-- stride -->
              {{ zDims }} {{ yDims }} {{ xDims }} 1 3 <!-- count -->
            </DataItem>
            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }} {{ tDims }} 3">
              {{ hdf5FileName }}:/{{ group }}/{{ vector_field }}
            </DataItem>
          </DataItem>
        </Attribute>
        <Attribute AttributeType="Scalar" Name="{{ vector_field }}_magnitude" Center="Node">
          <DataItem ItemType="Function" Function="SQRT($0*$0 + $1*$1 + $2*$2)" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }}">
            <DataItem ItemType="HyperSlab" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }}" Format="XML">
              <DataItem Dimensions="3 5" Format="XML">
                0 0 0 {{ curIdx }} 0 <!-- coordinate start -->
                1 1 1 1 1 <!-- stride -->
                {{ zDims }} {{ yDims }} {{ xDims }} 1 1 <!-- count -->
              </DataItem>
              <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }} {{ tDims }} 3">
                {{ hdf5FileName }}:/{{ group }}/{{ vector_field }}
              </DataItem>
            </DataItem>
            <DataItem ItemType="HyperSlab" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }}" Format="XML">
              <DataItem Dimensions="3 5" Format="XML">
                0 0 0 {{ curIdx }} 1 <!-- coordinate start -->
                1 1 1 1 1 <!-- stride -->
                {{ zDims }} {{ yDims }} {{ xDims }} 1 1 <!-- count -->
              </DataItem>
              <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }} {{ tDims }} 3">
                {{ hdf5FileName }}:/{{ group }}/{{ vector_field }}
              </DataItem>
            </DataItem>
            <DataItem ItemType="HyperSlab" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }}" Format="XML">
              <DataItem Dimensions="3 5" Format="XML">
                0 0 0 {{ curIdx }} 2 <!-- coordinate start -->
                1 1 1 1 1 <!-- stride -->
                {{ zDims }} {{ yDims }} {{ xDims }} 1 1 <!-- count -->
              </DataItem>
              <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{{ zDims }} {{ yDims }} {{ xDims }} {{ tDims }} 3">
                {{ hdf5FileName }}:/{{ group }}/{{ vector_field }}
              </DataItem>
            </DataItem>
          </DataItem>
        </Attribute>
        {%- endfor %}
      </Grid>
      {%- endfor %}
    </Grid>
  </Domain>
</Xdmf>
"""

if len(sys.argv) < 2:
    print("Usage: {} hdf5-file".format(sys.argv[0]))
    exit(1)

group = "flow-mri"
necesary_fiels = ('x_coordinates', 'y_coordinates', 'z_coordinates', 't_coordinates')
hdf5FileName = sys.argv[1]


hdf5file = h5py.File(hdf5FileName, 'r')
if len(hdf5file)==0:
    print("Missing group in file {}".format(hdf5FileName))
group =  list(hdf5file.keys())[0]
hdf5file = hdf5file[group]

for field in necesary_fiels:
    if not field in hdf5file:
        print("Missing field {} in file {}".format(field, hdf5FileName))
        exit(1)

scalar_fields = []
vector_fields = []
for field in hdf5file:
    if len(hdf5file[field].shape) >=4:
        for idx,coordinates in zip((0,1,2,3), ('z_coordinates', 'y_coordinates', 'x_coordinates', 't_coordinates')):
            if hdf5file[field].shape[idx] != hdf5file[coordinates].shape[0]:
                print("The field {} does not have the correct dimensions for {}. Skipping it in the xdmf file".format(field, coordinates))
                continue
    if len(hdf5file[field].shape) == 4:
        scalar_fields.append(field)
    if len(hdf5file[field].shape) == 5:
        if hdf5file[field].shape[4] != 3:
            print("The field {} has {} components. Expected would be 3 components for a vector field".format(field. hdf5file[field].shape[4]))
        else:
            vector_fields.append(field)

xdmfFile = os.path.splitext(hdf5FileName)[0]+'.xdmf'
xdmfFile = open(xdmfFile, 'w')

zDims = hdf5file['z_coordinates'].shape[0]
yDims = hdf5file['y_coordinates'].shape[0]
xDims = hdf5file['x_coordinates'].shape[0]
tDims = hdf5file['t_coordinates'].shape[0]

xdmfStr = jinja2.Template(xdmfTemplate).render(dict(xDims=xDims, yDims=yDims, zDims=zDims, tDims=tDims, scalar_fields=scalar_fields,
    vector_fields=vector_fields, hdf5FileName=os.path.basename(hdf5FileName), group=group, timeSlices=hdf5file['t_coordinates']))
xdmfFile.write(xdmfStr)

