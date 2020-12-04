<schema xmlns="http://purl.oclc.org/dsdl/schematron">

  <title>FlowMRI schema in Schematron</title>
  <!--
    Developed with XPath evaluation on `h5dump -H -x ...h5 | sed 's///g'` in IntelliJ IDEA and
    Can be tested in Python with

    from lxml import etree; from lxml.isoschematron import Schematron
    sch = Schematron(etree.parse('flow_CStest_Volunteer_R4_h5dump_H.sch'), store_report=True, error_finder=Schematron.ASSERTS_AND_REPORTS)
    sch.validate(etree.parse('flow_CStest_Volunteer_R4_h5dump_H.xml'))
    print(str(sch.validation_report).replace('\xa0',' '))
  -->

  <ns prefix="hdf5" uri="http://support.hdfgroup.org/HDF5/XML/schema/HDF5-File.xsd"/>

  <pattern>
    <title>Structural validation (presence of all datasets and attributes)</title>
    <rule context="//hdf5:Group[@Name='flow-mri']">
      <assert test="hdf5:Attribute[@Name='t_heart_cycle_period']"> t_heart_cycle_period attribute missing in <value-of select="name(.)" /> . </assert>
      <assert test="hdf5:Dataset[@Name='t_coordinates']"         > t_coordinates dataset missing in <value-of select="name(.)" /> . </assert>
      <assert test="hdf5:Dataset[@Name='x_coordinates']"         > x_coordinates dataset missing in <value-of select="name(.)" /> . </assert>
      <assert test="hdf5:Dataset[@Name='y_coordinates']"         > y_coordinates dataset missing in <value-of select="name(.)" /> . </assert>
      <assert test="hdf5:Dataset[@Name='z_coordinates']"         > z_coordinates dataset missing in <value-of select="name(.)" /> . </assert>
      <assert test="hdf5:Dataset[@Name='intensity'    ]"         > intensity dataset missing in <value-of select="name(.)" /> . </assert>
      <assert test="hdf5:Dataset[@Name='velocity_mean']"         > velocity_mean dataset missing in <value-of select="name(.)" /> . </assert>
      <assert test="hdf5:Dataset[@Name='velocity_cov' ]"         > velocity_cov dataset missing in <value-of select="name(.)" /> . </assert>
    </rule>
  </pattern>

  <pattern>
    <title>Dimensional validation (compatibility of spatiotemporal dimensions)</title>
    <rule context="//hdf5:Group[@Name='flow-mri']">
      <let name="t_heart_cycle_period_ds" value="hdf5:Attribute[@Name='t_heart_cycle_period']/hdf5:Dataspace"/>
      <let name="t_coordinates_ds" value="hdf5:Dataset[@Name='t_coordinates']/hdf5:Dataspace"/>
      <let name="x_coordinates_ds" value="hdf5:Dataset[@Name='x_coordinates']/hdf5:Dataspace"/>
      <let name="y_coordinates_ds" value="hdf5:Dataset[@Name='y_coordinates']/hdf5:Dataspace"/>
      <let name="z_coordinates_ds" value="hdf5:Dataset[@Name='z_coordinates']/hdf5:Dataspace"/>

      <let name="t_dim" value="$t_coordinates_ds/*/hdf5:Dimension/@DimSize"/>
      <let name="x_dim" value="$x_coordinates_ds/*/hdf5:Dimension/@DimSize"/>
      <let name="y_dim" value="$y_coordinates_ds/*/hdf5:Dimension/@DimSize"/>
      <let name="z_dim" value="$z_coordinates_ds/*/hdf5:Dimension/@DimSize"/>

      <let name="intensity_ds"     value="hdf5:Dataset[@Name='intensity']/hdf5:Dataspace"    />
      <let name="velocity_mean_ds" value="hdf5:Dataset[@Name='velocity_mean']/hdf5:Dataspace"/>
      <let name="velocity_cov_ds"  value="hdf5:Dataset[@Name='velocity_cov' ]/hdf5:Dataspace"/>

      <assert test="$t_heart_cycle_period_ds/hdf5:ScalarDataspace"> t_heart_cycle_period has invalid dataspace <value-of select="name($t_heart_cycle_period_ds/*)" /> . </assert>
      <assert test="$t_coordinates_ds/*/@Ndims=1"> t_coordinates number of dimensions <value-of select="$t_coordinates_ds/*/@Ndims" /> should be 1. </assert>
      <assert test="$x_coordinates_ds/*/@Ndims=1"> x_coordinates number of dimensions <value-of select="$x_coordinates_ds/*/@Ndims" /> should be 1. /> . </assert>
      <assert test="$y_coordinates_ds/*/@Ndims=1"> y_coordinates number of dimensions <value-of select="$y_coordinates_ds/*/@Ndims" /> should be 1. /> . </assert>
      <assert test="$z_coordinates_ds/*/@Ndims=1"> z_coordinates number of dimensions <value-of select="$z_coordinates_ds/*/@Ndims" /> should be 1. /> . </assert>

      <assert test="$intensity_ds/*/@Ndims=4"    > intensity number of dimensions = <value-of select="$intensity_ds/*/@Ndims" /> should be 4. /> . </assert>
      <assert test="$intensity_ds/*/hdf5:Dimension[1]/@DimSize=$z_dim"> intensity_dim[0] != z_dim: <value-of select="$intensity_ds/*/hdf5:Dimension[1]/@DimSize" /> != <value-of select="$z_dim" />. </assert>
      <assert test="$intensity_ds/*/hdf5:Dimension[2]/@DimSize=$y_dim"> intensity_dim[1] != y_dim: <value-of select="$intensity_ds/*/hdf5:Dimension[2]/@DimSize" /> != <value-of select="$y_dim" />. </assert>
      <assert test="$intensity_ds/*/hdf5:Dimension[3]/@DimSize=$x_dim"> intensity_dim[2] != x_dim: <value-of select="$intensity_ds/*/hdf5:Dimension[3]/@DimSize" /> != <value-of select="$x_dim" />. </assert>
      <assert test="$intensity_ds/*/hdf5:Dimension[4]/@DimSize=$t_dim"> intensity_dim[3] != t_dim: <value-of select="$intensity_ds/*/hdf5:Dimension[4]/@DimSize" /> != <value-of select="$t_dim" />. </assert>

      <assert test="$velocity_mean_ds/*/@Ndims=5"> velocity_mean number of dimensions = <value-of select="velocity_mean_ds/*/@Ndims" /> should be 5. /> . </assert>
      <assert test="$velocity_mean_ds/*/hdf5:Dimension[1]/@DimSize=$z_dim"> velocity_mean_dim[0] != z_dim: <value-of select="$velocity_mean_ds/*/hdf5:Dimension[1]/@DimSize" /> != <value-of select="$z_dim" />. </assert>
      <assert test="$velocity_mean_ds/*/hdf5:Dimension[2]/@DimSize=$y_dim"> velocity_mean_dim[1] != y_dim: <value-of select="$velocity_mean_ds/*/hdf5:Dimension[2]/@DimSize" /> != <value-of select="$y_dim" />. </assert>
      <assert test="$velocity_mean_ds/*/hdf5:Dimension[3]/@DimSize=$x_dim"> velocity_mean_dim[2] != x_dim: <value-of select="$velocity_mean_ds/*/hdf5:Dimension[3]/@DimSize" /> != <value-of select="$x_dim" />. </assert>
      <assert test="$velocity_mean_ds/*/hdf5:Dimension[4]/@DimSize=$t_dim"> velocity_mean_dim[3] != t_dim: <value-of select="$velocity_mean_ds/*/hdf5:Dimension[4]/@DimSize" /> != <value-of select="$t_dim" />. </assert>
      <assert test="$velocity_mean_ds/*/hdf5:Dimension[5]/@DimSize=3"     > velocity_mean_dim[4] != 3: <value-of select="$velocity_mean_ds/*/hdf5:Dimension[5]/@DimSize" /> != 3. </assert>

      <assert test="$velocity_cov_ds/*/@Ndims=6" > velocity_cov number of dimensions = <value-of select="$velocity_cov_ds/*/@Ndims" /> should be 6. /> . </assert>
      <assert test="$velocity_cov_ds/*/hdf5:Dimension[1]/@DimSize=$z_dim"> velocity_cov_dim[0] != z_dim: <value-of select="$velocity_cov_ds/*/hdf5:Dimension[1]/@DimSize" /> != <value-of select="$z_dim" />. </assert>
      <assert test="$velocity_cov_ds/*/hdf5:Dimension[2]/@DimSize=$y_dim"> velocity_cov_dim[1] != y_dim: <value-of select="$velocity_cov_ds/*/hdf5:Dimension[2]/@DimSize" /> != <value-of select="$y_dim" />. </assert>
      <assert test="$velocity_cov_ds/*/hdf5:Dimension[3]/@DimSize=$x_dim"> velocity_cov_dim[2] != x_dim: <value-of select="$velocity_cov_ds/*/hdf5:Dimension[3]/@DimSize" /> != <value-of select="$x_dim" />. </assert>
      <assert test="$velocity_cov_ds/*/hdf5:Dimension[4]/@DimSize=$t_dim"> velocity_cov_dim[3] != t_dim: <value-of select="$velocity_cov_ds/*/hdf5:Dimension[4]/@DimSize" /> != <value-of select="$t_dim" />. </assert>
      <assert test="$velocity_cov_ds/*/hdf5:Dimension[5]/@DimSize=3"     > velocity_cov_dim[4] != 3: <value-of select="$velocity_cov_ds/*/hdf5:Dimension[5]/@DimSize" /> != 3. </assert>
      <assert test="$velocity_cov_ds/*/hdf5:Dimension[6]/@DimSize=3"     > velocity_cov_dim[5] != 3: <value-of select="$velocity_cov_ds/*/hdf5:Dimension[6]/@DimSize" /> != 3. </assert>
    </rule>
  </pattern>

  <pattern>
    <title>Floating point format validation (should be double precision)</title>
    <rule context="//hdf5:Dataset/hdf5:DataType/hdf5:AtomicType/hdf5:FloatType">
      <assert test="@ByteOrder='LE' and @Size='8' and @SignBitLocation='63' and @ExponentBits='11' and @ExponentLocation='52' and @MantissaBits='52' and @MantissaLocation='0'"> Float type attributes not double-precision-compatible. /> . </assert>
    </rule>
  </pattern>

</schema>