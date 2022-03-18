# UGRIDOutPutSolver

## Documentation
TO DO

## Plotting results with QGIS

UGRID netcdf files can be visualized in QGIS using the [CrayFish plugin](https://www.lutraconsulting.co.uk/projects/crayfish/).  

Remarks:  

- To plot vectors the components should have the attribute "long_name" set to "u component of .." and "v component of .."; this can be added using ncatted:  

`
ncatted -a long_name,uobs 1,o,c,"u component of uobs" -a long_name,uobs 2,o,c,"v component of uobs" Case.nc
`

- Does not seems to like if there is variables of type "int".

## Plotting results with Paraview

There is a [UGRID reader for paraview](https://github.com/FeliciaBrisc/UGRID-Reader-for-ParaView). TO TEST.
