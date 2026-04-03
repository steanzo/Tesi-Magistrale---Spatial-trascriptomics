import numpy as np
import spatialdata as sd
import spatialdata_io as sio
import spatialdata_plot
from napari_spatialdata import Interactive

from dask_image.imread import imread
from xarray import DataArray
from spatialdata.models import Image2DModel
from spatialdata.models import PointsModel
from spatialdata.transformations.transformations import Identity

# ho fatto C_26_T_16, C_26_T_17, C_7_T_9
#

c3t5 = sd.read_zarr("dati/c26t16/c26t16_align1.zarr/")
c3t5


c3t5.table.to_df()["Foxo3"].max()





# rinominiamo aligned in global
rename_dict = {"aligned": "global"}
c3t5.rename_coordinate_systems(rename_dict)


# Se mi perdo C3T5 (Shapes)
# ccc = sd.read_zarr("dati/c3t5/c3t5.zarr/")
# c3t5.shapes["C3T5"] = ccc.shapes["C3T5"]
# c3t5


# c3t5["C3T5"]
# c3t5["C3T5_full_image"]




Interactive(c3t5)


c3t5.points 



from spatialdata.transformations import (
    align_elements_using_landmarks,
    get_transformation_between_landmarks
    )

affine = get_transformation_between_landmarks(
    references_coords=c3t5["Visium_landmarks"], moving_coords=c3t5["IF_landmarks"]
)
affine

c3t5





# allineo merge a full_image
affine = align_elements_using_landmarks(
    references_coords=c3t5["Visium_landmarks"],
    moving_coords=c3t5["IF_landmarks"],
    reference_element=c3t5["C3T5_full_image"],
    moving_element=c3t5["C3T5_merge"],
    reference_coordinate_system="global",
    moving_coordinate_system="global",
    new_coordinate_system="aligned"
)
affine


# allineo GFP a full_image
# affine = align_elements_using_landmarks(
#     references_coords=c3t5["Visium_landmarks"],
#     moving_coords=c3t5["IF_landmarks"],
#     reference_element=c3t5["C3T5_full_image"],
#     moving_element=c3t5["C3T5_GFP"],
#     reference_coordinate_system="global",
#     moving_coordinate_system="global",
#     new_coordinate_system="aligned"
# )
# affine



# allineo WGA a full_image
affine = align_elements_using_landmarks(
    references_coords=c3t5["Visium_landmarks"],
    moving_coords=c3t5["IF_landmarks"],
    reference_element=c3t5["C3T5_full_image"],
    moving_element=c3t5["C3T5_WGA"],
    reference_coordinate_system="global",
    moving_coordinate_system="global",
    new_coordinate_system="aligned"
)
affine

# inserisco lo C3T5 Shapes
affine = align_elements_using_landmarks(
    references_coords=c3t5["Visium_landmarks"],
    moving_coords=c3t5["Visium_landmarks"],
    reference_element=c3t5["C3T5_full_image"],
    moving_element=c3t5["C3T5"],
    reference_coordinate_system="global",
    moving_coordinate_system="global",
    new_coordinate_system="aligned"
)
affine

# allineo dapi a full_image
affine = align_elements_using_landmarks(
    references_coords=c3t5["Visium_landmarks"],
    moving_coords=c3t5["IF_landmarks"],
    reference_element=c3t5["C3T5_full_image"],
    moving_element=c3t5["C3T5_dapi"],
    reference_coordinate_system="global",
    moving_coordinate_system="global",
    new_coordinate_system="aligned"
)
affine


# allineo bungaro a full_image
affine = align_elements_using_landmarks(
    references_coords=c3t5["Visium_landmarks"],
    moving_coords=c3t5["IF_landmarks"],
    reference_element=c3t5["C3T5_full_image"],
    moving_element=c3t5["C3T5_bungaro"],
    reference_coordinate_system="global",
    moving_coordinate_system="global",
    new_coordinate_system="aligned"
)
affine




c3t5


sdata = c3t5.transform_to_coordinate_system("aligned")
sdata

Interactive(sdata)



sdata.write("dati/c3t5/c3t5_align1.zarr", overwrite=True)
# lust save in align1

