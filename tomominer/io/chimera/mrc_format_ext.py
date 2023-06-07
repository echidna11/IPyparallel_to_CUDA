

# -----------------------------------------------------------------------------
# Read MRC or CCP4 map file format electron microscope data.
# Byte swapping will be done if needed.
#

from VolumeData.mrc.mrc_format import MRC_Data 

class MRC_Data_ext(MRC_Data):

    def __init__(self, path, file_type):
        MRC_Data.__init__(self, path, file_type)


    # modified version of permute_matrix_to_xyz_axis_order(). The original version has a bug 
    def permute_matrix_to_xyz_axis_order(self, matrix):
        if self.ijk_to_crs == [0,1,2]:
            return matrix

        kji_to_src = [2-self.ijk_to_crs[2-a] for a in (0,1,2)]
        m = matrix.transpose(kji_to_src)

        return m

