from pyevtk.hl import imageToVTK
import os.path
from pyvista import examples

# define export path
export_folder = "./output"


def export_data(field, grid, frame_id):

    d = grid.d
    Ex = field.Ex_exp
    Ey = field.Ey_exp
    Ez = field.Ez_exp

    Bx = field.Bx_exp
    By = field.By_exp
    Bz = field.Bz_exp
    
    frame_path = os.path.join(export_folder, frame_id)
    
    imageToVTK(
        frame_path,
        spacing = (d,d,d),
        pointData = {
            "E_x": Ex,
            "E_y": Ey,
            "E_z": Ez,

            "B_x": Bx,
            "B_y": By,
            "B_z": Bz,
        }
    )