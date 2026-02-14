from pyevtk.hl import imageToVTK

def export_data(field, grid, path):

    d = grid.d
    Ex = field.Ex_exp
    Ey = field.Ey_exp
    Ez = field.Ez_exp

    Bx = field.Bx_exp
    By = field.By_exp
    Bz = field.Bz_exp

    imageToVTK(
        str("output/" + path),
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