set cartoon_transparency, 0
set ray_opaque_background, 1
set ray_shadows, 0
set valence, 0
set ray_trace_mode, 1
bg_color white

load PDB/7my8_fit.pdb, FP
spectrum b, blue white, minimum=0, maximum=80

hide all
show cartoon, FP
show sticks, resi 840+851 and (not name c+n+o)
util.cnc  resi 840+851
hide (hydro)

set_view (\
     0.071482144,   -0.236539513,   -0.968991101,\
     0.323477566,   -0.913463891,    0.246848419,\
    -0.943528533,   -0.331097096,    0.011221834,\
    -0.000026897,    0.000088301, -110.688896179,\
    28.818910599,   25.406803131,   26.836227417,\
   -17.465654373,  238.792404175,  -20.000000000 )
ray; png graph/structure_FP.png

