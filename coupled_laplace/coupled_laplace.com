# read data for different regions
gfx read node CoupledLaplace_1.part0.exnode node_offset 0
gfx read element CoupledLaplace_1.part0.exelem element_offset 0 node_offset 0
gfx define faces egroup Region1

gfx read node CoupledLaplace_2.part0.exnode node_offset 1000
gfx read element CoupledLaplace_2.part0.exelem element_offset 100 node_offset 1000
gfx define faces egroup Region2

gfx read node CoupledLaplace_Interface.part0.exnode node_offset 100000
gfx read element CoupledLaplace_Interface.part0.exelem element_offset 1000 node_offset 100000
gfx define faces egroup Interface

# create window
gfx create window 1
gfx modify window 1 background colour 1 1 1
gfx modify window 1 view interest_point 2.0,0.5,0.0 eye_point 2.0,0.5,10.0 up_vector 0.0,1.0,0.0

# modify spectrum
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range -1.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1

# color spheres with nodal values
gfx modify g_element Region1 node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 data Phi
gfx modify g_element Region2 node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 data Phi
gfx modify g_element Interface node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 data Lambda

# color element lines
gfx modify g_element Region1 lines select_on material default selected_material default_selected data Phi
gfx modify g_element Region2 lines select_on material default selected_material default_selected data Phi
gfx modify g_element Interface lines select_on material default selected_material default_selected data Lambda

# color element surfaces
gfx modify g_element Region1 surfaces select_on material default data Phi spectrum default selected_material default_selected render_shaded
gfx modify g_element Region2 surfaces select_on material default data Phi spectrum default selected_material default_selected render_shaded

# open scene editor
gfx edit scene
