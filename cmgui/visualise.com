#@exnodes=<./MainTime_1_*.part*.exnode>;
#@exelems=<./MainTime_1_*.part*.exelem>;
#foreach $filename (@exnodes) {
#    print "Reading $filename\n";
#    gfx read node "$filename";
#}
#foreach $filename (@exelems) {
#    print "Reading $filename\n";
#    gfx read elem "$filename";
#}

#gfx read node MainTime_1_0.part0.exnode
#gfx read elem MainTime_1_0.part0.exelem

gfx read node EMGExample_M.part0.exnode
gfx read elem EMGExample_M.part0.exelem

gfx define faces egroup Region1D
gfx create window 1
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
gfx modify spectrum default clear overwrite_colour
gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1
gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below banded number_of_bands 20 band_ratio 0.06 component 1
gfx modify g_element Region1D node_points glyph sphere general size "0.05*0.05*0.05" centre 0,0,0 
gfx modify g_element Region1D lines select_on material default selected_material default_selected
gfx modify g_element Region1D surfaces select_on material default data Phi spectrum default selected_material default_selected render_shaded

gfx draw axes
gfx edit scene
