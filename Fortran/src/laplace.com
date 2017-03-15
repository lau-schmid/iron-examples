########reading nodes and elements########

gfx read nodes GENERICEXAMPLE NODE_DATA_LAPLACE.part0.exnode

gfx read elements GENERICEXAMPLE ELEMENT_DATA_LAPLACE.part0.exelem

########defining faces########

gfx define faces

########defining depedent fields##########

gfx define field Temperature component DEPENDENT.1  ;


######## opening window to visualize results##########
gfx cre win 1 ;

gfx cre mat copper ambient 1 0.2 0 diffuse 0.6 0.3 0 specular 0.7 0.7 0.5 shininess 0.3

########creating nodes with labels##########
gfx modify g_element "/" general clear;
gfx modify g_element "/" node_points subgroup REGION coordinate GEOMETRY LOCAL glyph point general size "1*1*1" centre 0,0,0 font default label cmiss_number select_on material default selected_material default_selected;


########creating creating a spectrum to display temperature field ##########

gfx create spectrum temperature_spectrum autorange

########creating surfaces with Temperature field ##########

gfx modify g_element "/" surfaces coordinate GEOMETRY tessellation default LOCAL select_on material default data Temperature spectrum temperature_spectrum selected_material default_selected render_shaded;
gfx modify g_element "/" point coordinate GEOMETRY NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material copper selected_material copper;

gfx modify spectrum temperature_spectrum range 0 1


########displaying color bar on the visualization window. ##########

gfx create colour_bar spectrum temperature_spectrum

gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on material copper selected_material copper normalised_window_fit_left;

gfx modify window 1 image view_all

