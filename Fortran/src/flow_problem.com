########reading nodes and elements########

gfx read nodes GENERICEXAMPLE NODEDATA_FLUID_FLOW.part0.exnode

gfx read elements GENERICEXAMPLE ELEMENTDATA_FLUID_FLOW.part0.exelem

########defining faces########

gfx define faces

########defining depedent fields##########

gfx define field vel_x component DEPENDENTFIELD.1  ;

gfx define field vel_y component DEPENDENTFIELD.2  ;

gfx define field pressure component  DEPENDENTFIELD.3  ;

######## opening window to visualize results##########

gfx cre win 1 ;

gfx cre mat copper ambient 1 0.2 0 diffuse 0.6 0.3 0 specular 0.7 0.7 0.5 shininess 0.3

########creating surfaces with dependent fields ##########

gfx mod g_e REGION surfaces data vel_x spectrum default

########creating nodes with labels##########
#point label cmiss_number
gfx modify g_element REGION node_points glyph point  select_on material default selected_material default_selected;

########creating creating a spectrum to display dependent fields ##########

gfx create spectrum vel_spectrum autorange

gfx mod g_e REGION surfaces data vel_x spectrum vel_spectrum

gfx modify spectrum vel_spectrum autorange

########displaying color bar on the visualization window. ##########

gfx create colour_bar spectrum vel_spectrum

gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on material copper selected_material copper normalised_window_fit_left;


gfx modify window 1 image view_all



#gfx define field E2 add fields Geometry.x Dependent.1 ;

