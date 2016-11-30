

for ($i = 0 ; $i < 1000 ; $i=$i+10) #
  {
    gfx read node MainTime_1_$i.part0.exnode time $i
    gfx read node MainTime_M_2_$i.part0.exnode time $i node_offset 10000
  }

gfx read elem EMGExample_FE.part0.exelem generate
gfx read elem EMGExample_M.part0.exelem generate element_offset 10000 node_offset 10000 line_offset 10000 face_offset 10000

# define fields: pressure, deformed
gfx define field deformed_coordinates coordinate_system rectangular_cartesian composite DependentFE.1 DependentFE.2 DependentFE.3;
gfx define field pressure coordinate_system rectangular_cartesian composite DependentFE.4;
gfx define field velo coordinate_system rectangular_cartesian composite contraction_velocity.3;
gfx create spectrum VM;
gfx modify spectrum VM clear overwrite_colour;
gfx modify spectrum VM linear reverse range 0 1 rainbow colour_range 0 1 component 1;
gfx create spectrum ac;
gfx modify spectrum ac clear overwrite_colour;
gfx modify spectrum ac linear reverse range 1.6 3.7 rainbow colour_range 0 1 component 1;
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0 1 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx create spectrum shl;
gfx modify spectrum shl clear overwrite_colour;
gfx modify spectrum shl linear reverse range 1.2 1.21231 rainbow colour_range 0 1 component 1;
gfx create spectrum velo;
gfx modify spectrum velo clear overwrite_colour;
gfx modify spectrum velo linear reverse range -1 0.5 rainbow colour_range 0 1 component 1;
gfx create material black normal_mode ambient 0 0 0 diffuse 0 0 0 emission 0 0 0 specular 0.3 0.3 0.3 alpha 1 shininess 0.2;
gfx create material blue normal_mode ambient 0 0 0.5 diffuse 0 0 1 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material bone normal_mode ambient 0.7 0.7 0.6 diffuse 0.9 0.9 0.7 emission 0 0 0 specular 0.1 0.1 0.1 alpha 1 shininess 0.2;
gfx create material default normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material default_selected normal_mode ambient 1 0.2 0 diffuse 1 0.2 0 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material gold normal_mode ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 1 shininess 0.2;
gfx create material green normal_mode ambient 0 0.5 0 diffuse 0 1 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.1;
gfx create material muscle normal_mode ambient 0.4 0.14 0.11 diffuse 0.5 0.12 0.1 emission 0 0 0 specular 0.3 0.5 0.5 alpha 1 shininess 0.2;
gfx create material red normal_mode ambient 0.5 0 0 diffuse 1 0 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material silver normal_mode ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material tissue normal_mode ambient 0.9 0.7 0.5 diffuse 0.9 0.7 0.5 emission 0 0 0 specular 0.2 0.2 0.3 alpha 1 shininess 0.2;
gfx create material transparent_gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 0 shininess 0.2;
gfx create material white normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx modify g_element Region3D general clear circle_discretization 6 default_coordinate Geometry element_discretization "4*4*4" native_discretization none;
gfx modify g_element Region3D lines line_width 2 select_on material black selected_material default_selected;
gfx modify g_element Region3D lines coordinate deformed_coordinates select_on material gold selected_material default_selected;
gfx modify g_element Region1D general clear circle_discretization 6 default_coordinate GeometryM element_discretization "4*4*4" native_discretization none;
gfx modify g_element Region1D lines coordinate GeometryM3D select_on invisible material default selected_material default_selected;
gfx modify g_element Region1D node_points coordinate GeometryM3D glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default data Vm spectrum VM selected_material default;
gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes -x y eye_spacing 0.25 width 1222 height 1005;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 9.66769 13.5977 -8.1583 interest_point 3.6 1.5 0.75 up_vector -0.206884 -0.51072 -0.834485 view_angle 35.1398 near_clipping_plane 0.162027 far_clipping_plane 57.9031 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

gfx edit scene


