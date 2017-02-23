#Example_a1: Graphical element groups: viewing a cube
#
# In cmgui, every element group has a graphical rendition - a list of
# attributes such as lines, surfaces etc. that describes how it is to look.
# Combined, this 'graphical element group' monitors changes in its nodes and
# elements and automatically updates the graphics to reflect the changes. By
# default, every group is drawn using lines, ie. as a wireframe mesh, and
# later examples will show how other graphical representations can easily be
# created. In this example, a one-element cube is read in and displayed. It
# also shows how the graphics are updated in reponse to node changes, and
# shows some other useful graphics commands on the way.
#
# This example .com file will be loaded in a comfile window with 'All',
# 'Selected' and 'Close' buttons at the bottom. Lines beginning with a hash (#)
# are comments which are ignored. All other lines are commands, and all
# graphics commands begin with 'gfx'. Note that command tokens can be
# abbreviated, eg. 'cre' instead of 'create' - as long as the short form is
# not also short for any other token. Do not abbreviate commands too much
# otherwise they are not readable to other users or by scripts used to update
# comfiles when new commands with similar names are added.
#
# You should not run all commands at once since it is important to
# understand each step. Instead, perform the commands one at a time (in the
# order they appear) by double-clicking them. Alternatively, select a group of
# commands in the comfile window with the mouse and press the 'Selected'
# button. Between the commands are comments describing them, or instructions
# for carrying out interactive input, such as dragging the mouse.
#
#----------
#
# Read in the nodes and elements representing the cube.
gfx create region cube;
# Cmgui 2.8 users should replace the following with cube_0.4.xml in FieldML 0.4 format
gfx read region $example/cube.xml region cube;
gfx define faces egroup cube;
#
# Create a visualisation of the lines = edges of the cube
gfx modify g_element cube lines coordinate coordinates material default;
#
# Open a 3-D graphics window (named 1). You can also do this by selecting
#'3-D window' from the Tools menu on the Command Window.
gfx create window 1
#
# Create blue axes of 1.2 units length and draw them in the scene.
# The draw command will create a static graphic in the root region rendition.
gfx create material blue ambient 0.2 0.2 0.9 diffuse 0.2 0.2 0.9
gfx modify g_element "/" point  glyph axes general size "1.2*1.2*1.2" centre 0,0,0 font default select_on material blue selected_material default_selected;
#
# Now spend some time getting used to the mouse actions that change your view
# of the object. Press the left mouse button in the graphics window and drag it
# around to 'tumble' the object in 3-D. The middle and right mouse buttons can
# be used in the same way to translate and zoom the object. Try these, with and
# without the 'Perspective' button checked. You can always press the 'View_all'
# button to return to a view of the whole object, or alternatively type:
gfx modify window 1 image view_all
#
# To demonstrate moving nodes, first show node numbers. This involves adding
# node_points labelled with their cmiss_number to graphical element 'cube'.
# We will also give them their own colour.
gfx create material orange ambient 1 0.25 0 diffuse 1 0.25 0
gfx modify g_element cube node_points coordinate coordinates material orange label cmiss_number;
#
# Now bring up the node viewer so we can change view and change nodal
# coordinates.
#gfx create node_viewer
#
# To change the position of a node, type the node number in at the top of the
# dialog, press ENTER, choose the "coordinates" field, change the x, y, or z
# values, then click 'Apply'. Try changing the x coordinate of node 2 to 0.5.
# The cube will now be distorted in the graphics window.
#
# Now add blue surfaces to the distorted cube.
gfx create material bluey ambient 0 0.25 0.5 diffuse 0 0.4 1 specular 0.5 0.5 0.5 shininess 0.3;
gfx modify g_element cube surfaces coordinate coordinates mat bluey;
#
# You can also change the position of nodes by reading in nodes from a file.
# Since we have distorted the cube from its original position, just read in the
# original node file to demonstrate this effect.
if ($TESTING) {
	gfx read region $example/cube.xml region cube;
}
#
# If you want to continue entering commands, the following sets the command
# prompt to 'gfx ', saving you from typing it all the time. Most commands allow
# this capability.
gfx
#
#----------
#
# TIPS
#
# By studying the commands used in this example you will hopefully start to see
# some pattern to how they work. While many interactive editors are available
# to control parts of the program, text commands are also available to control
# most features so that command files such as this may be written. Many
# features are only available through text commands at this time.
#
# There is a very simple way to find out what commands you can use and what
# parameters they can take: using the ? and ?? tokens. If you enter 'gfx ?' on
# the command line of the command window, cmgui will list all the keywords you
# can possibly enter after gfx. Enter 'gfx create ?' for a further list of
# keywords. Now you can see where the 'gfx create axes' command used above
# fits in. Now enter 'gfx create axes ?'. The possible options are now
# placed in angled brackets indicating that you can supply any of the
# parameters. All of the options have default values, given in square brackets.
# You can put 'length 1.2' after the command stem, which is equivalent to
# 'lengths 1.2*1.2*1.2'. The 'material' option allows the colour of the axes
# to be controlled; you can list available materials using the command
# 'gfx list material'.
#
# The ?? help mode works much like ? except that it is recursive - it will list
# all the keywords that may follow the current one, and all that may follow
# them, and so on. Type gfx ?? to see all the graphics commands and parameters
# available - this will take a few seconds.
#
# Export the cube region in fieldml format.
if ($TESTING) {
	gfx write region cube_export.fieldml region cube;
}

