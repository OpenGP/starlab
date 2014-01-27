#--- This file compiles the full starlab. If you have compile problems, you 
#    should try starlab_mini.pro. Alternatively, you can comment out subdirs 
#    if something is not working correctly.
TEMPLATE = subdirs
CONFIG += ordered

#-------------------------------------------------------------------------------
#                               BASIC CHECKS
#-------------------------------------------------------------------------------
!exists(core):error(have you cloned starlab core in this folder?)
!exists(surfacemesh):error(have you cloned surfacemesh in this folder?)

#-------------------------------------------------------------------------------
#                               STARLAB CORE
#-------------------------------------------------------------------------------
SUBDIRS += core

#-------------------------------------------------------------------------------
#                           ADVANCED CORE PLUGINS
#-------------------------------------------------------------------------------
#SUBDIRS += core/plugins/gui_python
#SUBDIRS += plugins/project_io_starlab

#-------------------------------------------------------------------------------
#                               SURFACEMESH   
#-------------------------------------------------------------------------------
#--- Model Plugin
SUBDIRS += surfacemesh/surfacemesh #< DYNAMIC DATATYPE
#--- I/O Plugins
SUBDIRS += surfacemesh/surfacemesh_io_off
SUBDIRS += surfacemesh/surfacemesh_io_obj
#--- Rendering Plugins
SUBDIRS += surfacemesh/surfacemesh_render_verts
SUBDIRS += surfacemesh/surfacemesh_render_flat
SUBDIRS += surfacemesh/surfacemesh_render_smooth
SUBDIRS += surfacemesh/surfacemesh_render_wireframe
SUBDIRS += surfacemesh/surfacemesh_render_flatwire
SUBDIRS += surfacemesh/surfacemesh_render_transparent
#--- Filter Plugins [[ @TODO: FIX COMPILE ISSUES ]]
SUBDIRS += surfacemesh/surfacemesh_filter_normalize
SUBDIRS += surfacemesh/surfacemesh_filter_geoheat
#SUBDIRS += surfacemesh/surfacemesh_filter_ballpivoting
#SUBDIRS += surfacemesh/surfacemesh_filter_au_skeleton
SUBDIRS += surfacemesh/filter_depthscanner
#--- Mode Plugins
SUBDIRS += surfacemesh/surfacemesh_mode_info
SUBDIRS += surfacemesh/surfacemesh_mode_arapdeform
#--- Decorate Plugins [[ @TODO: INTERFACE NOT READY!! ]]
#SUBDIRS += surfacemesh/surfacemesh_decorate_normals
#SUBDIRS += surfacemesh/surfacemesh_decorate_selection

#-------------------------------------------------------------------------------
#                          EXAMPLES (CLOUD MODEL)
#-------------------------------------------------------------------------------
#SUBDIRS += core/example/cloud #< a *very* simple "cloud" model
#SUBDIRS += core/example/cloud_io_pts
#SUBDIRS += core/example/cloud_render_points
#SUBDIRS += core/example/cloud_filter_normalize
#SUBDIRS += core/example/cloud_mode_select
#SUBDIRS += core/example/cloud_decorate_selection

#-------------------------------------------------------------------------------
#                          EXAMPLES (TEMPLATES)
#-------------------------------------------------------------------------------
#SUBDIRS += core/example/example_mode_withwidget #< shows you how to create a mode plugin with widget
