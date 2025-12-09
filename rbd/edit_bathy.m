%
%  Main MATLAB script for viewing and editing .rdb files
%

% put subdirectory with the other mfiles in path

clear all

path(path,'./rbd_mfiles')

% get any default settings from defaults.m

rbd_defaults
close

% initialize stuff
ns=0;
nb=0;
nn=0;

format compact

% load first file

load_rbd
make_fig
plot_rbd
zoom

% use a main menu

edit_menu

close(gcf)
