%
%  Main menu script for viewing and editing .rdb files
%


k=1;

% leave program if gridding or quit is chosen

while k < 10;
k=menu('Choose an option','Add a Data File','RePlot','Delete Points',...
'Add Points','Get Depths','Select an Area',...
'Delete Selected Area ','Grid Selected Area','Add Spot Location','Quit');

if    k == 1
 load_rbd
 close(gcf)
 make_fig
 plot_rbd
 zoom on
elseif k == 2
 close(gcf)
 make_fig
 plot_rbd
 zoom on
elseif k == 3
 del_rbd
 close(gcf)
 make_fig
 plot_rbd
 zoom on
elseif k == 4
 add_rbd
 close(gcf)
 make_fig
 plot_rbd
 zoom on
elseif k == 5
 dep_rbd  
 close(gcf)
 make_fig
 plot_rbd
 zoom on
elseif k == 6
 zoom off
 sel_grd 
 zoom on
elseif k == 7
 del_area
elseif k == 8
 run_grid
elseif k == 9
 add_spot
elseif k == 10
 close(gcf)
end

end

write_rbd
clear all
