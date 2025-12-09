%
%  Main menu script for viewing and editing .rdb files
%

%set(0, "DefaultFigurePosition", f.Position)
%k=1;

% leave program if gridding or quit is chosen

%while k < 10;
UIControl_FontSize_bak = get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 16);
k=menu('CHOOSE AN OPTION','Top View','Northward View','Eastward View','3D View','Delete Points','UNdelete Points','Delete Area',...
'UNdelete Area','View Point Z(t)','Grid Valid Points',...
'Save QC struct','Change  MOP','Exit Editor');
set(0, 'DefaultUIControlFontSize', UIControl_FontSize_bak);

if    k == 1
 dview='T';
 PlotSAQC
 edit_menuSA
elseif k == 2
 dview='N';
 PlotSAQC
 edit_menuSA
elseif k == 3
 dview='E';
 PlotSAQC
 edit_menuSA
elseif k == 4
 dview='3';
 PlotSAQC
 edit_menuSA
elseif k == 5
 DeletePointsSA
elseif k == 6
 UNDeletePointsSA
elseif k == 7
 SelectAreaSA
elseif k == 8
 % add_spot
 UNSelectAreaSA
elseif k == 9
 % select point and plot time series
 SelectZtimeHistorySA
elseif k == 10
 % grid valid points
 PlotScatterGridSA
 edit_menuSA
elseif k == 11
 save(QCmatfile,'QC')
 fprintf('QC struct array saved to: %s\n',QCmatfile)
 edit_menuSA
elseif k == 12
 get_newSA
elseif k == 13
 close all   
end

%end

%write_rbd
%clear all
