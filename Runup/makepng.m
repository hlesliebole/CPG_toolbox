function makepng(pfile);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','-r300','-loose',pfile);
fprintf(1,'Output image written to: %s \n ', pfile);  