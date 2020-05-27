%% save figures in the file
files = dir;

for j = 1 : length(files)
    k = strfind(files(j).name,'.fig');
     if ~isempty(k)
         gcf = openfig(files(j).name);
         set(gcf,'Visible','off')
         saveas(gcf,files(j).name(1:end-4),'tiff')
         close(gcf)
     end
end