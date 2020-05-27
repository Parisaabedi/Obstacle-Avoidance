%%% This file extracts all the trials with changed decision in the middle!
subs = 1:7; subs = [subs 9];
for i = 5
    sheet = strcat('s',num2str(subs(i)));
    num = readtable('ChangedDecisions.xlsx','Sheet',sheet);
    filename = strcat('CompF_',sheet,'.mat');
    load(filename)
    first = num(:,{'FL','F0','FR'});
    for j = 1 : 3
        Comp{j,1} = [Comp{j,1} zeros(size(Comp{j,1},1),1)];
        cds = table2array(first(:,j));
        cds(isnan(cds)) = [];
        if ~isempty(cds)
            Comp{j,1}(cds,9) = 1;
        end
    end
    save(filename,'Comp');
%     filename = strcat('CompNF_',sheet,'.mat');
%     load(filename)
%     first = num(:,{'NFL','NF0','NFR'});
%     for j = 1 : 3
%         Comp{j,1} = [Comp{j,1} zeros(size(Comp{j,1},1),1)];
%         cds = table2array(first(:,j));
%         cds(isnan(cds)) = [];
%         if ~isempty(cds)
%             Comp{j,1}(cds,9) = 1;
%         end
%     end
%     save(filename,'Comp');
end