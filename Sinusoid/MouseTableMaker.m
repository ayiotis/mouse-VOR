%% Mouse Table Maker
%This script calls the MouseAnalysisSummary fucntion to compile a table
%with all the data for each frequency for each mouse.

%Written by Andrianna Ayiotis
%Last updated 06/25/2018
%Last run on 08/02/2018 
%% Make the table
mice = {'WUm252','WUm255','WUm276','WUm279','WUm282','WUm283','WUm284','WUm295','WUm296','WUm297'}';
type = {'double het','cko','cko','double het','cko','cko','single het','double het','double het','cko'}';
allmousetab = table();
cd Figures
for i = 1:length(mice)
   mousetab = MouseAnalysisSummary(mice{i});
   savefig([mice{i},'.fig'])
   saveas(gcf,[mice{i},'.jpg'])
   mouse = [cell2table([repmat(mice(i),[size(mousetab,1),1]),repmat(type(i),[size(mousetab,1),1])]),mousetab];
   allmousetab = [allmousetab;mouse];
   close;
end
cd ../
lab = [{'Mouse','Type'},mousetab.Properties.VariableNames];
allmousetab.Properties.VariableNames = lab;
save('MouseSummary.mat','allmousetab')