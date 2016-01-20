% run full inverse model
clear all;

dataFolder = '../data/'; % folder to store output
dataSubfolder = '/repaired/'; % where output txt files are stored
probingRange = {'both'}; % use single cell, multicell, or both probing datasetss
targetRange = [1,3]; % targets to evaluate (1=left, 3=bottom)
typeRange = [0, 1]; % type of perturbation (0=cell, 1=synape)
percentRange = [5,10]; % percentages of perturbation

learnDict = 0;

% folders
bmmFolderOriginal =  '../data/original';
bmmFolderPerturb =  '../data/perturbed_probing';

for iprobing = 1:length(probingRange)
    for itarget = 1:length(targetRange)
        for itype = 1:length(typeRange)
            for ipercent = 1:length(percentRange)
                probing = probingRange{iprobing};
                target_num = targetRange(itarget);
                perturbation_type = typeRange(itype);
                perturbation_percentage = percentRange(ipercent);
                save('iterparams.mat');
                if learnDict
                    Initial_Dictionary % This picks a single stim pattern as a starting point and starts to build the model (Dictionary).
                    if strcmp(probing, 'single') % read all data to build full dictionary
                        Read_All_Stim_Files
                    elseif strcmp(probing, 'multi')
                        Read_All_Stim_Files_Multi
                    elseif strcmp(probing, 'both')
                        Read_All_Stim_Files_Multi
                        Read_All_Stim_Files
                    end
                end
                Optimize_Stim % generate stim pattern
            end
            clear all; % to restore initial conditions and get ready for next iteration
            load('iterparams.mat');
        end
    end
end