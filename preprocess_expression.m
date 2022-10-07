%% Gene network data pre-processing script
% takes in gene expression data (single-cell RNA-seq)
% then pre-process via quantile normalization, log2 transform (optional based on user input)
% saves 3 data sets:
% list of sample names
% list of gene names
% matrix of gene expression (row=genes, cols=samples)

% Fair warning: because datasets come in different formats, user input may
% be required to correctly format the data for further downstream analysis.
% End goal is matrix of normalized gene expression and vectors of sample/gene names.

clear all; clc; close all
fprintf('Gene expression pre-processing script\nBy Kevin Murgas, Stony Brook University Dept. of Biomedical Informatics\n')

%% read in expression data

% basedir contains folder "Data" with raw expression data
basedir = pwd;

expNum = 75748; % enter GSE experiment ID to load
lookFor = strcat("GSE",num2str(expNum),"*"); % define string at start of data file names
filestruct = dir(fullfile(basedir,"Data",lookFor));
%lookFor = "GDSC_RNAexp.csv";
%filestruct = dir(fullfile(basedir,lookFor));
filestruct([filestruct.isdir])=[];
nFile=length(filestruct);

if ~nFile
    fprintf('\nNo files found! (for GSE %i)\n',expNum)
    return
else
    fprintf('\nFiles found (%i):\n',nFile)
    for i=1:nFile
        fprintf('%i: %s\n',i,filestruct(i).name)
    end
end
loadfiles=input("Please enter # of file data to load:\n");

fprintf('\nBegin loading experiment %i...\n',expNum)
temp={};
for i=loadfiles%1:nFile
    fprintf('Loading file %i/%i - %s...',i,nFile,['Data/' filestruct(i).name])
    % use tdfread(filename,delimiter) for tab delimited file
    temp{i} = readtable(fullfile(filestruct(i).folder,filestruct(i).name));
    
    % check if sample names all default, if so, reload
    if length(cell2mat(regexp(string(temp{i}.Properties.VariableNames),"Var")))==size(temp{i},2)
        fprintf('all variables default name, reloading headers...')
        fID = fopen(fullfile(filestruct(i).folder, filestruct(i).name));
        head = string(strsplit(fgetl(fID)));
        fclose(fID);
        temp{i}.Properties.VariableNames = head;
    end
    fprintf('done\n')
end

% pick one experiment loaded and begin parsing
% convert experiment from table to matrix + gene name vector
iExp = loadfiles;
expTable = temp{iExp};

% % 72056 Melanoma set only: construct sample names from first two rows, then skip first 3 rows
% sampName = strcat(string(expTable{1,2:end})','_',string(expTable{2,2:end})');
% expTable = expTable(4:end,:);

% get gene names from specified column
% user input, show first few columns and ask which is gene column
expTable(1:5,1:10)
iCol = input("Please enter # of column with gene names:\n");
geneName = string(table2array(expTable(:,iCol)));
geneForm = input("Please enter gene name form: 1=symbol, 2=ensemblID (ENSG000...), 3=EntrezID (7846)\n");
if geneForm==2 % ensembl only
    for i=1:length(geneName)
        dotind = regexp(geneName(i),'\.');
        if ~isempty(dotind)
            geneName{i}=geneName{i}(1:dotind-1);
        end
    end
end

% get expression data from remaining columns
iCol=input("Please enter # of leftmost column containing data:\n");

colName = string(expTable.Properties.VariableNames);
% convert cell 2 double
tableClass = string(varfun(@class,expTable,'OutputFormat','cell'));
if any(tableClass=="cell")
    fprintf("Converting data from cell to double...\n")
    ind = intersect(find(tableClass == "cell"),iCol:size(expTable,2));
    for i = 1:length(ind)
        expTable.(colName(ind(i))) = str2double(expTable.(colName(ind(i))));
    end
    fprintf('done\n')
end

% replace any missing values with 0
%TF = ismissing(expTable,{'' '.' 'NA' NaN -99});
TF = ismissing(expTable,{'NA' NaN});
if any(TF(:)) % if any missing values, replace with 0 and recast columns as double
    fprintf('Found missing values, replacing and recasting...\n')
    %[ri, ci] = ind2sub(size(expTable),find(TF));
    %expTable(ri,ci) = 0; % set to zero
    expTable(:,iCol:end) = fillmissing(expTable(:,iCol:end),'constant',{0});
end

fprintf('Extracting expression data:\n')
exprData = table2array(expTable(:,iCol:end));
exprData(1:5,1:6)
pause
doFilt = input("Filter low expressed genes? (0=no,1=yes)\n");
if doFilt
    fprintf('Filtering (>=10 total reads)...')
    ind = find(sum(exprData,2) >= 10); % minimum of 10 total reads
    exprData = exprData(ind,:);
    geneName = geneName(ind);
    fprintf('done\n')
end
doNorm = input("Do log2 normalization? (0=no,1=yes)\n");
if doNorm
    fprintf('Normalizing...')
    % normalize like this: quantile transform, with log2
    exprNorm = quantilenorm(exprData); % quantile normalization
    exprNorm(exprNorm<1)=1; % cut off values <1 so log2 is not negative
    exprNorm=log2(exprNorm+0.1); % and add 0.1 to all data for nonzero log
    fprintf('done\n')
else
    exprNorm=exprData;
end

% get sample names
sampName = colName(iCol:end);
nSamp = length(sampName);
fprintf('Data loaded and parsed\n')

% % For 72056 melanoma set, remove "unresolved" tumors (samp = XX_0)
% ind = extractAfter(sampName,3) ~= "0";
% sampName = sampName(ind);
% exprNorm = exprNorm(:,ind);

%% convert gene names to symbols (for matching with PPI gene names)
% for now just working with the HGNC symbols
% if necesary, will implement converting Ensembl and EntrezID to HGNC

%% save data
% save as geneName, exprNorm, sampName
filename = fullfile(basedir,"Data_proc",sprintf("GSE%i_%i_%s.mat",expNum,iExp,"fql2"));
%filename = fullfile(basedir,"GDSC_RNAproc_fql2.mat");
save(filename,'geneName','exprNorm','sampName')

fprintf('Saved data at %s\n',filename)
