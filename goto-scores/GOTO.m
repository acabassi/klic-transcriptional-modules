function [GOTObp GOTOmf GOTOcc] = GOTO(fileName)
clusters = importdata([fileName '.csv']);
load('orf2SGD')
clusterIDs    = clusters.data;
orfNames      = clusters.textdata;
[orfNames,IX] = sort(orfNames);
clusterIDs    = clusterIDs(IX);
nGenes        = length(clusterIDs);
matchingInds  = find(ismember(orf2SGD_orfs, orfNames));

if (~isempty(find(isnan(matchingInds),1)))
    disp('Error - could not find all ORFs in the dictionary')
end

matchingSGDIDs = orf2SGD_SGDIDs(matchingInds);

% We now need to find these SGDIDs in the GOTO matrix:

% Load the list of row labels for the GOTO matrix 
load('GOTO_MatrixRowLabels')

matchingRows = find(ismember(rowLabels, matchingSGDIDs));
% Find the order in which they appear (i.e. the non-sorted version of the
% above)
index  = getnameidx(rowLabels(matchingRows), matchingSGDIDs);

if (~isempty(find(isnan(matchingRows),1)))
    disp('Error - could not find all SGDIDs among the GOTO matrix row labels')
end



fid=fopen('gene_associationsgdgoaP.csv');
disp('Calculating GOTO (bp)...')
mycell = cell(nGenes,1);
count = 1;
for i = 1:max(matchingRows)
    %disp(i)
    abb = fgetl(fid);
    if(ismember(i, matchingRows))
        mycell{count} = abb;
        count = count+1;
    end
end

mymat = zeros(nGenes,6382);
for i = 1:nGenes
    parts = strread(mycell{i},'%f','delimiter',',')';
    mymat(i,:) = parts;
end

mymat2     = mymat(:,matchingRows);
mymat3     = tril(mymat2,-1)' + mymat2; 


myMat = mymat3(index,index);


uniqueIDs  = unique(clusterIDs);
nClusters  = length(uniqueIDs); 
savedDist  = zeros(1,nClusters);
genesPerCluster = zeros(1,nClusters);
for i = 1:length(uniqueIDs)
    currentID   = uniqueIDs(i);
    inds        = (clusterIDs == currentID);
    genesPerCluster(i) = sum(inds);
    if(sum(inds)>1)
        submat      = myMat(inds,inds);
        currentDist = abs(sum(sum(tril(submat,-1))));%mean(mean(submat));
        savedDist(i)= currentDist;
    else
        savedDist(i) = NaN;
    end
end

%nanmean(savedDist.*genesPerCluster/551)
innersu = (2./(genesPerCluster.*(genesPerCluster-1))).*savedDist;

weightNorm = sum(genesPerCluster(~isnan(innersu)));
weights    = genesPerCluster/weightNorm;

GOTObp = nansum(weights.*innersu);
disp(GOTObp)
fclose(fid);


fid=fopen('gene_associationsgdgoaF.csv');
disp('Calculating GOTO (mf)...')
mycell = cell(nGenes,1);
count = 1;
for i = 1:max(matchingRows)
    %disp(i)
    abb = fgetl(fid);
    if(ismember(i, matchingRows))
        mycell{count} = abb;
        count = count+1;
    end
end

mymat = zeros(nGenes,6382);
for i = 1:nGenes
    parts = strread(mycell{i},'%f','delimiter',',')';
    mymat(i,:) = parts;
end

mymat2     = mymat(:,matchingRows);
mymat3     = tril(mymat2,-1)' + mymat2; 


myMat = mymat3(index,index);


uniqueIDs  = unique(clusterIDs);
nClusters  = length(uniqueIDs); 
savedDist  = zeros(1,nClusters);
genesPerCluster = zeros(1,nClusters);
for i = 1:length(uniqueIDs)
    currentID   = uniqueIDs(i);
    inds        = (clusterIDs == currentID);
    genesPerCluster(i) = sum(inds);
    if(sum(inds)>1)
        submat      = myMat(inds,inds);
        currentDist = abs(sum(sum(tril(submat,-1))));%mean(mean(submat));
        savedDist(i)= currentDist;
    else
        savedDist(i) = NaN;
    end
end

%nanmean(savedDist.*genesPerCluster/551)
innersu = (2./(genesPerCluster.*(genesPerCluster-1))).*savedDist;

weightNorm = sum(genesPerCluster(~isnan(innersu)));
weights    = genesPerCluster/weightNorm;
GOTOmf = nansum(weights.*innersu);
disp(GOTOmf)
fclose(fid);


fid=fopen('gene_associationsgdgoaC.csv');
disp('Calculating GOTO (cc)...')
mycell = cell(nGenes,1);
count = 1;
for i = 1:max(matchingRows)
    %disp(i)
    abb = fgetl(fid);
    if(ismember(i, matchingRows))
        mycell{count} = abb;
        count = count+1;
    end
end

mymat = zeros(nGenes,6382);
for i = 1:nGenes
    parts = strread(mycell{i},'%f','delimiter',',')';
    mymat(i,:) = parts;
end

mymat2     = mymat(:,matchingRows);
mymat3     = tril(mymat2,-1)' + mymat2; 


myMat = mymat3(index,index);


uniqueIDs  = unique(clusterIDs);
nClusters  = length(uniqueIDs); 
savedDist  = zeros(1,nClusters);
genesPerCluster = zeros(1,nClusters);
for i = 1:length(uniqueIDs)
    currentID   = uniqueIDs(i);
    inds        = (clusterIDs == currentID);
    genesPerCluster(i) = sum(inds);
    if(sum(inds)>1)
        submat      = myMat(inds,inds);
        currentDist = abs(sum(sum(tril(submat,-1))));%mean(mean(submat));
        savedDist(i)= currentDist;
    else
        savedDist(i) = NaN;
    end
end

%nanmean(savedDist.*genesPerCluster/551)
innersu = (2./(genesPerCluster.*(genesPerCluster-1))).*savedDist;

weightNorm = sum(genesPerCluster(~isnan(innersu)));
weights    = genesPerCluster/weightNorm;
GOTOcc = nansum(weights.*innersu);
disp(GOTOcc)
fclose(fid);

csvwrite([fileName '_GOTO.csv'], [GOTObp GOTOmf GOTOcc]);
end