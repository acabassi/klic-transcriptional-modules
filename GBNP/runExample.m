clear
close all

dataType = 'Multinomial';

for i=1:200
    fileName = strcat('data/HarbisonDataSubsample', int2str(i));

    uniqueIdentifier = 1;
    nIterations      = 20;
    drawFigures      = true;

    hyperParameterSamplingFrequency = 1;
    verbose       = false;
    initialise    = true;
    gammaPrior    = [2 4];
    thinningFreq  = 5;

    fHandle = @ApproxDPCluster_Multinomial;

    fullFileName = strcat(fileName, '.csv');
    clustering = feval(fHandle, fullFileName, uniqueIdentifier, nIterations, dataType, ...
        drawFigures, hyperParameterSamplingFrequency, verbose, initialise, ...
        gammaPrior, thinningFreq, 10, 0);

    data = readtable(fullFileName);
    clusters = splitvars(table(data(:,1), array2table(clustering.'), 'VariableNames', {'Gene','Cluster'}));
    writetable(clusters, strcat(fileName, '_clusters.csv'), 'WriteVariableNames', 0)
end