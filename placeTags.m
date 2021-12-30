function placedData = placeTags(placementdir,placementfile,labData,seeding)
% perform QC checks on survey data, intersect labData and placementfile, and assign precision

%% file
% get taggedCSV 
placementfname = [placementdir,placementfile];
placementRaw = ConvCSV2Struct(placementfname,1);

% extract date from filename
numberstr = regexp(placementfile,'(\d+)','tokens');
placedData.date = datenum(numberstr{1},'yyyymmdd');

placedData.placementfname = placementfname;
%placedData.labfname = labData.labfname;

%% find placed tags in the lab data and create indexes
[Lia,Locb] = ismember(labData.tagnum,placementRaw.tagnum);
foundIdx = find(Lia);
[~,sortIdx] = sort(Locb(foundIdx)); 

% check for errors
nraw = length(placementRaw.tagnum);
nconfirm = length(sortIdx);
if nraw ~= nconfirm
    
    % check for double tags
    tagsort = sort(placementRaw.tagnum);
    doublediff = diff(tagsort);
    if any(doublediff==0)
        ndd = find(doublediff==0);
        doubletag = tagsort(ndd);
        for nd = 1:length(doubletag)
            disp(['double tagnum detected for tag ',num2str(doubletag(nd))]);
        end
    end
    % check for wrong numbers
    [a,~] = ismember(placementRaw.tagnum,labData.tagnum);
    wrongtag = placementRaw.tagnum(~a);
    for nw = 1:length(wrongtag)
        disp(['tagnum ',num2str(wrongtag(nw)),' is not found in the tracers characteristics file.']);
    end
    error('Error in survey data.  Please fix and rerun PITtrack');
end

%% transfer labData data to placedData
labDatanames = fieldnames(labData);
% initialize loop to list substructure fieldnames (e.g. Vel.x)
nntot = length(labDatanames);
for nn = 1:nntot
    eval(['placedData.',labDatanames{nn},' = labData.',labDatanames{nn},'(foundIdx);']);
end

%% transfer placementRaw data to placedData
placementRawnames = fieldnames(placementRaw);
% initialize loop to list substructure fieldnames (e.g. Vel.x)
nn2tot = length(placementRawnames);
eval(['placedData.nttot = length(placementRaw.',placementRawnames{1},');']);

for nn2 = 1:nn2tot
    chk = ~strcmp(placementRawnames{nn2},'tagnum');
    chk2 = ~strcmp(placementRawnames{nn2},'Month');
    chk3 = ~strcmp(placementRawnames{nn2},'Year');
    chk4 = ~strcmp(placementRawnames{nn2},'precision');
    
    if chk && chk2 && chk3 && chk4
        eval(['placedData.',placementRawnames{nn2},' = zeros(placedData.nttot,1);']);
        eval(['placedData.',placementRawnames{nn2},'(sortIdx) = placementRaw.',placementRawnames{nn2},''';']);
    end
end

%% assign precision
TF = isfield(placementRaw,'precision');
if seeding
    % on the day of seeding the error is assumed to be 5 cm
    placedData.precision = ones(placedData.nttot,1)*0.05;
elseif ~TF
    % specify precision
    % All are assumed to be 1 m error
    placedData.precision = ones(placedData.nttot,1)*1;
else
    % specify precision for Toronto studies
    % 'X' means tag was seen visually (error of 5 cm, 20 means the stick
    % antenna was used.  All others are assumed to be 1 m error
    precstr = {'X','20cm'};
    placedData.precision = ones(placedData.nttot,1);
    % loop to check each value in 
    for nrawi = 1:nraw
        Yprec = strcmpi(placementRaw.precision(nrawi),precstr);
        if Yprec(1)
            placedData.precision(nrawi) = 0.50; % note increased to 0.50 to match stick antenna based on what looked like errors close to zero
        elseif Yprec(2)
            placedData.precision(nrawi) = 0.50;
        end                    
    end
end
end

