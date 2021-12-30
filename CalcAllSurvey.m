function AllSurvey = CalcAllSurvey(labData,placedData)
% creates a set of matrices to store data for the desired tracer stones
nttot = length(labData.tagnum);
ndtot = length(placedData);

% NaN matrices for data read from files.  NaN is necessary so that
% positions are not read as 0 m
Easting = zeros(nttot,ndtot)/0; % storage for Easting;
Northing = zeros(nttot,ndtot)/0; % storage for Northing;
Tprecision = zeros(nttot,ndtot)/0; % storage for precision;
Sdist = zeros(nttot,ndtot)/0; % storage for centerline distance (was called ldist in previous code, updated to S to match Bradley 2012);
Ndist = zeros(nttot,ndtot)/0; % storage for lateral distance (was called hdist in previous coded);

%% determine whether particles were found on the survey dates
% found = 1 if the tracer was detected and located on the survey date
for nd = 1:ndtot
    % if there is no taglist, check if the tracerData tagnum are members of the labData tagnum
    % (done all at once per date using vector notation)
    [found(:,nd),idx] = ismember(labData.tagnum,placedData(nd).tagnum);
      
    fidx = find(idx); % find the row numbers from the index
    % save data from input file into the right locations (missing data will
    % still be NaN)
    Easting(fidx,nd) = placedData(nd).Easting(idx(fidx));
    Northing(fidx,nd) = placedData(nd).Northing(idx(fidx));
    Tprecision(fidx,nd) = placedData(nd).precision(idx(fidx));
%    Elevation(fidx,nd) = placedData(nd).elevation(idx(fidx));
    Ndist(fidx,nd) = placedData(nd).Ndist(idx(fidx));
    Sdist(fidx,nd) = placedData(nd).Sdist(idx(fidx));
end
% calculate mean travel lengths for found particles
meanL = zeros(1,ndtot);
for nd = 2:ndtot
    fimem = found(:,nd-1) & found(:,nd); %
    meanL(nd)=mean((Sdist(fimem,nd)-Sdist(fimem,nd-1)).^2+(Ndist(fimem,nd)-Ndist(fimem,nd-1)).^2).^0.5;
end

AllSurvey.found = found;
AllSurvey.Easting = Easting;
AllSurvey.Northing = Northing;
AllSurvey.Tprecision = Tprecision;
%AllSurvey.Elevation = Elevation;
AllSurvey.Ndist = Ndist;
AllSurvey.Sdist = Sdist;
AllSurvey.meanL = meanL;
