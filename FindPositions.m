function foundData = FindPositions(datenum,AllSurvey,goodate,cumul)
% determine found missing refound tags and infer positions

%% switch from structure to flat matrices
Easting = AllSurvey.Easting;
Northing = AllSurvey.Northing;
Tprecision = AllSurvey.Tprecision;
Ndist = AllSurvey.Ndist;
Sdist = AllSurvey.Sdist;
found = AllSurvey.found;
meanL = AllSurvey.meanL;

[nttot,ndtot] = size(Easting);

%% initialize last known position matrices based on the found particles (not
% found tracers will still be NaN)
Eastinglast = AllSurvey.Easting(:,datenum);
Northinglast = Northing(:,datenum);
Tprecisionlast = Tprecision(:,datenum);
Ndistlast = Ndist(:,datenum);
Sdistlast = Sdist(:,datenum);

%% create empty matrices to keep track of positions through time
%foundlater = zeros(nttot,1);
%missing = zeros(nttot,1); % missing = 1 if the tracer was not found but was found at a later date
%lost = zeros(nttot,1); % lost = 1 if the tracer was not found on this date or any date from this point
fl_inf = zeros(nttot,1); % inferred = 1 if the tracer was missing and its position could be inferred based on pre and post surveys that showed it had not moved
fl_likmov = zeros(nttot,1); % likelymov = 1 if the tracer was missing and it likely moved in this survey period
fl_likunmov = zeros(nttot,1); % likelyunmov = 1 if the tracer was missing and it likely did not move in this survey period
fl_ind = zeros(nttot,1); % indeterminite = 1 if the tracer was missing and its position could not be confirmed or a probability assigned
%removed = zeros(nttot,1); % removed = 1 after the date when the previously placed tracer was removed from the river ** not used yet but could be!

%% determine whether tracers were in or out on the day of interest
in = goodate & found(:,1) == 1; % placed = 1 after the date when the tracer was first placed in the river

%in = placed &~removed; % in = 1 if the tracer was in the river on the survey date
out = ~in; % out = 1 if the tracer was not in the river on the survey date

%% classify lost and missing particles
% missing is any particles in the river that were not found in the survey 
missing = in & ~found(:,datenum);
if datenum<ndtot
    lost = missing & ~any(found(:,datenum+1:end),2); % never found again
else
    lost = missing; % this is the last survey, so all particles that are missing are also considered lost (note that subsequent surveys can change classification of lost or missing for all dates)
end
foundlater = missing & ~lost;

%% subclassify foundlater particles and record last positions
if any(foundlater)
    % get the numbers of particles foundlater
    nflmem = find(foundlater); 
    nfltot = sum(foundlater);
           
    % for each foundlater particle
    for nfl = 1:nfltot
        % get index to missing data
        nfli = nflmem(nfl);
        % find previous position of foundlater particle
        previouspos = find(found(nfli,1:datenum-1),1,'last');
        % find future position of foundlater particle
        futurepos = datenum+find(found(nfli,datenum+1:end),1,'first');
            
        Sdistprevious = Sdist(nfli,previouspos);
        Ndistprevious = Ndist(nfli,previouspos);
        Tprecisionprevious = Tprecision(nfli,previouspos);
        
        % find future position of foundlater particle
        futurepos = datenum+find(found(nfli,datenum+1:end),1,'first');
        Sdistfuture = Sdist(nfli,futurepos);
        Ndistfuture = Ndist(nfli,futurepos);
        Tprecisionfuture = Tprecision(nfli,futurepos);

        % calculate distance of movement
        fl_length=((Sdistprevious- Sdistfuture)^2+(Ndistprevious- Ndistfuture)^2)^0.5;
        % find precision
        fl_precision = max([Tprecisionprevious Tprecisionfuture]);

        % if the travel length is less than the precision
        if fl_length<fl_precision
            % classify as fl_inf because it has not moved
            fl_inf(nfli) = 1;
            % use previous position as new position
            ndidx = previouspos;
%            ndidx = futurepos;
        % else see if there is one flood period in the missing period that is much bigger than the others 
        else
            % check if this is a cumulative assessment or not
            if ~cumul
                % sort the mean movments
                sortL = sort(meanL(previouspos+1:futurepos));
                meanLbefore = meanL(datenum);
            else
                % if it is a cumulative assessment, just isolate the max prior to and after the current date
                meanLbefore = max(meanL(previouspos+1:datenum));
                meanLafter = max(meanL(datenum:futurepos));
                sortL = sort([meanLbefore meanLafter]);
                
            end
            % check if there is a dominant period (1 period has mean travel length > 2x any of the others)
            chk = sortL(end)>2*sortL(end-1);
            if chk
                % check if the current period is the dominant period
                if meanLbefore == sortL(end)
                    % classify as fl_likelymov 
                    fl_likmov(nfli) = 1;
                    % use future position as the updated position
                    ndidx = futurepos;
                else
                    % the current period is not the dominant period, so particles are
                    % classify as fl_likelyunmove
                    fl_likunmov(nfli) = 1;
                    % use previous position as the updated position
                    ndidx = previouspos;
%                    ndidx = futurepos;
                end
            else
                % indeterminate position
                % classify as fl_ind
                fl_ind(nfli) = 1;
                % use previous position as the updated position
                ndidx = previouspos;
 %               ndidx = futurepos;
            end
        end
        % save future positions for all missing data
        Eastinglast(nfli) = Easting(nfli,ndidx);
        Northinglast(nfli) = Northing(nfli,ndidx);
        Sdistlast(nfli) = Sdist(nfli,ndidx);
        Ndistlast(nfli) = Ndist(nfli,ndidx);
        Tprecisionlast(nfli) = Tprecision(nfli,ndidx);
    end
end

%% record last positions of lost particles
if any(lost)
    lostmem = find(lost);
    nltot = length(lostmem);
    for nl = 1:nltot
        nli = lostmem(nl);
        % find previous positions of lost particles
        previouspos = find(found(nli,1:datenum-1),1,'last');
        if ~isempty(previouspos)
            % save all last known positions
            Eastinglast(nli) = Easting(nli,previouspos);
            Northinglast(nli) = Northing(nli,previouspos);
            Tprecisionlast(nli) = Tprecision(nli,previouspos);
            Sdistlast(nli) = Sdist(nli,previouspos);
            Ndistlast(nli) = Ndist(nli,previouspos);
        else
            
            disp(['tag ',num2str(Data.tagnum(nli)),' is lost but has no previous position data']);
        end
    end
end

%% save all data in foundData structure to return to main program
foundData.Eastinglast = Eastinglast;
foundData.Northinglast = Northinglast;
foundData.precisionlast = Tprecisionlast;
foundData.Sdistlast = Sdistlast;
foundData.Ndistlast = Ndistlast;
foundData.in = in;
foundData.out = out;
foundData.found = found(:,datenum);
foundData.missing = missing;
foundData.foundlater = foundlater;
foundData.lost = lost;
foundData.fl_inf = fl_inf;
foundData.fl_likmov = fl_likmov;
foundData.fl_likunmov = fl_likunmov;
foundData.fl_ind = fl_ind;

end
