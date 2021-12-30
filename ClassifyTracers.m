function clas = ClassifyTracers(found1,found2,classed,AllSurvey,dateend)
% classify tracers 

% subclass found intercept (fimem)  particles found on both the start and end dates
f1mem = found1.found & classed;
f2mem = found2.found & classed;

fimem = f1mem & f2mem; 
% subclass foundlater of particles that were found in 1, not in 2, and
% later refound
foundlatermem = found2.foundlater & classed;

% subclass found refound (refoundmem) particles that were missing on start date but found on end date
refoundmem = f2mem & ~f1mem; 

%% subclassify refound particles
nttot = length(found1.in);
rf_inf = zeros(nttot,1); % inferred = 1 if the tracer was refound and its travel length could be inferred based on pre and post surveys that showed it had not moved
rf_likmov = zeros(nttot,1); % likelymov = 1 if the tracer was refound and its travel length was assessed as likely in this survey period
rf_likunmov = zeros(nttot,1); % likelyunmov = 1 if the tracer was refound and its travel length was assessed as unlikely in this survey period
rf_ind = zeros(nttot,1); % indeterminite = 1 if the tracer was refound but its travel length could not be confirmed or a probability assigned
if any(refoundmem)
    % get the numbers of particles refound
    nrfmem = find(refoundmem); 
    nrftot = sum(refoundmem);
    % get meanL from AllSurvey
    meanL = AllSurvey.meanL;
           
    % for each refound particle
    for nrf = 1:nrftot
        % get index to missing data
        nrfi = nrfmem(nrf);
        
        % find previous position of refound particle
        previouspos = find(AllSurvey.found(nrfi,1:dateend-1),1,'last');
        Sdistprevious = AllSurvey.Sdist(nrfi,previouspos);
        Ndistprevious = AllSurvey.Ndist(nrfi,previouspos);
        Tprecisionprevious = AllSurvey.Tprecision(nrfi,previouspos);
        
        % find future position of refound particle
        futurepos = dateend; %particle was found
        Sdistfuture = AllSurvey.Sdist(nrfi,futurepos);
        Ndistfuture = AllSurvey.Ndist(nrfi,futurepos);
        Tprecisionfuture = AllSurvey.Tprecision(nrfi,futurepos);

        % calculate distance of movement
        rf_length=((Sdistprevious-Sdistfuture)^2+(Ndistprevious-Ndistfuture)^2)^0.5;
        % find precision
        rf_precision = max([Tprecisionprevious Tprecisionfuture]);

        % if the travel length is less than the precision
        if rf_length<rf_precision
            % classify as inferred because it has not moved
            rf_inf(nrfi) = 1;
        % else see if there is one flood period in the missing period that is much bigger than the others 
        else
            % sort the mean movments
            sortL = sort(meanL(previouspos+1:futurepos));
            % check if there is a dominant period (1 period has mean travel length > 2x any of the others)
            if sortL(end)>2*sortL(end-1)
                % check if the current period is the dominant period
                if meanL(dateend) == sortL(end)
                    % likely moved
                    rf_likmov(nrfi) = 1;
                else
                    % the current period is not the dominant period, so particles are
                    % likely unmoved
                    rf_likunmov(nrfi) = 1;
                end
            else
                % indeterminate position
                rf_ind(nrfi) = 1;
            end
        end
    end
end
%
% calculate event travel lengths
clas.travel_length = ((found2.Sdistlast-found1.Sdistlast).^2+(found2.Ndistlast-found1.Ndistlast).^2).^0.5;
precision = max([found1.precisionlast found2.precisionlast]')';

clas.placed = found1.in; % number per size class
clas.foundearlier = f1mem; % number found previous survey
clas.found = f2mem; % found per size class
clas.missing = found2.missing & classed; % missing per size class

clas.intersect = fimem; % found intersect per size class
clas.refound = refoundmem; % found refound per size class

clas.moved = clas.travel_length>=precision & clas.intersect; % intersect moved per size class
clas.unmoved = clas.travel_length<precision & clas.intersect; % intersect unmoved per size class

clas.lost = found2.lost & classed;% missing lost per size class
clas.foundlater = foundlatermem;% missing foundlater per size class

clas.fl_inf = found2.fl_inf & classed; % foundlater inferred per size class
clas.fl_likmov = found2.fl_likmov & classed; % foundlater likely moved per size class
clas.fl_likunmov = found2.fl_likunmov & classed; % foundlater likely unmoved per size class
clas.fl_ind = found2.fl_ind & classed; % foundlater indeterminate per size class

clas.rf_inf = rf_inf & classed; % refound inferred per size class
clas.rf_likmov = rf_likmov & classed; % refound likely moved per size class
clas.rf_likunmov = rf_likunmov & classed; % refound likely unmoved per size class
clas.rf_ind = rf_ind & classed;

clas.inf = clas.fl_inf|clas.rf_inf;
clas.lik = clas.fl_likmov|clas.fl_likunmov|clas.rf_likmov|clas.rf_likunmov;

clas.inf_mov = clas.travel_length>=precision & clas.inf;
clas.inf_unmov = clas.travel_length<precision & clas.inf;
clas.lik_mov = clas.travel_length>=precision & clas.lik;
clas.lik_unmov = clas.travel_length<precision & clas.lik;

clas.lost_mov = clas.travel_length>=precision & clas.lost; % lost previously moved
clas.lost_unmov = clas.travel_length<precision & clas.lost; % lost previously unmoved

end