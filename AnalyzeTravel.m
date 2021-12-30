function [nt,ni,pt,pi,mt,mi,plumet,plumei,mD50,clas] = AnalyzeTravel(AllSurvey,dates,Db,placeddatenum,datestart,dateend,D50,Db_cat,C)
% perform the analysis on tracer movements

%% confirm that all data meets size and placeddate criteria
% get the relevant size data for categories related to Dbcat
[~,~,bin_Db] = histcounts(Db,Db_cat);% get the bins for each tag
% find all data that was categorized by size
goo_Db = bin_Db>0;
% check placed date
goo_date = placeddatenum<=dates(datestart);
% intersect size and placement date criteria
classed = goo_Db&goo_date; 

%% find positions for start and end dates
if datestart == 1
    cumul = 1;% control parameter used in FindPositions for determination of likely moved/unmoved - difference when cumulative vs event based assessement being done
else
    cumul = 0;
end
% find positions at start date
found1 = FindPositions(datestart,AllSurvey,classed,cumul);
% find positions at end date
found2 = FindPositions(dateend,AllSurvey,classed,cumul);
% classify found missing etc
clas = ClassifyTracers(found1,found2,classed,AllSurvey,dateend);

%% count classes for all particle sizes
nt.placed = sum(clas.placed); % in river
nt.foundearlier = sum(clas.foundearlier); % found in previous survey
nt.found = sum(clas.found); % found
nt.missing = sum(clas.missing); % missing

nt.intersect = sum(clas.intersect); % intersect
nt.refound = sum(clas.refound); % refound

nt.moved = sum(clas.moved); % intersect moved
nt.unmoved = sum(clas.unmoved); % intersect not moved

nt.lost = sum(clas.lost); % missing lost
nt.foundlater = sum(clas.foundlater); % missing foundlater
nt.inf = sum(clas.inf);
nt.lik = sum(clas.lik);

nt.lost_mov = sum(clas.lost_mov); % lost previously moved
nt.lost_unmov = sum(clas.lost_unmov); % lost previously moved

nt.inf_mov = sum(clas.inf_mov);
nt.inf_unmov = sum(clas.inf_unmov);
nt.lik_mov = sum(clas.lik_mov);
nt.lik_unmov = sum(clas.lik_unmov);

nt.fl_inf = sum(clas.fl_inf); % foundlater inferred
nt.fl_likmov = sum(clas.fl_likmov); % foundlater indeterminate
nt.fl_likunmov = sum(clas.fl_likunmov); % foundlater indeterminate
nt.fl_ind = sum(clas.fl_ind); % foundlater indeterminate

nt.rf_inf = sum(clas.rf_inf); % refound inferred
nt.rf_likmov = sum(clas.rf_likmov); % refound indeterminate
nt.rf_likunmov = sum(clas.rf_likunmov); % refound indeterminate
nt.rf_ind = sum(clas.rf_ind); % refound indeterminate

%% count classes by size class
ni.placed = histcounts(Db(clas.placed),Db_cat)'; % number per size class
ni.foundearlier = histcounts(Db(clas.foundearlier),Db_cat)'; % number found previous survey
ni.found = histcounts(Db(clas.found),Db_cat)'; % found per size class
ni.missing = histcounts(Db(clas.missing),Db_cat)'; % missing per size class

ni.intersect = histcounts(Db(clas.intersect),Db_cat)'; % found intersect per size class
ni.refound = histcounts(Db(clas.refound),Db_cat)'; % found refound per size class

ni.moved = histcounts(Db(clas.moved),Db_cat)'; % intersect moved per size class
ni.unmoved = histcounts(Db(clas.unmoved),Db_cat)'; % intersect unmoved per size class

ni.lost = histcounts(Db(clas.lost),Db_cat)';% missing lost per size class
ni.foundlater = histcounts(Db(clas.foundlater),Db_cat)';% missing foundlater per size class
ni.inf = histcounts(Db(clas.inf),Db_cat)';
ni.lik = histcounts(Db(clas.lik),Db_cat)';

ni.lost_mov = histcounts(Db(clas.lost_mov),Db_cat)'; % lost previously moved
ni.lost_unmov = histcounts(Db(clas.lost_unmov),Db_cat)'; % lost previously unmoved

ni.inf_mov = histcounts(Db(clas.inf_mov),Db_cat)';
ni.inf_unmov = histcounts(Db(clas.inf_unmov),Db_cat)';
ni.lik_mov = histcounts(Db(clas.lik_mov),Db_cat)';
ni.lik_unmov = histcounts(Db(clas.lik_unmov),Db_cat)';

ni.fl_inf = histcounts(Db(clas.fl_inf),Db_cat)'; % foundlater inferred per size class
ni.fl_likmov = histcounts(Db(clas.fl_likmov),Db_cat)'; % foundlater likely moved per size class
ni.fl_likunmov = histcounts(Db(clas.fl_likunmov),Db_cat)'; % foundlater likely unmoved per size class
ni.fl_ind = histcounts(Db(clas.fl_ind),Db_cat)'; % foundlater indeterminate per size class

ni.rf_inf = histcounts(Db(clas.rf_inf),Db_cat)'; % refound inferred per size class
ni.rf_likmov = histcounts(Db(clas.rf_likmov),Db_cat)'; % refound likely moved per size class
ni.rf_likunmov = histcounts(Db(clas.rf_likunmov),Db_cat)'; % refound likely unmoved per size class
ni.rf_ind = histcounts(Db(clas.rf_ind),Db_cat)'; % refound indeterminate per size class


%% set data to be analysed based on assumptions about indeterminate and lost data
% note that C is set in PITtrack_default.m. The value is noted in the output file name

% case 1 - based on fimem (particles found on both days)
if C == 1
    nTmov = nt.moved;
    nTtot = nt.intersect;
    nImov = ni.moved;
    nItot = ni.intersect;
    gooTLmem = clas.moved; % for travel length calculations
    gooLJmem = clas.intersect; % for plume (lajeunesse) calculations

elseif C == 2
    % add in found later inferred (refound are already included
    nTmov = nt.moved + nt.inf_mov;
    nTtot = nt.intersect + nt.inf;
    nImov = ni.moved + ni.inf_mov;
    nItot = ni.intersect + ni.inf;
    gooTLmem = clas.moved | clas.inf_mov;
    gooLJmem = clas.intersect | clas.inf;

elseif C == 3
    % add in likely
    nTmov = nt.moved + nt.inf_mov + nt.lik_mov;
    nTtot = nt.intersect + nt.inf + nt.lik;
    nImov = ni.moved + ni.inf_mov + ni.lik_mov;
    nItot = ni.intersect + ni.inf + ni.lik;
    gooTLmem = clas.moved | clas.inf_mov | clas.lik_mov;
    gooLJmem = clas.intersect | clas.inf | clas.lik;

elseif C == 4
    % add in likely
    nTmov = nt.moved + nt.inf_mov + nt.lik_mov + nt.lost_mov;
    nTtot = nt.intersect + nt.inf + nt.lik + nt.lost;
    nImov = ni.moved + ni.inf_mov + ni.lik_mov + ni.lost_mov;
    nItot = ni.intersect + ni.inf + ni.lik + ni.lost;
    gooTLmem = clas.moved | clas.inf_mov | clas.lik_mov | clas.lost_mov;
    gooLJmem = clas.intersect | clas.inf | clas.lik | clas.lost;

end

travel_length_fix = zeros(length(clas.travel_length),1)/0; %make all NaN
travel_length_fix(gooLJmem) = 0; %change all that were used for plume statistics to 0
travel_length_fix(gooTLmem) = clas.travel_length(gooTLmem); %change all that were used for movement statistics to travel length value
clas.travel_length_fix = travel_length_fix;


%% all together
% Calculate mobility
pt = CalcMobility(nTmov,nTtot);

% Calc Travel Length statistics
mt = CalcTravelLength(clas.travel_length(gooTLmem),Db(gooTLmem));

% Calc Plume data
% parameters for Lajeunesse 2018 are based on the plume statistics rather than
% travel lengths
plumet = CalcPlume(found1,found2,gooLJmem);

%% by size
% Calculate mobility
pi = CalcMobility(nImov,nItot); % works with vector format, so all three sizes calculated at once

% loop for each size - Travel Length and Lajeunesse require different particles to be isolated
nDbtot = length(Db_cat)-1;
for nDb = 1:nDbtot
    % Calc Travel Length statistics
    gooS = bin_Db==nDb & gooTLmem;
    mi(nDb) = CalcTravelLength(clas.travel_length(gooS),Db(gooS));
    % Calc Lajeunesse data
    gooL = bin_Db==nDb & gooLJmem;
    plumei(nDb) = CalcPlume(found1,found2,gooL);
end

%% save info about the D50
nD50 = find(histcounts(D50,Db_cat));
mD50 = mi(nD50); % geometric mean

end



function p = CalcMobility(nmov,ntot)
% calculate mobility and upper and lower confidence intervals

zse = 1.96; % 95% confidence

% mobility
p.pmean = (nmov./ntot);
% calculate standard error
pse = (p.pmean.*(1-p.pmean)./ntot).^0.5;
% upper and lower bounds
p.pupper = p.pmean+zse*pse;
p.plower = p.pmean-zse*pse;

end

function m = CalcTravelLength(travel_length,Db)

if ~isempty(travel_length)
    mn = length(travel_length);
    m.mmean = mean(travel_length);
    mstd = std(travel_length);
    skew = skewness(travel_length);
    % calculate standard error
    pse = mstd/sqrt(mn);
    zse = 1.96; % 95% confidence
    m.mupper = m.mmean+zse*pse;
    m.mlower = m.mmean-zse*pse;
    % calculate geometric mean
    m.mgmean = 10^(mean(log10(travel_length)));
    gstd = std(log10(travel_length));
    gpse = gstd/sqrt(m.mn);
    m.mgupper = 10^(log10(m.mgmean)+zse*gpse);
    m.mglower = 10^(log10(m.mgmean)-zse*gpse);
    % calculate median
    m.mmedian = median(travel_length);
    n_lower = m.mn/2-zse*sqrt(m.mn/2);
    tlsort = sort(travel_length);
    if n_lower<1
        n_lower = 1;
    end
    m.mmlower = tlsort(round(n_lower));
    n_upper = m.mn/2+zse*sqrt(m.mn/2);
    if n_upper>m.mn
        n_upper = m.mn;
    end
    m.mmupper = tlsort(round(n_upper));

    % calculate mean size
    m.Dbmean = 2^(mean(log2(Db)));
    m.Dbmedian = median(Db);
else
    m.mn = 0;
    m.mmean = 0;
    m.mstd = 0;
    m.skew = 0;
    m.mupper = 0;
    m.mlower = 0;
    % calculate geometric mean
    m.mgmean = 0;
    m.mgupper = 0;
    m.mglower = 0;
    % calculate median
    m.mmedian = 0;
    m.mmlower = 0;
    m.mmupper = 0;
    % calculate mean size
    m.Dbmean = 0;
    m.Dbmedian = 0;
end
end

function plume = CalcPlume(found1,found2,goomem)

plume.Svardiff = var(found2.Sdistlast(goomem))-var(found1.Sdistlast(goomem));
plume.Sskew = skewness(found2.Sdistlast(goomem));

end