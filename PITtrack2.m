function PITtrack2(site,C)
% program for analysing repeated surveys of RFID (or perhaps magnetic) coarse particle tracers
% written in support of MacVicar and Papangelakis (submission to ESPL) and
% Papangelakis, MacVicar, Montakhab, and Ashmore (submission to GRL)

%% get needed information about site from PITtrack_default
def_data = PITtrack_default(site); % load default parameters
% simplify how default params are written
D50 = def_data.D50;
szcat = def_data.szcat;
fdir = def_data.fdir;

%% movement data assumption
def_data.C = C;% assumption for how to analyse movement data
% C = 1 - movement calculated based on fimem (intersect particles found on both days)
% only
% C = 2 - movement stats include inferred (missing particles that did not move in a later survey)
% C = 3 - movement stats include inferred and likely tracers 
% likelymov (missing particles that moved in a later survey, but one survey period had a flood for which the mean travel distance of all found particles was more than twice that of any of the other flood in the uncertain period) and
% likelyunmov (missing particles that moved in a later survey and for which the current survey period was not the dominant flood)
% C = 4 - movement stats include inferred, likely, and lost
% lostmov (never found again but had previously moved) and
% lostunmov (never found again but had never moved from the starting position up until the point that they were lost)
% note that this control parameter C is used in AnalyzeTravel.m


%% get taglist if included in def_data
% taglist specifies the tag numbers to use in the analysis.  All others are
% ignored
if isfield(def_data,'taglist')
    taglistfname = [def_data.fdir,def_data.taglist];
    taglist = ConvCSV2Struct(taglistfname,1);
else
    taglist.tagnum = [];
end
% save list of tags to process to labData
taglist = taglist.tagnum;


%% get tag characteristics (i.e. the labData) file
labData = ConvCSV2Struct([fdir,def_data.fchar],1);
if isfield(labData,'placeddate')
    labData.placeddatenum = datenum(labData.placeddate,'yyyymmdd');
else
    nttot = length(labData.tagnum);
    labData.placeddatenum = zeros(nttot,1);
end

%% get tagged files
fnames = dir([fdir,'*tagged*.csv']);
nftot = length(fnames);

%% Place tags
for nf = 1:nftot
    placementfile = fnames(nf).name;
    % create placedData structure
    if nf == 1
        seeding = 1;
    else
        seeding = 0;
    end
     
    placedData(nf) = placeTags(fdir,placementfile,labData,seeding);
end

% Add Sdist if not included in the tagged.csv data
if ~isfield(placedData,'Sdist')
    thalname = def_data.thalname;
    for nf = 1:nftot
        % calculate Sdist and Ndist
        [Sdist,Ndist] = CalcSdist([fdir,thalname],placedData(nf).Easting,placedData(nf).Northing);
        placedData(nf).Sdist = Sdist;
        placedData(nf).Ndist = Ndist;
    end
end

% use taglist to cut labData to only the requested particles
if ~isempty(taglist)
    % shorten the labdata to only the requested tags
    [~,labidx] = ismember(taglist,labData.tagnum);
    labDatanames = fieldnames(labData);
    labDatai = labData;
    %initialize loop to list substructure fieldnames (e.g. Vel.x)
    nntot = length(labDatanames);
    for nn = 1:nntot
        eval(['labData.',labDatanames{nn},' = labDatai.',labDatanames{nn},'(labidx);']);
    end
    
end

%% calculate AllSurvey matrices
% sort by date
[~,idx] = sort([placedData.date]);
placedData = placedData(idx);
% get AllSurvey matrices (easier for calculations than the structure file format of placedData)
AllSurvey = CalcAllSurvey(labData,placedData);

%% Analyze Mobility and travel distances
% nt are the total number of tracers found moved, etc
% ni is the number of tracers found moved, etc per size class
% pt is the mobility
% pi is the mobility of each size class
% mt is the total travel length data
% mi is the travel length data per size classe
% plumet is the descriptors of the total plume variance and skewness data 
% plumei is the descriptors of the  plume variance and skewness data per size class
% mD50 is the travel length data of the median particle size
% clas are the classification matrices

% for each survey period
for nf = 2:nftot
    % event based results
     [nt_event(nf-1),ni_event{nf-1},pt_event(nf-1),pi_event{nf-1},mt_event(nf-1),mi_event{nf-1},plumet_event(nf-1),plumei_event{nf-1},mD50_event(nf-1),clas_event(nf-1)] = AnalyzeTravel(AllSurvey,[placedData.date],labData.db_mm,labData.placeddatenum,nf-1,nf,D50,szcat,def_data.C);
    % cumulative results since initial placement
     [nt_cumul(nf-1),ni_cumul{nf-1},pt_cumul(nf-1),pi_cumul{nf-1},mt_cumul(nf-1),mi_cumul{nf-1},plumet_cumul(nf-1),plumei_cumul{nf-1},mD50_cumul(nf-1),clas_cumul(nf-1)] = AnalyzeTravel(AllSurvey,[placedData.date],labData.db_mm,labData.placeddatenum,1,nf,D50,szcat,def_data.C);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output

%% AllSurvey Data
% outputs series of matrices with rows by tagnum and columns by date to
% show 'found' classification and position data (Easting, Northing, Sdist, and Ndist
if def_data.writeAllSurvey
    foundmem = AllSurvey.found(:,1);
    % found table
    Tfound = makePITTrackDataTable(AllSurvey.found(foundmem,:),labData.tagnum(foundmem),labData.db_mm(foundmem),[placedData.date]);
    writetable(Tfound,[fdir,'AllSurvey_found.csv'],'Delimiter',',');

    TEasting = makePITTrackDataTable(AllSurvey.Easting(foundmem,:),labData.tagnum(foundmem),labData.db_mm(foundmem),[placedData.date]);
    writetable(TEasting,[fdir,'AllSurvey_Easting.csv'],'Delimiter',',');

    TNorthing = makePITTrackDataTable(AllSurvey.Northing(foundmem,:),labData.tagnum(foundmem),labData.db_mm(foundmem),[placedData.date]);
    writetable(TNorthing,[fdir,'AllSurvey_Northing.csv'],'Delimiter',',');

    TSdist = makePITTrackDataTable(AllSurvey.Sdist(foundmem,:),labData.tagnum(foundmem),labData.db_mm(foundmem),[placedData.date]);
    writetable(TSdist,[fdir,'AllSurvey_Sdist.csv'],'Delimiter',',');

end

%% Travel Length Data
%save([fdir,'lengthclassed.mat'],'class_event','class_cumul')
% travel length data by event
Travel_Length_event = makePITTrackDataTable([clas_event.travel_length_fix],labData.tagnum,labData.db_mm,[placedData(2:end).date]);
writetable(Travel_Length_event,[fdir,'Travel_distance_event_C',num2str(def_data.C),'.csv'],'Delimiter',',');
% travel length data cumulative
Travel_Length_cumul = makePITTrackDataTable([clas_cumul.travel_length_fix],labData.tagnum,labData.db_mm,[placedData(2:end).date]);
writetable(Travel_Length_cumul,[fdir,'Travel_distance_cumul_C',num2str(def_data.C),'.csv'],'Delimiter',',');

%% Summary data for event-based analysis 

% headers
% date table data
date.start = datestr([placedData(1:nftot-1).date]);
date.end = datestr([placedData(2:nftot).date]);
dtt = struct2table(date);
% size table data
sizeD.small = szcat(1:end-1)';
sizeD.big = szcat(2:end)';
stt = struct2table(sizeD);

% totals for tracer size classes grouped together all dates
ptt = struct2table(pt_event);
ntt = struct2table(nt_event);
mtt = struct2table(mt_event);
plumet = struct2table(plumet_event);
tout = [dtt ntt ptt mtt plumet];
writetable(tout,[fdir,'Summary_AllSizes_event_C',num2str(def_data.C),'.csv'],'Delimiter',',');

% totals for tracer size classes one file per date, all size classes 
% for nf = 1:nftot-1
%     fname_i = [datestr(placedData(nf).date,29),'to',datestr(placedData(nf+1).date,29)];
%     iout = [stt struct2table(ni_event{nf}) struct2table(pi_event{nf}) struct2table(mi_event{nf}) struct2table(plumei_event{nf})];
%     writetable(iout,[fdir,'Output_event_C',num2str(def_data.C),fname_i,'.csv'],'Delimiter',',');
%     
% end

% totals for tracer size classes one file per size class all dates
if length(szcat)>2
    for nDb = 1:length(szcat)-1
        fname_d = ['Summary_',num2str(round(szcat(nDb)),'%d'),' - ',num2str(round(szcat(nDb+1)),'%d'),' mm_event_C',num2str(def_data.C)];
        for nf = 1:nftot-1
            iout = [struct2table(ni_event{nf}) struct2table(pi_event{nf}) struct2table(mi_event{nf}) struct2table(plumei_event{nf})];
            if nf == 1
                dout = iout(nDb,:);
            else
                dout = [dout;iout(nDb,:)];
            end
        end
        writetable([dtt dout],[fdir,fname_d,'.csv'],'Delimiter',',');
    end
end

%% Summary data for cumulative analysis 

% totals for tracer size classes grouped together all dates
ptt = struct2table(pt_cumul);
ntt = struct2table(nt_cumul);
mtt = struct2table(mt_cumul);
plumet = struct2table(plumet_cumul);
tout = [dtt ntt ptt mtt plumet];
writetable(tout,[fdir,'Summary_AllSizes_cumul_C',num2str(def_data.C),'','.csv'],'Delimiter',',');

% totals for tracer size classes one file per date, all size classes 
% for nf = 1:nftot-1
%     fname_i = [datestr(placedData(nf).date,29),'to',datestr(placedData(nf+1).date,29)];
%     iout = [stt struct2table(ni_cumul{nf}) struct2table(pi_cumul{nf}) struct2table(mi_cumul{nf}) struct2table(plumei_cumul{nf})];
%     writetable(iout,[fdir,'Output_C',num2str(def_data.C),'_cumul_',fname_i,'.csv'],'Delimiter',',');
%     
% end

% totals for tracer size classes one file per size class all dates
if length(szcat)>2
    for nDb = 1:length(szcat)-1
        fname_d = ['Summary_',num2str(round(szcat(nDb)),'%d'),' - ',num2str(round(szcat(nDb+1)),'%d'),' mm_cumul_C',num2str(def_data.C)];
        for nf = 1:nftot-1
            iout = [struct2table(ni_cumul{nf}) struct2table(pi_cumul{nf}) struct2table(mi_cumul{nf}) struct2table(plumei_cumul{nf})];
            if nf == 1
                dout = iout(nDb,:);
            else
                dout = [dout;iout(nDb,:)];
            end
        end
        writetable([dtt dout],[fdir,fname_d,'.csv'],'Delimiter',',');
    end
end

function T = makePITTrackDataTable(data,tagnum,db_mm,date)
% subprogram to create a table with the date, tagnum, and particle size information included

nftot = length(date);
datee = datestr(date(1),'mmm_dd_yyyy');
eval([datee,' = data(:,1);']);
evalstr = ['tagnum,db_mm,',datee];
for n = 2:nftot
    datee = datestr(date(n),'mmm_dd_yyyy');
    eval([datee,' = data(:,',num2str(n),');']);
    evalstr = [evalstr,',',datee];
end
eval(['T = table(',evalstr,');']);