function def_data = PITtrack_default(site)

%% write data control
def_data.writeAllSurvey = 1;
def_data.writeTravelLength = 1;

%% site specific data
% site 1 = Brad = 0; % Halfmoon Creek (Bradley 2017)
% site 2 = Ganny = 0; % Wilket Creek (Papangelakis et al 2019)

if site == 1
    def_data.D50 = 37; % mm
    def_data.szcat = 2.^[5:1:8];
    def_data.fdir = 'Ganatsekaigon\';
    def_data.fchar = 'Tag Characteristics 2.csv';
    def_data.thalname = 'GCrk_Thalweg_rev.csv';
    def_data.taglist = 'DS_Ganny_tagnum.csv';
end
