function def_data = PITtrack_default(site)

%% write data control
def_data.writeAllSurvey = 1;
def_data.writeTravelLength = 1;

%% site specific data
% site 1 = Brad = 0; % Halfmoon Creek (Bradley 2017)
% site 2 = Ganny = 0; % Wilket Creek (Papangelakis et al 2019)

if site == 1
    def_data.D50 = 37; % mm
    def_data.szcat = 2.^[5:1:8]; % size categories in -phi values (log base 2 classification system)
    def_data.fdir = 'Ganatsekaigon\'; % directory path
    def_data.fchar = 'Tag Characteristics 2.csv'; % name of file that contains tracer characteristics such as sizes
    def_data.thalname = 'GCrk_Thalweg_rev.csv'; % name of file that has the thalweg data
    def_data.taglist = 'DS_Ganny_tagnum.csv'; % list of tag numbers to include in analysis (may be a subset of all tags in river)
end
