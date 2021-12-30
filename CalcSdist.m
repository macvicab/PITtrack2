function [Sdist,Ndist] = CalcSdist(thalname,E,N)
% calculates streamwise and normal distance from thalweg or channel
% centerline

thalXY = ConvCSV2Struct(thalname,1);


% calculate thalweg distances
thalweg = CalcThalweg(thalXY.Easting,thalXY.Northing,thalXY.elevation);
[Sdist,Ndist] = CalcLdistance(E,N,thalweg.Emore,thalweg.Nmore,thalweg.ldistmore);

end

function [ldist,hdist] = CalcLdistance(E,N,Ethalweg,Nthalweg,Ldistance)
% create zero matrices for the positions

nttot = length(E);

ldist = zeros(nttot,1);
hdist = zeros(nttot,1);

for nt = 1:nttot
    % calculate l and h distances based on the thalweg distance
    [~,hdall] = cart2pol(E(nt)-Ethalweg,N(nt)-Nthalweg);
    % find the thalweg position closest to the cartesian coordinate
    [hdist(nt),thalnum] = min(hdall);
    % save the ldistance
    ldist(nt) = Ldistance(thalnum);
end

end

function thalweg = CalcThalweg(X,Y,Z)
xi = ~isnan(X)&~isnan(Y)&~isnan(Z);

X = X(xi);
Y = Y(xi);
Z = Z(xi);

% interpolate thalweg to a finer grid
int = 0.1;

% calculate ldistance
nthaltot = length(X);
thalweg.ldistance = zeros(1,nthaltot);

[~,rdist] = cart2pol(X(2:end)-X(1:end-1),Y(2:end)-Y(1:end-1));
thalweg.Easting = X;
thalweg.Northing = Y;
thalweg.ldistance(2:end) = cumsum(rdist);

% increase resolution of thalweg data
thalweg.ldistmore = 0:int:thalweg.ldistance(end);

thalweg.Nmore = interp1(thalweg.ldistance,Y,thalweg.ldistmore);
thalweg.Emore = interp1(thalweg.ldistance,X,thalweg.ldistmore);
thalweg.Zmore = interp1(thalweg.ldistance,Z,thalweg.ldistmore);

end