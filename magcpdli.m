%%
%%Method proposed by Li et al. (2013) 
%%
%By 
%Carrillo-de la Cruz, J. L. : juanluiscc9@comunidad.unam.mx
%Gomez, Dario : dariogr9@outlook.com
%Velazques, Pablo : pablovelsan04@gmail.com
%%%%%%%%%%%%%%%%%%%% Description %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INVERSION 
     % 1- Import data [GEOtiff (*.tif) or radially averaged amplitude spectrum (*.txt)] and flight height (km)
     % example : defractal('synt_200.tif',0) or defractal('RAAS.txt',0)
     % 2- Choose the wavenumber range for to the top of the magnetic source 
     % 3- Choose the wavenumber range for to the centroid of the magnetic source
     % NOTE: The .txt are files with two columns: wavenumber (1/km) and
     % radially averaged amplitude spectrum,
     % optinally, the condifence intervals can be add in a third column.
     % View Outputs : Depth to the top, centroid and bottom of magnetic sources, 
     %                and their uncertainty.
     % View plots : Corrected radially averaged amplitude spectrum and corrected scaled spectrum 
     %              forward modeling and the table with results.
     %   Select yes (y) or no (n) for the calculation of Heat flow and temperature at depth modeling 
     %
     %
     % Optional modeling: The temperature model needs the parameters of 
     %                       Zb (km) is the Curie Point depth 
     %                       Ts (°C) the mean annual surface temperature 
     %                       Ao (microW/m3) is heat production at the surface
     %                       b (km) is the empirical heat-production depth-distribution parameter
     %                       k (W/mK) is thermal conductivity
     %
function magcpdli(varargin)
switch nargin
    case 1
        data = varargin{1};
        Zc = 0;
    case 2
        data = varargin{1};
        Zc = varargin{2};
    otherwise
        error('Unexpected input')
end
%%calculation of DBMS/Curie depth
if fileid(data) == 1
    b=data_aero(data);
    dt=detrend_2d(b); % Polynomial Trend Surface Estimation (degree = 1)
    a2=double(raps(dt)); % Radially average power spectrum
elseif fileid(data) == 0
    reply = input('Does your data have confidence intervals? (y/n): ','s');
    clc
    if strcmp(reply,'n')
        a1=load(data);
        a2=[a1 zeros(length(a1),1)];
        cfint = 0;
    elseif strcmp(reply,'y')
        a2=load(data);
        cfint=1;
    else
        disp('Error: Unexpected input')
    end
else
    disp('Error: Unexpected input')
end
Ans=1;
while Ans==1  
%clear all;
clc;
close all;
%%
% fractal exponent 
b2=input('Enter the fractal exponent: ');
clc;
%%
%Calling the function basedepth
[w2,pw2,ps2,ci1]=basedepth_scaling_defractal_interactive(a2,b2); 
%%
%Ploting figures to select data used in slope computation
f1=figure;
plot(w2,pw2,'k.-');
ylabel('ln(A(k))','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
xlabel('k ($\frac{1}{km}$)','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
ylim([floor(min(pw2)) ceil(max(pw2)+1)]);
[~,xs2,ys2] = selectdata('sel','br','BrushShape','rect','BrushSize',0.01);
close(f1)
p2=polyfit(xs2,ys2,1); %slope Zt and intercept
wh2=min(xs2)-0.005:0.001:max(xs2)+0.005; %wavenumber to graphs
kk2=polyval(p2,wh2); %values for files
Zt=(-1)*(p2(1)/(2*pi));
sd2=sqrt(sum((ys2-polyval(p2,xs2)).^2)/(length(xs2))); %standart deviation
Zt_error= sd2/((2*pi)*(max(xs2)-min(xs2)));
%%
f2=figure;
plot(w2,ps2,'k.-');
ylabel('ln($\frac{A(k)}{k}$)','Interpreter', ...
    'latex','FontSize',12,'FontName','Times New Roman')
xlabel('k ($\frac{1}{km}$)','Interpreter','latex', ...
    'FontSize',12,'FontName','Times New Roman')
[~,xs,ys] = selectdata('sel','br','BrushShape','rect','BrushSize',0.01);
close(f2)
pl=polyfit(xs,ys,1); %slope Z0 and intercept
w1=min(xs)-0.005:0.001:max(xs)+0.005; %wavenumber to graphs
kk=polyval(pl,w1); %values for files
Z0=(-1)*(pl(1)/(2*pi));
sd1=sqrt(sum((ys-polyval(pl,xs)).^2)/(length(xs))); %standart deviation
p1_err= sd1/((2*pi)*(max(xs)-min(xs)));
Zb=2*Z0-Zt; %depth to the bottom of the magnetic source
Zb_err=sqrt((4*(p1_err))^2+(Zt_error^2)); %DBMS uncertainity (Martos et al., 2018)
%% Forward model
g2=a2(1:size(a2),2); %Amplitude spectrum
sx=w2(1:find(w2==xs2(length(ys2)))); sy=g2(1:find(w2==xs2(length(ys2))));
[Mx,I]=max(sx); 
x0=Mx + 2*pi*sx(I)*Zt + 0.5*(b2-1)*log(2*pi*sx(I));
func = @(x0,x)x0 - ((2*pi*x)*Zt)- log( (2*pi*x).^(0.5*(b2-1)) )...
    + log((1-exp(-2*pi*x*(Zb-Zt)))); %Synthetic amplitud Spectrum from Li et al. (2013)
[xc,resnorm] = lsqcurvefit(func,x0,sx,sy);
mf=sqrt(resnorm/length(sx)); %misfit
wn=linspace(w2(1),w2(end));
fm1=func(xc,wn); %Synthetic Spectrum
fit_fm=mf; %Misfit
%%
%ploting figures 
figure('units','centimeters','position',[1 1 25 25]);
ax1 = subplot(2,2,1);
hold(ax1,'on')
if (fileid(data) == 1 || cfint == 1)
    fill([w2', fliplr(w2')],[pw2'+ci1',fliplr(pw2'-ci1')],...
        [0.6 0.6 0.6],'edgecolor','k','Facealpha',0.3,'LineStyle','none');
end
plot(w2,pw2,'g-');
plot(wh2,kk2,'-r','Linewidth',1,'Markersize',12);
box on
ylabel('ln(A(k))','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
xlabel('k ($\frac{1}{km}$)','Interpreter','latex','FontSize',12,'FontName','Times New Roman');
xlim(ax1,[0 w2(end)])
ylim([floor(min(pw2)) ceil(max(pw2))])
hold(ax1,'off')
%
ax2 = subplot(2,2,2);
hold(ax2,'on')
plot(w2,ps2,'g-'); xlim([0 max(w2)]); box on
plot(w1,kk,'-r','Linewidth',0.8,'Markersize',12);
ylabel('ln($\frac{A(k)}{k}$)', ...
    'Interpreter','latex','FontSize',12,'FontName','Times New Roman')
xlabel('k ($\frac{1}{km}$)','Interpreter', ...
    'latex','FontSize',12,'FontName','Times New Roman')
hold(ax2,'off')
%
ax3 = subplot(2,2,3);
hold(ax3,'on')
if (fileid(data) == 1 || cfint == 1)
    p3=fill([w2', fliplr(w2')],[a2(1:size(a2),2)'+ci1',...
        fliplr(a2(1:size(a2),2)'-ci1')],[0.6 0.6 0.6],...
        'edgecolor','k','Facealpha',0.3, 'LineStyle','none');
end
p1=plot(w2,a2(1:size(a2),2),'b-','Linewidth',1.5,'Markersize',12); %Observed Spectrum
p2=plot(wn,fm1,'-r','Linewidth',1,'Markersize',12); %Synthetic Spectrum
box on
ylabel('ln(A(k))','Interpreter','latex','FontSize',12, ...
    'FontName','Times New Roman')
xlabel('k ($\frac{1}{km}$)','Interpreter','latex', ...
    'FontSize',12,'FontName','Times New Roman');
xlim(ax3,[0 0.2])
if (fileid(data) == 1 || cfint == 1)
    legend([p1 p2 p3],{'Observed spectrum',...
        ['Synthetic spectrum (\beta = ' num2str(b2) ')'],'95% confidence interval'})
elseif fileid(data) == 0
    legend([p1 p2],{'Observed spectrum',...
        ['Synthetic spectrum (\beta = ' num2str(b2) ')']})
end
hold(ax3,'off')
%
subplot(2,2,4)
rectangle('Position',[0,0,1,1],'FaceColor','w')
str = {'Values',['\beta = ' num2str(b2)],...
    ['Z_t = ' num2str(Zt-Zc,'%#.2f') ' km ^+_- ' num2str(Zt_error,'%.2f') ' km'],...
    ['Z_0 = ' num2str(Z0-Zc,'%#.2f') ' km ^+_- ' num2str(p1_err,'%.2f') ' km'],...
    ['Z_b = ' num2str(Zb-Zc,'%#.2f') ' km ^+_- ' num2str(Zb_err,'%.2f') ' km'],...
    ['Misfit = ' num2str(fit_fm,'%#.4f')]};
text(0.05,0.5,str,'FontSize',14); axis off
%%
fprintf(['Fractal exponent: %.2f\n' 'Zt = %.4f km\n' 'Z0 = %.4f km\n' 'Zb = %.4f km\n'],b2,Zt-Zc,Z0-Zc,Zb-Zc); %Depths of magnetic body
fprintf(['Zt error = %.4f km\n' 'Z0 error = %.4f km\n' 'Zb error = %.4f km\n' 'Misfit = %.4f\n'],Zt_error,p1_err,Zb_err,fit_fm); %Errors
%% Yes/No rerun
reply = input('Would you like to run the program again? (y/n): ','s');
if strcmp(reply,'n')
    clc
    fprintf(['Fractal exponent: %.2f\n' 'Zt = %.4f km\n' 'Z0 = %.4f km\n' 'Zb = %.4f km\n'],b2,Zt-Zc,Z0-Zc,Zb-Zc); %Depths of magnetic body
    fprintf(['Zt error = %.4f km\n' 'Z0 error = %.4f km\n' 'Zb error = %.4f km\n' 'Misfit = %.4f\n'],Zt_error,p1_err,Zb_err,fit_fm); %Errors
    p=input('\nDo you want to model the temperature profile? (y/n): ','s');
    if strcmp(p,'y')
        geotherms(Zb)
    elseif strcmp(p,'n')
        break
    else
        disp('Error: Unexpected input')
        break
    end
    break
elseif strcmp(reply,'y')
    Ans=1;
else
    disp('Error: Unexpected input')
    break
end
end
%%
function x=fileid(data)
% This program identifies the input file extension. If the files have one 
% of these extensions '.csv', '.dat', '.txt', '.SPC', the return value is 0.
% If the file has an extension of '.tif', the return value is 1.
[~,~,ext] = fileparts(data);
if (strcmpi(strtrim(ext),'.csv') || strcmpi(strtrim(ext),'.txt') || strcmpi(strtrim(ext),'.dat') || strcmpi(strtrim(ext),'.SPC'))
    x=0;
elseif (strcmpi(strtrim(ext),'.tif'))
    x=1;
else
    disp('Error: Unexpected file')
end
end
%%
%function to import the aeromagnetic data (GEOtif)
function [b]=data_aero(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read
%  Auto-generated by MATLAB on 23-May-2019 17:59:12
% Import the file
rawData1 = importdata(fileToRead1);
% For some simple files (such as a CSV or JPEG files), IMPORTDATA might
% return a simple array.  If so, generate a structure so that the output
% matches that from the Import Wizard.
[~,name] = fileparts(fileToRead1);
newData1.(matlab.lang.makeValidName(name)) = rawData1;
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
[b]=rawData1;
end
%%
%%
%function to Polynomial detrend Surface Estimation (degree = 1)
function Z_f = detrend_2d(Z)
%This function is written by Munther Gdeisat-The General Engineering
%  Research Institute (GERI) at Liverpool John Moores University.
% This program is written on 9th October 2011
%This function is the 2D equivalent of detrend function in Matlab
%  Z_f = DETREND(Z) removes the best plane fit trend from the
%     data in the 2D array Z and returns the residual in the 2D array Z_f
if size(Z,2) < 2
    disp('Z must be a 2D array')
    return
end 
M = size(Z,2);
N = size(Z,1);
[X,Y] = meshgrid(1:M,1:N);
%Make the 2D data as 1D vector
Xcolv = X(:); % Make X a column vector
Ycolv = Y(:); % Make Y a column vector
Zcolv = Z(:); % Make Z a column vector
Const = ones(size(Xcolv)); % Vector of ones for constant term
% find the coeffcients of the best plane fit
Coefficients = [Xcolv Ycolv Const]\Zcolv; % Find the coefficients
XCoeff = Coefficients(1); % X coefficient
YCoeff = Coefficients(2); % X coefficient
CCoeff = Coefficients(3); % constant term
% detrend the data
Z_p = XCoeff * X + YCoeff * Y + CCoeff;
Z_f = Z - Z_p;
end
%%
%function to calculate radially average power spectrum 
function [cv3]= raps(img)
% function raps(img)
%
% Computes and plots radially an average
% spectrum of magnetic data in a GEOtif image
% of the amplitude spectrum
% Based and modified  on (C) E. Ruzanski, RCG, 2009 and 
% Tom Ridsdill-Smith  March 2000
%% Process image size information
[N, M] = size(img);
%% Compute spectrum
imgf = fftshift(fft2(img));
imgfp = (abs(imgf)).^2;
%% Adjust PSD size
dimDiff = abs(N-M);
dimMax = max(N,M);
% Make square
if N > M                                                                    % More rows than columns
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(N,dimDiff/2) imgfp NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
    else                                                                    % Odd difference
        imgfp = [NaN(N,floor(dimDiff/2)) imgfp NaN(N,floor(dimDiff/2)+1)];
    end
elseif N < M                                                                % More columns than rows
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(dimDiff/2,M); imgfp; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
    else
        imgfp = [NaN(floor(dimDiff/2),M); imgfp; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
    end
end
%% Compute radially average amplitude spectrum
[X, Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
[~, rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
rho = round(rho);
i = cell(floor(dimMax/2) + 1, 1);
for r = 0:floor(dimMax/2)
    i{r + 1} = find(rho == r);
end
Pf = zeros(1, floor(dimMax/2)+1);
CF = zeros(1, floor(dimMax/2)+1);
for r = 0:floor(dimMax/2)
    Pf(1, r + 1) = nanmean( log( sqrt( imgfp( i{r+1} ) ) ) ); %amplitde spectra
    CF(1, r + 1) = 1.96 * ( nanstd( log( sqrt( imgfp( i{r+1} ) ) ) ) ... 
                   / sqrt( length( sqrt( imgfp( i{r+1} ) ) ) ) ); %confidence interval
end
%% calculate the 1D angular frequencies
% (in radians) for a signal of length n
%using the fcoef fuction
w=fcoef1(dimMax);
w=w(1:floor(dimMax/2)+1);
w1=w';
Pf2=Pf';
CF2=CF'; %confidence interval 
cv2=[w1(2:size(w1)),Pf2(2:size(Pf2)),CF2(2:size(CF2))];
[cv3]=cv2(cv2(:,3) <= 1, :);
end
%%
function [w]=fcoef1(n)
% function [w]=fcoef1(n)
%
% Calculates the 1D Fourier angular frequencies
% for a signal of length n
% in (1/K_unit)
% % Tom Ridsdill-Smith   March 2000
%
if (rem(n,2)==0)
   % Even signal length
   w=[(0:(n/2)) (((-n/2)+1):-1)].*(1/n);
   
else
   % Odd signal length
   w=[(0:((n-1)/2)) ((-(n-1)/2):-1)].*(1/(n-1));
   
end
end
%%
%function to calculate the temperature profie at depth
function geotherms(Zb)
% Geotherms calculates the heat flow at the surface using ecuation from Martos  et
% al., 2017 and plot the profile temperatura using ecuation from Ravat et
% al,. 2016. 
% geo(Zb,Ts,Ao,b)
% Zb (km) is the Curie Point depth 
% Ts (ï¿½C) the mean annual surface temperature 
% Ao (microW/m3) is heat production at the surface
% b (km) is the empirical heat-production depth-distribution parameter
% k (W/mK) is thermal conductivity
%%
Zb=1000*Zb; %Zb (m)
%%
c=1;
while c==1
Ts=20; b=10; Ao=2.5; Tc = 580;
fprintf(['\nDefault Values \nSurface temperature (Ts): %.2f °C\n' ...
    'Curie temperature (Tc): %.2f °C\n' ...
    'Surface heat production (Ao) = %.2f microW/m3\n' ...
    'Sacaling length of surface radiogenic heat production (b) = %.2f km\n'],Ts,Tc,Ao,b);
P=input('\nDo you want to use these default values? (y/n): ','s');
if strcmp(P,'y')
    V=[Ts,Tc,Ao,b];
elseif strcmp(P,'n')
    V=input('Enter values [Ts,Tc,Ao,b]: ');
else
    disp('Error: Unexpected input')
end
Ts=V(1); Tc=V(2); A=V(3)*1e-6; b=V(4)*1e3;
k=input('Enter thermal conductivity: '); %Thermal conductivity
qs = k*(Tc-Ts)/Zb + A*b - ((A*b^2)/Zb)*(1-exp(-Zb/b)); %Martos et al., 2017 %heat flow at the surface (W/m2)
fprintf('K = %.2f W/(m°C) and Qs = %.4f mW/m2\n',k,qs*1e3); %Surface Heat Flow
zb=Zb/1000; dlim=ceil(zb/5)*5;
z=0:1000*dlim;
Tcrust=Ts + z*qs/k + A*b*(b-z)/k - A*b^2*exp(-z/b)/k; %Ravat et al., 2016
tm=max(Tcrust); Tlim=ceil(tm/100)*100;
figure('units','centimeters','position',[6 6 10 20]);
hold on
plot(Tcrust,z'/1000)
plot([0, Tlim], [zb, zb],'--k','LineWidth',1.5) %Curie Point Depth
plot([Tc, Tc], [0, dlim],'--r','LineWidth',1.5) %Curie temperature
axis ij
hold off
legend({['K = ' num2str(k,'%#.2g') ' W/(mK)'],...
    [' CPD = ' num2str(zb,'%.2f') ' km'],...
    [' CT = ' num2str(Tc) ' °C']})
title('Steady-state one-dimensional geotherm'); xlabel('Temperature (°C)'); ylabel('Depth (km)'); grid on; box on; ylim([0,dlim]); xlim([0,Tlim]);
%% Yes/No rerun
reply = input('\nWould you want to make a new model? (y/n): ','s');
if strcmp(reply,'n')
    break
elseif strcmp(reply,'y')
    c=1;
else
    disp('Error: Unexpected input')
    break
end
end
end
%%
%Calculation of Zt and Z0
function[w,pw,ps,ci]=basedepth_scaling_defractal_interactive(x,b1) %fractal correction used by li 2013
w=x(1:size(x),1); % wavenumber (1/K_unit)
ps=(x(1:size(x),2)- log( (2*pi).*w )) + log((( (2*pi) .*w ).^( 0.5*(b1-1) ) ));% for Z0
pw=x(1:size(x),2)+log(( (2*pi) .*w ).^(0.5*(b1-1)));% for Zt
ci=x(1:size(x),3);
end 
%%
function [pointslist,xselect,yselect] = selectdata(varargin)
% selectdata: graphical selection of data points on a plot using the mouse
% usage: pointslist = selectdata         % uses all default options
% usage: pointslist = selectdata(prop1,val1,prop2,val2,...)
%
% SELECTDATA allows the user to select data points on a given plot
% using the mouse, in a variety of modes. 'Lasso' mode allows the
% selection by a user directed lasso around the points. 'Brush' mode
% selects points as you brush over them with the mouse. 'Rect' mode
% draws a rectangle, selecting any points inside the rectangle.
% 'Closest' mode looks for a aingle mouse click, finding the closest
% point to the mouse.
%
% Returned is a list of the point indexes selected, additionally you
% can specify that the points be deleted from the plot.
%
% Arguments: (input)
%  The input arguments to SELECTDATA are all property/value pairs.
%  (See PARSE_PV_PAIRS for more details.)
%  http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9082&objectType=FILE
% 
%  Property names and character values can be shortened as long as
%  the shortening is unambiguous. Capitalization is ignored.
%
%  Valid Properties: 'Action', 'Axes', 'BrushSize', 'Identify',
%           'Ignore' , 'Pointer', 'SelectionMode'
%
%  'Action' - {'list', 'delete'}
%           'delete' causes the selected points to be deleted from
%           the figure after their selection is final.
%
%           DEFAULT VALUE: 'List'
%
%  'Axes' - denotes the axes from which to select points
%
%           DEFAULT VALUE: gca
%  
%  'BrushShape' - {'rect', 'circle'}
%
%           DEFAULT VALUE: 'circle'
%
%           Sets the shape of the brush to be used. Both brush
%           shapes are relative to the figure axes, so a nominally
%           "circular" brush is actually elliptical if the axis
%           units/lengths are unequal.
%
%           Only used when SelectionMode is 'brush'.
%
%  'BrushSize' - Only used when SelectionMode is 'brush'.
%
%           DEFAULT VALUE: 0.05
%
%           The default value will specify a rectangular brush
%           that has dimensions of 5% of the current axes. Note
%           that on a set of square axes, the brush will always
%           look square, even if the axes have very different
%           units.
%
%           0 < brushsize <= 0.25
%
%  'Identify' - {'on', 'off'}
%           Causes the selected points to be temporarily over-plotted
%           with a new filled red marker ('o') as they are selected.
%
%           Points selected with a lasso, or rect may be deselected
%           as the tool is modified. Brush selections are cumulative.
%
%           DEFAULT VALUE: 'on'
%
%  'Ignore' - a data handle, or []
%
%           A list of data handles to be ignored. This allows you to
%           use selectdata on only some sets of points, while others
%           in the same figure are ignored. This is a useful option
%           when you may have plotted some data points but also a
%           curve fit through your data. You can then cause the plotted
%           curve to be ignored by selectdata.
%
%           DEFAULT VALUE: []
%
%  'Label' - {'off', 'on'}
%
%           Causes text labels with their (x,y) coordinates to appear
%           next to each point selected.
%
%           Beware that selecting large numbers of points and creating
%           and displaying the label for them can be time consuming.
%           This option is a great one for single point selection,
%           but I have seen system-related problems when rapidly
%           selecting & deselecting large numbers of points with the
%           rect tool.
%
%           DEFAULT VALUE: 'off'
%
%  'Pointer' - {'crosshair' | 'fullcrosshair' | 'arrow' | 'ibeam' |
%           'watch' | 'topl' | 'topr' | 'botl' | 'botr' | 'left' |
%           'top' | 'right' | 'bottom' | 'circle' | 'cross' | 'fleur' |
%           'hand' }
%
%           Changes the cursor pointer while selection is active.
%           After selection terminates, the figure pointer is
%           restored to its old setting.
%
%           DEFAULT VALUE: 'crosshair'
%
%  'Return' - {'selected' | 'unselected' }
%
%           Selection of data points, perhaps if used to indicate
%           outliers in data, would normally return that set of 
%           points selected. But some users may prefer to see the
%           list of points returned to be those points which were
%           NOT selected.
%
%           DEFAULT VALUE: 'selected'
%
%  'SelectionMode' - {'Lasso', 'Brush', 'Rect', 'Closest'}
%
%           DEFAULT VALUE: 'Lasso'
%
%           If 'Brush' is used, then the brush will be a rectangular
%           brush, with a size as defined by the 'BrushSize' property. 
%           The brush will be centered at the current mouse coordinates.
%           Brush size will be a fraction of the axis limits for the
%           current axes. Click inside the axes, then drag the "brush"
%           around. Any points the brush crosses over will be selected.
%
%           If 'Lasso' is chosen, then click inside the axes to define
%           one end of the lasso, then drag with the mouse still down,
%           causing the mouse to define a general curvilinear region.
%           The polygon will close automatically when the mouse button
%           is released. A point is "inside" the lasso if inpolygon
%           identifies it as so. BEWARE of convoluted lassos that
%           intersect themselves.
%
%           If 'Rect' is chosen, then click inside the axes to define
%           one corner of the rect, dragging to specify the opposite
%           corner, just as rbbox would do.
%
%           If 'Closest' was chosen, then a single mouse click in the
%           figure window is used, then that point which is closest in
%           in Euclidean distance (in window units) is chosen. You
%           can move the cursor around (don't release the mouse until
%           you are done) and the currently selected point will be
%           highlighted.
%
%  'Verify' - { 'off' | 'on' }
%
%           If set to 'on', this causes a dialog box to pop up after
%           selection. The user can then acccept the selection, redo
%           it, or cancel out, causing no points to be selected.
%
%           Note, if cancel is chosen from the dialog, and 'return'
%           was specify to return those points 'unselected', then ALL
%           the points will actually be returned.
%
%           DEFAULT VALUE: 'off'
%
%
% Note: other properties are available for use, but I've chosen to
% leave them semi-hidden because they don't seem terribly useful to
% most users. These properties allow you to specify the colors of the
% selection tool itself, the colors of the selected point markers,
% the transparency of the selection tool, the marker itself, etc.
% Default values for these parameters are:
%
%           FlagMarker    = 'o'
%           FlagColor     = 'r'
%           Fill          = 'on'
%           FillColor     = 'y'
%           FillEdgeColor = 'b'
%           FillTrans     = 0.5
%           MaxBrush      = 0.25
%           RemoveTool    = 'on'
%           RemoveFlagged = 'on'
%           RemoveLabels  = 'on'
%
% Further documentation on these parameters can be found by editting
% selectdata.m itself.
%
%
% Arguments: (output)
%  pointslist - list of points that were selected. If only one
%           dataset was found, then points list is a simple vector.
%           If multiple sets of points were found in the axes, then
%           points list will be a cell array.
%
%           NOTE: Each set of points is peeled off the stach in the
%           matlab stores it.
%
%  xselect - array (or cell array) containing the x coordinates of
%           those points identified in the selection
%
%  yselect - array (or cell array) containing the y coordinates of
%           those points identified in the selection
%
%
% Example:
%  Plot a set of points, then select some with the mouse
%  using a rectangular brush, return the indices of those
%  points selected. (Note: user interaction on the plot
%  will be necessary.)
%
%    x = 0:.1:1;
%    y = x.^2;
%    plot(x,y,'o')
%    pl = selectdata('selectionmode','brush');
%  
% Example:
%  Select a single point with the mouse, then delete that
%  selected point from the plot.
%
%    pl = selectdata('selectionmode','closest','action','delete');
%
% Example:
%  Select some points using a rect(rbbox) tool, also return
%  the (x,y) coordinates from multiple curves plotted. Use
%  shortened versions of the properties and values.
%  
%    plot(rand(5,2),rand(5,2),'o')
%    [pl,xs,ys] = selectdata('sel','r');
%  
% Example:
%  Plot a curve and some data points on one plot, select
%  some points from the data plotted, but ignore the
%  smooth curve, even if the lasso passes over it.
%    x = 0:.01:1;
%    y = exp(x);
%    ynoisy = y + randn(size(y))/2;
%    h1 = plot(x,y,'-');
%    hold on
%    h2 = plot(x,ynoisy,'o');
%    [pl,xs,ys] = selectdata('sel','lasso','ignore',h1);
%
% See also:
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 2/19/07
% defaults for the parameters
params.Axes = gca;
params.SelectionMode = 'lasso';
params.Action = 'list';
params.BrushShape = 'circle';
params.BrushSize = .05;
params.Identify = 'on';
params.Ignore = [];
params.Pointer = 'cross';
params.Return = 'selected';
params.Verify = 'off';
params.Label = 'off';
% Undocumented options, also unchecked for validity.
% These parameters control the marker used to identify
% those points which are currently selected.
% FlagMarker must be a valid plot marker, FlagColor
% must be a valid color spec. The default values are...
params.FlagMarker = 'o';
params.FlagColor = 'r';
% More (unchecked) options that are yours to fiddle with
% (or not.) These control the fill color to be applied
% to the interior of the lasso, rect, and brush. Also
% controlled are the degree of transparency to be applied.
params.Fill = 'on';
params.FillColor = 'y';
params.FillEdgeColor = 'b';
params.FillTrans = 0.5; % must be in the interval [0,1]
% The maximum relative brushsize allowed is also
% controlled here, just in case I ever wanted to allow
% the brush to be a bit larger.
params.MaxBrush = 0.25;
% After selection has been accomplished, the flagged points
% and the selection tool itself are normally deleted. But
% in some circumstances it may be useful to not delete
% those objects. 'off' will cause these objects to remain.
params.RemoveTool = 'on';
params.RemoveFlagged = 'on';
params.RemoveLabels = 'on';
% save this to restore later
axissize = axis;
% process any Property/value pairs.
if nargin>0
  params = parse_pv_pairs(params,varargin);
end
% check the supplied parameters for validity
params = check_params(params);
% bring the focus to the figure containing the
% designated axes
fighandle = get(params.Axes,'parent');
figure(fighandle)
% get the current figure pointer, so we can
% restore it later on
oldpointer = get(fighandle,'Pointer');
% extract xdata and ydata from the specified axes
% get the children of the axes
hc = get(params.Axes,'children');
% are any of the data handles to be ignored?
if ~isempty(params.Ignore)
  hc = setdiff(hc,params.Ignore);
end
% strip out xdata and ydata
xdata = get(hc,'xdata');
ydata = get(hc,'ydata');
% if we must highlight the points as they are
% selected, then for efficiency we need to know
% how many we may expect.
flaghandle = [];
if ~iscell(xdata)
  xdata = xdata(:);
  ydata = ydata(:);
  % total number of data points
  npoints = length(xdata);
else
  for i = 1:length(xdata)
    xdata{i} = xdata{i}(:);
    ydata{i} = ydata{i}(:);
  end
  % total number of data points
  npoints = cellfun('length',xdata);
end
% for textlabels if desired
texthandles = [];
% set up a while loop in case we need to verify
% satisfaction
selectionflag = true;
while selectionflag
  % and the total number currently selected
  nsel = 0;
  
  % what SelectionMode was specified
  switch params.SelectionMode
    case 'closest'
      % Find the single closest point (in Euclidean distance)
      % to the mouse click
      
      % set the figure pointer
      if ~isempty(params.Pointer)
        set(fighandle,'Pointer',params.Pointer)
      end
      
      % mouse click?
      waitforbuttonpress;
      
      % Button Motion and Button Up
      set(fighandle,'WindowButtonMotionFcn',@CPmotion);
      set(fighandle,'WindowButtonUpFcn',@selectdone);
      
      % dx, dy to scale the distance
      dx = (axissize(2) - axissize(1));
      dy = (axissize(4) - axissize(3));
      
      % current closest point is
      cc = get(gca,'CurrentPoint');
      xy = cc(1,1:2);
      
      % what point was closest?
      [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy);
      nsel = 1;
      
      % identify any points?
      if strcmp(params.Identify,'on')
        flagpoints
      end
      
      % label them?
      if strcmp(params.Label,'on')
        maketextlabels
      end
      
      % selecthandle is not needed for this mode op operation
      selecthandle = [];
      
      % wait until the mouse button was released
      uiwait
      
      % ....
      
      % resume.
      
      % all we need to do here is restore the figure pointer
      % if we changed it before
      if ~isempty(params.Pointer)
        set(fighandle,'Pointer',oldpointer)
      end
      
    case 'rect'
      % Selection rect as a polygon
      % mouse click?
      waitforbuttonpress;
      % set the figure pointer
      if ~isempty(params.Pointer)
        set(fighandle,'Pointer',params.Pointer)
      end
      % button down detected
      cc = get(gca,'CurrentPoint');
      rectxy1 = cc(1,1:2);
      rectxy2 = rectxy1 + eps(rectxy1);
      % make a polygon of the box, initially of nil area
      xv = [rectxy1(1), rectxy2(1), rectxy2(1), rectxy1(1), rectxy1(1)];
      yv = [rectxy1(2), rectxy1(2), rectxy2(2), rectxy2(2), rectxy1(2)];
      % no points should been selected
      [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
      % and plot it, filled or not
      hold on
      if strcmp(params.Fill,'on')
        % filled
        selecthandle = fill(xv,yv,params.FillColor);
        set(selecthandle,'facealpha',params.FillTrans, ...
          'linestyle','--','edgecolor',params.FillEdgeColor)
      else
        % unfilled
        selecthandle = plot(xv,yv,'r:');
      end
      % we can undo the hold now
      hold off
      % Button Motion and Button Up
      set(fighandle,'WindowButtonMotionFcn',@rectmotion);
      set(fighandle,'WindowButtonUpFcn',@selectdone);
      % wait until the selection is done
      uiwait
      % ....
      % resume.
      % The rect already is a polygon, stored in (xv,yv)
    case 'lasso'
      % Selection lasso as a polygon
      % mouse click?
      waitforbuttonpress;
      % set the figure pointer
      if ~isempty(params.Pointer)
        set(fighandle,'Pointer',params.Pointer)
      end
      % button down detected
      cc = get(gca,'CurrentPoint');
      xlasso = cc(1,1);
      ylasso = cc(1,2);
      % form the polygon
      xv = xlasso;
      yv = ylasso;
      % and plot it, filled or not
      hold on
      if strcmp(params.Fill,'on')
        % filled
        selecthandle = fill(xv,yv,params.FillColor);
        set(selecthandle,'facealpha',params.FillTrans, ...
          'linestyle','--','edgecolor',params.FillEdgeColor)
      else
        % unfilled
        selecthandle = plot(xv,yv,'r:');
      end
      % we can undo the hold now
      hold off
      % Button Motion and Button Up
      set(fighandle,'WindowButtonMotionFcn',@lassomotion);
      set(fighandle,'WindowButtonUpFcn',@selectdone);
      % wait until the selection is done
      uiwait
      % ....
      % resume.
      % The lasso already is a polygon, stored in (xv,yv)
    case 'brush'
      % paint over the data, with a rectangular brush
      % mouse click?
      waitforbuttonpress;
      % set the figure pointer
      if ~isempty(params.Pointer)
        set(fighandle,'Pointer',params.Pointer)
      end
      % button down detected
      bc = get(gca,'CurrentPoint');
      brushcenter = bc(1,1:2);
      % dx, dy for the brush
      bdx = params.BrushSize*(axissize(2) - axissize(1));
      bdy = params.BrushSize*(axissize(4) - axissize(3));
      if strcmpi(params.BrushShape,'rect')
        % make the brush polygon as a fixed size rectangle
        % that we can slide around
        xv = brushcenter(1) + [-1, 1, 1, -1, -1]*bdx/2;
        yv = brushcenter(2) + [-1, -1, 1, 1, -1]*bdy/2;
      else
        % a circle was specified
        theta = linspace(0,2*pi,100);
        xv = cos(theta)*bdx/2 + brushcenter(1);
        yv = sin(theta)*bdy/2 + brushcenter(2);
      end
      % draw the initial brush polygon, filled or not
      hold on
      if strcmp(params.Fill,'on')
        % filled
        selecthandle = fill(xv,yv,params.FillColor);
        set(selecthandle,'facealpha',params.FillTrans, ...
          'linestyle','-','edgecolor',params.FillEdgeColor)
      else
        % unfilled
        selecthandle = plot(xv,yv,'r:');
      end
      hold off
      % have any points been selected?
      [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
      % identify any points?
      if strcmp(params.Identify,'on') && (nsel>0)
        flagpoints
      end
      
      % label them?
      if strcmp(params.Label,'on') && (nsel>0)
        maketextlabels
      end
      
      % Button Motion and Button Up
      set(fighandle,'WindowButtonMotionFcn',@brushmotion);
      set(fighandle,'WindowButtonUpFcn',@selectdone);
      % wait until the selection is done
      uiwait
      % ....
      % resume.
  end
  
  % verify?
  if strcmpi(params.Verify,'on')
    % pop up a dialog
    ButtonName = questdlg( ...
      'Are you satisfied with the points selected?','???', ...
      'Yes','Redo Selection','Cancel Selection','Yes');
    
    switch ButtonName
      case 'Yes'
        % we can drop through
        selectionflag = false;
        
      case 'Cancel Selection'
        % drop out, with nothing selected
        if ~iscell(xdata)
          pointslist = [];
          xselect = [];
          yselect = [];
        else
          for i = 1:numel(xdata);
            pointslist{i} = [];
            xselect{i} = [];
            yselect{i} = [];
          end
        end
        
        % we can drop through
        selectionflag = false;
        
      case 'Redo Selection'
        % or try again. The while loop will cycle
        % until happy or canceled
        
    end
  
  else
    % no verification was requested, so we want to
    % terminate the while loop after only one pass through.
    selectionflag = false;
  end
  
end
% pointslist and xselect, yselect are already complete.
% Do we delete the selected points?
if strcmpi(params.Action,'delete')
  if ~iscell(xdata)
    % only one set, so xdata and ydata are not cell arrays
    
    % Which points from the data fall in the selection polygon?
    xdata(pointslist) = [];
    ydata(pointslist) = [];
    
    % drop those points from the plot
    set(hc,'xdata',xdata,'ydata',ydata)
  else
    % it was a cell array, so there were multiple sets.
    for i = 1:numel(xdata);
      
      xdata{i}(pointslist{i}) = [];
      ydata{i}(pointslist{i}) = [];
      
      % drop those points from the plot
      set(hc(i),'xdata',xdata{i},'ydata',ydata{i})
    end
  end
end
% was 'return' set to be the selected list or the unselected one?
% Do nothing if 'selected', we are already done.
if strcmpi(params.Return,'unselected')
  if ~iscell(xdata)
    % only one set, so xdata and ydata are not cell arrays
    pointslist = setdiff((1:npoints)',pointslist);
    xselect = xdata(pointslist);
    yselect = ydata(pointslist);
  else
    % it was a cell array, so there were multiple sets.
    for i = 1:numel(xdata);
      pointslist{i} = setdiff((1:npoints(i))',pointslist{i});
      
      xselect{i} = xdata{i}(pointslist{i});
      yselect{i} = ydata{i}(pointslist{i});
    end
  end
end
% =====================================================
%     begin nested functions
% =====================================================
function brushmotion(src,evnt) %#ok
  % nested function for motion of the brush
  
  % get the new mouse position
  mousenew = get(params.Axes,'CurrentPoint');
  mousenew = mousenew(1,1:2);
  
  % make sure the axes are fixed
  axis(axissize)
  
  % how far did it move
  brushoffset = mousenew - brushcenter;
  brushcenter = mousenew;
  
  xv = xv + brushoffset(1);
  yv = yv + brushoffset(2);
  
  % did we brush over any new points
  [pl,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
  
  % for a brush, we need to append any selected points to
  % the already selected list
  if ~iscell(xdata)
    % only one set, so xdata and ydata are not cell arrays
    pointslist = union(pointslist,pl);
    xselect = xdata(pointslist);
    yselect = ydata(pointslist);
    nsel = length(pointslist);
  else
    % it was a cell array, so there were multiple sets.
    for j = 1:numel(pointslist);
      pointslist{j} = union(pointslist{j},pl{j});  %#ok
      pointslist{j} = pointslist{j}(:); %#ok
      
      if ~isempty(pointslist{j})
        xselect{j} = xdata{j}(pointslist{j});
        yselect{j} = ydata{j}(pointslist{j});
      end
    end
    
    % total of points selected
    nsel = sum(cellfun('length',pointslist));
  end
  
  % identify any points?
  if strcmp(params.Identify,'on')
    flagpoints
  end
  
  % label them?
  if strcmp(params.Label,'on')
    maketextlabels
  end
  
  % replot the brush in its new position
  set(selecthandle,'xdata',xv,'ydata',yv)
  
end
% =====================================================
% =====================================================
function CPmotion(src,evnt) %#ok
  % nested function to select the closest point
  
  % get the new mouse position
  mousenew = get(params.Axes,'CurrentPoint');
  xy = mousenew(1,1:2);
  
  % make sure the axes stay fixed
  axis(axissize)
  
  % what point was closest?
  [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy);
  nsel = 1;
    
  % identify any points?
  if strcmp(params.Identify,'on')
    flagpoints
  end
  
  % label them?
  if strcmp(params.Label,'on')
    maketextlabels
  end
  
end
% =====================================================
% =====================================================
function rectmotion(src,evnt) %#ok
  % nested function for expansion or contraction of the rect
  
  % get the new mouse position
  mousenew = get(params.Axes,'CurrentPoint');
  rectxy2 = mousenew(1,1:2);
  
  % make sure the axes are fixed
  axis(axissize)
  
  % update the rect polygon of the box, changing the second corner
  xv = [rectxy1(1), rectxy2(1), rectxy2(1), rectxy1(1), rectxy1(1)];
  yv = [rectxy1(2), rectxy1(2), rectxy2(2), rectxy2(2), rectxy1(2)];
    
  % did we brush over any new points
  [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
  
  % identify any points?
  if strcmp(params.Identify,'on')
    flagpoints
  end
  % label them?
  if strcmp(params.Label,'on')
    maketextlabels
  end
    
  % replot the rect in its new position
  set(selecthandle,'xdata',xv,'ydata',yv)
  
end
% =====================================================
% =====================================================
function lassomotion(src,evnt) %#ok
  % nested function for expansion of the lasso
  
  % get the new mouse position
  mousenew = get(params.Axes,'CurrentPoint');
  mousenew = mousenew(1,1:2);
  % append the new point on the end of the last lasso
  xlasso = [xlasso,mousenew(1,1)];
  ylasso = [ylasso,mousenew(1,2)];
  
  % and close it to form the polygon
  xv = [xlasso,xlasso(1)];
  yv = [ylasso,ylasso(1)];
  % replot the newly extended lasso
  set(selecthandle,'xdata',xv,'ydata',yv)
  
  % did we enclose any new points?
  [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
  
  % identify any points?
  if strcmp(params.Identify,'on')
    flagpoints
  end
  % label them?
  if strcmp(params.Label,'on')
    maketextlabels
  end
    
  % make sure the axes are maintained in size
  axis(axissize)
  
end
% =====================================================
% =====================================================
function selectdone(src,evnt) %#ok
  % nested function for mouse up
  
  % do we remove the selection tool?
  if strcmpi(params.RemoveTool,'on')
    % delete the selection object from the plot
    delete(selecthandle)
    selecthandle = [];
  end
  
  % do we remove the flagged points?
  if strcmpi(params.RemoveFlagged,'on')
    % also remove the flagged/plotted points
    if ~isempty(flaghandle)
      delete(flaghandle)
      flaghandle = [];
    end
  end
  
  % do we remove the flagged points?
  if strcmpi(params.RemoveLabels,'on')
    % also remove the labels
    delete(texthandles)
    texthandles = [];
  end
  
  % cancel the WindowButtonFcn's that we had set
  set(fighandle,'WindowButtonMotionFcn',[]);
  set(fighandle,'WindowButtonUpFcn',[]);
  
  % restore the figure pointer to its original setting
  if ~isempty(params.Pointer)
    set(fighandle,'Pointer',oldpointer)
  end
  
  % and resume execution, back in the mainline
  uiresume
end
% =====================================================
% =====================================================
function flagpoints
  % nested function for flagging the selected points
  
  % Are these the first points flagged? If so,
  % we need to plot them and set the marker, etc.
  if isempty(flaghandle) && (nsel > 0)
    % hold the figure, so we can add the flagged points
    hold on
    
    if ~iscell(xselect)
      flaghandle = plot(xselect,yselect,params.FlagMarker);
      set(flaghandle,'Color',params.FlagColor,'MarkerFaceColor',params.FlagColor)
    else
      flaghandle = plot(vertcat(xselect{:}),vertcat(yselect{:}),params.FlagMarker);
      set(flaghandle,'Color',params.FlagColor,'MarkerFaceColor',params.FlagColor)
    end
    
    % now release the hold
    hold off
  elseif ~isempty(flaghandle)
    % otherwise, we just need to update xdata and ydata
    
    if nsel == 0
      set(flaghandle,'xdata',[],'ydata',[]);
      
    elseif ~iscell(xselect)
      set(flaghandle,'xdata',xselect,'ydata',yselect);
    else
      set(flaghandle,'xdata',vertcat(xselect{:}),'ydata',vertcat(yselect{:}));
    end
  end
end
% =====================================================
% =====================================================
function maketextlabels
  % nested function for generation of text labels
  % over each point selected
  
  % We need to remove the last set of text labels
  delete(texthandles)
  
  % creat a new set of handles 
  if ~iscell(xselect)
    xtext = xselect;
    ytext = yselect;
  else
    xtext = vertcat(xselect{:});
    ytext = vertcat(yselect{:});
  end
  
  % anything selected?
  if ~isempty(xtext)
    textlabels = cell(1,nsel);
    for L = 1:nsel
      textlabels{L} = ['(',num2str(xtext(L)),',',num2str(ytext(L)),')'];
    end
    texthandles = text(xtext,ytext,textlabels);
  end
end
% =====================================================
end % mainline end
% ================================================
%               end main function
% ================================================
% ================================================
%                  subfunctions
% ================================================
function [pl,xsel,ysel,nsel] = testpoly(xv,yv,xdata,ydata)
% checks which points are inside the given polygon
% was there more than one set of points found in the plot?
if ~iscell(xdata)
  % only one set, so xdata and ydata are not cell arrays
  
  % Which points from the data fall in the selection polygon?
  pl = find(inpolygon(xdata,ydata,xv,yv));
  nsel = length(pl);
  
  xsel = xdata(pl);
  ysel = ydata(pl);
else
  % it was a cell array, so there were multiple sets.
  pl = cell(size(xdata));
  xsel = pl;
  ysel = pl;
  nsel = 0;
  for i = 1:numel(xdata);
    pl{i} = find(inpolygon(xdata{i},ydata{i},xv,yv));
    nsel = nsel + length(pl{i});
    
    if ~isempty(pl{i})
      xsel{i} = xdata{i}(pl{i});
      ysel{i} = ydata{i}(pl{i});
    end
    
  end
end
end % subfunction end
% ================================================
%                  subfunction
% ================================================
function [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy)
% find the single closest point to xy, in scaled units
if ~iscell(xdata)
  % just one set of points to consider
  D = sqrt(((xdata - xy(1))/dx).^2 + ((ydata - xy(2))/dy).^2);
  [junk,pointslist] = min(D(:)); %#ok
  xselect = xdata(pointslist);
  yselect = ydata(pointslist);
else
  % there is more than one set of points
  Dmin = inf;
  pointslist = cell(size(xdata));
  for i = 1:numel(xdata);
    D = sqrt(((xdata{i} - xy(1))/dx).^2 + ((ydata{i} - xy(2))/dy).^2);
    [mind,ind] = min(D(:)); %#ok
    
    if mind < Dmin
      % searching for the closest point
      Dmin = mind;
      
      pointslist = cell(size(xdata));
      xselect = cell(size(xdata));
      yselect = cell(size(xdata));
      
      pointslist{i} = ind;
      xselect{i} = xdata{i}(ind);
      yselect{i} = ydata{i}(ind);
    end
  end
end
end % subfunction end
% ================================================
%                  subfunction
% ================================================
% ============================================
% subfunction - check_params
% ============================================
function par = check_params(par)
% check the parameters for acceptability
%
% Defaults
%  Axes = gca;
%  SelectionMode = 'lasso';
%  Action = 'list';
%  BrushSize = .05;
% Axes == gca by default
if isempty(par.Axes)
  par.Axes = gca;
else
  if ~ishandle(par.Axes)
    error 'Axes must be the handle to a valid set of axes.'
  end
end
% SelectionMode == 'brush' by default
if isempty(par.SelectionMode)
  par.SelectionMode = 'brush';
else
  valid = {'rect', 'brush', 'lasso', 'closest'};
  if ~ischar(par.SelectionMode)
    error 'Invalid Style: Must be character'
  end
  
  ind = strmatch(lower(par.SelectionMode),valid); %#ok
  if isempty(ind) || (length(ind)>1)
    error(['Invalid SelectionMode: ',par.SelectionMode])
  end
  par.SelectionMode = valid{ind};
end
% BrushShape == 'circle' by default
if isempty(par.BrushShape)
  par.BrushShape = 'circle';
else
  valid = {'rect', 'circle'};
  if ~ischar(par.BrushShape)
    error 'Invalid Style: Must be character'
  end
  
  ind = strmatch(lower(par.BrushShape),valid); %#ok
  if isempty(ind) || (length(ind)>1)
    error(['Invalid SelectionMode: ',par.BrushShape])
  end
  par.BrushShape = valid{ind};
end
% Action == 'list' by default
if isempty(par.Action)
  par.Action = 'list';
else
  valid = {'list', 'delete'};
  if ~ischar(par.Action)
    error 'Invalid Action: Must be character'
  end
  
  ind = strmatch(lower(par.Action),valid); %#ok
  if isempty(ind) || (length(ind)>1)
    error(['Invalid Action: ',par.Action])
  end
  par.Action = valid{ind};
end
% Pointer == 'crosshair' by default, but
% if empty, will not change the pointer.
if ~isempty(par.Pointer)
  valid = {'crosshair', 'fullcrosshair', 'arrow', 'ibeam', ...
    'watch', 'topl', 'topr', 'botl', 'botr', 'left', 'top', ...
    'right', 'bottom', 'circle', 'cross', 'fleur', ...
    'custom', 'hand'};
  
  if ~ischar(par.Pointer)
    error 'Invalid Pointer: Must be character'
  end
  
  ind = strmatch(lower(par.Pointer),valid,'exact');
  if isempty(ind)
    ind = strmatch(lower(par.Pointer),valid); %#ok
    if isempty(ind) || (length(ind)>1)
      error(['Invalid Pointer: ',par.Pointer])
    end
  end
  par.Pointer = valid{ind};
end
% Identify == 'on' by default
if isempty(par.Identify)
  par.Identify = 'on';
else
  valid = {'on', 'off'};
  if ~ischar(par.Identify)
    error 'Value for Identify is invalid: Must be character'
  end
  
  ind = strmatch(lower(par.Identify),valid); %#ok
  if isempty(ind) || (length(ind)>1)
    error(['Invalid Action: ',par.Identify])
  end
  par.Identify = valid{ind};
end
% Label == 'off' by default
if isempty(par.Label)
  par.Label = 'off';
else
  valid = {'on', 'off'};
  if ~ischar(par.Label)
    error 'Value for Label is invalid: Must be character'
  end
  
  ind = strmatch(lower(par.Label),valid); %#ok
  if isempty(ind) || (length(ind)>1)
    error(['Invalid Action: ',par.Label])
  end
  par.Label = valid{ind};
end
% Return == 'selected' by default
if isempty(par.Return)
  par.Return = 'selected';
else
  valid = {'selected', 'unselected'};
  if ~ischar(par.Return)
    error 'Value for Return is invalid: Must be character'
  end
  
  ind = strmatch(lower(par.Return),valid); %#ok
  if isempty(ind) || (length(ind)>1)
    error(['Invalid Action: ',par.Return])
  end
  par.Return = valid{ind};
end
% Verify == 'off' by default
if isempty(par.Verify)
  par.Verify = 'off';
else
  valid = {'on', 'off'};
  if ~ischar(par.Verify)
    error 'Value for Verify is invalid: Must be character'
  end
  
  ind = strmatch(lower(par.Verify),valid); %#ok
  if isempty(ind) || (length(ind)>1)
    error(['Invalid Action: ',par.Verify])
  end
  par.Verify = valid{ind};
end
% BrushSize == 0.05 by default
if isempty(par.BrushSize)
  par.BrushSize = 0.05;
else
  if (length(par.BrushSize)>1) || (par.BrushSize<=0) || (par.BrushSize>par.MaxBrush) 
    error 'Brushsize must be scalar, and 0 < BrushSize <= 0.25'
  end
end
% Ignore == [] by default
if ~isempty(par.Ignore) && any(~ishandle(par.Ignore))
  error 'Ignore must be empty, or a data handle'
end
end % check_params
% ============================================
% Included subfunction - parse_pv_pairs
% ============================================
function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we
% have four parameters that we wish to use optionally in
% the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be
% supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params = 
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity'
% was truncated as supplied. Also note that the order the pairs were
% supplied was arbitrary.
npv = length(pv_pairs);
n = npv/2;
if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  return
end
if ~isstruct(params)
  error 'No structure for defaults was supplied'
end
% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
  p_i = lower(pv_pairs{2*i-1});
  v_i = pv_pairs{2*i};
  
  ind = strmatch(p_i,lpropnames,'exact');
  if isempty(ind)
    ind = find(strncmp(p_i,lpropnames,length(p_i)));
    if isempty(ind)
      error(['No matching property found for: ',pv_pairs{2*i-1}])
    elseif length(ind)>1
      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
    end
  end
  p_i = propnames{ind};
  
  % override the corresponding default in params
  params = setfield(params,p_i,v_i); %#ok
  
end
end % parse_pv_pairs
end
