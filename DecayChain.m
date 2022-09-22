function [A, N, varargout]=DecayChain(numN,Th,N0,tunit,ts,te,tn,YNact,varargin)

% --- Update 22 SEP 2022 ---
% To probably run, this file MUST be named "DecayChain.m"
% No exception
% This file MUST be in the same directory as the main body of the code.
% See the main body of the code for my complaints about this function
% It works fine but I could have made things much neater by using multiple functions
% Instead this function does about everything
% And isn't commented very well
% --- Update Concluded ---
  
QOIreturn=[0]; % In case function provides no output. varargout will be
% set equal to QOIreturn at end of code.
%% QOI-related
% may eliminate QOI sections in future. Simply adjusting ending time should
% be sufficient.

YNQOI=varargin{1}; % QOI Y/N?
if nargin==13 % QOI + error
    QOInuclide=varargin{2};
    QOItime=varargin{3};
    YNQOIerror=varargin{4};
    QOIanalytic=varargin{5};
elseif nargin==11 % Just QOI
    QOInuclide=varargin{2};
    QOItime=varargin{3};
    YNQOIerror=0; % don't attempt to calculate QOI error
else % No QOI
    YNQOIerror=0; % don't calculate QOI error (there's no QOI)
end

%% Time Conversions

if tunit =='s'; % s->s
    tconv=1;
elseif tunit =='m'; % m->s
    tconv=60;
elseif tunit =='h'; % h->s
    tconv=3600;
elseif tunit =='d'; % d->s
    tconv=3600*24;
elseif tunit =='a'; % a->s
    tconv=3600*24*365;
else
    tconv=1; % if unit given doesn't match, will assume seconds.
    fprintf('Specified time unit invalid. Assuming times are given in seconds.');
end

%% Forward Euler Coefficient Matrix
digits(10); % precision.

L = vpa(log(2)./Th); % decay constants. vpa ensures precision.
dt = (te-ts)/tn; % time step length
% The for loop + if statement below will create our forward Euler coefficient
% matrix.
C=zeros(numN,numN); % initializing the coefficient matrix
for i=1:numN;
    for j=1:numN;
        if i==j;
            C(i,j)=1-L(j)*dt;
        elseif j==i-1
            C(i,j)=L(j)*dt;
        else
            C(i,j)=0;
        end
    end
end

%% Activity Plotting - initialization
% I will denote variables/secions pertaining to plotting activity with
% commenting P! after the line of code.
if YNact==1; % Do they want to plot activity?
    tnplot = 200; % P! activity plot & .txt file will only utilize this many points
    % may make this an input variable of the function in the future.
    dtplot = (te-ts)/tnplot; % P! time step
    tplot=ts; % P! intial time
    
    P=zeros(numN,tnplot);% P! this matrix will denote the activity measurements.
    P(:,1)=N0;
    index=1;
end
% P! End of activity plotting section

%% Forward Euler Method Execution
N=N0; % step #1

% if-else statement dependent on whether or not we will be plotting
% activity.
if YNact==1; % Do they want to plot activity?
    for p=ts:dt:te;
        N=C*N;
        % P! this section pertains to plotting activity
        if p>=dtplot+tplot;
            index=index+1; % Indexing for plot coordinates
            for z=1:numN;
                P(z,index)=N(z); % logging activity of each value, for plotting. Num of nuclides.
            end
            tplot=dtplot+tplot;
        end
        % P! end of activity plotting loop
        % QOI:
        if YNQOI==1; % might come up with a better way to code this in, eventually.
            if p==QOItime; % If you wish to find some quantity of interest (activity).
                QOI=L(QOInuclide).*N(QOInuclide)/tconv; % nuclides -> activity
                QOIreturn(1)=dt*tconv; % function will return these values (timestep length in s)
                QOIreturn(2)=QOI; % function will return this value (activity in Bq)
            end
        end
        % P! end of Euler w/plotting loop
    end
else % No plot. Bypassing the plotting loop should speed up the program.
    for p=ts:dt:te;
        N=C*N;
        % QOI
        if YNQOI==1;
            if p==QOItime; % If you wish to find some quantity of interest (activity).
                QOI=L(QOInuclide).*N(QOInuclide)/tconv;
                QOIreturn(1)=dt*tconv; % function will return these values (timestep length in s)
                QOIreturn(2)=QOI; % function will return this value (activity in Bq)
            end
        end
    end
    
end



%% QOI error return
if YNQOIerror==1;
    QOIerror=QOI-QOIanalytic;
    QOIreturn(3)=QOIerror; % function will return these values
end

%% Activity Plotting
% Activity calculation, in case it's relevant.
A=L.*N; % decay constants multiplied by respective number of each nuclide
% conversion factor.
A=A/tconv; % activity in Bq.

% P! Creating the plot.
% Creating an array of the corresponding time values for each point
if YNact==1;
    timeP=zeros(1,tnplot);
    timeP(1,1)=ts;
    for i=1:tnplot % Either years or seconds.
        if tunit=='a'
        timeP(1,i+1) = (ts + i*dtplot);
        else
        timeP(1,i+1) = (ts + i*dtplot)/tconv;
        end
    end
    
    % converting number densities to activity
    Ac=L.*P/tconv;
    
    % Writing file
    fileID=fopen('DecayChain.txt','w');
    fprintf(fileID,'Number Densities of nuclides versus Time.\n');
    fprintf(fileID,'Note: The time step size used to generate the values for this txt file\nis distinct from that used in main file''s computations. That is, they can adjusted separately.\n\n');
    if tunit=='a' % years
    fprintf(fileID,'time (a)    ');
    else % seconds
    fprintf(fileID,'time (s)    ');
    end
    for i=1:numN
        fprintf(fileID,'Nuclide %d  ||  ',i);
    end
   fprintf(fileID,'\n');
    Ptranspose=P';
    Ptimetranspose=timeP';
    for j=1:tnplot;
    fprintf(fileID,'%d   ',Ptimetranspose(j));
    for i=1:numN
        fprintf(fileID,'   %d',Ptranspose(j,i));
    end
    fprintf(fileID,'\n');
    end
    %End writing file
    
    % Plotting
    figure
    hold on
    for i=1:numN;
        plot(timeP,Ac(i,:));
    end
    % Labeling Plot
    
    if tunit=='a' % years
    title('Activity (Bq) vs time (years)');
    xlabel('time (years)');
    else % seconds
    title('Activity (Bq) vs time (seconds)');    
    xlabel('time (seconds)')
    end
    ylabel('Activity (Bq)');
    legend
end

%% Return
varargout={QOIreturn}; % return these variables
%% Stability check
    maxL=max(L);
    globalerror=maxL*dt; % on the order of, not precisely or exactly
    if globalerror>=2;
        fprintf('Warning. Numerical approximation(s) may be unstable.\n');
    end
    %approximate error
end
