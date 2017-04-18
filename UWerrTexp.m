function [value,dvalue,ddvalue,tauint,dtauint,qval,wopt, gammaFbb, drho] = ...
  UWerrTexp(Data,Parm,Nrep,Name,Quantity,varargin)
%UWERRTEXP Autocorrelation-analysis of MC time-series.
% [value,dvalue] = UWerrTexp(Data) Returns the average of the first column of 
%  Data and the error of the average (taking into account autocorrelations)
%
% following the Gamma-method in
% ``Monte Carlo errors with less errors''
% by Ulli Wolff, hep-lat/0306017
%
% and the TauExp bias correction
%   following the method in
%  ``Critical slowing down and error analysis in lattice QCD simulations''
%           arxiv: 1009.5228
% 
% 
%--------------------------------------------------------------------------
%  Based on : Ulli Wolff,   Nov. 2006, UWerr Version V6
%--------------------------------------------------------------------------
%  Francesco Virotta, Jan. 2011, Version V4
%--------------------------------------------------------------------------
%
% detailed description of the call:
% [value,dvalue,ddvalue,tauint,dtauint,qval,wopt,gammaFbb,drho] = ...
%    UWerrTexp(Data,Parm,Nrep,Name,Quantity,P1,P2,....)
%
% omitted or set to [] input values get default value indicated in [D=..]
%
% Data     -- (N x Nalpha) matrix of measured (equilibrium!) data 
%             N = total number of measurements
%             Nalpha = number of (primary) observables
%             if your data are in a different format, consider 
%             the MATLAB commands `reshape' and `permute'
%
% Parm     -- A vector of (up to 6) values:
%             Parm(1) : guess for the ratio S of tau/tauint 
%             if set = 0, absence of autocorrelations is *assumed*
%             Parm(2) : Tauexp used to add a tail and give an upper bound to the
%             error. If set = 0, performs only 'standard' UWerr analysis
%             Parm(3) : Nsigma, determines where to attach the tail (larger
%             values attach the tail at earlier times, the tail *wont* be
%             attached past the first negative value of the autocorrelation
%             function)
%             Parm(4) : Wsmall, check for fast decay of gamma (if wopt < Wsmall
%             the amplitude for the tail is given by max(rho(wopt+1),2*drho(wopt+1)), 
%             only if this results in a smaller tauexp than the one at wopt).
%             This amplitude is always used if rho(wopt)<0
%             Parm(5) : MDR, Molecular dynamics times R (factor used for
%             rescaling plots)
%             Parm(6) : If ~=0 outputs both upper bound and lower bound
%             values for error, tauint and wopt
%             [D = [1.5 0 3 5 1 0]]
%
% Nrep     -- vector [N1 N2 N3 ...] specifying a breakup of the N rows
%             of Data into replica of length N1,N2,.. (N=sum(Nrep)!)  [D=[N]]
%             The replica distribution is histogrammed and a Q-value 
%             (=probability to find this much or more scatter) is given 
%             for R >= 2 replica; if you have one history you may set Nrep to 
%             artificially split it into >=2 replica to see their distribution                
%
% Name     -- if string: name of observable for titles of generated plots
%             if not string: all plots are supressed [D='NoName']
%
% Quantity -- either:
%             scalar function handle (@functionname) for the derived 
%             observable F; it has to operate on a row-vector of length 
%             Nalpha as first argument; optional parameters P1,P2,... are 
%             passed on to this function as 2nd, 3rd, .... argument
%          -- or:
%             integer between 1 and Nalpha to analyze primary observable [D=1]
%--------------------------------------------------------------------------
% output::
%  value    -- estimate for F
%  dvalue   -- its error (incl. autocorrelation effects). If param(6)=1,
%              returns a vector of 2 values: [lowerbound upperbound]
%  ddvalue  -- statistical error of the error (only lowerbound)
%  tauint   -- integrated autocorrelation time: If param(6)=1,
%              returns a vector of 2 values: [lowerbound upperbound]
%  dtauint  -- statistical error of tauint. If param(6)=1,
%              returns a vector of 2 values: [lowerbound upperbound], upperbound
%              assumes no error on tauexp.
%  qval     -- Q-value of the replica distribution if R >= 2
%              (goodness of fit to a constant)
%  wopt     -- returns the numerical value of the summation window. If
%              param(6)=1, returns a vector of 2 values: [lowerbound upperbound]
%  gammaFbb -- Autocorrelation function (only up to 2*wopt)
%  drho     -- Error on the normalized autocorrelation function

%--------------------------------------------------------------------------
%

mlock  % this locks the function in memory and preserves persistent variables 
       % between calls; munlock('UWerrTexp');clear('UWerrTexp') required 
       % before a clear command 
persistent f % figure handles for plot windows


% analyze input arguments:
[N,Nalpha]=size(Data); % number of measurements and observables
Dparm = [1.5 0 3 5 1 0];

if nargin < 5 || isempty(Quantity), Quantity=1;    end
if nargin < 4 || isempty(Name),     Name='NoName'; end
if nargin < 3 || isempty(Nrep),     Nrep=N;      end
if nargin < 2 || isempty(Parm),     Parm=Dparm;      end

L = length(Parm); DL = length(Dparm);
if L>DL, error('UWerrTexp:IncorrectArgin','Incorrect value of Parm');end

Parm(L+1:DL) = Dparm(L+1:end);

Stau  = Parm(1);
Texp  = Parm(2);
Nsigma= Parm(3);
Wsmall= Parm(4);
MDR   = Parm(5);
BOTH  = (Parm(6)~=0);

if any(Nrep ~= round(Nrep)) || any(Nrep < 1) || sum(Nrep) ~= N
  error('UWerrTexp:IncorrectArgin','inconsistent N,Nrep')
end

if isnumeric(Quantity)
   if round(Quantity) ~= Quantity || Quantity < 1 || Quantity > Nalpha
       error('UWerrTexp:IncorrectArgin','illegal numeric value of Quantity')
   end
   primary = 1; % logical variable
   primaryindex = Quantity;
else
   primary = 0;
   evalQ = @(x) feval(Quantity,x,varargin{:});
end
Nrep = Nrep(:);   % make column
R    = length(Nrep); % number of Replica

%--------------------------------------------------------------------------
% means of primaries:

abb = mean(Data,1);    % total mean, primary obs.
abr = zeros(R,Nalpha); % replicum-mean, primary obs.
i0=1;
for r=1:R
  i1 = i0-1+Nrep(r);
  abr(r,:) = mean(Data(i0:i1,:),1);
  i0 = i0+Nrep(r);
end

% means of derived (or primary, depending on Quantity):
if primary
  Fbb = abb(primaryindex);
  Fbr = abr(:,primaryindex);	
else
  Fbb = evalQ(abb);    % total mean, derived obs.
  Fbr = zeros(R,1);      % replica-means, derived obs.
  for i=1:R
    Fbr(i)=evalQ(abr(i,:));
  end
end
Fb=sum(Fbr.*Nrep)/N;  % weighed mean of replica means

%--------------------------------------------------------------------------
% form the gradient of f and project fluctuations:

if primary
  delpro = Data(:,primaryindex)-abb(primaryindex);
else
  delpro = project(Data,evalQ);
end

%--------------------------------------------------------------------------
% Computation of Gamma, automatic windowing
% note: W=1,2,3,.. in Matlab-vectorindex <-> W=0,1,2,.. in the paper!
% sick case:
Fvar = mean(delpro.^2);
if Fvar == 0
  error('UWerrTexp:NoFluctuations','No fluctuations!')
end

% compute Gamma(t), find the optimal W
if Stau == 0  % no autocorrelations assumed
  wopt = 0;
  tmax = 0;
else
  wopt = 1;
  tmax = floor(min(Nrep)/2);    % do not take t larger than this
end

gammaFbb=computeGamma(delpro,Nrep,tmax);

if wopt % compute wopt
  wopt = findWopt(gammaFbb,tmax , Stau, N);
  tmax = min(tmax,2*wopt);
end

gammaFbb = gammaFbb(1:tmax+1);

CFbbopt = gammaFbb(1) + 2*sum(gammaFbb(2:wopt+1));   % first estimate
if CFbbopt <= 0
  error('UWerrTexp:GammaPathological','Gamma pathological: estimated error^2 < 0')
end

gammaFbb = gammaFbb+CFbbopt/N;                     % bias in Gamma corrected
CFbbopt  = gammaFbb(1) + 2*sum(gammaFbb(2:wopt+1));% refined estimate
sigmaF   = sqrt(CFbbopt/N);                        % error of F
rho      = gammaFbb/gammaFbb(1);                   % normalized autocorr.
rhosum   = cumsum(rho)-0.5;                        % tau_int(t)

%----------------------------------------------------------------------------

% bias cancellation for the mean value

if R >= 2
  bF = (Fb-Fbb)/(R-1);
  Fbb= Fbb - bF;
  if abs(bF) > sigmaF/4
    warning('UWerrTexp:BiasCancel','A %.1f sigma bias of the mean has been cancelled.',bF/sigmaF)
  end
  Fbr = Fbr - bF*N./Nrep;
  Fb  = Fb  - bF*R;
end

%----------------------------------------------------------------------------

% answers to be returned: (besides gammaFbb)
value   = Fbb;
dvalue  = sigmaF;
ddvalue = dvalue*sqrt((wopt+0.5)/N);
tauint  = rhosum(wopt+1);
dtauint = tauint*2*sqrt((wopt-tauint+0.5)/N);
drho    = computeDrho(rho,wopt,N);

if Texp
  [tauintu wupper] = computeUpperBound(rho,drho,Texp,Nsigma,Wsmall);
  CFbboptu   = gammaFbb(1)*2*tauintu;  % upper estimate
  IDX = 1 + BOTH;
  dvalue(IDX)  = sqrt(CFbboptu/N);
  taudraw = [tauint tauintu];
  tauint(IDX)  = tauintu;
  dtauint(IDX) = sqrt(dtauint(1)^2+(drho(wupper+1)*Texp)^2); % assumes error on texp=0
  wdraw = [wopt wupper];
  wopt(IDX)    = wupper;
else
  taudraw = tauint;
  wdraw = wopt;
end

% Q-value for replica distribution if R >= 2:
if R >= 2
  chisq = sum((Fbr-Fb).^2.*Nrep)/CFbbopt;
  qval  = 1-gammainc(chisq/2,(R-1)/2);
else
  qval  = [];
end

%----------------------------------------------------------------------------
% Plotting:
doplot = [ Stau ~= 0, primary ~=0, R >= 2];
% make plots, if demanded:
if ~ ischar(Name)   % no plots
  for n=1:3
    if ~isempty(f) && ishandle(f(n)), close(f(n)); end
  end
  f = [];
  return
elseif isempty(f)
  for n=1:3
    if doplot(n), f(n) = figure; else f(n)=NaN; end
  end
elseif any(~ishandle(f))
  for n=1:3
    if doplot(n) && ~ishandle(f(n)), f(n) = figure;  end
  end
end

if ismatlab % || ~(ismatlab || verLessThan('???'))
  scrsz=get(0,'screensize'); % to place plots on the screen
  position = [0.01*scrsz(3) 0.45*scrsz(4) 0.45*scrsz(3) 0.45*scrsz(4);
              0.51*scrsz(3) 0.52*scrsz(4) 0.45*scrsz(3) 0.35*scrsz(4);
              0.51*scrsz(3) 0.03*scrsz(4) 0.45*scrsz(3) 0.35*scrsz(4)];
  for n=1:3
    if doplot(n)
      set(f(n),'DefaultAxesFontSize',16);
      set(f(n),'Position',position(n,:));
    end
  end
end

% autocorrelation:

if doplot(1)
  figure(f(1));
  clf
  if ismatlab || ~(ismatlab || verLessThan('3.2'))
    set(f(1),'Name','Autocorrelation');
  end
  % plot of GammaF(t)/GammaF(0) = rhoF(t)
  % put the right dimensions to everything
  
  time    = MDR*(0:tmax)';
  w       = wdraw*MDR;
  texp    = Texp*MDR;
  rhosum  = rhosum*MDR;
  drhosum = 2*rhosum.*sqrt( time/(N*MDR) );
  tau     = taudraw*MDR;
  
  subplot(1,2,1)
  errorbar(time,rho,drho);
  v=[0 time(end) min(min(rho-drho)*1.2,0) 1]; axis(v); 
  hold on;
  
  plot([0 time(end)],[0 0],'--');     % zero line
  
  if Texp
    plot([w(1) w(1)], v(3:4),'r--');  % Bar at woptl
    plot([w(2) w(2)], v(3:4),'r--');  % Bar at woptu 
  else
    plot([w w], v(3:4),'r-');         % Bar at wopt
  end
  title(['normalized autocorrelation of ' Name])
  ylabel('\rho')
  hold off;
  
  %----------------------------------------------------------------
  % plot of rhosum
  subplot(1,2,2)
  % tauint vs. W:
  errorbar(time, rhosum, drhosum); 
  v = [0 time(end) min(rhosum-drhosum)*0.8 max(rhosum+drhosum)*1.2];
  
  if Texp
    rhosumupper  = rhosum + max(rho,2*drho)*texp;
    drhosumupper = sqrt(drhosum.^2+(drho*texp).^2);    
    v(4) = max( v(4), max(rhosumupper(2:end)+drhosumupper(2:end))*1.1);
  end
  axis(v);
  hold on
  if Texp
    errorbar(time,rhosumupper,drhosumupper);
    plot([0 time(end)],[tau(1) tau(1)],'r--')  % hor. line at estimate
    plot([0 time(end)],[tau(2) tau(2)],'r--')  % hor. line at estimate
    plot([w(1) w(1)],v(3:4),'r--');            % Bar at wopt
    plot([w(2) w(2)],v(3:4),'r--');            % Bar at wopt
  else
    plot([0 time(end)],[tau tau],'r--')        % hor. line at estimate
    plot([w w],v(3:4),'r-');                   % Bar at wopt
  end
  title(['\tau_{int} with statistical errors of ' Name]);
  ylabel('\tau_{int}')
  xlabel('W')
  hold off
end

% histogram of data, if primary:

if doplot(2)
  figure(f(2));
  clf
  if ismatlab || ~(ismatlab || verLessThan('3.2'))
    set(f(2),'Name','Histogram');
  end
  hist(Data(:,primaryindex),20);
  title(['Distribution of all data for ' Name]);
end

% histogram of replica values if R >= 2:
if doplot(3)
  figure(f(3));
  clf
  p=(Fbr-Fb)./(sigmaF*sqrt(N./Nrep-1));  
  hist(p);
  title(['replica distribution (mean=0,var=1) for ' Name ' >> Q = ' num2str(qval,'%.2g') ' <<'])
end
drawnow
if ~any(ishandle(f)), f=[]; end
end % UWerrTexp
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

function g = autoCF(data)
% Autocorrelation function of data array. 

n   = size(data,1);
siz = 2^nextpow2(3*n);

fdata  = fft(data,siz);
afdata = fdata.*conj(fdata);

g=real(ifft(afdata));
g=g(1:n,:);
end % autoCF
%----------------------------------------------------------------------------

function gamma = computeGamma(delpro, Nr, Tmax)

gamma = zeros(Tmax+1,1);

i0=1;
for r = 1:length(Nr)
  i1=i0-1+Nr(r);
  tmp = autoCF(delpro(i0:i1));
  gamma = gamma + tmp(1:(Tmax+1));
  i0=i0+Nr(r);
end
gamma = gamma(1:Tmax+1)./(length(delpro) - r*(0:Tmax)');
end % computeGamma
%----------------------------------------------------------------------------

function delpro = project(Data,evalQ)

[N, Nalpha]=size(Data);
fgrad = zeros(Nalpha,1);

h  =  std(Data,1,1)/sqrt(N)/2; % spread for num.derivative=sigma_i/sqrt(N)/2
v  =  mean(Data,1);
dv = v;

for a = 1:Nalpha
  if h(a) == 0     % Data for this observable do not fluctuate
    fgrad(a)    = 0;
  else
    dv(a)    = v(a) + h(a);
    fgrad(a) = evalQ(dv);
    dv(a)    = v(a) - h(a);
    fgrad(a) = fgrad(a)-evalQ(dv);
    dv(a)    = v(a);
    fgrad(a) = fgrad(a)/(2*h(a));
  end
end

 % projected deviations (on total statistics)
delpro = Data*fgrad - v*fgrad;
end % project
%--------------------------------------------------------------------------

function w=findWopt(Gamma, Tmax, S,N)

rint = cumsum(Gamma(2:(Tmax+1))/Gamma(1));
tauW = S./(log(abs((rint+1)./rint)));
tauW(rint<=0) = eps;
t = (1:Tmax)';
gW = exp(-t./tauW)-tauW./sqrt(t*N);
w = find(gW<=0,1);
if isempty(w)
  warning('UWerrTexp:FailedCond','Windowing condition failed up to W = %i',Tmax)
  w = Tmax;
end
end% findWopt
%----------------------------------------------------------------------------

function drho = computeDrho(rho, wopt, N)
% changed subleading terms (prev. versions!) here:
% contruct errors acc. to hep-lat/0409106 eq. (E.11)
% pad zeros to simplify summation:

L = length(rho);
paddedsize = 2*L + wopt + 1;
rho((L+1):paddedsize) = 0;
drho = zeros(L,1);
T = 0:L-1;
BEG = max(1,T-wopt);

for t=T; %0:L-1
  k = BEG(t+1):t+wopt;                          % summation range
  drho(t+1)=sum((rho(k+t+1)+rho(abs(k-t)+1)-2*rho(t+1)*rho(k+1)).^2);
end

drho = sqrt(drho/N);
end%computeDrho

function [tau w] = computeUpperBound(Rho, DRho, Texp, Nsigma, Wsmall)
persistent warn
if isempty(warn), warn=1; end

Wmax = length(Rho) - 1;
rhosum = cumsum(Rho)-0.5;

t = find(Rho - Nsigma*DRho<0,1); % t is a matlab vector index ( 1...L)
                                 % stop at first negative occurence
                                 
if isempty(t)
  warning('UWerrTexp:WtooSmall',['Could not meet %d sigma criterion.'...
    ' Setting Wupper = %i (Wmax)'],Nsigma,Wmax);
  t = Wmax + 1;
else
  try % try to get the window closest to Nsigma from zero
    myfun   = @(T,N) abs(abs(Rho(T)) - N*DRho(T));
    [M IDX] = min(myfun([t t-1],Nsigma)); % compare t and t-1, select the smallest
    t = t - (IDX-1);
  catch ERR
    display(ERR.message)
    t=1;
  end
end

if t<=1, t=2; end % always get a window of size at least one

if t <= Wsmall || Rho(t) < 0, % small window case  
  t1=t; t2=t+1;
  tmp1 = rhosum(t1)+max(Rho(t1),2*DRho(t1))*Texp;
  tmp2 = rhosum(t2)+max(Rho(t2),2*DRho(t2))*Texp;
  
  [tauintu IDX] = min([tmp1 tmp2]);
  
  if IDX==2,
    warning('UWerrTexp:Changingwopt2','Estimating tau_int^u at %d, instead of wopt(2) = %d',t,t-1);
    if warn
      warning('OFF','UWerrTexp:Changingwopt2');
      display(' warnings for UWerrTexp:Changingwopt2 are now  turned ''OFF''')
      display(' to turn them back on use:')
      display(' warning(''ON'',''UWerrTexp:Changingwopt2'')')
      warn=0;
    end
    t = t+1;
  end
else
  tauintu = rhosum(t)+Rho(t)*Texp;
end

tau  = tauintu;
w = t-1; % w is in physical notation (0...L-1)
end

function val = ismatlab()
val = ~isequal(license(),'GNU General Public License');
end

function result = verLessThan(verstr)

verParts = getParts(verstr);
curParts = getParts(version);
result = (sign(curParts - verParts) * [1; .1; .01]) < 0;
end

function parts = getParts(V)
parts = sscanf(V, '%d.%d.%d')';
if length(parts) < 3
  parts(3) = 0; % zero-fills to 3 elements
end
end