function testcontactpatches(tfun, tpfi, tpei, dpf, dpe, ppo)
%testcontactpatches(tfun, tpi, dpe, ppo) runs tests on foot and environment
%patches of different types
%
%   tfun must be the handle of a function expecting two arguments, the foot and
%   environment patch to test, and returning a boolean indicating test success.
%
%   tpfi and tpei is a vector of indices to into an internally generated array
%   of test (foot and environment) patches.  Indices tpfi=1:2 correspond to
%   a foot test patch of each type pr, sc in order.  Indices tpei = 1:15
%   correspond to a test patch of each type ee, he, yr, oc, pc, pe, pr, pq, sc,
%   cr in order.
%
%   dpf/dpe indicates whether to draw the test patches.  If dpf/dpe=1, each test
%   patch will be plotted in its own axes before the corresponding call to tfun.
%   The axes will remain current with hold on so that tfun may add graphics, if
%   desired. dpe=2 sets up axes but does not actually draw the patches.  dpe=0
%   disables graphics.
%
%   ppo is a cell array of name, value pairs to pass to patchplot().
%
%   testpatches(tfun, tpfi, tpei, dpf, dpe) uses default ppo
%
%   testpatches(tfun, tpfi, tpei, dpf) defaults dpe=1
%
%   testpatches(tfun, tpfi, tpei) defaults dpf=1 and dpe=1
%   
%   testpatches(tfun, tpfi) defaults tpei=1:15
%
%   testpatches(tfun) defaults tpfi=1:2 and tpei=1:15
%
%   testpatches() defaults tfun=@(pf,pe)(1)
%
% Copyright (C) 2016- Dimitrios Kanoulas

% adjust these to control which plots are drawn
pf_first = 1; pf_last = 2;
pe_first = 1; pe_last = 15;

maxnc = 4; % maximum number of subplot columns
dss = 0.4; % default sample size

defppo = {'aw',3,'as',1.5,'bw',2,'gw',1};

av = 0; % draw axes
dt = 0; % draw title

% end of configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 1); tfun = @(pf,pe)(1); end
if (nargin < 2); tpfi = pf_first:pf_last; end
if (nargin < 3); tpei = pe_first:pe_last; end
if (nargin < 4); dpf = 0; end
if (nargin < 5); dpe = 1; end
if (nargin < 6); ppo = defppo; end

% 2 primary test foot patches
ipf=1; % planar foot patch
pf(ipf).name = 'plane (aa rect)';
pf(ipf).s = 'p'; pf(ipf).b = 'r'; pf(ipf).d = [1.2 0.5]; pf(ipf).ss = 2*dss;

ipf=ipf+1; % spherical foot patch
pf(ipf).name = 'concave sphere';
pf(ipf).s = 's'; pf(ipf).b = 'c'; pf(ipf).k = 1.3; pf(ipf).d = 1/pf(ipf).k;
pf(ipf).ss = dss/4; pf(ipf).gd = 2;

% 15 primary test environment patches
r = [0 0 pi];
i=1;
pe(i).name = 'convex elliptic paraboloid';
pe(i).s = 'e'; pe(i).b = 'e'; pe(i).k = [-1 -2]; pe(i).d = [1.5 1];
pe(i).ss = 0.5*dss; pe(i).gd = 2;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'concave elliptic paraboloid';
pe(i).s = 'e'; pe(i).b = 'e'; pe(i).k = [1 2]; pe(i).d = [1.5 1];
pe(i).ss = 0.5*dss; pe(i).gd = 2;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'hyperbolic paraboloid';
pe(i).s = 'h'; pe(i).b = 'e'; pe(i).k = [0.5 -1]; pe(i).d = [1.5 1];
pe(i).ss = 0.5*dss; pe(i).gd = 2;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'convex cylindric paraboloid';
pe(i).s = 'y'; pe(i).b = 'r'; pe(i).k = -2; pe(i).d = [2 0.8];
pe(i).ss = 0.5*dss; pe(i).gd = 4;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'concave cylindric paraboloid';
pe(i).s = 'y'; pe(i).b = 'r'; pe(i).k = 2; pe(i).d = [2 0.8];
pe(i).ss = 0.5*dss; pe(i).gd = 4;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'convex circular paraboloid';
pe(i).s = 'o'; pe(i).b = 'c'; pe(i).k = -1; pe(i).d = 1.8; pe(i).ss = 0.5*dss;
pe(i).gd = 4;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'concave circular paraboloid';
pe(i).s = 'o'; pe(i).b = 'c'; pe(i).k = 1; pe(i).d = 1.8; pe(i).ss = 0.5*dss;
pe(i).gd = 4;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'plane (circle)';
pe(i).s = 'p'; pe(i).b = 'c'; pe(i).d = 2; pe(i).ss = dss;
pe(i).gd = 2;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'plane (ellipse)';
pe(i).s = 'p'; pe(i).b = 'e'; pe(i).d = [2.5 1.7]; pe(i).ss = dss;
pe(i).gd = 2;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'plane (aa rect)';
pe(i).s = 'p'; pe(i).b = 'r'; pe(i).d = [2 1.2]; pe(i).ss = 2*dss;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'plane (convex quad)';
pe(i).s = 'p'; pe(i).b = 'q'; pe(i).d = [2.1 1.7 1.8 2.2 pi/4]; pe(i).ss = 2*dss;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'convex sphere';
pe(i).s = 's'; pe(i).b = 'c'; pe(i).k = -1/2; pe(i).d = 1.8; pe(i).ss = dss;
pe(i).gd = 2;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'concave sphere';
pe(i).s = 's'; pe(i).b = 'c'; pe(i).k = 1/2; pe(i).d = 1.8; pe(i).ss = dss;
pe(i).gd = 2;
pe(i).c = [1 2 1]; pe(i).r = r;

i=i+1;
pe(i).name = 'convex circular cylinder';
pe(i).s = 'c'; pe(i).b = 'r'; pe(i).k = -1/2; pe(i).d = [2 1.8];
pe(i).ss = 0.5*dss; pe(i).gd = 4;
pe(i).c = [1 2 1]; pe(i).r = [0 0 pi/2];

i=i+1;
pe(i).name = 'concave circular cylinder';
pe(i).s = 'c'; pe(i).b = 'r'; pe(i).k = 1/2; pe(i).d = [2 1.8];
pe(i).ss = 0.5*dss; pe(i).gd = 4;
pe(i).c = [1 2 1]; pe(i).r = [0 0 pi/2];

n = length(tpei); % number of test patches
nr = ceil(n/maxnc); % num axes rows
if (nr>1); nc = maxnc; else; nc = n; end % num axes cols

% run tests
t = 1; tpfis = length(tpfi);
for j = tpfi
  for i = tpei
    tpe(i) = patchchk(pe(i),'elaborate','all');
    if (dpf>0 || dpe>0)
      subplot(nr,tpfis*nc,t);
      if (dpf==1)
        patchplot(pf(i), ppo{:},'dg', 1, 'gr',0, 'bw',2, 'aw',3, 'as', 1,...
                  'fa', 0.9);
        hold('on');
      end
      if (dpe==1)
        patchplot(tpe(i), ppo{:},'dg', 1, 'gr',0, 'bw',2, 'aw',3, 'as', 1,...
                  'fa', 0.9);
        hold('on');
      end
      if (dt); title(tpe(i).name); end
      if (~av); set(gca(),'Visible','off'); end
    end

    if (~tfun(pf(j),tpe(i))); error('test %d failed (test patch %d)', t, i); end
    
    camzoom(2); t = t+1;
  end
end

end