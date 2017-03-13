function [ss, hss, ssm, pp, hpp] = demoautopatchcontact(asopts,varargin)
% demoautopatchcontact() runs demospl() in auto segmentation mode with newrock.pcd
% and runs a contact analysis between the foot patches and the fitted
% environment ones.
%
%   Optional arg asopts selects autoseg options, default 5, see demospl.m
%
%   Additional optional args are passed on to demospl, e.g. use 'dp',0 to
%   disable drawing patches for timing.
%
%   Run example: [ss, hss, ssm, pp, hpp] = demoautopatchcontact();
%
% Copyright (C) 2016- Dimitrios Kanoulas

if (nargin<1); asopts = 7; end

fn = datafn('iros2017/terrain4.pcd');

sco = 7; spo = 9; spmax = 3000;

opts = {'sc',1,'scopts',sco,'fit',0,'pfopts',4,...
        'sp',0,'av',0,'dt',0,'df',1,'spmax',spmax,'spopts',spo,...
        'smopts',3,'as',1,'asopts',asopts};

fprintf('automatic segmentation: fitting paraboloids with ellipse boundary\n');

[ss, hss, ssm, pp, hpp] = demospl(fn,opts{:});

% Foot patches
r = 2; d = 1.8; c = 1.2/sqrt(2); dss = 0.4;
p.name = 'plane (aa rect)';
p.s = 'p'; p.b = 'r'; p.d = [0.08 0.05]; p.ss = 2*dss;

tot_res = 0; cgp = 0;
for i=1:size(pp{1},2)
  pe = pp{1}(i); pe = patchchk(pe,'gb',1); p.c=pe.c; p.r = pe.r;
  if (pe.residual<0.003)
    cgp = cgp+1;
  else
    tot_res = tot_res + pe.residual;
  end
  if (contactpatchcheck(p,pe,varargin{:}))
    [pc,phi,cp,~] = patchcontact(p,pe,'da',0,varargin{:});
  end
end

fprintf('Avg residual for contact patches=%g\n',tot_res/(size(pp{1},2)-cgp));
end