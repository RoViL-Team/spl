function [ss, hss, ssm, pp, hpp] = demomanpatchcontact(pcd,varargin)
% demomanpatchcontact() runs demospl(), contactpatchcheck(), and patchcontact() in
% manual segmentation mode.
%
% Copyright (C) 2016- Dimitrios Kanoulas

fn = datafn(pcd); sco = 3; % grid organized, enables mesh and bfs

spo = 9; spmax = 3000;
opts = {'sc',1,'scopts',sco,'fit',0,'pfopts',4,...
        'sp',0,'av',0,'dt',0,'spmax',spmax,'spopts',spo,...
        'smopts',3,'ms',1};

fprintf('manseg contact analysis: fitting paraboloids with ellipse boundary\n');

[ss, hss, ssm, pp, hpp] = demospl(fn,opts{:});

% Foot patches
dx = 0.2; dy = 0.1; dss = 0.4;
p.name = 'plane (aa rect)';
p.s = 'p'; p.b = 'r'; p.d = [dx dy]; p.ss = 2*dss;

for i=1:size(pp{1},2)
  pe = pp{1}(i); p.c=pe.c; p.r = pe.r;
  if (contactpatchcheck(p,pe,varargin{:}))
    [pc,~,~,~] = patchcontact(p,pe,'da',0,varargin{:});
  end
end
end