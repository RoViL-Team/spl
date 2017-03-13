function [ic, sec] = contactpatchcheck(pf, pe, varargin)
% ok = contactpatchcheck(p) checks whether the foot patch pf and the environment
% patch pe can be in contact in the local 2D coordinate system of the pe, i.e. 
% if of is contained into pe's bounadry in their nominal pose.
%
%   INPUTS
%
%   The parameters pf and pe are (scalar) Matlab structs describing the patches.
%   See the documentation for the function patchchk() for details.
%
%   The contact evaluation is described in the paper "Reasoning About Foot
%   Contact using Curved Surface Patches" by D. Kanoulas, P. Kryczka,
%   N. Tsagarakis and Marsette Vona.  In particular, the evaluation takes place
%   in the local xy-plane of the environment patch and checks whether the 
%   borders of the two patches are not intersecting and pf is in pe's boundary.
%
%   OUTPUTS
%
%   The output IC (in contact) represents whether the patches' boundaries are
%   not intersecting and takes values 1 (not intersecting) and 0 (intersecting).
%
%   The output SEC (seconds) provides the total process time.
%
%   OPTIONS
%
%   A list of (name,value) pairs may be given to set options.  Unrecognized 
%   names cause warnings.  If a name is given more than once the last-given 
%   (rightmost) value takes precedence.
%
%   Option 'dbg' (default 0): may have any nonnegative value. dbg=0 disables
%   debug.  Larger values enable successively more verbose debug.
%
% Copyright (C) 2016 Dimitrios Kanoulas

tstart = tic(); % timing

% process (name,value) options
nva = length(varargin);
if (mod(nva,2)~=0)
  warning('expect even number of varargs (name,value pairs)');
  nva = nva-1;
end
nopt = nva/2;

% option defaults
dbg = 0;

for i=1:nopt
  n = varargin{2*i-1}; v = varargin{2*i};
  if (ischar(n))
    switch (n)
      case 'dbg'; dbg = v;
      otherwise; warning('unexpected optional arg %s',n);
    end
  else
    warning('non-string, expected name');
  end
end

% determine parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pf = patchchk(pf,'gb',1); % update foot patch boundary fields
pe = patchchk(pe,'gb',1); % update env patch boundary fields

% read or infer surface and boundary type
if (isfield(pf,'b')); btf = pf.b; end
if (isfield(pe,'b')); bte = pe.b; end

% check contacts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('r',btf) % if foot patch is rect
  if (dbg); fprintf('Foot patch is a rect.\n'); end

  rx = pf.d(1); ry = pf.d(2); % dx and dy radii
  
  % if inside the boundary
  if ((pe.bl(rx,ry)<=0) && (pe.bl(rx,-ry)<=0) &&...
      (pe.bl(-rx,-ry)<=0) && (pe.bl(-rx,ry)<=0))
    ic = 1; % pf boundary is inside pe boundary
  else
    ic = 0; % pf boundary is outside pe boundary
  end
end

if strcmp('c',btf) % if foot patch is circle
if (dbg); fprintf('Foot patch is a circle.\n'); end

  r = pf.d(1); % radius
  
  switch (bte)
    case 'c' % circle
      if (r<pe.d(1)); ic = 1; else; ic = 0; end
    case {'e','r'} % ellipse, axis aligned rectangle
      if (r<min(pe.d)); ic = 1; else; ic = 0; end
    case 'q' % convex quadrilateral
      % quad verts v 4x2 in CCW order starting in quadrant 1
      g = pe.d(5); p = [g; pi-g; pi+g; -g];
      re = [pe.d(1); pe.d(2); pe.d(3); pe.d(4)];
      ve = [cos(p).*re, sin(p).*re];
      
      % dist from origin to line segment
      d1 = abs(det([ve(2,:)-ve(1,:);-ve(1,:)]))/norm(ve(2,:)-ve(1,:));
      d2 = abs(det([ve(3,:)-ve(2,:);-ve(2,:)]))/norm(ve(3,:)-ve(2,:));
      d3 = abs(det([ve(4,:)-ve(3,:);-ve(3,:)]))/norm(ve(4,:)-ve(3,:));
      d4 = abs(det([ve(1,:)-ve(4,:);-ve(4,:)]))/norm(ve(1,:)-ve(4,:));
      if (r<=min([d1, d2, d3, d4])); ic = 1; else; ic = 0; end 
    otherwise
      if (cki); iee('unknown bounds type'); end
  end
end

% generate outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dbg); fprintf('ic=%g\n',ic); end

t = toc(tstart);
if (nargout>1); sec = t; end

end