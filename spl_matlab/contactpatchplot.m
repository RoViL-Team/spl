function h = contactpatchplot(phi,pc,pe,pc_a,cp_a,phi_steps,varargin)
%h = contactpatchplot(phi,pc,pe,pc_a,cp_a,phi_steps) draws a set of patch
%contact circle arcs in counter-clockwise (ccw) representation.
%
%   The parameter PHI is a Nx2 matrix containing the set of N rotational angles
%   in the local xy-local patch frame for which the contacts are valid.
%
%   The parameter PC is the contact patch.
%
%   The parameter PE is the environment patch.
%
%   The parameter PC_A is a cell array of all contact patches.
%
%   The parameter CP_A is a cell array with all contact points.
%
%   The plot is made in the current axes of the current figure, if any,
%   else a new figure window is opened.  The return value h is the handle
%   of a new Matlab hgtransform group containing the generated graphics.
%
%   If called with hold off then the circle plot replaces the current
%   contents of the current axes, if any, and the axes are reconfigured for
%   the circle.  If called with hold on, the circle is added to the current
%   axes, which are not reconfigured.
%
%   OPTIONS
%
%   A list of name,value pairs may be given to set options. Unrecognized
%   names cause warnings. If a name is given more than once the last-given
%   (rightmost) value takes precedence.
%
%   Option 'fa' (default 0.5): face alpha, a scalar in the interval [0,1].
%
%   Option 'dc' (default 1): whether to draw the rotation angle circle.
%
%   Option 'drot' (default 1): whether to draw the negative angles in the
%   circle in red.
%
%   Option 'dcp' (default 1): whether to draw the contact patch.
%
%   Option 'dcpts' (default 1): whether to draw the contact points.
%
%   Option 'da' (default 1): whether to draw contact patch axis (see patchplot).
%
%   Option 'dbg' (default 0): may have any nonnegative value.  dbg=0 disables
%   debug.  Larger values enable successively more verbose debug.  dbg=1 enables
%   the graphics pause.  dbg=2 enables the gaphics pause and user input to
%   continue.
%
% Copyright (C) 2016- Dimitrios Kanoulas

% process name,value options
nva = length(varargin);
if (mod(nva,2)~=0)
  warning('expect even number of varargs (name,value pairs)');
  nva = nva-1;
end
nopt = nva/2;

% option defaults
fa = 0.5; dc = 1; drot = 1; dcp = 1; dcpts = 1; da = 1; dbg = 0;

for i=1:nopt
  n = varargin{2*i-1}; v = varargin{2*i};
  if (ischar(n))
    switch (n)
      case 'fa'; fa = v; case 'dc'; dc = v; case 'drot'; drot = v;
      case 'dcp'; dcp = v; case 'dcpts'; dcpts = v; case 'da'; da = v;
      case 'dbg'; dbg = v;
      otherwise; warning('unexpected optional arg %s',n);
    end
  else
    warning('non-string, expected name');
  end
end

% Determine parameters
phi_size = size(phi);
if (length(pc.d)==1); ry = pc.d(1); else; ry = pc.d(2); end
rx = pc.d(1); rf = sqrt(rx^2+ry^2); % pc radius
cx = pc.c(1); cy = pc.c(2); cz = pc.c(3); cr = pe.r; r = rf;

% Hold the figure
washold = ishold(); hold('on'); if (~washold); hold('off'); end
h = []; hp = []; hpts = []; % fig handles

% Plot the circle
if (drot && dc)
  t = linspace(0,2*pi,100); plot3(r*sin(t)+cx,r*cos(t)+cy,cz+ones(1,100)*0.01);
  %rotate(t,cr,norm(cr)); set(t,'Rotate',90);
end

% Plot the invalid angles
if (drot)
  h = plot_arc(0.0,phi(1,1),cx,cy,cz,cr,r,'r',fa);
  for i=1:phi_size(1)-1
    h = plot_arc(phi(i,2),phi(i+1,1),cx,cy,cz,cr,r,'r',fa);
  end
  h = plot_arc(phi(end,2),2*pi,cx,cy,cz,cr,r,'r',fa);

  % Plot the valid angles
  for i=1:phi_size(1)
    h = plot_arc(phi(i,1),phi(i,2),cx,cy,cz,cr,r,'g',fa);
  end
end

% Plot the contact patch
if (dbg==2); pause('on'); end
for i=1:size(phi_steps,2)
  tic_t = tic();
  if (dcp) % draw the contact patch
    if (ishandle(hp)); delete(hp); end
    pc_c = pc_a{i}; hold('on');
    hp = patchplot(pc_c,'gr',0,'gc',[0 0 1],'fc','g','fa',0.2,'da',da);
  end
  
  if (dcpts) % draw the contact points
    if (ishandle(hpts)); delete(hpts); end
    cp_c = cp_a{i};
    hpts = plot3(cp_c(:,1),cp_c(:,2),cp_c(:,3),'o','MarkerSize',10,...
                 'MarkerFaceColor',[1 0 0]);
  end
  if (~dbgpause(toc(tic_t),'new rotation\n')); return; end
end

  function h = plot_arc(ai,ae,cx,cy,cz,cr,r,color,fa)
  % h = plot_arc(ai,ae,cx,cy,r) plots a circular arc
  % ai is starting angle of the arc in radians
  % ai is ending angle of the arc in radians 
  % (cx,cy) is the center of the circle
  % r is the radius
  % color is the color of filling
  % fa is the alpha value of the color
    
  t = linspace(ai,ae); x = r*cos(t) + cx; y = r*sin(t) + cy;
  x = [x cx x(1)]; y = [y cy y(1)];
  h = fill3(x,y,cz+ones(size(x))*0.1,color,'FaceAlpha',fa);
  rotate(h,cr,rad2deg(norm(cr)),[cx cy cz]);
  end

  function o = dbgpause(t,varargin)
  % o = dbgpause(t) pause a for-loop run for t sec, used for drwaing multiple 
  % contact patches while rotating the foot patch ccw. 
  o = 1;
  if (dbg==2)
    es = input('Rotate foot patch (q <enter> quits, <enter> continues):','s');
    if (strcmpi(es,'q')); o = 0; else; evalin('caller',es); end
  elseif (dbg==1)
    pause((10/length(phi_steps))-t);
  end
  end

end