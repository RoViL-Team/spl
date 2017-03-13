function [pc,cp,ct,sec] = contactpatchpts(pf,pe,phi,varargin)
% [pc,cp,ct,sec] = contactpatchpts(pf,pe) computes the contact points between
% the foot and the environment patches.
%
%   INPUTS
%
%   The parameters pf and pe are (scalar) Matlab structs describing the foot and
%   environment patches respectively.  See the documentation for the function
%   patchchk() for details.
%
%   The parameter phi (in rad) is the ccw rotational angle aroudn the local
%   z-axis for the foot patch pf. 
%
%   OUTPUTS
%
%   The output pc is the output foot contact patch after the translations and
%   rotations of the foot patch pf to contact the environment patch pe.
%
%   The output cp is an Ncp-by-3 matrix with the set of contact points between
%   the foot patch pf and the environment patch pe at the particular
%   pose configuration.
%
%   The output ct is the type of relation between the contact points cp.  If ct 
%   equals 0, then the set of points are separate points.  If ct equals 1, then
%   the contact points are connected with contact lines between them.  If ct
%   equals 2, then the contacts points form a polygon.
%
%   The output sec is the total time complexity
%
%   OPTIONS
%
%   A list of (name,value) pairs may be given to set options.  Unrecognized 
%   names cause warnings.  If a name is given more than once the last-given 
%   (rightmost) value takes precedence.
%
%   Option 'dbg' (default 0): may have any nonnegative value.  dbg=0 disables
%   debug.  Larger values enable successively more verbose debug.
%
% Copyright (C) 2016- Dimitrios Kanoulas

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

for iopt=1:nopt
  n = varargin{2*iopt-1}; v = varargin{2*iopt};
  if (ischar(n))
    switch (n)
      case 'dbg'; dbg = v;
      otherwise; warning('unexpected optional arg %s',n);
    end
  else
    warning('non-string, expected name');
  end
end

% init the contact patch pc as the foot patch and the outputs: note that the
% contact patch is rotated ccw by phi around the local z-axis
pc = pf; pc.c = pe.c; pc.r = rlog(rotz(radtodeg(-phi))); cp = []; ct = 0;

% read or infer surface and boundary type
if (isfield(pf,'b')); btf = pf.b; end;
if (isfield(pe,'s')); ste = pe.s; end;

% Optimizaton function for the hyperbolic paraboloid
function [c,ceq] = boundaryfun(x,ve)
  lx = x(1); ly = x(2);
  bl = blquad([ve(1,1) ve(1,2);
               ve(2,1) ve(2,2);
               ve(3,1) ve(3,2);
               ve(4,1) ve(4,2)]);
  
  c = bl(lx,ly);
  ceq = []; %bl(lx,ly);
end

if strcmp('r',btf) % if foot patch is a planar rect
  if (dbg); fprintf('Foot patch is a planar rect.\n'); end;

  % vertices of the foot rect in the local frame
  rx = pf.d(1); ry = pf.d(2); % dx and dy radii size
  ve1 = [ rx -ry 0]; ve2 = [ rx  ry 0]; % upper right/left
  ve3 = [-rx  ry 0]; ve4 = [-rx -ry 0]; % lower left/right
  vel = [ve1; ve2; ve3; ve4]; % local
  ve = [ve1; ve2; ve3; ve4];

  % rotate locally the foot patch wrt the environment patch
  [ve(:,1),ve(:,2),ve(:,3)] = xform3 (ve(:,1),ve(:,2),ve(:,3),pc.r,pc.c,0,1);

  % update the z-component of the projection of the rotated pf vertices on pe and 
  % find the max z-value
  ve(:,3) = pe.zl(ve(:,1),ve(:,2)); [max_ve_val,max_ve_ind]= max(ve(:,3));
  max_ve_allind = find(abs(ve(:,3) - max_ve_val) <= 10^-4);
  
  % Find the exact points of contact (in local frame first)
  switch (ste) % environment patch type
    case 'e' %elliptic paraboloid
      if (all(pe.k<0)) % convex
        if (dbg); fprintf('PE: convex elliptic paraboloid.\n'); end;
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      else
        if (dbg); fprintf('PE: concave elliptic paraboloid.\n'); end;
        pc.c = [0 0 ve(max_ve_ind,3)]; %the max (in local frame) indexed vertex
        cpz(1:length(max_ve_allind),1) = pc.c(3);
        cp = [vel(max_ve_allind,1), vel(max_ve_allind,2), cpz];
        ct = 0;
      end
    case 'h' %hyperbolic paraboloid
      if (dbg); fprintf('PE: hyperbolic paraboloid.\n'); end;
      ts_opt = tic(); % timing
      
      % fmincon solution (inside the foot patch boundaries)
      fun = @(x)-pe.zl(x(1),x(2));
      x0 = [pe.d(1) 0]; % local not rotated env patch max point 
      A = []; b = []; Aeq=[]; beq=[];
      lb = [min(ve(:,1)) min(ve(:,2))]; ub = [max(ve(:,1)) max(ve(:,2))];
      options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',10,...
                                      'Display', 'off');
      [x_opt,~] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)boundaryfun(x,ve),...
                          options);
       
      pc.c = [0 0 pe.zl(x_opt(1),x_opt(2))];
      [x_opt(1),x_opt(2),~] = xform3 (x_opt(1),x_opt(2),0,pc.r,pc.c,1,1);
      cp = [x_opt(1) x_opt(2) pc.c(3); -x_opt(1) -x_opt(2) pc.c(3);];
      
      te_opt = toc(ts_opt);
      if (dbg); fprintf('Time for contact search=%g\n',te_opt); end

    case 'y' % cylindric paraboloid
      if (dbg); fprintf('PE: cylindric paraboloid.\n'); end;
      if (pe.k<0) % convex
        if (dbg); fprintf('PE: hyperbolic paraboloid (convex y).\n'); end;
        pc.c = [0 0 0];
        theta = atan(ry/rx);
        
        phi_m = mod(phi,pi); % the analysis for the first two quarters
        
        % the intersection between the rect and the xy-axes with
        % trigonometry

        if (phi_m>theta && phi_m<(pi-theta))
          cp = [ ry/tan(phi_m)  ry pc.c(3);
                -ry/tan(phi_m) -ry pc.c(3)];
        else
          cp = [ rx  rx*tan(phi_m) pc.c(3);
                -rx -rx*tan(phi_m) pc.c(3)];
        end
      else
        pc.c = [0 0 ve(max_ve_ind,3)];
        cpz(1:length(max_ve_allind),1) = pc.c(3);
        cp = [vel(max_ve_allind,1), vel(max_ve_allind,2), cpz];
        ct = 0;
      end   
    case 'o' % circular paraboloid
      if (dbg); fprintf('PE: circular paraboloid.\n'); end;
      if (pe.k<0) % convex
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      else
        pc.c = [0 0 ve(max_ve_ind,3)];
        cp = [ve1(1), ve1(2), pc.c(3); ve2(1), ve2(2), pc.c(3);
              ve3(1), ve3(2), pc.c(3); ve4(1), ve4(2), pc.c(3);];
        ct = 0;
      end
    case 'p' % plane
      if (dbg); fprintf('PE: plane.\n'); end;
      pc.c = [0 0 0];
      cp = [ve1; ve2; ve3; ve4]; ct = 2;
    case 's' % sphere
      if (dbg); fprintf('PE: sphere.\n'); end;
      if (pe.k<0) % convex
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      else
        pc.c = [0 0 ve(max_ve_ind,3)];
        cp = [ve1(1), ve1(2), pc.c(3); ve2(1), ve2(2), pc.c(3);
              ve3(1), ve3(2), pc.c(3); ve4(1), ve4(2), pc.c(3);];
        ct = 0;
      end
    case 'c' % circular cylinder
      if (dbg); fprintf('PE: circular cylinder.\n'); end;
      if (pe.k<0) % convex
        pc.c = [0 0 0];
        theta = atan(ry/rx);
        
        phi_m = mod(phi,pi); % the analysis for the first two quarters
        
        % the intersection between the rect and the xy-axes with
        % trigonometry
        if (phi_m>theta && phi_m<(pi-theta))
          cp = [ ry/tan(phi_m)  ry pc.c(3);
                -ry/tan(phi_m) -ry pc.c(3)];
        else
          cp = [ rx  rx*tan(phi_m) pc.c(3);
                -rx -rx*tan(phi_m) pc.c(3)];
        end

        ct = 0;
      else
        pc.c = [0 0 ve(max_ve_ind,3)];
        cpz(1:length(max_ve_allind),1) = pc.c(3);
        cp = [vel(max_ve_allind,1), vel(max_ve_allind,2), cpz];
        ct = 0;
      end
    otherwise
      iee('unknown surface type');
  end
end

if strcmp('c',btf) % if foot patch is circle
  if (dbg); fprintf('Foot patch is a circular sphere.\n'); end
  % foot patch radius
  rf = 1/pf.k(1); phi = 0; %radius and angle due to rotation symmetry
  
  switch (ste)
    case 'e' %elliptic paraboloid
      if (all(pe.k<0)) % convex
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      else % concave
        pc.c = [0 0 0];
        if (pe.k(1)>pe.k(2)); kx = pe.k(1); else; kx = pe.k(2); end
        k = pf.k(1);
        if (k<kx)
          pc.c(3) = ((k-kx)^2)/(2*k*k*kx);
          
          % solve to x and z to find the contact between parabola and circle
          af = rf+pc.c(3);
          zf = real(roots([1; (2/kx)-(2*af); (af^2)-(rf^2)]));
          xf = sqrt(2*abs(zf(1))/kx);
          
          if (pe.k(1)>pe.k(2))
            cp = [xf(1) 0 zf(1); -xf(1) 0 zf(1)];
          else
            cp = [0 xf(1) zf(1); 0 -xf(1) zf(1)];
          end
          ct = 0;
        else
          cp = [pc.c]; ct = 0;
        end
      end
    case 'h' %hyperbolic paraboloid
      pc.c = [0 0 0];
      if (pe.k(1)>pe.k(2)); kx = pe.k(1); else; kx = pe.k(2); end
      k = pf.k(1);
      if (k<kx); pc.c(3) = ((k-kx)^2)/(2*k*k*kx); end
      
      % solve to x and z to find the contact between parabola and circle
      af = rf+pc.c(3);
      zf = real(roots([1; (2/kx)-(2*af); (af^2)-(rf^2)]));
      xf = sqrt(2*abs(zf(1))/kx);
          
      if (pe.k(1)>pe.k(2))
      	cp = [xf(1) 0 zf(1); -xf(1) 0 zf(1)];
      else
        cp = [0 xf(1) zf(1); 0 -xf(1) zf(1)];
      end
      ct = 0;
    case {'y','o'} %cylindric/circular paraboloid
      if (pe.k<0) % convex
        if (dbg); fprintf('PE: convex cylindric/circular paraboloid.\n'); end;
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      else
        if (dbg); fprintf('PE: concave cylindric/circular paraboloid.\n'); end;
        pc.c = [0 0 0];
        kx = pe.k(1); k = pf.k(1);
        if (k<kx); pc.c(3) = ((k-kx)^2)/(2*k*k*kx); end
        
        % solve to x and z to find the contact between parabola and circle
        af = rf+pc.c(3);
        zf = real(roots([1; (2/kx)-(2*af); (af^2)-(rf^2)]));
        xf = sqrt(2*abs(zf(1))/kx);
        cp = [0 xf(1) zf(1); 0 -xf(1) zf(1)];
        ct = 0;
      end
    case 'p'
      pc.c = [0 0 0];
      cp = [pc.c]; ct = 0;
    case 'c' % circular cylinder
      if (pe.k<0) % convex
        if (dbg); fprintf('PE: convex circular cylinder.\n'); end;
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      else
        if (dbg); fprintf('PE: concave circular cylinder.\n'); end;
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      end
    case 's' % sphere
      if (pe.k<0) % convex
        if (dbg); fprintf('PE: convex sphere.\n'); end;
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      else
        if (dbg); fprintf('PE: concave sphere.\n'); end;
        pc.c = [0 0 0];
        cp = [pc.c]; ct = 0;
      end
    otherwise
      fprintf('unknown surface type\n');
  end
end

% transform the contact patch to the world frame
pc.r = rlog(rotz(radtodeg(-phi)));
[pc.c(1),pc.c(2),pc.c(3)] = xform3 (pc.c(1),pc.c(2),pc.c(3),pe.r,pe.c,0,0);
pc.r = rlog((rexp(pe.r))*rexp(pc.r));

% transform the contact points to the world frame
[cp(:,1),cp(:,2),cp(:,3)] = xform3 (cp(:,1),cp(:,2),cp(:,3),pc.r,pe.c,0,0);

t = toc(tstart);
if (dbg); fprintf('Time for contact search=%g\n',t); end
if (nargout>3); sec = t; end

end

function bl = blquad(v)
% generate implicit boundary function for quadrilateral
% v are 4 quad verts in CCW order (4x2)

% edge direction vectors
d = [v(2,:)-v(1,:); v(3,:)-v(2,:); v(4,:)-v(3,:); v(1,:)-v(4,:)];
l = [norm(d(1,:)); norm(d(2,:)); norm(d(3,:)); norm(d(4,:))];
vi = (l>eps); ns = sum(vi);

% TBD for now don't get crazy with degenerate cases
if (ns<3); bl = @(x,y)(sqrt(x.*x+y.*y)); return; end

% make s and d are the start points and unit direction vectors of each
% non-degenerate side
s = zeros(ns,2); dd = zeros(ns,2); j = 1;
for i=1:4
  if (vi(i)); s(j,:) = v(i,:); dd(j,:) = d(i,:)/l(i); j = j+1; end;
end
d = dd;

n = zeros(ns,2); % outward pointing unit normals
for i=1:ns; nn = cross([d(i,:),0],[0,0,1]); n(i,:) = nn(1,1:2); end

c = zeros(ns,1); % perp dist to origin
for i=1:ns; c(i) = -n(i,:)*s(i,:)'; end

  function o = f(x,y)
  sz = size(x); nd = sz(1)*sz(2);
  x = x(:)'; y = y(:)'; % always row vectors
  d = zeros(ns,nd); % distance from each pt (cols) to each line (rows)
  for i=1:ns; d(i,:) = x*n(i,1)+y*n(i,2)+c(i); end
  o = zeros(1,nd);
  % for pts inside, return closest negative distance
  ni = (d<0); vi = all(ni); dd = max(d); o(vi) = dd(vi);
  % for pts outside, return closest positive distance
  d(ni) = inf; dd = min(d); vi = ~vi; o(vi) = dd(vi);
  o = reshape(o,sz);
  end

bl = @f;

end % blquad