function [pc,phi,cp,ct,sec] = patchcontact(pf,pe,varargin)
% [pc,phi,cp,ct] = patchcontact(pf,pe) computes the contacts between the foot
% and the environment patch.
%
%   INPUTS
%
%   The parameters pf and pe are (scalar) Matlab structs describing the foot and
%   environment patches respectively.  See the documentation for the function
%   patchchk() for details.
%
%   OUTPUTS
%
%   The output PC is the output foot contact patch after the translations and
%   rotations of the foot patch pf to contact the environment patch pe.
%
%   The output PHI is the set of the ccw angles (in rad) around the z local axis
%   of the foot patch pf, that are valid for contact.  The 0 rad angle is the 
%   one of the patch's local axis (x-up, y-left, z-out).  The phi Na-by-2 matrix
%   contains the starting (first column) and the ending (second column) angles 
%   in radians.
%
%   The output CP is an Ncp-by-3 matrix with the set of contact points between
%   the foot patch pf and the environment patch pe at the particular pose 
%   configuration.
%
%   The output CT is the type of relation between the contact points cp.  If ct 
%   equals 0, then the set of points are separate points.  If ct equals 1, then
%   the contact points are conected with contact lines between them.  If ct
%   equals 2, then the contacts points form a polygon.
%
%   OPTIONS
%
%   A list of (name,value) pairs may be given to set options.  Unrecognized 
%   names cause warnings.  If a name is given more than once the last-given 
%   (rightmost) value takes precedence.
%
%   Option 'dorot' (default 0): whether to rotate the foot patch.
%
%   Option 'drot' (default 0): whether to draw the rotation angle circle.
%
%   Option 'dcp' (default 1): whether to draw the contact patch.
%
%   Option 'dcpts' (default 1): whether to draw the contact points.
%
%   Option 'da' (default 1): whether to draw contact patch axis (see patchplot).
%
%   Option 'phis' (default 0.1): the step angle for the rotation around the
%   z-axis for the foot patch wrt the environment patch.
%
%   Option 'dbg' (default 0): may have any nonnegative value.  dbg=0 disables
%   debug.  dbg=2 enables the gaphics pause.  Larger values enable successively
%   more verbose debug.
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
dorot = 0; drot = 0; dcp = 1; dcpts = 1; phis = 0.01; da = 1; dbg = 0;
%TBD replace this with option selection

for iopt=1:nopt
  n = varargin{2*iopt-1}; v = varargin{2*iopt};
  if (ischar(n))
    switch (n)
      case 'dorot'; dorot = v; case 'drot'; drot = v; case 'dcp'; dcp = v;
      case 'dcpts'; dcpts = v; case 'phis'; phis = v; case 'da'; da = v;
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

% Find the possible rotations of the patch around its z-axis
[phi] = contactpatchrot(pe,pf);

% Pick a set of angles
phi_steps = [];
if (dorot)
  for i=1:size(phi,1); phi_steps = [phi_steps, phi(i,1):phis:phi(i,2)]; end
else
  phi_steps = phi(1,1);
end

% Find the exact contact patch for the particular set of phi values
pc_a = cell(size(phi_steps,2),1); cp_a = cell(size(phi_steps,2),1);
for i=1:size(phi_steps,2)
  [pc,cp,ct,sec] = contactpatchpts(pf,pe,phi_steps(i));
  pc_a{i}=pc; cp_a{i}=cp;
end

t = toc(tstart);
if (dbg); fprintf('Time for contact search=%g\n',t/size(phi_steps,2)); end
if (nargout>4); sec = t; end

% Plot the contact
contactpatchplot(phi,pc,pe,pc_a,cp_a,phi_steps,'drot',drot,'dcp',dcp,...
                 'dcpts',dcpts,'da',da,'dbg',dbg);
end